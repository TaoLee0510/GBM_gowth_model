# install.packages(c("dplyr","readr","purrr","stringr","ggplot2","mgcv"), Ncpus = 2)
library(dplyr)
library(readr)
library(purrr)
library(stringr)
library(ggplot2)
library(mgcv)

# Helper: read all events CSV files under a simulation's events directory
read_events_dir <- function(dir_events) {
  files <- list.files(dir_events, pattern = "_events\\.csv$", full.names = TRUE)
  if (length(files) == 0) return(NULL)
  df <- files %>% map_dfr(\(f) readr::read_csv(f, show_col_types = FALSE) %>% mutate(.file = basename(f)))
  # Expected columns include at least:
  # step, id, Label, g_alloc, div_event, death_event, risk_div,
  # death_resource, death_random, death_timeout
  df
}

# Core fitting function: Fit P(g) and D(g) from per-hour event logs
#
# Arguments:
#   events_df     : data.frame with events (see columns above)
#   exclude_daily : if TRUE, drop hours where step %% 24 == 0 (daily random-death hours)
#   bins          : number of quantile bins for the binning estimator
#   fit_gam       : if TRUE, also fit GAM curves for P(g) and D(g)
#   fd, rt        : optional references for plotting net rate; pass NA to skip
#   out_prefix    : prefix for output files (CSV/PNG)
#   death_mode    : which death type to treat as D(g) signal.
#                   one of c("resource_only","all","random_only","timeout_only").
#                   Default = "resource_only".
fit_PD_from_events <- function(events_df,
                               exclude_daily = TRUE,
                               bins = 30,
                               fit_gam = TRUE,
                               fd = NA_real_,
                               rt = NA_real_,
                               out_prefix = "PD_fit",
                               death_mode = c("resource_only","all","random_only","timeout_only")) {
  death_mode <- match.arg(death_mode)
  stopifnot(all(c("step","Label","g_alloc","div_event","death_event","risk_div",
                  "death_resource","death_random","death_timeout") %in% names(events_df)))

  # Base table
  df <- events_df %>%
    mutate(
      Label = as.integer(Label),
      g = as.numeric(g_alloc),
      risk_div = as.integer(risk_div),
      death_event = as.integer(death_event),
      death_resource = as.integer(death_resource),
      death_random = as.integer(death_random),
      death_timeout = as.integer(death_timeout)
    ) %>%
    filter(is.finite(g), g >= 0)

  # For random-only death analysis we must NOT drop daily ticks
  local_exclude_daily <- if (death_mode == "random_only") FALSE else isTRUE(exclude_daily)
  if (local_exclude_daily) {
    df <- df %>% filter((step %% 24L) != 0L)
  }

  # ---------- (A) Binning estimator ----------
  # Helper to choose the death indicator according to death_mode
  death_indicator <- function(.df) {
    if (death_mode == "resource_only") {
      as.integer(.df$death_event == 1 & .df$death_resource == 1)
    } else if (death_mode == "random_only") {
      as.integer(.df$death_event == 1 & .df$death_random == 1)
    } else if (death_mode == "timeout_only") {
      as.integer(.df$death_event == 1 & .df$death_timeout == 1)
    } else {
      # "all"
      as.integer(.df$death_event == 1)
    }
  }

  make_bin_table_label <- function(subdf, label_tag) {
    if (nrow(subdf) == 0) return(tibble())
    cutp <- unique(quantile(subdf$g, probs = seq(0, 1, length.out = bins + 1), na.rm = TRUE))
    if (length(cutp) < 3) {
      n_risk <- sum(subdf$risk_div == 1, na.rm = TRUE)
      d_ind  <- death_indicator(subdf)
      return(tibble(
        Label = first(subdf$Label),
        g_mid = median(subdf$g, na.rm = TRUE),
        P_hat = ifelse(n_risk > 0, mean(subdf$div_event[subdf$risk_div == 1] == 1, na.rm = TRUE), NA_real_),
        D_hat = mean(d_ind, na.rm = TRUE),
        n = nrow(subdf),
        label_name = label_tag
      ))
    }
    subdf %>%
      mutate(g_bin = cut(g, breaks = cutp, include.lowest = TRUE, right = TRUE),
             d_ind = death_indicator(cur_data())) %>%
      group_by(g_bin) %>%
      summarise(
        Label = first(Label),
        g_mid = median(g, na.rm = TRUE),
        n = n(),
        n_risk = sum(risk_div == 1, na.rm = TRUE),
        P_hat = ifelse(n_risk > 0, sum(div_event == 1 & risk_div == 1, na.rm = TRUE)/n_risk, NA_real_),
        D_hat = sum(d_ind == 1, na.rm = TRUE)/n,
        .groups = "drop"
      ) %>%
      mutate(label_name = label_tag)
  }

  bins_tumor  <- df %>% filter(Label == 1L) %>% make_bin_table_label("tumor")
  bins_normal <- df %>% filter(Label == 0L) %>% make_bin_table_label("normal")
  bins_all <- bind_rows(bins_tumor, bins_normal)

  if (nrow(bins_all) > 0) {
    readr::write_csv(bins_all, paste0(out_prefix, "_bins.csv"))
  }

  # ---------- (B) GAM smoothing (optional) ----------
  gam_models <- list()
  pred_grid  <- NULL
  if (fit_gam) {
    fit_gam_div <- function(subdf) {
      sdf <- subdf %>% filter(risk_div == 1)
      if (nrow(sdf) < 100) return(NULL)
      gs <- sdf$g
      lo <- quantile(gs, 0.01, na.rm = TRUE); hi <- quantile(gs, 0.99, na.rm = TRUE)
      sdf <- sdf %>% filter(g >= lo, g <= hi)
      response <- sdf$div_event
      if (length(unique(response)) < 2) return(NULL)
      mgcv::gam(response ~ s(g, k = 10), data = sdf, family = binomial(link = "logit"), method = "REML")
    }
    fit_gam_death <- function(subdf) {
      sdf <- subdf
      if (death_mode == "random_only") {
        # Random deaths are generated on daily ticks; keep those hours
        sdf <- sdf %>% filter((step %% 24L) == 0L)
      }
      if (nrow(sdf) < 100) return(NULL)
      gs <- sdf$g
      lo <- quantile(gs, 0.01, na.rm = TRUE); hi <- quantile(gs, 0.99, na.rm = TRUE)
      sdf <- sdf %>% filter(g >= lo, g <= hi)
      # Choose death response by mode
      response <- death_indicator(sdf)
      if (length(unique(response)) < 2) return(NULL)
      mgcv::gam(response ~ s(g, k = 10), data = sdf, family = binomial(link = "logit"), method = "REML")
    }

    df_t <- df %>% filter(Label == 1L)
    df_n <- df %>% filter(Label == 0L)

    gam_models$P_tumor  <- fit_gam_div(df_t)
    gam_models$D_tumor  <- fit_gam_death(df_t)
    gam_models$P_normal <- fit_gam_div(df_n)
    gam_models$D_normal <- fit_gam_death(df_n)

    g_grid <- seq(max(1e-6, quantile(df$g, 0.001, na.rm = TRUE)),
                  quantile(df$g, 0.999, na.rm = TRUE), length.out = 200)
    pred <- function(mod) if (is.null(mod)) rep(NA_real_, length(g_grid)) else as.numeric(predict(mod, newdata = data.frame(g = g_grid), type = "response"))

    pred_grid <- tibble(
      g = g_grid,
      P_tumor  = pred(gam_models$P_tumor),
      D_tumor  = pred(gam_models$D_tumor),
      P_normal = pred(gam_models$P_normal),
      D_normal = pred(gam_models$D_normal)
    )
    readr::write_csv(pred_grid, paste0(out_prefix, "_gam_curves.csv"))
  }

  # ---------- (C) Visualization ----------
  plot_one <- function(label_tag, bins_tbl, P_col, D_col, ggrid_tbl, P_pred, D_pred) {
    p <- ggplot() +
      geom_point(data = bins_tbl, aes(x = g_mid, y = !!sym(P_col)), color = "#377eb8", alpha = 0.8, size = 2) +
      geom_point(data = bins_tbl, aes(x = g_mid, y = !!sym(D_col)), color = "#e41a1c", alpha = 0.8, size = 2) +
      labs(x = "g (per-hour glucose allocation)", y = "rate per hour",
           title = paste0("P(g) & D(g) - ", label_tag)) +
      theme_bw()
    if (!is.null(ggrid_tbl)) {
      if (all(is.finite(ggrid_tbl[[P_pred]]))) p <- p + geom_line(data = ggrid_tbl, aes(x = g, y = !!sym(P_pred)), color = "#377eb8", linewidth = 1)
      if (all(is.finite(ggrid_tbl[[D_pred]]))) p <- p + geom_line(data = ggrid_tbl, aes(x = g, y = !!sym(D_pred)), color = "#e41a1c", linewidth = 1)
    }
    p
  }

  if (nrow(bins_tumor) > 0) {
    p_tum <- plot_one(
      "Tumor",
      bins_tumor, "P_hat", "D_hat",
      if (!is.null(pred_grid)) pred_grid %>% select(g, P_tumor, D_tumor) else NULL,
      "P_tumor", "D_tumor"
    )
    ggsave(paste0(out_prefix, "_tumor_PD.png"), p_tum, width = 6, height = 4, dpi = 150)
  }
  if (nrow(bins_normal) > 0) {
    p_nor <- plot_one(
      "Normal",
      bins_normal, "P_hat", "D_hat",
      if (!is.null(pred_grid)) pred_grid %>% select(g, P_normal, D_normal) else NULL,
      "P_normal", "D_normal"
    )
    ggsave(paste0(out_prefix, "_normal_PD.png"), p_nor, width = 6, height = 4, dpi = 150)
  }

  # (Optional) Net fitness comparison plot using GAM curves (if fd & rt provided)
  if (is.finite(fd) && is.finite(rt) && !is.null(pred_grid)) {
    comp_tbl <- pred_grid %>%
      transmute(
        g,
        net_tumor  = P_tumor - D_tumor,
        net_normal = P_normal - D_normal,
        ref_scale  = rt
      )
    p_net <- ggplot(comp_tbl, aes(x = g)) +
      geom_hline(yintercept = 0, linetype = 2, color = "grey50") +
      geom_line(aes(y = net_tumor,  color = "Tumor  P-D"), linewidth = 1) +
      geom_line(aes(y = net_normal, color = "Normal P-D"), linewidth = 1) +
      scale_color_manual(values = c("Tumor  P-D"="#984ea3","Normal P-D"="#4daf4a")) +
      labs(x = "g", y = "net growth rate per hour", color = "",
           title = "Net growth = P(g) - D(g) (GAM)") +
      theme_bw()
    ggsave(paste0(out_prefix, "_net_PD.png"), p_net, width = 6, height = 4, dpi = 150)
  }

  list(bins = bins_all, gam = gam_models, curves = pred_grid)
}

# Batch: traverse Results/ and fit each supply level separately
fit_PD_from_results <- function(root,
                                exclude_daily = TRUE,
                                bins = 30,
                                fit_gam = TRUE,
                                fd = NA_real_, rt = NA_real_,
                                death_mode = c("resource_only","all","random_only","timeout_only")) {
  death_mode <- match.arg(death_mode)
  supply_dirs <- list.dirs(root, full.names = TRUE, recursive = FALSE)
  if (length(supply_dirs) == 0) stop("No subdirectories found in Results/")
  out <- list()

  for (sup in supply_dirs) {
    sim_dirs <- list.dirs(sup, full.names = TRUE, recursive = FALSE)
    ev_all <- NULL
    for (sd in sim_dirs) {
      evd <- file.path(sd, "events")
      if (dir.exists(evd)) {
        df <- read_events_dir(evd)
        if (!is.null(df)) {
          df$supply_label <- basename(sup)
          df$sim_id <- basename(sd)
          ev_all <- bind_rows(ev_all, df)
        }
      }
    }
    if (is.null(ev_all)) next

    prefix <- file.path(sup, paste0("PD_fit_", basename(sup)))
    dir.create(sup, showWarnings = FALSE, recursive = TRUE)
    message("Fitting P(g), D(g) under supply: ", basename(sup),
            " | death_mode=", death_mode,
            " | exclude_daily=", if (death_mode=="random_only") FALSE else exclude_daily)

    out[[basename(sup)]] <- fit_PD_from_events(
      events_df = ev_all,
      exclude_daily = exclude_daily,
      bins = bins,
      fit_gam = fit_gam,
      fd = fd, rt = rt,
      out_prefix = prefix,
      death_mode = death_mode
    )
  }

  invisible(out)
}

summarize_PD_all <- function(results_root,
                             exclude_daily = TRUE,
                             bins = 30,
                             fit_gam = TRUE,
                             fd = NA_real_, rt = NA_real_,
                             death_mode = c("resource_only","all","random_only","timeout_only"),
                             g_points = 200,
                             overwrite_fit = FALSE) {
  death_mode <- match.arg(death_mode)

  suppressPackageStartupMessages({
    library(dplyr); library(readr); library(purrr)
    library(stringr); library(ggplot2)
  })

  # ---------- helpers ----------
  read_events_dir <- function(dir_events) {
    files <- list.files(dir_events, pattern = "_events\\.csv$", full.names = TRUE)
    if (length(files) == 0) return(NULL)
    files %>% map_dfr(\(f) readr::read_csv(f, show_col_types = FALSE) %>% mutate(.file = basename(f)))
  }

  read_supply_events <- function(supply_dir) {
    sim_dirs <- list.dirs(supply_dir, full.names = TRUE, recursive = FALSE)
    ev_all <- NULL
    for (sd in sim_dirs) {
      evd <- file.path(sd, "events")
      if (dir.exists(evd)) {
        df <- read_events_dir(evd)
        if (!is.null(df)) {
          df$supply_label <- basename(supply_dir)
          df$sim_id <- basename(sd)
          ev_all <- bind_rows(ev_all, df)
        }
      }
    }
    ev_all
  }

  find_latest <- function(supply_dir, what = c("bins","curves")) {
    what <- match.arg(what)
    patt <- if (what == "bins") "^PD_fit_.*_bins\\.csv$" else "^PD_fit_.*_gam_curves\\.csv$"
    cand <- list.files(supply_dir, pattern = patt, full.names = TRUE)
    if (length(cand)==0) return(NULL)
    cand[which.max(file.info(cand)$mtime)]
  }

  # ---------- 1) Fit each supply (ALL merged) ----------
  supply_dirs <- list.dirs(results_root, full.names = TRUE, recursive = FALSE)
  if (length(supply_dirs) == 0) stop("No supply subdirectories under: ", results_root)

  for (sup in supply_dirs) {
    bins_ok <- !is.null(find_latest(sup, "bins"))
    curv_ok <- !is.null(find_latest(sup, "curves"))
    if (overwrite_fit || !(bins_ok && curv_ok)) {
      ev <- read_supply_events(sup)
      if (is.null(ev)) next
      prefix <- file.path(sup, paste0("PD_fit_", basename(sup), "_ALL"))
      message("[FIT] supply=", basename(sup), " | death_mode=", death_mode)
      fit_PD_from_events(
        events_df    = ev,
        exclude_daily = exclude_daily,
        bins         = bins,
        fit_gam      = fit_gam,
        fd = fd, rt = rt,
        out_prefix   = prefix,
        death_mode   = death_mode
      )
    }
  }

  # ---------- 2) Summarize bins / curves ----------
  bins_all <- purrr::map_dfr(supply_dirs, function(sup) {
    f <- find_latest(sup, "bins")
    if (is.null(f)) return(NULL)
    readr::read_csv(f, show_col_types = FALSE) %>% mutate(supply = basename(sup))
  })
  curv_all_raw <- purrr::map_dfr(supply_dirs, function(sup) {
    f <- find_latest(sup, "curves")
    if (is.null(f)) return(NULL)
    readr::read_csv(f, show_col_types = FALSE) %>% mutate(supply = basename(sup))
  })

  # Write summarized CSV
  if (nrow(bins_all) > 0) readr::write_csv(bins_all, file.path(results_root, "ALL_bins.csv"))
  if (nrow(curv_all_raw) > 0) readr::write_csv(curv_all_raw, file.path(results_root, "ALL_curves_raw.csv"))

  # Standardize g grid for overlay comparison
  curves_all <- NULL
  if (nrow(curv_all_raw) > 0) {
    g_min <- quantile(curv_all_raw$g, 0.01, na.rm = TRUE)
    g_max <- quantile(curv_all_raw$g, 0.99, na.rm = TRUE)
    g_grid <- seq(g_min, g_max, length.out = g_points)
    interp <- function(x, y, xout) {
      keep <- is.finite(x) & is.finite(y)
      x <- x[keep]; y <- y[keep]
      if (length(unique(x)) < 2) return(rep(NA_real_, length(xout)))
      approx(x = x, y = y, xout = xout, ties = "ordered")$y
    }
    curves_all <- purrr::map_dfr(split(curv_all_raw, curv_all_raw$supply), function(df) {
      tibble(
        supply = unique(df$supply),
        g = g_grid,
        P_tumor  = interp(df$g, df$P_tumor,  g_grid),
        D_tumor  = interp(df$g, df$D_tumor,  g_grid),
        P_normal = interp(df$g, df$P_normal, g_grid),
        D_normal = interp(df$g, df$D_normal, g_grid)
      )
    })
    readr::write_csv(curves_all, file.path(results_root, "ALL_curves_interp.csv"))
  }

  # Merge manifest metadata (if present)
  manifest_path <- file.path(results_root, "manifest.csv")
  manifest <- if (file.exists(manifest_path)) {
    readr::read_csv(manifest_path, show_col_types = FALSE) %>%
      mutate(supply = supply_label) %>%
      select(supply, supply_mode, supply_value, supply_ratio)
  } else NULL

  if (!is.null(manifest) && nrow(bins_all) > 0) {
    bins_all <- bins_all %>% left_join(manifest, by = "supply")
    readr::write_csv(bins_all, file.path(results_root, "ALL_bins_w_manifest.csv"))
  }
  if (!is.null(manifest) && !is.null(curves_all)) {
    curves_all <- curves_all %>% left_join(manifest, by = "supply")
    readr::write_csv(curves_all, file.path(results_root, "ALL_curves_interp_w_manifest.csv"))
  }

  # ---------- 3) Plotting ----------
  out_png <- function(name) file.path(results_root, name)

  # 3.1 Overlay binned lines (Tumor / Normal)
  if (nrow(bins_all) > 0) {
    bins_t <- bins_all %>% filter(Label == 1L) %>% group_by(supply) %>% arrange(g_mid, .by_group = TRUE) %>% ungroup()
    bins_n <- bins_all %>% filter(Label == 0L) %>% group_by(supply) %>% arrange(g_mid, .by_group = TRUE) %>% ungroup()

    p_bins_t <- ggplot(bins_t, aes(x = g_mid, color = supply, group = supply)) +
      geom_line(aes(y = P_hat)) +
      geom_line(aes(y = D_hat), linetype = 2) +
      labs(x = "g", y = "rate per hour", color = "supply", title = "Tumor: P(g) (solid) & D(g) (dashed) — binned") +
      theme_bw()
    ggsave(out_png("ALL_bins_tumor_lines.png"), p_bins_t, width = 7, height = 5, dpi = 150)

    p_bins_n <- ggplot(bins_n, aes(x = g_mid, color = supply, group = supply)) +
      geom_line(aes(y = P_hat)) +
      geom_line(aes(y = D_hat), linetype = 2) +
      labs(x = "g", y = "rate per hour", color = "supply", title = "Normal: P(g) (solid) & D(g) (dashed) — binned") +
      theme_bw()
    ggsave(out_png("ALL_bins_normal_lines.png"), p_bins_n, width = 7, height = 5, dpi = 150)
  }

  # 3.2 Overlay GAM curves (Tumor / Normal)
  if (!is.null(curves_all) && nrow(curves_all) > 0) {
    p_curv_t <- ggplot(curves_all, aes(x = g, color = supply)) +
      geom_line(aes(y = P_tumor)) +
      geom_line(aes(y = D_tumor), linetype = 2) +
      labs(x = "g", y = "rate per hour", color = "supply",
           title = "Tumor: P(g) (solid) & D(g) (dashed) — GAM") +
      theme_bw()
    ggsave(out_png("ALL_curves_tumor_overlay.png"), p_curv_t, width = 7, height = 5, dpi = 150)

    p_curv_n <- ggplot(curves_all, aes(x = g, color = supply)) +
      geom_line(aes(y = P_normal)) +
      geom_line(aes(y = D_normal), linetype = 2) +
      labs(x = "g", y = "rate per hour", color = "supply",
           title = "Normal: P(g) (solid) & D(g) (dashed) — GAM") +
      theme_bw()
    ggsave(out_png("ALL_curves_normal_overlay.png"), p_curv_n, width = 7, height = 5, dpi = 150)

    # 3.3 Net growth (P-D) overlay (Tumor / Normal)
    net_tbl <- curves_all %>%
      transmute(supply, g,
                net_tumor  = P_tumor - D_tumor,
                net_normal = P_normal - D_normal,
                supply_ratio = supply_ratio %||% NA_real_)
    p_net <- ggplot(net_tbl, aes(x = g, color = supply)) +
      geom_line(aes(y = net_tumor)) +
      geom_line(aes(y = net_normal), linetype = 2) +
      labs(x = "g", y = "net rate per hour", color = "supply",
           title = "Net growth = P(g) - D(g) (Tumor solid / Normal dashed)") +
      theme_bw()
    ggsave(out_png("ALL_curves_net_PD_overlay.png"), p_net, width = 7, height = 5, dpi = 150)
  }

  message("Done. Outputs are in: ", results_root,
          "\n - ALL_bins.csv / ALL_bins_w_manifest.csv",
          "\n - ALL_curves_raw.csv / ALL_curves_interp.csv / ALL_curves_interp_w_manifest.csv",
          "\n - ALL_*_overlay.png / ALL_*_lines.png")
}
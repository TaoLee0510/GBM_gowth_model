# ======================= Diversity & Time-Series Analysis =======================
# Read the "daily snapshot (every 24 steps)" csv and output:
# 1) Diversity (Shannon/Simpson/Entropy) over time (merging multiple simulations):
#    For each supply, show: boxplot at each time point + mean trend line
# 2) G distribution over time for a single simulation: violin plot (one plot per simulation)
# 3) Karyotype proportions over time for a single simulation: stacked area plot ("stream" plot; one plot per simulation)

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(stringr); library(purrr)
  library(tidyr); library(ggplot2); library(forcats); library(arrow); library(ragg)
})

# ---- Helpers ----
.parse_day <- function(fname) {
  m <- stringr::str_match(basename(fname), "_day(\\d+)\\.(parquet|csv)$")
  as.integer(m[,2])
}

.kt_cols <- function() paste0("K", 1:22)

.kt_str <- function(df_kt_cols) {
  # df_kt_cols: data.frame with columns K1..K22
  as.character(do.call(paste, c(df_kt_cols, list(sep = "."))))
}

# List daily snapshot files (prefer Parquet over legacy CSV)
.list_daily_files <- function(dir_path) {
  fpq <- list.files(dir_path, pattern = "_day\\d+\\.parquet$", full.names = TRUE)
  if (length(fpq) > 0) return(fpq[order(.parse_day(fpq))])
  fcsv <- list.files(dir_path, pattern = "_day\\d+\\.csv$", full.names = TRUE)
  fcsv[order(.parse_day(fcsv))]
}

# Read a daily snapshot regardless of format
.read_daily_file <- function(f) {
  if (grepl("\\.parquet$", f, ignore.case = TRUE)) {
    as_tibble(arrow::read_parquet(f))
  } else {
    readr::read_csv(f, show_col_types = FALSE, progress = FALSE, col_types = cols())
  }
}

# Fixed color map for karyotypes (from karyotype library CSV)
.read_karyolib <- function(path) {
  stopifnot(file.exists(path))
  kdf <- read.csv(path, stringsAsFactors = FALSE)
  if (ncol(kdf) >= 2) names(kdf)[1:2] <- c("karyotype","fitness")
  kdf
}


.diversity_from_counts <- function(counts) {
  # counts: numeric vector of counts per karyotype
  n <- sum(counts)
  if (n <= 0 || length(counts) == 0)
    return(tibble(H = NA_real_, Simpson = NA_real_, Entropy = NA_real_, Richness = 0L))
  p <- counts / n
  H <- -sum(p * log(p))       # Shannon (nat log)
  D <- 1 - sum(p * p)         # Simpson diversity (1 - lambda)
  tibble(H = H, Simpson = D, Entropy = H, Richness = length(counts))
}

#
# ---- (1) Compute diversity over time for a single simulation (one csv directory) ----
compute_diversity_for_csv_dir <- function(csv_dir) {
  files <- .list_daily_files(csv_dir)
  if (length(files) == 0) return(tibble())

  purrr::map_dfr(files, function(f) {
    day <- .parse_day(f)
    df  <- .read_daily_file(f)
    if (!all(.kt_cols() %in% names(df))) return(tibble())
    kt  <- .kt_str(df[, .kt_cols(), drop = FALSE])
    counts <- as.integer(table(kt))
    mets <- .diversity_from_counts(counts)
    tibble(day = day, H = mets$H, Simpson = mets$Simpson, Entropy = mets$Entropy, Richness = mets$Richness)
  })
}

#
# ---- (1a) Aggregate diversity over all simulations under Results/ (group by supply) ----
collect_diversity_across_results <- function(results_root) {
  supply_dirs <- list.dirs(results_root, full.names = TRUE, recursive = FALSE)
  purrr::map_dfr(supply_dirs, function(sup) {
    sim_dirs <- list.dirs(sup, full.names = TRUE, recursive = FALSE)
    purrr::map_dfr(sim_dirs, function(sim) {
      csv_dir <- file.path(sim, "csv")
      if (!dir.exists(csv_dir)) return(tibble())
      dd <- compute_diversity_for_csv_dir(csv_dir)
      if (nrow(dd) == 0) return(tibble())
      dd %>% mutate(supply = basename(sup), sim_id = basename(sim))
    })
  })
}

#
# ---- (1b-1) Boxplot + mean line (faceted by supply) ----
plot_diversity_boxplots <- function(div_all, out_path = NULL) {
    if (nrow(div_all) == 0) return(invisible(NULL))
   # Only keep Entropy and bin by 30 days
  long <- div_all %>%
    select(day, supply, Entropy) %>%
    mutate(
      day_bin = floor(day / 30),
      day_bin_label = paste0(day_bin * 30, "-", day_bin * 30 + 29)
    ) %>%
    rename(value = Entropy) %>%
    mutate(
      # Ensure chronological order of bins on x-axis within each supply panel
      day_bin_label = factor(day_bin_label,
                             levels = unique(day_bin_label[order(day_bin)]))
    )

  p <- ggplot(long, aes(x = day_bin_label, y = value, group = day_bin_label)) +
    geom_boxplot(outlier.size = 0.7, width = 0.6, alpha = 0.85) +
    stat_summary(fun = mean, geom = "line", aes(group = 1), linewidth = 1) +
      facet_grid(supply ~ ., scales = "free_y") +
      coord_cartesian(ylim = c(0, 1.5)) +
    labs(x = "Day bin (30 days)", y = "Entropy (Shannon)",
         title = "Entropy over time across replicates (30-day bins)") +
    theme_bw()+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  if (!is.null(out_path)) ggsave(out_path, p, width = 12, height = 8, device = cairo_pdf)
  p
}
# ---- (1b-2) Boxplot + mean line for Simpson index (faceted by supply, 30-day bins) ----
plot_diversity_boxplots_simpson <- function(div_all, out_path = NULL) {
  if (nrow(div_all) == 0) return(invisible(NULL))
  # Keep Simpson and bin by 30 days
  long <- div_all %>%
    select(day, supply, Simpson) %>%
    mutate(
      day_bin = floor(day / 30),
      day_bin_label = paste0(day_bin * 30, "-", day_bin * 30 + 29)
    ) %>%
    rename(value = Simpson) %>%
    mutate(
      day_bin_label = factor(day_bin_label,
                             levels = unique(day_bin_label[order(day_bin)]))
    )

  p <- ggplot(long, aes(x = day_bin_label, y = value, group = day_bin_label)) +
    geom_boxplot(outlier.size = 0.7, width = 0.6, alpha = 0.85) +
    stat_summary(fun = mean, geom = "line", aes(group = 1), linewidth = 1) +
    facet_grid(supply ~ ., scales = "free_y") +
    coord_cartesian(ylim = c(0, 1.5)) +
    labs(x = "Day bin (30 days)", y = "Simpson diversity (1 - λ)",
         title = "Simpson diversity over time across replicates (30-day bins)") +
    theme_bw()+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  if (!is.null(out_path)) ggsave(out_path, p, width = 12, height = 8, device = cairo_pdf)
  p
}
#
# ---- (2) G distribution over time for a single simulation: violin plot ----
plot_G_violin_for_sim <- function(csv_dir, out_path, sample_cells = 5000) {
  files <- .list_daily_files(csv_dir)
  if (length(files) == 0) return(invisible(NULL))
  g_tbl <- purrr::map_dfr(files, function(f) {
    day <- .parse_day(f)
    df <- .read_daily_file(f)
    if (!all(c("G","Label") %in% names(df))) return(tibble())
    df <- dplyr::filter(df, Label == 1L)
    if (nrow(df) == 0) return(tibble())
    if (nrow(df) > sample_cells) df <- df[sample.int(nrow(df), sample_cells), , drop = FALSE]
    tibble(day = day, G = df$G)
  })
  if (nrow(g_tbl) == 0) return(invisible(NULL))

  # --- NEW: aggregate into 30-day bins ---
  g_tbl <- g_tbl %>%
    mutate(
      day_bin = floor(day / 30),
      day_bin_label = paste0(day_bin * 30, "-", day_bin * 30 + 29)
    ) %>%
    mutate(
      # sort by time
      day_bin_label = factor(day_bin_label, levels = unique(day_bin_label))
    )

  # Compute per-bin means for annotation
  means_tbl <- g_tbl %>%
    group_by(day_bin_label) %>%
    summarise(mean_G = mean(G, na.rm = TRUE), .groups = "drop")

  p <- ggplot(g_tbl, aes(x = day_bin_label, y = G, group = day_bin_label)) +
    # violin as the main distribution shape
    geom_violin(scale = "width", trim = FALSE) +
    # overlay a narrow boxplot on top of the violin
    geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.7) +
    # mark the mean with a point
    geom_point(data = means_tbl, aes(x = day_bin_label, y = mean_G),
               shape = 23, size = 2.2, fill = "white", inherit.aes = FALSE) +
    # annotate the numeric mean (rounded)
    #geom_text(data = means_tbl, aes(x = day_bin_label, y = mean_G, label = sprintf("%.2f", mean_G)),
    #         vjust = -0.9, size = 2.8, inherit.aes = FALSE) +
    labs(x = "Day bin (30 days)", y = "G",
         title = "Distribution of G over time (30-day bins, per simulation)") +
    theme_bw()+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(out_path, p, width = 10, height = 5)
  p
}

# ---- (2b) Tumor/normal scatter per day with fixed karyotype colors ----
plot_karyotype_colored_scatter_for_sim <- function(csv_dir, out_dir) {
  files <- .list_daily_files(csv_dir)
  if (length(files) == 0) return(invisible(NULL))
  files <- files[order(.parse_day(files))]
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  # Build karyotype color palette
  all_kt_levels <- compute_all_karyotypes_for_sim(csv_dir)
  kcols <- .build_karyotype_palette_distinct(all_kt_levels)
  col_norm <- "#D3D3D3"  # light gray for normal

  for (f in files) {
    day <- .parse_day(f)
    df <- .read_daily_file(f)
    if (!all(c("X","Y","Label", .kt_cols()) %in% names(df))) next
    df$Label <- as.integer(df$Label) # Ensure Label is integer (0/1)

    # karyotype string per cell (coerce to integer to avoid '2.0' formatting)
    kt_df <- df[, .kt_cols(), drop = FALSE]
    kt_df[] <- lapply(kt_df, function(x) as.integer(round(as.numeric(x))))
    kt  <- .kt_str(kt_df)

    # Extend palette to include any new tumor karyotypes appearing today
    lab_is_tumor <- (df$Label == 1L)
      lab_is_tumor[is.na(lab_is_tumor)] <- FALSE

    kt_tumor_today <- unique(kt[lab_is_tumor])
    kcols <- .ensure_palette_levels(kcols, kt_tumor_today)

    kt_chr <- as.character(kt)
    col_vec <- rep(col_norm, length(kt_chr))
    col_vec[lab_is_tumor] <- kcols[kt_chr[lab_is_tumor]]

    na_t <- is.na(col_vec) & lab_is_tumor
    if (any(na_t)) col_vec[na_t] <- "#000000" 

    # Bounds
    xmax <- max(df$X, na.rm = TRUE); ymax <- max(df$Y, na.rm = TRUE)

    # Save PNG under png_karyo
    out_png <- file.path(out_dir, sprintf("Cells_karyo_day%03d.png", day))
    png(filename = out_png, width = 1000, height = 1000, units = "px")
    op <- par(mar = c(0,0,0,0), pty = "s")
    plot(NA, xlim = c(1, xmax), ylim = c(1, ymax), axes = FALSE, xlab = "", ylab = "", asp = 1)
    points(df$X, df$Y, pch = 16, cex = 1.5, col = col_vec)
    mtext(sprintf("Tumor karyotypes (colored) / normals (gray)  Day:%d", day),
          side = 3, line = -1.5, adj = 0, cex = 1.2)
    par(op); dev.off()
  }
  invisible(TRUE)
}

#
# ---- (3) Karyotype proportions over time for a single simulation: stacked area plot (merge all but Top-N as Other) ----
compute_karyotype_proportions <- function(csv_dir, top_n = 12) {
  files <- .list_daily_files(csv_dir)
  if (length(files) == 0) return(tibble())
  files <- files[order(.parse_day(files))]

  by_day <- purrr::map_dfr(files, function(f) {
    day <- .parse_day(f)
    df  <- .read_daily_file(f)
    if (!all(c(.kt_cols(), "Label") %in% names(df))) return(tibble())
    df <- df %>% filter(Label == 1L)  # tumor only
    if (nrow(df) == 0) return(tibble())
    kt  <- .kt_str(df[, .kt_cols(), drop = FALSE])
    tibble(day = day, kt = kt)
  })
  if (nrow(by_day) == 0) return(tibble())

  # Select the Top-N most frequent karyotypes across the entire simulation
  top_levels <- by_day %>% count(kt, sort = TRUE) %>% head(top_n) %>% pull(kt)

  by_day %>%
    mutate(kt2 = if_else(kt %in% top_levels, kt, "Other")) %>%
    count(day, kt2, name = "n") %>%
    group_by(day) %>% mutate(prop = n / sum(n)) %>% ungroup()
}

plot_karyotype_stream_for_sim <- function(csv_dir, out_path) {
  prop_tbl <- compute_karyotype_union99_props(csv_dir)
  if (nrow(prop_tbl) == 0) return(invisible(NULL))

  # sort karyotypes by numeric order (e.g., "1.2.3" < "1.2.10")
  prop_tbl <- prop_tbl %>%
    mutate(kt = factor(kt, levels = .kt_numeric_levels(levels(kt))))

  p <- ggplot(prop_tbl, aes(x = day, y = prop, fill = kt)) +
    geom_area(position = "fill", alpha = 0.95) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(x = "Day", y = "Proportion",
         title = "Tumor-only karyotype composition over time (union of daily Top-k covering ≥99%)") +
    theme_bw() +
    theme(legend.position = "none")

  # Use a fixed palette for karyotypes
  levs <- levels(prop_tbl$kt)
  all_kt_levels <- compute_all_karyotypes_for_sim(csv_dir)
  pal_all <- .build_karyotype_palette_distinct(all_kt_levels)
  pal_all <- .ensure_palette_levels(pal_all, levs)   
  pal_use <- pal_all[levs]
  if (any(is.na(pal_use))) {
    miss_n <- sum(is.na(pal_use))
    pal_use[is.na(pal_use)] <- grDevices::hcl(
      h = seq(0, 360, length.out = miss_n + 1)[-1], c = 80, l = 60
    )
  }
  p <- p + scale_fill_manual(values = pal_use, drop = FALSE)
  ggsave(out_path, p, width = 10, height = 6)
  p
}

# ---- (3b) Karyotype counts over time (tumor-only), symmetric funnel plot ----
.kt_numeric_levels <- function(kts) {
  sp <- strsplit(kts, "\\.")
  mat <- do.call(rbind, lapply(sp, function(v) as.integer(v)))
  ord <- do.call(order, as.data.frame(mat))
  unique(kts[ord])
}

# gather karyotype and sort by numeric order
compute_all_karyotypes_for_sim <- function(csv_dir) {
  files <- .list_daily_files(csv_dir)
  if (length(files) == 0) return(character(0))
  files <- files[order(.parse_day(files))]
  levs <- purrr::map(files, function(f) {
    df <- .read_daily_file(f)
    if (!all(c(.kt_cols(), "Label") %in% names(df))) return(character(0))
    df$Label <- as.integer(df$Label) # Ensure Label is integer (0/1)
    df <- dplyr::filter(df, Label == 1L)
    if (nrow(df) == 0) return(character(0))
    kt_df <- df[, .kt_cols(), drop = FALSE]
    kt_df[] <- lapply(kt_df, function(x) as.integer(round(as.numeric(x))))
    .kt_str(kt_df) |> unique()
  }) |> unlist(use.names = FALSE) |> unique()
  if (length(levs) > 0) levs <- .kt_numeric_levels(levs)
  levs
}

# HCL colors for distinct karyotypes
.build_karyotype_palette_distinct <- function(levels) {
# Return a set of high-contrast, visually friendly HEX colors while avoiding near-black/near-white.
# Strategy:
# 1) For a small number of classes, prefer the Okabe–Ito colorblind-safe palette.
# 2) For larger class counts, sample evenly on the HCL color wheel and alternate L/C to increase separability.
# 3) Filter out near-black/near-white; if still insufficient, augment with a phase-shifted HCL wheel.

  levs <- unique(as.character(levels))
  levs <- levs[!is.na(levs) & nzchar(levs)]
  n <- length(levs)
  if (n == 0) return(setNames(character(0), character(0)))
  if (n == 1) return(setNames("#1f77b4", levs))

  # 1) Okabe–Ito (colorblind-safe) + a few extensions
  okabe <- c(
    "#E69F00", "#56B4E9", "#009E73", "#F0E442",
    "#0072B2", "#D55E00", "#CC79A7", "#999999"
  )
  ext8 <- c("#A6761D", "#1B9E77", "#D95F02", "#7570B3")
  base_pal <- c(okabe, ext8)
  if (n <= length(base_pal)) {
    cols <- base_pal[seq_len(n)]
    names(cols) <- levs
    return(cols)
  }

  # 2) Many classes: HCL hue wheel with alternating L/C
  make_hcl_wheel <- function(k) {
    h <- seq(0, 360, length.out = k + 1)[-1]
    L_pattern <- c(70, 60, 80)  
    C_pattern <- c(70, 60, 65)  
    L <- rep(L_pattern, length.out = k)
    C <- rep(C_pattern, length.out = k)
    grDevices::hcl(h = h, c = C, l = L)
  }

  # Build an oversized candidate pool + phase-shifted wheel to reduce adjacency similarity
  k_pool <- max(n, 72L)
  wheel1 <- make_hcl_wheel(k_pool)
  wheel2 <- make_hcl_wheel(k_pool)
  wheel2 <- wheel2[c(((k_pool %/% 2) + 1):k_pool, 1:(k_pool %/% 2))]
  pool <- unique(c(base_pal, as.vector(rbind(wheel1, wheel2))))

  # Filter near-black / near-white colors
  is_bad <- function(hex) {
    rgb <- grDevices::col2rgb(hex) / 255
    lum <- 0.2126 * rgb[1,] + 0.7152 * rgb[2,] + 0.0722 * rgb[3,]
    lum < 0.22 | lum > 0.92
  }
  pool <- pool[!is_bad(pool)]

  # If still insufficient, relax slightly and add more candidates
  if (length(pool) < n) {
    more <- make_hcl_wheel(n * 3)
    more <- more[!is_bad(more)]
    pool <- unique(c(pool, more))
  }

  cols <- pool[seq_len(n)]
  names(cols) <- levs
  cols
}


.ensure_palette_levels <- function(pal, new_levels) {
  new_levels <- unique(as.character(new_levels))
  miss <- setdiff(new_levels, names(pal))
  if (length(miss) == 0) return(pal)
  c(pal, .build_karyotype_palette_distinct(miss))
}

compute_karyotype_counts_union99 <- function(csv_dir) {
  files <- .list_daily_files(csv_dir)
  if (length(files) == 0) return(tibble())
  files <- files[order(.parse_day(files))]

  # Collect counts per day for tumor-only
  by_day <- purrr::map_dfr(files, function(f) {
    day <- .parse_day(f)
    df <- .read_daily_file(f)
    if (!all(c(.kt_cols(), "Label") %in% names(df))) return(tibble())
    df$Label <- as.integer(df$Label) 
    df <- dplyr::filter(df, Label == 1L)
    if (nrow(df) == 0) return(tibble())
    kt  <- .kt_str(df[, .kt_cols(), drop = FALSE])
    tibble(day = day, kt = kt)
  })
  if (nrow(by_day) == 0) return(tibble())

  # Per-day select minimal Top-k whose cumulative proportion reaches ≥ 99%, then take union across days
  daily_top <- by_day %>%
    count(day, kt, name = "n") %>%
    group_by(day) %>%
    mutate(prop = n / sum(n)) %>%
    arrange(desc(prop), .by_group = TRUE) %>%
    mutate(cumprop = cumsum(prop)) %>%
    # keep all rows strictly below 0.99, plus the first row that crosses 0.99
    filter(cumprop < 0.99 | (cumprop >= 0.99 & dplyr::lag(cumprop, default = 0) < 0.99)) %>%
    ungroup()
  sel_levels <- unique(daily_top$kt)
  if (length(sel_levels) == 0) return(tibble())

  # Numeric-lexicographic order of karyotype strings
  sel_levels <- .kt_numeric_levels(sel_levels)

  # Build counts table restricted to selected union levels
  counts_tbl <- by_day %>%
    filter(kt %in% sel_levels) %>%
    count(day, kt, name = "count") %>%
    group_by(day) %>%
    mutate(total = sum(count)) %>%
    ungroup() %>%
    mutate(kt = factor(kt, levels = sel_levels))
  counts_tbl
}

compute_karyotype_union99_props <- function(csv_dir) {
  counts_tbl <- compute_karyotype_counts_union99(csv_dir)
  if (nrow(counts_tbl) == 0) return(tibble())
  counts_tbl %>%
    group_by(day) %>%
    mutate(prop = count / total) %>%
    ungroup()
}

plot_karyotype_counts_funnel_for_sim <- function(csv_dir, out_path) {
  counts_tbl <- compute_karyotype_counts_union99(csv_dir)
  if (nrow(counts_tbl) == 0) return(invisible(NULL))

  # For each day, compute symmetric bounds around zero for each ordered karyotype band
  counts_ordered <- counts_tbl %>% arrange(day, kt)
  bands <- counts_ordered %>%
    group_by(day) %>%
    mutate(
      cum_prev = dplyr::lag(cumsum(count), default = 0),
      center   = total / 2,
      ymin     = cum_prev - center,
      ymax     = cum_prev + count - center
    ) %>%
    ungroup()

  p <- ggplot(bands, aes(x = day, ymin = ymin, ymax = ymax, fill = kt)) +
    geom_ribbon(alpha = 0.95, color = NA) +
    labs(x = "Day", y = "Count (centered)", fill = "Karyotype",
         title = "Tumor-only karyotype counts over time (symmetric funnel of selected union)") +
    theme_bw() +
    theme(legend.position = "none")
    levs <- levels(bands$kt)
    all_kt_levels <- compute_all_karyotypes_for_sim(csv_dir)
    pal_all <- .build_karyotype_palette_distinct(all_kt_levels)
    pal_all <- .ensure_palette_levels(pal_all, levs)   # ensure palette covers plot levels
    pal_use <- pal_all[levs]
    if (any(is.na(pal_use))) {
      miss_n <- sum(is.na(pal_use))
      pal_use[is.na(pal_use)] <- grDevices::hcl(
        h = seq(0, 360, length.out = miss_n + 1)[-1], c = 80, l = 60
      )
    }
  p <- p + scale_fill_manual(values = pal_use, drop = FALSE)
  ggsave(out_path, p, width = 10, height = 6)
  p
}


#
# ---- Orchestrators ----
# 1) Merge all supplies/simulations: diversity boxplot + mean line
run_diversity_analysis_all <- function(base_output) {
  results_root <- file.path(base_output, "Results")
  div_all <- collect_diversity_across_results(results_root)
  out_pdf <- file.path(results_root, "diversity_boxplots_by_supply_Entroy.pdf")
  plot_diversity_boxplots(div_all, out_pdf)
  out_pdf_simpson <- file.path(results_root, "diversity_boxplots_by_supply_Simpson.pdf")
  plot_diversity_boxplots_simpson(div_all, out_pdf_simpson)
  invisible(div_all)
}

#
# 2) Single simulation plots: G violin plot & karyotype stacked area plot (one plot per simulation)
run_per_sim_plots <- function(base_output, sample_cells = 5000) {
  results_root <- file.path(base_output, "Results")
  supply_dirs <- list.dirs(results_root, full.names = TRUE, recursive = FALSE)
  purrr::walk(supply_dirs, function(sup) {
    sim_dirs <- list.dirs(sup, full.names = TRUE, recursive = FALSE)
    purrr::walk(sim_dirs, function(sim) {
      csv_dir <- file.path(sim, "csv")
      if (!dir.exists(csv_dir)) return(NULL)
      out_dir <- file.path(sim, "analysis")
      dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

      # (2) G violin
      plot_G_violin_for_sim(csv_dir, file.path(out_dir, "G_violin_over_time.pdf"),
                            sample_cells = sample_cells)
      # (2b) Daily scatter with fixed karyotype colors → png_karyo
      out_karyo_png <- file.path(sim, "png_karyo")
      plot_karyotype_colored_scatter_for_sim(csv_dir, out_karyo_png)
      # (3a) Karyotype stacked area using union of daily Top-k covering ≥99% and fixed colors
      plot_karyotype_stream_for_sim(csv_dir,
      file.path(out_dir, "Karyotype_stream_union99.pdf"))
      # (3b) Karyotype counts symmetric funnel (tumor-only, union of >99% share per day)
      plot_karyotype_counts_funnel_for_sim(csv_dir,
        file.path(out_dir, "Karyotype_counts_funnel_union99.pdf"))
    })
  })
  invisible(TRUE)
}


# ---- High-level orchestrator (single-call entrypoint) ----
# Runs both analyses by default:
#  - Across-replicate diversity boxplots with mean line
#  - Per-simulation G violin and karyotype stacked area plots
run_diversity_timeseries <- function(base_output,
                                     sample_cells = 5000,
                                     do_diversity = TRUE,
                                     do_per_sim  = TRUE) {
  stopifnot(is.character(base_output), length(base_output) == 1)
  if (isTRUE(do_diversity)) {
    run_diversity_analysis_all(base_output)
  }
  if (isTRUE(do_per_sim)) {
  run_per_sim_plots(base_output,
                    sample_cells = sample_cells)
  }
  invisible(TRUE)
}

# ---- Example usage (uncomment to run directly) ----
# base_output <- "/Users/4482173/Documents/Project/GBM_Model"
# # 1) Merge multiple simulations: Diversity boxplot (faceted by supply)
# div_all <- run_diversity_analysis_all(base_output)
# # 2) Per-simulation: G violin plot + karyotype proportion stacked plot
# run_per_sim_plots(base_output, sample_cells = 5000)
# # Or manually for a single simulation csv directory:
# csv_dir <- file.path(base_output, "Results", "0.25", "Sim_004", "csv")
# compute_diversity_for_csv_dir(csv_dir)
# plot_G_violin_for_sim(csv_dir, "G_violin_over_time.pdf")
# plot_karyotype_stream_for_sim(csv_dir, "Karyotype_stream_top12.pdf")
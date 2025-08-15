# ======================= Diversity & Time-Series Analysis =======================
# Read the "daily snapshot (every 24 steps)" csv and output:
# 1) Diversity (Shannon/Simpson/Entropy) over time (merging multiple simulations):
#    For each supply, show: boxplot at each time point + mean trend line
# 2) G distribution over time for a single simulation: violin plot (one plot per simulation)
# 3) Karyotype proportions over time for a single simulation: stacked area plot ("stream" plot; one plot per simulation)

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(stringr); library(purrr)
  library(tidyr); library(ggplot2); library(forcats)
})

# ---- Helpers ----
.parse_day <- function(fname) {
  m <- stringr::str_match(basename(fname), "_day(\\d+)\\.csv$")
  as.integer(m[,2])
}

.kt_cols <- function() paste0("K", 1:22)

.kt_str <- function(df_kt_cols) {
  # df_kt_cols: data.frame with columns K1..K22
  do.call(paste, c(df_kt_cols, list(sep = ".")))
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
  files <- list.files(csv_dir, pattern = "_day\\d+\\.csv$", full.names = TRUE)
  if (length(files) == 0) return(tibble())
  files <- files[order(.parse_day(files))]

  purrr::map_dfr(files, function(f) {
    day <- .parse_day(f)
    df  <- readr::read_csv(f, show_col_types = FALSE, progress = FALSE, col_types = cols())
    if (!all(.kt_cols() %in% names(df)))
      return(tibble(day = day, H = NA_real_, Simpson = NA_real_, Entropy = NA_real_, Richness = NA_integer_))
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
# ---- (1b) Boxplot + mean line (faceted by supply) ----
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
  geom_boxplot(outlier_size = 0.7, width = 0.6, alpha = 0.85) +
  stat_summary(fun = mean, geom = "line", aes(group = 1), linewidth = 1) +
  facet_grid(. ~ supply, scales = "free_y") +
  labs(x = "Day bin (30 days)", y = "Entropy (Shannon)",
       title = "Entropy over time across replicates (30-day bins)") +
  theme_bw()
if (!is.null(out_path)) ggsave(out_path, p, width = 12, height = 8, device = cairo_pdf)
p
}

#
# ---- (2) G distribution over time for a single simulation: violin plot ----
plot_G_violin_for_sim <- function(csv_dir, out_path, sample_cells = 5000) {
  files <- list.files(csv_dir, pattern = "_day\\d+\\.csv$", full.names = TRUE)
  if (length(files) == 0) return(invisible(NULL))
  files <- files[order(.parse_day(files))]

  g_tbl <- purrr::map_dfr(files, function(f) {
    day <- .parse_day(f)
    df <- readr::read_csv(f, show_col_types = FALSE, progress = FALSE,
                          col_types = cols(G = col_double()))
    if (!"G" %in% names(df)) return(tibble())
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

  p <- ggplot(g_tbl, aes(x = day_bin_label, y = G, group = day_bin_label)) +
    geom_violin(scale = "width", trim = FALSE) +
    stat_summary(fun = median, geom = "point", size = 1.2) +
    labs(x = "Day bin (30 days)", y = "G",
         title = "Distribution of G over time (30-day bins, per simulation)") +
    theme_bw()
  ggsave(out_path, p, width = 10, height = 5)
  p
}

#
# ---- (3) Karyotype proportions over time for a single simulation: stacked area plot (merge all but Top-N as Other) ----
compute_karyotype_proportions <- function(csv_dir, top_n = 12) {
  files <- list.files(csv_dir, pattern = "_day\\d+\\.csv$", full.names = TRUE)
  if (length(files) == 0) return(tibble())
  files <- files[order(.parse_day(files))]

  by_day <- purrr::map_dfr(files, function(f) {
    day <- .parse_day(f)
    df  <- readr::read_csv(f, show_col_types = FALSE, progress = FALSE, col_types = cols())
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

plot_karyotype_stream_for_sim <- function(csv_dir, out_path, top_n = 12) {
  prop_tbl <- compute_karyotype_proportions(csv_dir, top_n = top_n)
  if (nrow(prop_tbl) == 0) return(invisible(NULL))

  prop_tbl <- prop_tbl %>% mutate(kt2 = forcats::fct_reorder2(kt2, day, prop, .fun = last))
  p <- ggplot(prop_tbl, aes(x = day, y = prop, fill = kt2)) +
    geom_area(position = "fill", alpha = 0.95) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(x = "Day", y = "Proportion", fill = "Karyotype",
        title = sprintf("Tumor-only karyotype composition over time (Top-%d)", length(unique(prop_tbl$kt2)))) +
    theme_bw()
  ggsave(out_path, p, width = 10, height = 6)
  p
}

#
# ---- Orchestrators ----
# 1) Merge all supplies/simulations: diversity boxplot + mean line
run_diversity_analysis_all <- function(base_output) {
  results_root <- file.path(base_output, "Results")
  div_all <- collect_diversity_across_results(results_root)
  out_pdf <- file.path(results_root, "diversity_boxplots_by_supply.pdf")
  plot_diversity_boxplots(div_all, out_pdf)
  invisible(div_all)
}

#
# 2) Single simulation plots: G violin plot & karyotype stacked area plot (one plot per simulation)
run_per_sim_plots <- function(base_output, top_n_karyotypes = 12, sample_cells = 5000) {
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
      # (3) Karyotype stacked area
      plot_karyotype_stream_for_sim(csv_dir,
        file.path(out_dir, sprintf("Karyotype_stream_top%d.pdf", top_n_karyotypes)),
        top_n = top_n_karyotypes)
    })
  })
  invisible(TRUE)
}


# ---- High-level orchestrator (single-call entrypoint) ----
# Runs both analyses by default:
#  - Across-replicate diversity boxplots with mean line
#  - Per-simulation G violin and karyotype stacked area plots
run_diversity_timeseries <- function(base_output,
                                     top_n_karyotypes = 12,
                                     sample_cells = 5000,
                                     do_diversity = TRUE,
                                     do_per_sim  = TRUE) {
  stopifnot(is.character(base_output), length(base_output) == 1)
  if (isTRUE(do_diversity)) {
    run_diversity_analysis_all(base_output)
  }
  if (isTRUE(do_per_sim)) {
    run_per_sim_plots(base_output,
                      top_n_karyotypes = top_n_karyotypes,
                      sample_cells = sample_cells)
  }
  invisible(TRUE)
}

# ---- Example usage (uncomment to run directly) ----
# base_output <- "/Users/4482173/Documents/Project/GBM_Model"
# # 1) Merge multiple simulations: Diversity boxplot (faceted by supply)
# div_all <- run_diversity_analysis_all(base_output)
# # 2) Per-simulation: G violin plot + karyotype proportion stacked plot
# run_per_sim_plots(base_output, top_n_karyotypes = 12, sample_cells = 5000)
# # Or manually for a single simulation csv directory:
# csv_dir <- file.path(base_output, "Results", "0.25", "Sim_004", "csv")
# compute_diversity_for_csv_dir(csv_dir)
# plot_G_violin_for_sim(csv_dir, "G_violin_over_time.pdf")
# plot_karyotype_stream_for_sim(csv_dir, "Karyotype_stream_top12.pdf")
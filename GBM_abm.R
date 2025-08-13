# GBM Agent-Based Model Simulation on a 2-D Lattice
# -------------------------------------------------
# This script implements an agent-based model (ABM) for glioblastoma (GBM)
# on an N×N grid. It now uses:
#   • Stable IDs in the grid (no row-index remapping on death).
#   • Bulk growth-rate computation via calc_rates_all_cpp (Rcpp).

library(yaml)
library(dplyr)
library(tibble)
library(stringr)
library(Rcpp)
library(future.apply)


# Avoid serializing external pointers and reloading on workers
if (isFALSE(getOption("GBM_ABM_CPP_LOADED", FALSE))) {
  Rcpp::sourceCpp("src/abm_kernels.cpp", rebuild = FALSE, verbose = FALSE)
  options(GBM_ABM_CPP_LOADED = TRUE)
}

# Small helper: null coalescing operator
"%||%" <- function(a,b) if (is.null(a)) b else a
# ------------------------------- I/O ------------------------------------

read_config <- function(path) {
  cfg_raw <- yaml::read_yaml(path)
  cfg_raw$N  <- as.numeric(cfg_raw$Nd)
  cfg_raw$Ctd  <- as.numeric(cfg_raw$Ctd)
  cfg_raw$Rt   <- as.numeric(cfg_raw$Rt)
  cfg_raw$Rn   <- as.numeric(cfg_raw$Rn)
  cfg_raw$N_division_possibility <- as.numeric(cfg_raw$N_division_possibility)
  cfg_raw$Cg  <- as.numeric(cfg_raw$Cg)
  cfg_raw$m   <- as.numeric(cfg_raw$m)
  cfg_raw$N_death_rate <- as.numeric(cfg_raw$N_death_rate)
  cfg_raw$T_death_rate <- as.numeric(cfg_raw$T_death_rate)
  cfg_raw$MSR <- as.numeric(cfg_raw$MSR)
  cfg_raw$radius <- as.numeric(cfg_raw$radius)
  cfg_raw$CooldownFactor <- as.numeric(cfg_raw$CooldownFactor %||% 5)
  # Parse total_supply as vector: accept numeric, YAML seq, or comma-separated string
  ts_raw <- cfg_raw$total_supply
  if (is.null(ts_raw)) {
    ts_vec <- cfg_raw$N * cfg_raw$N
  } else if (is.numeric(ts_raw) && length(ts_raw) > 1) {
    ts_vec <- as.numeric(ts_raw)
  } else if (is.character(ts_raw) && length(ts_raw) == 1 && grepl(",", ts_raw)) {
    parts <- strsplit(ts_raw, ",")[[1]]
    ts_vec <- as.numeric(trimws(parts))
  } else {
    ts_vec <- as.numeric(ts_raw)
  }
  ts_vec <- ts_vec * (cfg_raw$N * cfg_raw$N)
  cfg_raw$total_supply_vec <- ts_vec
  # For single-run functions, default to the first value
  cfg_raw$total_supply <- as.numeric(ts_vec[1])
  cfg_raw$DTr <- as.numeric(cfg_raw$DTr)
  cfg_raw$DNr <- as.numeric(cfg_raw$DNr)
  if (is.na(cfg_raw$total_supply)) cfg_raw$total_supply <- cfg_raw$N * cfg_raw$N
  if (is.na(cfg_raw$DTr)) cfg_raw$DTr <- 0.2
  if (is.na(cfg_raw$DNr)) cfg_raw$DNr <- 0.2
  # Parse supply_mode (default "proportional", allow "equal")
  mode <- tolower(as.character(cfg_raw$supply_mode %||% "proportional"))
  if (!mode %in% c("proportional","equal")) mode <- "proportional"
  cfg_raw$supply_mode <- mode
  # Replicates from YAML (default 1)
  cfg_raw$replicates <- as.integer(cfg_raw$replicates %||% 1L)
  if (is.na(cfg_raw$replicates) || cfg_raw$replicates < 1L) cfg_raw$replicates <- 1L
  # Simulation horizon (steps / hours)
  cfg_raw$Time <- as.integer(cfg_raw$Time %||% 2400L)  # default 2400 hours (~100 days)
  if (is.na(cfg_raw$Time) || cfg_raw$Time < 1L) cfg_raw$Time <- 240L
  return(cfg_raw)
}

read_karyolib <- function(path) {
  kdf <- read.csv(path, stringsAsFactors = FALSE)
  names(kdf)[1:2] <- c("karyotype", "fitness")
  kdf
}

kt_vec2str <- function(v) paste(v, collapse = ".")
fitness_of <- function(karyolib, v) {
  s <- kt_vec2str(v)
  f <- karyolib$fitness[match(s, karyolib$karyotype)]
  if (is.na(f)) NA else f
}

# -------------------------- Initialisation ------------------------------

init_cells <- function(cfg, karyolib) {
  N   <- cfg$N
  m   <- cfg$m
  n_g <- N * N
  grid <- matrix(NA_integer_, nrow = N, ncol = N)  # stores stable IDs

  diploid <- rep(2L, 22)
  fit_dip <- fitness_of(karyolib, diploid)

  # row constructor
  make_row <- function(id, x, y, kt_vec, label, r_base, f, mr) {
    tibble(
      id = id,
      X = x, Y = y,
      !!!setNames(as.list(kt_vec), paste0("K", 1:22)),
      f  = f,
      r  = r_base,
      G  = cfg$Cg,            # tumour overwritten after f set
      Label = label,          # 0=normal, 1=tumour
      DivisionTime = NA_real_,
      Time = 0,
      Status = 1L,
      Mr  = mr,
      TimeToDivide = 0,
      division = ifelse(label == 1L, 1L, 0L),
      Cooldown = 0
    )
  }

  speed_ratio <- cfg$Ctd
  next_id <- 1L

  # seed tumour at centre
  centre <- c(ceiling(N/2), ceiling(N/2))
  Cells  <- make_row(next_id, centre[1], centre[2], diploid, 1L, cfg$Rt, fit_dip, cfg$Mtr)
  grid[centre[1], centre[2]] <- next_id
  Cells$G[1] <- speed_ratio * fit_dip

  # seed normals
  n_norm <- floor(n_g * m)
  if (n_norm > 0) {
    empty  <- which(is.na(grid), arr.ind = TRUE)
    sel    <- empty[sample(seq_len(nrow(empty)), min(n_norm, nrow(empty))), , drop = FALSE]
    for (i in seq_len(nrow(sel))) {
      next_id <- next_id + 1L
      row <- make_row(next_id, sel[i,1], sel[i,2], diploid, 0L, 1, 1, cfg$Mnr)
      row$r <- runif(1, 0.95*cfg$Rn, 1.05*cfg$Rn)
      Cells <- bind_rows(Cells, row)
      grid[sel[i,1], sel[i,2]] <- next_id
    }
  }

  # Build per-id arrays (length = max_id); we will refresh per step
  max_id <- max(Cells$id)
  list(Cells = Cells, grid = grid, max_id = max_id)
}

# ----------------------- Neighbourhood & Rates --------------------------

# (kept for occasional single queries; uses per-id maps inside)
get_neighbour_counts <- function(Cells, grid, idx, radius, label_by_id, alive_by_id) {
  res <- neighbour_counts_xy(grid, label_by_id, alive_by_id, Cells$X[idx], Cells$Y[idx], radius)
  list(Nt = as.integer(res$Nt), Nn = as.integer(res$Nn))
}

# Build per-id maps quickly
build_id_maps <- function(Cells) {
  max_id <- max(Cells$id)
  x_by_id <- integer(max_id); y_by_id <- integer(max_id)
  label_by_id <- integer(max_id); label_by_id[] <- -1L
  alive_by_id <- integer(max_id)
  f_by_id <- numeric(max_id)

  alive_rows <- which(Cells$Status == 1L)
  ids <- Cells$id[alive_rows]
  x_by_id[ids] <- Cells$X[alive_rows]
  y_by_id[ids] <- Cells$Y[alive_rows]
  label_by_id[ids] <- as.integer(Cells$Label[alive_rows])
  alive_by_id[ids] <- 1L
  f_by_id[ids] <- Cells$f[alive_rows]

  list(max_id = max_id,
       active_ids = ids,
       x_by_id = x_by_id,
       y_by_id = y_by_id,
       label_by_id = label_by_id,
       alive_by_id = alive_by_id,
       f_by_id = f_by_id)
}

# ------------------------------- Division -------------------------------

divide_cell <- function(Cells, grid, idx, cfg, karyolib, id_state) {
  if (Cells$Status[idx] == 0L) return(list(Cells = Cells, grid = grid, id_state = id_state))

  # Mother info
  mother_id <- Cells$id[idx]
  pos <- c(Cells$X[idx], Cells$Y[idx])

  # Allocate new stable id
  new_id <- id_state$max_id + 1L
  # Ask C++ to pick a free neighbour and place the new id into grid
  dp <- divide_place_cpp(grid, pos[1], pos[2], new_id, cfg$N)
  if (!isTRUE(dp$success)) {
    return(list(Cells = Cells, grid = grid, id_state = id_state))
  }
  grid <- dp$grid
  target_x <- as.integer(dp$x); target_y <- as.integer(dp$y)

  # Prepare daughters by copying mother row
  mother    <- Cells[idx, ]
  daughters <- bind_rows(mother, mother)
  daughters$Time <- 0L

  # MS events through C++
  mother_kt <- as.integer(mother[paste0("K", 1:22)])
  ms_res <- apply_ms_events_cpp(
    mother_kt = mother_kt,
    label     = as.integer(Cells$Label[idx]),
    MSR       = as.numeric(cfg$MSR),
    karyo_lib_str = karyolib$karyotype,
    karyo_lib_fit = karyolib$fitness,
    max_attempts  = 200L
  )
  kt1 <- as.integer(ms_res$kt1); kt2 <- as.integer(ms_res$kt2)
  f1  <- as.numeric(ms_res$f1);  f2  <- as.numeric(ms_res$f2)

  speed_ratio <- cfg$Ctd
  for (j in 1:22) {
    daughters[1, paste0("K", j)] <- kt1[j]
    daughters[2, paste0("K", j)] <- kt2[j]
  }
  daughters$f[1] <- f1; daughters$f[2] <- f2
  daughters$G[1] <- speed_ratio * f1
  daughters$G[2] <- speed_ratio * f2

  # Overwrite mother in-place (keeps mother_id)
  Cells[idx, ] <- daughters[1, ]

  # Append daughter #2 as a new row with new stable id & coordinates
  d2 <- daughters[2, ]
  d2$id <- new_id
  d2$X  <- target_x
  d2$Y  <- target_y
  # Normal cells: set division=0 after division
  if (d2$Label == 0L) d2$division <- 0L
  Cells <- bind_rows(Cells, d2)

  # Update id state
  id_state$max_id <- new_id

  # DivisionTime will be set after bulk r update
  list(Cells = Cells, grid = grid, id_state = id_state)
}

# -------------------------------- Death ---------------------------------

death_update <- function(Cells, grid, cfg, step) {
  # (1) Resource threshold death
  total_supply <- if (!is.null(cfg$total_supply)) cfg$total_supply else (cfg$N * cfg$N)
  total_G      <- sum(Cells$G[Cells$Status == 1L])
  if (total_G <= 0) return(list(Cells = Cells, grid = grid))

  avg_supply <- total_supply / total_G
  death_thr  <- 0.2 * avg_supply
  dead_idx   <- which(Cells$Status == 1L & Cells$G < death_thr)

  # (2) Random daily death
  if (step %% 24L == 0L) {
    u <- runif(nrow(Cells))
    rand_dead <- which(Cells$Status == 1L &
                         ((Cells$Label == 1L & u <= cfg$T_death_rate) |
                          (Cells$Label == 0L & u <= cfg$N_death_rate)))
    dead_idx <- unique(c(dead_idx, rand_dead))
  }

  # (3) Division-timeout death
  div_delay_dead <- which(Cells$Status == 1L &
                            Cells$division == 1L &
                            abs(Cells$TimeToDivide) > 5 * Cells$DivisionTime)
  dead_idx <- unique(c(dead_idx, div_delay_dead))

  # Apply deaths: mark status=0 and clear grid positions
  if (length(dead_idx) > 0) {
    for (i in dead_idx) {
      if (!is.na(Cells$X[i]) && !is.na(Cells$Y[i])) {
        grid[Cells$X[i], Cells$Y[i]] <- NA_integer_
      }
    }
    Cells$Status[dead_idx] <- 0L
  }
  list(Cells = Cells, grid = grid)
}

# ------------------------------- Migration ------------------------------

migrate_cells <- function(Cells, grid, cfg) {
  alive_rows <- which(Cells$Status == 1L)
  for (i in alive_rows) {
    if (runif(1) < Cells$Mr[i]) {
      res_cpp <- migrate_one_cpp(grid, Cells$X[i], Cells$Y[i], Cells$id[i], cfg$N)
      if (isTRUE(res_cpp$success)) {
        grid <- res_cpp$grid
        Cells$X[i] <- as.integer(res_cpp$x)
        Cells$Y[i] <- as.integer(res_cpp$y)
      }
    }
  }
  list(Cells = Cells, grid = grid)
}

# ----------------------- Normal division activation ---------------------

activate_normal_division <- function(Cells, grid, cfg, maps) {
  norm_idxs <- which(Cells$Status == 1L & Cells$Label == 0L & Cells$division == 0L)
  if (length(norm_idxs) == 0) return(Cells)
  for (i in norm_idxs) {
    nb_counts <- neighbour_counts_xy(
      grid,
      maps$label_by_id,
      maps$alive_by_id,
      Cells$X[i], Cells$Y[i], cfg$radius
    )
    Nt <- as.integer(nb_counts$Nt)
    if (Nt > 0 && runif(1) <= cfg$N_division_possibility) {
      Cells$division[i]     <- 1L
      Cells$Time[i]         <- 0L
      Cells$TimeToDivide[i] <- 0L
    } else if (Nt == 0 && runif(1) <= (cfg$N_division_possibility / 10)) {
      Cells$division[i]     <- 1L
      Cells$Time[i]         <- 0L
      Cells$TimeToDivide[i] <- 0L
    }
  }
  Cells
}

# --------------------------- Main simulation ----------------------------

run_simulation <- function(cfg_file, karyolib_file, outputdir, prefix = "Cells") {
  cfg <- read_config(cfg_file)
  karyolib <- read_karyolib(karyolib_file)
  cfg$Fd <- karyolib$fitness[match(kt_vec2str(rep(2L, 22)), karyolib$karyotype)]
  cfg$Tg <- 1.2 * cfg$Cg
  run_simulation_cfg(cfg, karyolib, outputdir, prefix)
}

# Helper function: run simulation with parsed cfg and karyolib directly
run_simulation_cfg <- function(cfg, karyolib, outputdir, prefix = "Cells", global_id = NA_integer_) {
  init        <- init_cells(cfg, karyolib)
  Cells       <- init$Cells
  Cells$global_id <- as.integer(global_id)
  grid        <- init$grid
  id_state    <- list(max_id = init$max_id)

  # Buffer for hourly events
  events_buffer <- list()

  # initial per-id maps + bulk r
  maps <- build_id_maps(Cells)
  r_active <- calc_rates_all_cpp(grid,
                                 maps$active_ids,
                                 maps$x_by_id, maps$y_by_id,
                                 maps$label_by_id, maps$f_by_id, maps$alive_by_id,
                                 cfg$Rt, cfg$Rn, cfg$Cg, cfg$Fd, cfg$Tg,
                                 cfg$radius)
  # assign back by id
  if (length(maps$active_ids) > 0) {
    idx_alive <- which(Cells$Status == 1L)
    pos <- match(Cells$id[idx_alive], maps$active_ids)
    Cells$r[idx_alive] <- r_active[pos]
  }
  Cells$DivisionTime <- ifelse(Cells$Status == 1L, 1 / pmax(Cells$r, 1e-12), NA_real_)
  Cells$TimeToDivide <- ifelse(Cells$Status == 1L, Cells$DivisionTime - Cells$Time, NA_real_)

  step        <- 0L
  day_index   <- 0L

  # day 0 snapshot
  dir.create(file.path(outputdir, "csv"), showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(outputdir, "png"), showWarnings = FALSE, recursive = TRUE)
  write.csv(Cells, sprintf(paste0(outputdir, "/csv/%s_day%03d.csv"), prefix, day_index), row.names = FALSE)
  cat("Day", day_index, ":", sum(Cells$Status==1L), "living cells saved\n")

  # Build K matrix once (will be replaced by C++ returns after each hour)
  k_cols <- paste0("K", 1:22)
  K_mat <- as.matrix(Cells[, k_cols])

  max_steps <- as.integer(cfg$Time)
  if (is.na(max_steps) || max_steps < 1L) max_steps <- 240L
  while (step < max_steps) {
    step <- step + 1L

    # Call one full hourly step in C++
    res <- hourly_step_core_cpp(
      grid = grid,
      id   = Cells$id,
      X    = Cells$X,
      Y    = Cells$Y,
      Label= as.integer(Cells$Label),
      Status= Cells$Status,
      division = Cells$division,
      Time  = Cells$Time,
      DivisionTime = Cells$DivisionTime,
      TimeToDivide = Cells$TimeToDivide,
      Cooldown = Cells$Cooldown,
      r = Cells$r,
      f = Cells$f,
      G = Cells$G,
      Mr = Cells$Mr,
      K  = K_mat,
      Rt = cfg$Rt, Rn = cfg$Rn, Cg = cfg$Cg, Fd = cfg$Fd, Tg = cfg$Tg,
      radius = cfg$radius, N_division_possibility = cfg$N_division_possibility,
      T_death_rate = cfg$T_death_rate, N_death_rate = cfg$N_death_rate,
      total_supply = cfg$total_supply, DTr = cfg$DTr, DNr = cfg$DNr, supply_mode = cfg$supply_mode, MSR = cfg$MSR, Ctd = cfg$Ctd, CooldownFactor = cfg$CooldownFactor,
      step = step, N = cfg$N, max_id = id_state$max_id,
      karyo_lib_str = karyolib$karyotype,
      karyo_lib_fit = karyolib$fitness
    )

    # Rebuild Cells from the C++ return to avoid size-mismatch issues
    grid <- res$grid
    Cells <- tibble(
      id = as.integer(res$id),
      X  = as.integer(res$X),
      Y  = as.integer(res$Y),
      f  = as.numeric(res$f),
      r  = as.numeric(res$r),
      G  = as.numeric(res$G),
      Label = as.integer(res$Label),
      DivisionTime = as.numeric(res$DivisionTime),
      Time = as.integer(res$Time),
      Status = as.integer(res$Status),
      Mr = as.numeric(res$Mr),
      TimeToDivide = as.numeric(res$TimeToDivide),
      division = as.integer(res$division),
      Cooldown = as.numeric(res$Cooldown),
      g_alloc = as.numeric(res$g_alloc),
      div_event = as.integer(res$div_event),
      death_event = as.integer(res$death_event),
      risk_div     = as.integer(res$risk_div),
      death_resource = as.integer(res$death_resource),
      death_random   = as.integer(res$death_random),
      death_timeout  = as.integer(res$death_timeout)
    )
    Cells$global_id <- as.integer(global_id)

    # Write K columns back to the data frame
    K_mat <- res$K
    for (j in 1:22) Cells[[paste0("K", j)]] <- K_mat[, j]

    id_state$max_id <- res$max_id

    # Push per-hour events to buffer
    events_buffer[[length(events_buffer) + 1L]] <- tibble::tibble(
      step = step,
      id   = Cells$id,
      Label = Cells$Label,
      g_alloc = Cells$g_alloc,
      div_event = Cells$div_event,
      death_event = Cells$death_event,
      risk_div = Cells$risk_div,
      death_resource = Cells$death_resource,
      death_random   = Cells$death_random,
      death_timeout  = Cells$death_timeout
    )

    # Daily snapshot
    if (step %% 24L == 0L) {
      day_index <- day_index + 1L

      # --- living-only snapshot ---
      live_idx   <- which(Cells$Status == 1L)
      Cells_live <- Cells[live_idx, , drop = FALSE]

      csv_path <- sprintf(paste0(outputdir, "/csv/%s_day%03d.csv"), prefix, day_index)
      png_path <- sprintf(paste0(outputdir, "/png/%s_day%03d.png"), prefix, day_index)

      # Write CSV with living cells only
      write.csv(Cells_live, csv_path, row.names = FALSE)

      # Plot PNG (living cells only), no axes; Label=1 red, Label=0 green
      png(filename = png_path, width = 1000, height = 1000, units = "px")
      op <- par(mar = c(0, 0, 0, 0), oma = c(0, 0, 2, 0), pty = "s")
      plot(NA, xlim = c(1, cfg$N), ylim = c(1, cfg$N), axes = FALSE, xlab = "", ylab = "", asp = 1)
      if (nrow(Cells_live) > 0) {
        idx_norm <- which(Cells_live$Label == 0L)
        idx_tum  <- which(Cells_live$Label == 1L)
        if (length(idx_norm)) points(Cells_live$X[idx_norm], Cells_live$Y[idx_norm], pch = 16, cex = 1.5, col = "green")
        if (length(idx_tum))  points(Cells_live$X[idx_tum],  Cells_live$Y[idx_tum],  pch = 16, cex = 1.5, col = "red")
      }
      # Build title "total_supply: XX  Day:XX" with left alignment; use supply ratio for XX
      title_txt <- sprintf("total supply factor: %s  Day:%d",
                           format(signif(cfg$total_supply / (cfg$N * cfg$N), 6), trim = TRUE),
                           day_index)
      mtext(title_txt, side = 3, line = 0.5, adj = 0, outer = TRUE, cex = 1.2)
      par(op); dev.off()

      # Flush hourly events collected for the past day
      if (length(events_buffer) > 0) {
        ev_dir <- file.path(outputdir, "events")
        dir.create(ev_dir, showWarnings = FALSE, recursive = TRUE)
        ev_tbl <- dplyr::bind_rows(events_buffer)
        ev_path <- sprintf(paste0(ev_dir, "/%s_day%03d_events.csv"), prefix, day_index)
        write.csv(ev_tbl, ev_path, row.names = FALSE)
        events_buffer <- list()
      }

      # Log message
      cat("Day", day_index, ":", nrow(Cells_live), "living cells saved\n")

      # Compact in-memory state: drop dead rows from Cells and K_mat
      if (length(live_idx) < nrow(Cells)) {
        Cells <- Cells[live_idx, , drop = FALSE]
        K_mat <- K_mat[live_idx, , drop = FALSE]
      }
    }

  }
  cat("Reached configured Time =", cfg$Time, "hours. Simulation finished (living cells:", sum(Cells$Status==1L), ").\n")
  invisible(Cells)
}

# ----------------------- Batch runner from YAML -----------------------
run_ABM_simulation <- function(cfg_file, karyolib_file, base_output, workers = max(1, parallel::detectCores() - 1)) {
  cfg0 <- read_config(cfg_file)
  karyolib <- read_karyolib(karyolib_file)
  cfg0$Fd <- karyolib$fitness[match(kt_vec2str(rep(2L, 22)), karyolib$karyotype)]
  cfg0$Tg <- 1.2 * cfg0$Cg

  supplies <- cfg0$total_supply_vec
  reps <- cfg0$replicates

  # Prepare all jobs (supply x replicate) with simulation index reset per supply
  jobs <- list()
  for (s in supplies) {
    supply_label <- format(signif(s / (cfg0$N * cfg0$N), 6), trim = TRUE)
    # reset counter per supply
    for (r in seq_len(reps)) {
      out_dir <- file.path(base_output, "Results", supply_label, sprintf("Sim_%03d", r))
      cfg_i <- cfg0
      cfg_i$total_supply <- s
      jobs[[length(jobs)+1L]] <- list(
        cfg = cfg_i,
        out_dir = out_dir,
        prefix = sprintf("Cells_supply_%s_rep%03d", gsub("[^0-9.]", "", supply_label), r),
        cfg_path = cfg_file,
        karyolib_path = karyolib_file,
        global_id = as.integer(r)  # within-supply replicate id (001..reps)
      )
    }
  }

  # Create a manifest to verify mapping (supply x replicate -> folder/prefix)
  results_root <- file.path(base_output, "Results")
  dir.create(results_root, recursive = TRUE, showWarnings = FALSE)
  if (length(jobs) > 0) {
    manifest <- tibble::tibble(
      supply_mode  = vapply(jobs, function(j) j$cfg$supply_mode, character(1)),
      supply_value = vapply(jobs, function(j) j$cfg$total_supply, numeric(1)),
      supply_ratio = vapply(jobs, function(j) j$cfg$total_supply / (cfg0$N * cfg0$N), numeric(1)),
      supply_label = vapply(jobs, function(j) basename(dirname(j$out_dir)), character(1)),
      CooldownFactor = vapply(jobs, function(j) j$cfg$CooldownFactor, numeric(1)),
      replicate    = vapply(jobs, function(j) as.integer(j$global_id), integer(1)),
      global_id    = vapply(jobs, function(j) as.integer(j$global_id), integer(1)),
      out_dir      = vapply(jobs, function(j) j$out_dir, character(1)),
      prefix       = vapply(jobs, function(j) j$prefix, character(1))
    )
    write.csv(manifest, file.path(results_root, "manifest.csv"), row.names = FALSE)
  }

  # Run in parallel, avoiding global serialization and reloading code/kernels in each worker
  old_plan <- future::plan()
  on.exit({ try(future::plan(old_plan), silent = TRUE) }, add = TRUE)
  future::plan(future::multisession, workers = workers)

  future_lapply(
    jobs,
    FUN = function(job) {
      # Resolve project root from cfg path
      proj_dir <- dirname(job$cfg_path)
      # Load compiled kernels in the worker (no rebuild)
      Rcpp::sourceCpp(file.path(proj_dir, "src", "abm_kernels.cpp"), rebuild = FALSE, verbose = FALSE)
      # Source this script to get all R helpers into the worker session
      suppressMessages(source(file.path(proj_dir, "GBM_abm.R"), local = TRUE, chdir = TRUE))
      # Read karyotype library in the worker
      karyolib <- read_karyolib(job$karyolib_path)
      # Prepare cfg for this job (Fd, Tg depend on karyolib)
      cfg_i <- job$cfg
      cfg_i$Fd <- karyolib$fitness[match(kt_vec2str(rep(2L, 22)), karyolib$karyotype)]
      cfg_i$Tg <- 1.2 * cfg_i$Cg

      # Log which supply/replicate this worker is running
      ratio_str <- format(signif(cfg_i$total_supply / (cfg_i$N * cfg_i$N), 6), trim = TRUE)
      cat(sprintf("[RUN] supply_ratio=%s, replicate=%s -> %s\n",
                  ratio_str,
                  sprintf("%03d", as.integer(job$global_id)),
                  job$out_dir))

      # Ensure output directory exists and drop a small job_info file for auditing
      dir.create(job$out_dir, recursive = TRUE, showWarnings = FALSE)
      write.csv(
        data.frame(
          supply_mode  = cfg_i$supply_mode,
          supply_value = cfg_i$total_supply,
          supply_ratio = as.numeric(cfg_i$total_supply) / (cfg_i$N * cfg_i$N),
          supply_label = basename(dirname(job$out_dir)),
          CooldownFactor = cfg_i$CooldownFactor,
          replicate    = as.integer(job$global_id),
          global_id    = as.integer(job$global_id),
          out_dir      = job$out_dir,
          prefix       = job$prefix
        ),
        file = file.path(job$out_dir, "job_info.csv"), row.names = FALSE
      )

      # Run single simulation
      run_simulation_cfg(cfg_i, karyolib, job$out_dir, prefix = job$prefix, global_id = job$global_id)
      invisible(TRUE)
    },
    future.globals = FALSE,
    future.packages = c("Rcpp","yaml","dplyr","tibble","future.apply"),
    future.seed = TRUE
  )
}

# GBM Agent‑Based Model Simulation in a 2‑D lattice
# -------------------------------------------------
# This script implements an agent‑based model (ABM) for the growth of
# glioblastoma (GBM) tumours in two distinct brain regions (e.g., frontal
# and temporal lobes).  Each grid cell hosts at most one biological cell
# (tumour or normal).  Key processes implemented:
#   • Initial placement of one tumour cell at the centre and a user‑defined
#     fraction of normal cells elsewhere.
#   • Cell‑specific growth‑rate calculation that incorporates karyotype‑
#     dependent fitness and Lotka–Volterra–type competition with neighbours.
#   • Stochastic cell division with optional chromosome mis‑segregation.
#   • Glucose‑limited death.
#   • Hourly time steps (Δt = 1 h); grids are written to CSV once per day.
#
# External dependencies ---------------------------------------------------
# install.packages(c("yaml","dplyr","tibble","stringr","purrr"))

library(yaml)
library(dplyr)
library(tibble)
library(stringr)

# -------------------------------------------------------------------------
#                               I/O HELPERS
# -------------------------------------------------------------------------

read_config <- function(path) {
  cfg_raw <- yaml::read_yaml(path)
  cfg_raw$N <- as.numeric(cfg_raw$Nd)
  cfg_raw$Ctd  <- as.numeric(cfg_raw$Ctd)
  cfg_raw$Rt  <- as.numeric(cfg_raw$Rt)
  cfg_raw$Cg  <- as.numeric(cfg_raw$Cg)
  cfg_raw$m   <- as.numeric(cfg_raw$m)
  cfg_raw$MSR <- as.numeric(cfg_raw$MSR)
  return(cfg_raw)
}

read_karyolib <- function(path) {
  # Expected: first column =karyotype string (e.g. "2.2.2…"), second = fitness
  kdf <- read.csv(path, stringsAsFactors = FALSE)
  names(kdf)[1:2] <- c("karyotype", "fitness")
  return(kdf)
}

kt_vec2str <- function(v) paste(v, collapse = ".")
kt_str2vec <- function(s) as.integer(str_split(s, "\\.")[[1]])

fitness_of <- function(karyolib, v) {
  s <- kt_vec2str(v)
  f <- karyolib$fitness[match(s, karyolib$karyotype)]
  if (is.na(f)) NA else f
}

# -------------------------------------------------------------------------
#                           INITIALISATION
# -------------------------------------------------------------------------

init_cells <- function(cfg, karyolib) {
  N   <- cfg$N                           # grid length
  m   <- cfg$m                           # initial normal‑cell density (0‑1)
  n_g <- N * N                           # total lattice sites
  grid <- matrix(NA_integer_, nrow = N, ncol = N)  # stores row index of Cells
  
  diploid <- rep(2L, 22)
  fit_dip <- fitness_of(karyolib, diploid)
  
  # Helper to build a cell row -------------------------------------------
  make_row <- function(x, y, kt_vec, label, r_base,f, mr) {
    tibble(
      X = x, Y = y,
      !!!setNames(as.list(kt_vec), paste0("K", 1:22)),
      f  = f,
      r  = r_base,                      # will be updated later
      G  = cfg$Cg,
      Label = label,                   # 0 = normal, 1 = tumour
      DivisionTime = (1 / r_base)/f,      # hours
      Time = 0,                        # hours since last division
      Status = 1L,                      # 1 = alive, 0 = dead
      Mr  = mr,
      TimeToDivide = 0
    )
  }
  speed_ratio <- cfg$Ctd
  # 1. Single tumour cell in the centre ----------------------------------
  centre <- c(ceiling(N / 2), ceiling(N / 2))
  Cells  <- make_row(centre[1], centre[2], diploid, 1L, cfg$Rt,fit_dip,cfg$Mtr)
  grid[centre[1], centre[2]] <- 1L
  Cells$G<-speed_ratio * fit_dip
  Cells$TimeToDivide<-Cells$DivisionTime-Cells$Time
  
  # 2. Random normal cells ------------------------------------------------
  n_norm <- floor(n_g * m)
  empty  <- which(is.na(grid), arr.ind = TRUE)
  sel    <- empty[sample(seq_len(nrow(empty)), n_norm), , drop = FALSE]
  for (i in seq_len(n_norm)) {
    row <- make_row(sel[i, 1], sel[i, 2], diploid, 0L,
                    1 ,1, cfg$Mnr)
    row$TimeToDivide<-NA
    row$r<-NA
    row$DivisionTime<-NA
    Cells <- bind_rows(Cells, row)
    grid[row$X, row$Y] <- nrow(Cells)
  }
  list(Cells = Cells, grid = grid)
}

# -------------------------------------------------------------------------
#                   LOCAL NEIGHBOURHOOD QUANTITIES
# -------------------------------------------------------------------------

get_neighbour_counts <- function(Cells, grid, idx, radius = 2) {
  N    <- nrow(grid)
  pos  <- which(grid == idx, arr.ind = TRUE)[1, ]
  if (length(pos) == 0) return(list(Nt = 0, Nn = 0))
  x    <- pos[1]; y <- pos[2]
  xmin <- max(1, x - radius)
  xmax <- min(N, x + radius)
  ymin <- max(1, y - radius)
  ymax <- min(N, y + radius)
  
  nbr_idx <- as.vector(grid[xmin:xmax, ymin:ymax])
  nbr_idx <- nbr_idx[!is.na(nbr_idx) & nbr_idx != idx]
  Nt <- sum(Cells$Label[nbr_idx] == 1L)
  Nn <- sum(Cells$Label[nbr_idx] == 0L)
  list(Nt = Nt, Nn = Nn)
}

# -------------------------------------------------------------------------
#                         GROWTH‑RATE FUNCTION
# -------------------------------------------------------------------------

calc_rate <- function(Cells, grid, idx, cfg) {
  N          <- cfg$N * cfg$N           # carrying capacity proxy
  counts     <- get_neighbour_counts(Cells, grid, idx)
  Nt <- counts$Nt; Nn <- counts$Nn
  
  a <- cfg$Cg / (((Cells$f[idx]) / cfg$Fd) * cfg$Tg)
  b <- 1 / a
  
  if (Cells$Label[idx] == 1L) {
    # tumour
    r <- (cfg$Rt * (Cells$f[idx] / cfg$Fd)) * (1 - ((Nt + a * Nn) / N))
  } else {
    # normal
    r <- cfg$Rt * (1 - ((b * Nt + Nn) / N))
  }
  if(r <=0)# mark negative r as dead
  {
    Cells$Status[idx] == 0L
  }
  return(r)
}

# -------------------------------------------------------------------------
#                         CELL DIVISION LOGIC
# -------------------------------------------------------------------------

divide_cell <- function(Cells, grid, idx, cfg, karyolib) {
  # 0. Skip dead cells
  if (Cells$Status[idx] == 0L) {
    return(list(Cells = Cells, grid = grid))
  }
  
  # 1. Find mother cell position and its free von Neumann neighbors
  pos <- which(grid == idx, arr.ind = TRUE)[1, ]
  nbr <- rbind(
    c(pos[1] - 1, pos[2]),  # up
    c(pos[1] + 1, pos[2]),  # down
    c(pos[1],     pos[2] - 1),  # left
    c(pos[1],     pos[2] + 1)   # right
  )
  valid    <- nbr[,1] >= 1 & nbr[,1] <= cfg$N &
    nbr[,2] >= 1 & nbr[,2] <= cfg$N
  valid_nbr <- nbr[valid, , drop = FALSE]
  free_pos  <- valid_nbr[is.na(grid[cbind(valid_nbr[,1], valid_nbr[,2])]), , drop = FALSE]
  if (nrow(free_pos) == 0) {
    return(list(Cells = Cells, grid = grid))
  }
  
  # 2. Copy mother into two daughters and reset their Time
  mother    <- Cells[idx, ]
  daughters <- bind_rows(mother, mother)
  daughters$Time <- 0L
  
  # 3. Assign karyotype, fitness and resource for each daughter
    speed_ratio <- cfg$Ctd
    for (i in 1:2) {
      repeat {
        kt    <- as.integer(mother[paste0("K", 1:22)])
        chr   <- sample(1:22, 1)
        delta <- sample(c(-1,1), 1)
        kt[chr] <- pmin(pmax(kt[chr] + delta, 0), 10)
        f_new   <- fitness_of(karyolib, kt)
        if (!is.na(f_new)) break
      }
      # write validated karyotype and fitness
      for (j in 1:22) {
        daughters[i, paste0("K", j)] <- kt[j]
      }
      daughters$f[i] <- f_new
      daughters$G[i] <- speed_ratio * f_new
  }
  
  # 4. Place daughters into Cells & grid
  # overwrite mother with first daughter
  Cells[idx, ] <- daughters[1, ]
  
  # pick a random free spot for second daughter
  target        <- free_pos[sample(nrow(free_pos), 1), ]
  daughters$X[2] <- target[1]
  daughters$Y[2] <- target[2]
  
  # append second daughter to Cells and update grid
  Cells       <- bind_rows(Cells, daughters[2, ])
  new_idx     <- nrow(Cells)
  grid[target[1], target[2]] <- new_idx
  
  # 5. Recompute growth rate r and DivisionTime for both daughters
  Cells$DivisionTime[idx]    <- 1 / Cells$r[idx]
  Cells$DivisionTime[new_idx] <- 1 / Cells$r[new_idx]
  
  return(list(Cells = Cells, grid = grid))
}
# -------------------------------------------------------------------------
#                             DEATH CHECK
# -------------------------------------------------------------------------

death_update <- function(Cells, grid, cfg) {
  # 1. Compute total supply and average per-unit-G supply
  N            <- cfg$N
  total_supply <- 1 * N * N
  total_G      <- sum(Cells$G)
  
  # Avoid division by zero
  if (total_G <= 0) {
    warning("Total G is zero or negative; skipping death update.")
    return(list(Cells = Cells, grid = grid))
  }
  
  avg_supply   <- total_supply / total_G
  death_thr    <- 0.2 * avg_supply
  
  # 2. Identify cells with G below threshold
  dead_idx <- which(Cells$G < death_thr)
  
  if (length(dead_idx) > 0) {
    # Mark them dead
    Cells$Status[dead_idx] <- 0L
    
    # 3. Remove dead rows from Cells
    keep_idx <- which(Cells$Status == 1L)
    Cells    <- Cells[keep_idx, , drop = FALSE]
    
    # 4. Update grid: clear dead positions, then reassign surviving IDs
    grid[ grid %in% dead_idx ] <- NA_integer_
    
    # Because Cells rows have been reindexed, map old IDs → new IDs
    # keep_idx[i] was old ID for new row i
    for (new_i in seq_along(keep_idx)) {
      old_i <- keep_idx[new_i]
      grid[ grid == old_i ] <- new_i
    }
  }
  
  list(Cells = Cells, grid = grid)
}




# -------------------------------------------------------------------------
#                             MIGRATION
# -------------------------------------------------------------------------

migrate_cells <- function(Cells, grid, cfg) {
  N <- cfg$N
  for (i in seq_len(nrow(Cells))) {
    if (Cells$Status[i] != 1L) next
    # decide whether to move
    if (runif(1) < Cells$Mr[i]) {
      # find von Neumann neighbours
      pos <- c(Cells$X[i], Cells$Y[i])
      nbr <- rbind(
        pos + c(-1,  0),
        pos + c( 1,  0),
        pos + c( 0, -1),
        pos + c( 0,  1)
      )
      # filter valid coords
      valid <- nbr[,1] >= 1 & nbr[,1] <= N &
        nbr[,2] >= 1 & nbr[,2] <= N
      nbr   <- nbr[valid, , drop = FALSE]
      # pick only empty spots
      empties <- nbr[ is.na(grid[cbind(nbr[,1], nbr[,2])]), , drop = FALSE]
      if (nrow(empties) == 0) next
      # choose one at random
      tgt <- empties[sample(nrow(empties), 1), ]
      # update grid
      grid[pos[1], pos[2]]          <- NA_integer_
      grid[tgt[1], tgt[2]]          <- i
      # update cell coords
      Cells$X[i] <- tgt[1]
      Cells$Y[i] <- tgt[2]
    }
  }
  list(Cells = Cells, grid = grid)
}


# -------------------------------------------------------------------------
#                          MAIN SIMULATION LOOP
# -------------------------------------------------------------------------

run_simulation <- function(cfg_file, karyolib_file, prefix = "Cells") {
  cfg         <- read_config(cfg_file)
  karyolib    <- read_karyolib(karyolib_file)
  cfg$Fd      <- karyolib$fitness[match(kt_vec2str(rep(2L, 22)), karyolib$karyotype)]
  cfg$Tg      <- 1.2 * cfg$Cg
  
  init        <- init_cells(cfg, karyolib)
  Cells       <- init$Cells
  grid        <- init$grid
  
  step        <- 0L
  day_index   <- 0L
  write.csv(Cells, sprintf("./Results/%s_day%03d.csv", prefix, day_index), row.names = FALSE)
  cat("Day", day_index, ":", nrow(Cells), "living cells saved\n")
  repeat {
    step <- step + 1L
    
    # advance internal clocks -------------------------------------------
    Cells$Time <- Cells$Time + 1
    Cells$TimeToDivide<-Cells$DivisionTime-Cells$Time
    # identify cells ready to divide using TimeToDivide <= 0 and Status == 1
    to_divide_idx <- which(Cells$Status == 1L & Cells$TimeToDivide <= 0 & Cells$Label ==1L)
    for (idx in to_divide_idx) {
      res <- divide_cell(Cells, grid, idx, cfg, karyolib)
      Cells <- res$Cells; grid <- res$grid
    }
    for (idx in 1:nrow(Cells))
    {
        Cells$r[idx] <- calc_rate(Cells, grid, idx, cfg)
    }
    # death --------------------------------------------------------------
    res <- death_update(Cells, grid, cfg)
    Cells <- res$Cells; grid <- res$grid
    # migration------------------------------------------------------------
    mig <- migrate_cells(Cells, grid, cfg)
    Cells <- mig$Cells
    grid  <- mig$grid
    
    # daily snapshot -----------------------------------------------------
    if (step %% 24L == 0L) {
      day_index <- day_index + 1L
      write.csv(Cells, sprintf("./Results/%s_day%03d.csv", prefix, day_index), row.names = FALSE)
      cat("Day", day_index, ":", nrow(Cells), "living cells saved\n")
    }
    
    # termination: grid full --------------------------------------------
    if (all(!is.na(grid))) {
      cat("Grid saturated (", nrow(Cells), "cells). Stopping.\n", sep = "")
      break
    }
  }
  invisible(Cells)
}

# -------------------------------------------------------------------------
#                            EXAMPLE CALL
# -------------------------------------------------------------------------

# Uncomment to run with your own paths ------------------------------------
# result <- run_simulation("config.yaml", "karyolib.csv")

result <- run_simulation('/Volumes/Protable Disk/HIP-IMO/config.yaml', '/Volumes/Protable Disk/HIP-IMO/TEMPORAL_landscape_summary.csv')

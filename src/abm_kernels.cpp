// src/abm_kernels.cpp
// Rcpp kernels for performance-critical parts of the GBM ABM:
//   (1) Moore-neighbour queries and free-spot discovery
//   (2) Divide placement (choose a free neighbour cell and write to grid)
//   (3) Single-cell migration to a free neighbour cell
//   (4) Chromosome mis-segregation (MS) events per division with
//       Poisson-distributed event count, distinct chromosomes, random
//       direction, per-chromosome sum conservation, and karyotype-library validation.
//   (5) Bulk growth-rate evaluation for all alive cells (calc_rates_all_cpp)
//       using a stable-id grid and per-id attribute vectors.

#include <Rcpp.h>
#include <unordered_map>
#include <sstream>
#include <vector>
#include <algorithm>
#include <random>
using namespace Rcpp;

inline bool isNAint(int v) {
  return Rcpp::IntegerVector::is_na(v);
}

// ---------------------------------------------------------------------
// (1) Return free Moore neighbours as a flat integer vector
//     [x1, y1, x2, y2, ...] using 1-based coordinates (R style).
// ---------------------------------------------------------------------
// [[Rcpp::export]]
IntegerVector moore_free_neighbors(const IntegerMatrix& grid, int x, int y) {
  const int N = grid.nrow();
  std::vector<int> out; out.reserve(16);

  for (int dx = -1; dx <= 1; ++dx) {
    for (int dy = -1; dy <= 1; ++dy) {
      if (dx == 0 && dy == 0) continue;
      int nx = x + dx;
      int ny = y + dy;
      if (nx >= 1 && nx <= N && ny >= 1 && ny <= N) {
        int val = grid(nx - 1, ny - 1); // C++ is 0-based
        if (isNAint(val)) {
          out.push_back(nx);
          out.push_back(ny);
        }
      }
    }
  }
  return Rcpp::wrap(out);
}

// ---------------------------------------------------------------------
// (2a) Count tumour/normal neighbours (Nt/Nn) within a Moore window
//      by row-id lookup (legacy). Self is excluded.
//      grid holds IDs (stable), or NA for empty.
//      label_by_id: length >= max_id+1, 0=normal, 1=tumour, -1=unused.
//      alive_by_id: 0/1 mask, length >= max_id+1.
// ---------------------------------------------------------------------
// [[Rcpp::export]]
List neighbour_counts_idx(const IntegerMatrix& grid,
                          const IntegerVector& label_by_id,
                          const IntegerVector& alive_by_id,
                          int idx,
                          int radius) {
  const int N = grid.nrow();

  // Locate coordinates of id=idx by scanning grid (slow path; kept for compatibility)
  int x = -1, y = -1;
  for (int i = 0; i < N && x < 0; ++i) {
    for (int j = 0; j < N; ++j) {
      if (!isNAint(grid(i, j)) && grid(i, j) == idx) {
        x = i + 1; y = j + 1; break;
      }
    }
  }
  if (x < 0) {
    return List::create(_["Nt"]=0, _["Nn"]=0);
  }

  int xmin = std::max(1, x - radius);
  int xmax = std::min(N, x + radius);
  int ymin = std::max(1, y - radius);
  int ymax = std::min(N, y + radius);

  int Nt = 0, Nn = 0;
  for (int i = xmin; i <= xmax; ++i) {
    for (int j = ymin; j <= ymax; ++j) {
      if (i == x && j == y) continue; // exclude self
      int id = grid(i - 1, j - 1);
      if (!isNAint(id)) {
        if (alive_by_id[id] == 1) {
          int lab = label_by_id[id];
          if (lab == 1) ++Nt; else if (lab == 0) ++Nn;
        }
      }
    }
  }
  return List::create(_["Nt"]=Nt, _["Nn"]=Nn);
}

// ---------------------------------------------------------------------
// (2b) Count neighbours when caller provides (x,y). Faster.
// [[Rcpp::export]]
List neighbour_counts_xy(const IntegerMatrix& grid,
                         const IntegerVector& label_by_id,
                         const IntegerVector& alive_by_id,
                         int x,
                         int y,
                         int radius) {
  const int N = grid.nrow();
  int xmin = std::max(1, x - radius);
  int xmax = std::min(N, x + radius);
  int ymin = std::max(1, y - radius);
  int ymax = std::min(N, y + radius);
  int Nt = 0, Nn = 0;
  for (int i = xmin; i <= xmax; ++i) {
    for (int j = ymin; j <= ymax; ++j) {
      if (i == x && j == y) continue; // exclude self
      int id = grid(i - 1, j - 1);
      if (!isNAint(id)) {
        if (alive_by_id[id] == 1) {
          int lab = label_by_id[id];
          if (lab == 1) ++Nt; else if (lab == 0) ++Nn;
        }
      }
    }
  }
  return List::create(_["Nt"]=Nt, _["Nn"]=Nn);
}

// ---------------------------------------------------------------------
// (3) Divide placement: choose a free Moore neighbour for the daughter
//     and write its ID into the grid. Mother stays in place.
// ---------------------------------------------------------------------
// [[Rcpp::export]]
List divide_place_cpp(IntegerMatrix grid,
                      int mother_x, int mother_y,
                      int new_id,
                      int N) {
  std::vector<std::pair<int,int>> freepos;
  freepos.reserve(8);
  for (int dx = -1; dx <= 1; ++dx) {
    for (int dy = -1; dy <= 1; ++dy) {
      if (dx == 0 && dy == 0) continue;
      int nx = mother_x + dx;
      int ny = mother_y + dy;
      if (nx >= 1 && nx <= N && ny >= 1 && ny <= N) {
        int v = grid(nx - 1, ny - 1);
        if (isNAint(v)) freepos.emplace_back(nx, ny);
      }
    }
  }
  if (freepos.empty()) {
    return List::create(_["success"]=false,
                        _["grid"]=grid,
                        _["x"]=NA_INTEGER,
                        _["y"]=NA_INTEGER);
  }
  int k = (int) std::floor(R::unif_rand() * freepos.size());
  if (k == (int)freepos.size()) k = (int)freepos.size() - 1;
  int tx = freepos[k].first;
  int ty = freepos[k].second;

  grid(tx - 1, ty - 1) = new_id;

  return List::create(_["success"]=true,
                      _["grid"]=grid,
                      _["x"]=tx,
                      _["y"]=ty);
}

// ---------------------------------------------------------------------
// (4) Single-cell migration: move from (x,y) to a random free Moore
//     neighbour if available; update grid and return new coordinates.
//     row_id is now the stable ID.
// ---------------------------------------------------------------------
// [[Rcpp::export]]
List migrate_one_cpp(IntegerMatrix grid, int x, int y, int id, int N) {
  std::vector<std::pair<int,int>> freepos;
  freepos.reserve(8);
  for (int dx = -1; dx <= 1; ++dx) {
    for (int dy = -1; dy <= 1; ++dy) {
      if (dx == 0 && dy == 0) continue;
      int nx = x + dx;
      int ny = y + dy;
      if (nx >= 1 && nx <= N && ny >= 1 && ny <= N) {
        int v = grid(nx - 1, ny - 1);
        if (isNAint(v)) freepos.emplace_back(nx, ny);
      }
    }
  }
  if (freepos.empty()) {
    return List::create(_["success"]=false, _["grid"]=grid, _["x"]=x, _["y"]=y);
  }
  int k = (int) std::floor(R::unif_rand() * freepos.size());
  if (k == (int)freepos.size()) k = (int)freepos.size() - 1;
  int tx = freepos[k].first;
  int ty = freepos[k].second;

  grid(x - 1, y - 1) = NA_INTEGER;
  grid(tx - 1, ty - 1) = id;

  return List::create(_["success"]=true, _["grid"]=grid, _["x"]=tx, _["y"]=ty);
}

// ============================ MS event kernel ===========================

static inline std::string kt_to_str(const IntegerVector& kt) {
  std::ostringstream oss;
  for (int i = 0; i < kt.size(); ++i) {
    if (i) oss << ".";
    oss << kt[i];
  }
  return oss.str();
}

static inline double lookup_fitness(const std::unordered_map<std::string,double>& fmap,
                                    const IntegerVector& kt) {
  auto it = fmap.find(kt_to_str(kt));
  if (it == fmap.end()) return NA_REAL;
  return it->second;
}

// ---------------------------------------------------------------------
// Apply MS events for a single division (tumour cells only):
//   * ms_count ~ Poisson(MSR)
//   * select ms_count DISTINCT chromosomes (if available) with bounds
//     allowing +/-1 (mother in [1,9])
//   * per-chromosome random direction: (c-1,c+1) or (c+1,c-1)
//   * sum conservation per chromosome; both daughters must exist in the
//     provided karyotype library; otherwise retry up to max_attempts.
//   * on failure, fallback to no-MS (both daughters == mother).
// ---------------------------------------------------------------------
// [[Rcpp::export]]
List apply_ms_events_cpp(IntegerVector mother_kt,
                         int label,
                         double MSR,
                         CharacterVector karyo_lib_str,
                         NumericVector  karyo_lib_fit,
                         int max_attempts = 200) {
  // Build karyotype->fitness map
  std::unordered_map<std::string,double> fmap;
  fmap.reserve(karyo_lib_str.size());
  for (int i = 0; i < karyo_lib_str.size(); ++i) {
    fmap[ Rcpp::as<std::string>(karyo_lib_str[i]) ] = (double)karyo_lib_fit[i];
  }

  // Default: no MS (daughters identical to mother)
  IntegerVector kt1 = clone(mother_kt);
  IntegerVector kt2 = clone(mother_kt);
  double f1 = lookup_fitness(fmap, kt1);
  double f2 = lookup_fitness(fmap, kt2);

  // Only tumour cells undergo MS
  if (label != 1) {
    return List::create(_["kt1"]=kt1, _["kt2"]=kt2, _["f1"]=f1, _["f2"]=f2);
  }

  // Draw number of MS events
  int ms_count = R::rpois(MSR);
  if (ms_count <= 0) {
    return List::create(_["kt1"]=kt1, _["kt2"]=kt2, _["f1"]=f1, _["f2"]=f2);
  }

  // Valid chromosomes must allow -1/+1 (mother in [1,9])
  std::vector<int> valid_chr;
  valid_chr.reserve(mother_kt.size());
  for (int i = 0; i < mother_kt.size(); ++i) {
    int c = mother_kt[i];
    if (c >= 1 && c <= 9) valid_chr.push_back(i);
  }
  if (valid_chr.empty()) {
    return List::create(_["kt1"]=kt1, _["kt2"]=kt2, _["f1"]=f1, _["f2"]=f2);
  }

  ms_count = std::min<int>(ms_count, (int)valid_chr.size());

  // Try until both daughters exist in library or attempts exhausted
  for (int att = 0; att < max_attempts; ++att) {
    // Sample ms_count distinct chromosomes without replacement
    IntegerVector pool = IntegerVector(valid_chr.begin(), valid_chr.end());
    IntegerVector idx = Rcpp::sample(pool, ms_count, false);

    IntegerVector t1 = clone(mother_kt);
    IntegerVector t2 = clone(mother_kt);

    // Per-chromosome random direction, sum conservation
    for (int k = 0; k < idx.size(); ++k) {
      int i_chr = idx[k];
      if (unif_rand() < 0.5) {
        t1[i_chr] = mother_kt[i_chr] - 1;
        t2[i_chr] = mother_kt[i_chr] + 1;
      } else {
        t1[i_chr] = mother_kt[i_chr] + 1;
        t2[i_chr] = mother_kt[i_chr] - 1;
      }
    }

    double nf1 = lookup_fitness(fmap, t1);
    double nf2 = lookup_fitness(fmap, t2);
    if (!Rcpp::NumericVector::is_na(nf1) && !Rcpp::NumericVector::is_na(nf2)) {
      kt1 = t1; kt2 = t2; f1 = nf1; f2 = nf2;
      break;
    }
  }

  return List::create(_["kt1"]=kt1, _["kt2"]=kt2, _["f1"]=f1, _["f2"]=f2);
}

// ---------------------------------------------------------------------
// (5) Bulk growth-rate evaluation for all alive cells.
//     Inputs are stable-id based: per-id x,y,label,f,alive,
//     and the grid stores IDs (or NA).
//     Returns a numeric vector of r for each active id in the same order.
// ---------------------------------------------------------------------
// [[Rcpp::export]]
NumericVector calc_rates_all_cpp(const IntegerMatrix& grid,
                                 const IntegerVector& active_ids,
                                 const IntegerVector& x_by_id,
                                 const IntegerVector& y_by_id,
                                 const IntegerVector& label_by_id,
                                 const NumericVector& f_by_id,
                                 const IntegerVector& alive_by_id,
                                 double Rt, double Rn, double Cg, double Fd, double Tg,
                                 int radius) {
  const int N = grid.nrow();
  NumericVector r_out(active_ids.size());

  for (int idx = 0; idx < active_ids.size(); ++idx) {
    int id = active_ids[idx];
    if (alive_by_id[id] != 1) { r_out[idx] = 0.0; continue; }

    int lab = label_by_id[id];
    double f  = f_by_id[id];
    int x = x_by_id[id];
    int y = y_by_id[id];

    // neighbourhood bounds
    int xmin = std::max(1, x - radius);
    int xmax = std::min(N, x + radius);
    int ymin = std::max(1, y - radius);
    int ymax = std::min(N, y + radius);
    int capacity = (xmax - xmin + 1) * (ymax - ymin + 1);

    // neighbour counts
    int Nt = 0, Nn = 0;
    for (int i = xmin; i <= xmax; ++i) {
      for (int j = ymin; j <= ymax; ++j) {
        if (i == x && j == y) continue;
        int nid = grid(i - 1, j - 1);
        if (!isNAint(nid) && alive_by_id[nid] == 1) {
          int nlab = label_by_id[nid];
          if (nlab == 1) ++Nt; else if (nlab == 0) ++Nn;
        }
      }
    }

    // LV-style rate
    double r = 0.0;
    if (lab == 1) {
      // tumour
      double a = Cg / ((f / Fd) * Tg);
      r = (Rt * (f / Fd)) * (1.0 - ((Nt + a * Nn) / (double)capacity));
    } else {
      // normal
      double a = Cg / ((f / Fd) * Tg); // f/Fd should be 1 for normals
      double b = 1.0 / a;
      r = Rn * (1.0 - ((b * Nt + Nn) / (double)capacity));
    }
    if (r < 0.0) r = 0.0;
    r_out[idx] = r;
  }

  return r_out;
}

// ---------------------------------------------------------------------
// (6) One full hourly step executed in C++ to minimize R<->C++ overhead.
//     Operates on row-oriented vectors with a stable-id grid.
//     Notes:
//       - `id` is stable and stored in `grid` cells (or NA for empty)
//       - We maintain an id->row index map each call for O(1) lookups
//       - Appending a daughter resizes all column vectors and K matrix
//       - DivisionTime is recomputed from r after bulk update
// ---------------------------------------------------------------------
// [[Rcpp::export]]
List hourly_step_core_cpp(IntegerMatrix grid,
                     IntegerVector id,
                     IntegerVector X,
                     IntegerVector Y,
                     IntegerVector Label,
                     IntegerVector Status,
                     IntegerVector division,
                     IntegerVector Time,
                     NumericVector DivisionTime,
                     NumericVector TimeToDivide,
                     NumericVector Cooldown,
                     NumericVector r,
                     NumericVector f,
                     NumericVector G,
                     NumericVector Mr,
                     IntegerMatrix K,
                     double Rt, double Rn, double Cg, double Fd, double Tg,
                     int radius, double N_division_possibility,
                     double T_death_rate, double N_death_rate,
                     double total_supply, double DTr, double DNr, std::string supply_mode, double MSR, double Ctd, double CooldownFactor,
                     int step, int N, int max_id,
                     CharacterVector karyo_lib_str,
                     NumericVector  karyo_lib_fit) {
  int n = id.size();

  // Per-cell hourly glucose allocation and event flags
  NumericVector g_alloc(n, 0.0);
  IntegerVector div_event(n, 0);
  IntegerVector death_event(n, 0);
  IntegerVector risk_div(n, 0);        // 1 if cell is in the division risk set this hour
  IntegerVector death_resource(n, 0);  // 1 if a resource-based death occurred this hour
  IntegerVector death_random(n, 0);    // 1 if a random (daily) death occurred this hour
  IntegerVector death_timeout(n, 0);   // 1 if a timeout death occurred this hour

  // Build id->row map (size max_id+1)
  IntegerVector id2row(max_id + 1, NA_INTEGER);
  for (int i = 0; i < n; ++i) {
    if (i < id.size() && !IntegerVector::is_na(id[i])) {
      int myid = id[i];
      if (myid >= 0 && myid <= max_id) id2row[myid] = i + 1; // 1-based row index
    }
  }

  // Helper lambda: neighbour counts for a given row index
  auto nb_counts_row = [&](int row)->std::pair<int,int>{
    int x = X[row];
    int y = Y[row];
    int xmin = std::max(1, x - radius);
    int xmax = std::min(N, x + radius);
    int ymin = std::max(1, y - radius);
    int ymax = std::min(N, y + radius);
    int Nt = 0, Nn = 0;
    for (int i = xmin; i <= xmax; ++i) {
      for (int j = ymin; j <= ymax; ++j) {
        if (i == x && j == y) continue;
        int nid = grid(i - 1, j - 1);
        if (!IntegerVector::is_na(nid)) {
          int rj = id2row[nid];
          if (!IntegerVector::is_na(rj)) {
            int jrow = rj - 1; // back to 0-based
            if (Status[jrow] == 1) {
              if (Label[jrow] == 1) ++Nt; else if (Label[jrow] == 0) ++Nn;
            }
          }
        }
      }
    }
    return {Nt, Nn};
  };

  // (A) Activate normal division where applicable
  for (int irow = 0; irow < n; ++irow) {
    if (Status[irow] != 1) continue;
    if (Label[irow] != 0) continue; // only normal
    if (division[irow] != 0) continue;
    if (Cooldown[irow] > 0.0) continue; // still cooling down
    auto nb = nb_counts_row(irow);
    int Nt = nb.first;
    // Only activate if tumour neighbours exist; no baseline activation when Nt==0
    if (Nt > 0 && unif_rand() <= N_division_possibility) {
      division[irow] = 1;
      Time[irow] = 0;
      TimeToDivide[irow] = 0.0;
    }
  }

  // Mark division risk set at the start of the hour (before clocks advance)
  for (int irow = 0; irow < n; ++irow) {
    if (Status[irow] == 1 && division[irow] == 1) risk_div[irow] = 1; else risk_div[irow] = 0;
  }

  // (B) Advance clocks by 1 hour for all alive cells
  for (int irow = 0; irow < n; ++irow) {
    if (Status[irow] == 1) {
      Time[irow] += 1;
      if (Cooldown[irow] > 0.0) {
        Cooldown[irow] = std::max(0.0, Cooldown[irow] - 1.0);
      }
    }
  }

  // (C) Division attempts: collect ready rows and shuffle
  std::vector<int> ready;
  ready.reserve(n);
  for (int irow = 0; irow < n; ++irow) {
    if (Status[irow] != 1) continue;
    if (division[irow] != 1) continue;
    if (Time[irow] < DivisionTime[irow]) continue;
    ready.push_back(irow);
  }
  std::random_device rd;
  std::mt19937 gen(rd());
  std::shuffle(ready.begin(), ready.end(), gen);

  // Utility: apply MS and update mother & daughter karyotypes/fitness
  auto apply_ms = [&](int mother_row, IntegerVector& kt1, IntegerVector& kt2, double& f1, double& f2){
    IntegerVector mother_kt(22);
    for (int c = 0; c < 22; ++c) mother_kt[c] = K(mother_row, c);
    List ms_res = apply_ms_events_cpp(mother_kt, (int)Label[mother_row], MSR,
                                      karyo_lib_str, karyo_lib_fit, 200);
    kt1 = as<IntegerVector>(ms_res["kt1"]);
    kt2 = as<IntegerVector>(ms_res["kt2"]);
    f1  = as<double>(ms_res["f1"]);
    f2  = as<double>(ms_res["f2"]);
  };

  // (C.1) Iterate over ready rows
  for (size_t u = 0; u < ready.size(); ++u) {
    int irow = ready[u];
    if (Status[irow] != 1) continue; // may have died earlier in loop

    int mx = X[irow], my = Y[irow];

    int new_id = max_id + 1;
    List dp = divide_place_cpp(grid, mx, my, new_id, N);
    if (!as<bool>(dp["success"])) {
      continue; // no space
    }
    grid = as<IntegerMatrix>(dp["grid"]);
    int dx = as<int>(dp["x"]);
    int dy = as<int>(dp["y"]);

    // Mark a division event for the mother row in this hour
    div_event[irow] = 1;

    // MS and fitness for both daughters
    IntegerVector kt1, kt2;
    double f1v = f[irow], f2v = f[irow];
    apply_ms(irow, kt1, kt2, f1v, f2v);

    // Determine cell type
    bool is_tumor = (Label[irow] == 1);

    // Update mother in place -> daughter #1
    for (int c = 0; c < 22; ++c) K(irow, c) = kt1[c];
    double cd = 0.0;
    if (is_tumor) {
      // Tumor uses MS-derived fitness and tumor glucose rule
      f[irow]  = f1v;
      G[irow]  = Ctd * f1v;
      // Tumor cooldown remains 0
      Cooldown[irow] = 0.0;
    } else {
      // Normal keeps baseline fitness and uses normal glucose Cg
      // (ignore MS-derived fitness for normals)
      f[irow]  = f[irow];  // keep baseline (typically 1.0)
      G[irow]  = Cg;
      division[irow] = 0;  // normal mother leaves divisible state after division
      // Cooldown = CooldownFactor × cycle time; if CooldownFactor==0, disable cooldown
      if (CooldownFactor > 0.0) {
        cd = CooldownFactor * DivisionTime[irow];
      } else {
        cd = 0.0;
      }
      Cooldown[irow] = cd;
    }
    Time[irow] = 0;

    // Append daughter #2 as a new row
    int new_n = n + 1;
    id.push_back(new_id);
    X.push_back(dx); Y.push_back(dy);
    Label.push_back(Label[irow]);
    Status.push_back(1);
    division.push_back( is_tumor ? 1 : 0 );
    Time.push_back(0);
    DivisionTime.push_back(0.0);
    TimeToDivide.push_back(0.0);
    r.push_back(0.0);
    if (is_tumor) {
      // Tumor daughter uses MS-derived fitness and tumor glucose rule
      f.push_back(f2v);
      G.push_back(Ctd * f2v);
      Cooldown.push_back(0.0);
    } else {
      // Normal daughter keeps baseline fitness and uses normal glucose Cg
      f.push_back(f[irow]);
      G.push_back(Cg);
      // For normals, daughter inherits cooldown window
      Cooldown.push_back(cd);
    }
    Mr.push_back(Mr[irow]);

    // Extend per-hour allocation & event flags to keep vector sizes consistent
    g_alloc.push_back(0.0);
    div_event.push_back(0);
    death_event.push_back(0);
    risk_div.push_back(0);
    death_resource.push_back(0);
    death_random.push_back(0);
    death_timeout.push_back(0);

    // Resize K (n+1 x 22)
    IntegerMatrix Knew(new_n, 22);
    for (int rr = 0; rr < n; ++rr) for (int c = 0; c < 22; ++c) Knew(rr, c) = K(rr, c);
    for (int c = 0; c < 22; ++c) Knew(n, c) = kt2[c];
    K = Knew; // assign

    // Update id2row for appended row
    id2row.push_back(new_n); // one-based row index for new_id at tail

    // Update loop state
    n = new_n;
    max_id = new_id;
  }

  // (D) Bulk r update for all alive rows
  // Build per-id arrays for calc based on current row state
  // active_ids: all rows with Status==1 and with valid id
  std::vector<int> active_vec; active_vec.reserve(n);
  IntegerVector x_by_id(max_id + 1), y_by_id(max_id + 1), label_by_id(max_id + 1, -1), alive_by_id(max_id + 1);
  NumericVector f_by_id(max_id + 1);
  for (int irow = 0; irow < n; ++irow) {
    int myid = id[irow];
    if (IntegerVector::is_na(myid)) continue;
    if (myid < 0 || myid > max_id) continue;
    id2row[myid] = irow + 1;
    if (Status[irow] == 1) {
      active_vec.push_back(myid);
      x_by_id[myid] = X[irow];
      y_by_id[myid] = Y[irow];
      label_by_id[myid] = Label[irow];
      alive_by_id[myid] = 1;
      f_by_id[myid]     = f[irow];
    }
  }
  IntegerVector active_ids = wrap(active_vec);
  NumericVector r_active = calc_rates_all_cpp(grid, active_ids,
                                              x_by_id, y_by_id,
                                              label_by_id, f_by_id, alive_by_id,
                                              Rt, Rn, Cg, Fd, Tg, radius);
  // assign back to rows
  for (int k = 0; k < active_ids.size(); ++k) {
    int aid = active_ids[k];
    int row1 = id2row[aid];
    if (!IntegerVector::is_na(row1)) {
      int rr = row1 - 1;
      r[rr] = r_active[k];
    }
  }
  for (int irow = 0; irow < n; ++irow) if (Status[irow] == 1) {
    double rv = r[irow]; if (rv < 1e-12) rv = 1e-12;
    DivisionTime[irow] = 1.0 / rv;
    TimeToDivide[irow] = DivisionTime[irow] - (double)Time[irow];
  }

  // (E) Death updates
  // First compute per-cell hourly glucose allocation g_alloc[i]
  double total_G = 0.0;
  int    alive_cnt = 0;
  for (int irow = 0; irow < n; ++irow) if (Status[irow] == 1) { total_G += G[irow]; ++alive_cnt; }
  // Lowercase supply_mode for case-insensitive comparison
  std::string mode_lc = supply_mode;
  std::transform(mode_lc.begin(), mode_lc.end(), mode_lc.begin(), ::tolower);
  if (alive_cnt > 0) {
    if (mode_lc == "proportional") {
      // --- Rank-based culling under proportional mode ---
      // Each alive cell has a "need" = frac_i * G[i].
      // Tie-break policy: random tie-break among equal-need cells using a shuffled order before a stable sort.
      // If total need exceeds supply, we cull the largest-need cells first
      // (ties broken at random) until the remaining total need fits the supply.
      // This enforces that high-demand cells are more vulnerable under scarcity.

      // 1) Collect alive indices and compute per-cell need
      std::vector<int> alive_idx; alive_idx.reserve(alive_cnt);
      std::vector<double> need;   need.reserve(alive_cnt);
      double total_need = 0.0;
      for (int irow = 0; irow < n; ++irow) if (Status[irow] == 1) {
        const double frac = (Label[irow] == 1) ? DTr : DNr;
        const double nd = frac * G[irow];
        alive_idx.push_back(irow);
        need.push_back(nd);
        total_need += nd;
      }

      // 2) If feasible, nobody dies
      if (total_need <= total_supply) {
        // Keep g_alloc as already computed (proportional to G/total_G)
        // If total_G > 0, assign proportional allocation; else split equally
        if (total_G > 0) {
          for (int irow = 0; irow < n; ++irow) if (Status[irow] == 1) {
            g_alloc[irow] = total_supply * (G[irow] / total_G);
          }
        } else {
          double base = total_supply / static_cast<double>(alive_cnt);
          for (int irow = 0; irow < n; ++irow) if (Status[irow] == 1) g_alloc[irow] = base;
        }
        // Nothing else to do here.
      } else {
        // 3) Sort by need ASC with random tie-break to avoid systematic bias
        //    We shuffle first, then stable-sort by need.
        {
          std::random_device rd2;
          std::mt19937 gen2(rd2());
          std::shuffle(alive_idx.begin(), alive_idx.end(), gen2);
        }
        // Build pair (idx, need) using the shuffled order
        std::vector<std::pair<int,double>> pairs; pairs.reserve(alive_idx.size());
        for (size_t k = 0; k < alive_idx.size(); ++k) {
          int irow = alive_idx[k];
          const double frac = (Label[irow] == 1) ? DTr : DNr;
          pairs.emplace_back(irow, frac * G[irow]);
        }
        std::stable_sort(pairs.begin(), pairs.end(),
                         [](const std::pair<int,double>& a, const std::pair<int,double>& b){
                           return a.second < b.second; // ascending by need
                         });

        // 4) Keep smallest needs until capacity; mark the rest as dead (resource)
        double cum_need = 0.0;
        std::vector<char> keep(n, 0);
        for (size_t k = 0; k < pairs.size(); ++k) {
          int irow = pairs[k].first;
          double nd = pairs[k].second;
          if (cum_need + nd <= total_supply) {
            keep[irow] = 1;
            cum_need += nd;
          } else {
            keep[irow] = 0;
          }
        }

        // 5) Apply deaths for those not kept
        for (size_t k = 0; k < pairs.size(); ++k) {
          int irow = pairs[k].first;
          if (keep[irow] == 0) {
            int gx = X[irow], gy = Y[irow];
            if (gx >= 1 && gx <= N && gy >= 1 && gy <= N) {
              int gid = grid(gx - 1, gy - 1);
              if (!IntegerVector::is_na(gid) && gid == id[irow]) grid(gx - 1, gy - 1) = NA_INTEGER;
            }
            Status[irow] = 0;
            death_event[irow] = 1;
            death_resource[irow] = 1;
            death_random[irow] = 0;
            death_timeout[irow] = 0;
            g_alloc[irow] = 0.0; // no allocation for culled rows
          }
        }

        // 6) Optionally recompute g_alloc for survivors as proportional among survivors
        double total_G_surv = 0.0;
        for (int irow = 0; irow < n; ++irow) if (Status[irow] == 1) total_G_surv += G[irow];
        if (total_G_surv > 0) {
          for (int irow = 0; irow < n; ++irow) if (Status[irow] == 1) {
            g_alloc[irow] = total_supply * (G[irow] / total_G_surv);
          }
        } else {
          // no survivors; nothing to reassign
        }
      }
    } else {
      // equal mode
      double base = total_supply / static_cast<double>(alive_cnt);
      for (int irow = 0; irow < n; ++irow) if (Status[irow] == 1) g_alloc[irow] = base;
    }
  }
  // Then perform resource-based death checks using per-cell allocation.
  // Proportional mode: compare g_alloc[i] to (DTr/DNr) * G[i] (cell-specific).
  // Equal mode: compare equal base allocation to (DTr/DNr) * G[i].
  if (alive_cnt > 0) {
    if (mode_lc == "proportional") {
      // (resource-based deaths already handled in rank-based culling above)
    } else {
      const double base = total_supply / static_cast<double>(alive_cnt); // absolute supply per cell
      for (int irow = 0; irow < n; ++irow) if (Status[irow] == 1) {
        const double frac = (Label[irow] == 1) ? DTr : DNr;
        if (base < G[irow] * frac) {
          int gx = X[irow], gy = Y[irow];
          if (gx >= 1 && gx <= N && gy >= 1 && gy <= N) {
            int gid = grid(gx - 1, gy - 1);
            if (!IntegerVector::is_na(gid) && gid == id[irow]) grid(gx - 1, gy - 1) = NA_INTEGER;
          }
          Status[irow] = 0;
          death_event[irow] = 1;
          death_resource[irow] = 1; // mark resource-based death
          death_random[irow] = 0;
          death_timeout[irow] = 0;
        }
      }
    }
  }
  // Random daily
  if (step % 24 == 0) {
    for (int irow = 0; irow < n; ++irow) if (Status[irow] == 1) {
      double p = (Label[irow] == 1) ? T_death_rate : N_death_rate;
      if (unif_rand() <= p) {
        int gx = X[irow], gy = Y[irow];
        if (gx >= 1 && gx <= N && gy >= 1 && gy <= N) {
          int gid = grid(gx - 1, gy - 1);
          if (!IntegerVector::is_na(gid) && gid == id[irow]) grid(gx - 1, gy - 1) = NA_INTEGER;
        }
        Status[irow] = 0;
        death_event[irow] = 1;
        death_random[irow] = 1;
        death_timeout[irow] = 0;
        death_resource[irow] = 0;
      }
    }
  }
  // Timeout death (10x cycle as per your earlier rule; adjust factor if needed)
  for (int irow = 0; irow < n; ++irow) if (Status[irow] == 1 && division[irow] == 1) {
    if (std::abs(TimeToDivide[irow]) > 10.0 * DivisionTime[irow]) {
      int gx = X[irow], gy = Y[irow];
      if (gx >= 1 && gx <= N && gy >= 1 && gy <= N) {
        int gid = grid(gx - 1, gy - 1);
        if (!IntegerVector::is_na(gid) && gid == id[irow]) grid(gx - 1, gy - 1) = NA_INTEGER;
      }
      Status[irow] = 0;
      death_event[irow] = 1;
      death_random[irow] = 0;
      death_timeout[irow] = 1;
      death_resource[irow] = 0;
    }
  }

  // (F) Migration
  for (int irow = 0; irow < n; ++irow) if (Status[irow] == 1) {
    if (unif_rand() <= Mr[irow]) {
      List mv = migrate_one_cpp(grid, X[irow], Y[irow], id[irow], N);
      if (as<bool>(mv["success"])) {
        grid = as<IntegerMatrix>(mv["grid"]);
        X[irow] = as<int>(mv["x"]);
        Y[irow] = as<int>(mv["y"]);
      }
    }
  }

  // Return updated state
  return List::create(
    _["grid"] = grid,
    _["id"] = id,
    _["X"] = X,
    _["Y"] = Y,
    _["Label"] = Label,
    _["Status"] = Status,
    _["division"] = division,
    _["Time"] = Time,
    _["DivisionTime"] = DivisionTime,
    _["TimeToDivide"] = TimeToDivide,
    _["Cooldown"] = Cooldown,
    _["r"] = r,
    _["f"] = f,
    _["G"] = G,
    _["Mr"] = Mr,
    _["K"] = K,
    _["max_id"] = max_id,
    _["g_alloc"] = g_alloc,
    _["div_event"] = div_event,
    _["death_event"] = death_event,
    _["risk_div"] = risk_div,
    _["death_resource"] = death_resource,
    _["death_random"] = death_random,
    _["death_timeout"] = death_timeout
  );
}
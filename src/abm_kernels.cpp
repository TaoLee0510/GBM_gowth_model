// src/abm_kernels.cpp
#include <Rcpp.h>
using namespace Rcpp;

inline bool isNAint(int v) {
  return Rcpp::IntegerVector::is_na(v);
}

// [[Rcpp::export]]
IntegerVector moore_free_neighbors(const IntegerMatrix& grid, int x, int y) {
  // Return concatenated vector c(x1,y1,x2,y2,...) of free Moore neighbors
  // x,y are 1-based (R-style)
  const int N = grid.nrow();
  IntegerVector res; res.reserve(16);

  for (int dx = -1; dx <= 1; ++dx) {
    for (int dy = -1; dy <= 1; ++dy) {
      if (dx == 0 && dy == 0) continue;
      int nx = x + dx;
      int ny = y + dy;
      if (nx >= 1 && nx <= N && ny >= 1 && ny <= N) {
        int val = grid(nx - 1, ny - 1); // 0-based C++ index
        if (isNAint(val)) {
          res.push_back(nx);
          res.push_back(ny);
        }
      }
    }
  }
  return res;
}

// [[Rcpp::export]]
List neighbour_counts_idx(const IntegerMatrix& grid,
                          const IntegerVector& label_by_row,
                          int idx,
                          int radius) {
  // Count Nt, Nn within a Moore window of given radius around the cell with row id=idx
  // grid stores row indices (1..nrow(Cells)) or NA for empty
  // label_by_row is length nrow(Cells), 0=normal, 1=tumour
  // idx is a row id (1-based) existing in grid
  const int N = grid.nrow();

  // find position of idx in grid
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
        int lab = label_by_row[id - 1]; // id -> row -> label
        if (lab == 1) ++Nt;
        else if (lab == 0) ++Nn;
      }
    }
  }
  return List::create(_["Nt"]=Nt, _["Nn"]=Nn);
}

// [[Rcpp::export]]
List divide_place_cpp(IntegerMatrix grid,
                      int mother_x, int mother_y,
                      int new_row_index,
                      int N) {
  // Choose a free Moore neighbor for the daughter, write grid(target)=new_row_index
  // Return: list(success, grid, target_x, target_y)
  // mother stays in place; caller will overwrite mother's row in Cells with daughter1.

  // collect free neighbors
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
  // random pick
  int k = (int) std::floor(R::unif_rand() * freepos.size());
  if (k == (int)freepos.size()) k = (int)freepos.size() - 1;
  int tx = freepos[k].first;
  int ty = freepos[k].second;

  // write daughter index into target cell
  grid(tx - 1, ty - 1) = new_row_index;

  return List::create(_["success"]=true,
                      _["grid"]=grid,
                      _["x"]=tx,
                      _["y"]=ty);
}

// [[Rcpp::export]]
List migrate_one_cpp(IntegerMatrix grid, int x, int y, int row_id, int N) {
  // Try to move a cell (row_id) from (x,y) to a random free Moore neighbor.
  // Return: list(success, grid, x, y) with updated position if moved.
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

  // move: clear old, set new
  grid(x - 1, y - 1) = NA_INTEGER;
  grid(tx - 1, ty - 1) = row_id;

  return List::create(_["success"]=true, _["grid"]=grid, _["x"]=tx, _["y"]=ty);
}
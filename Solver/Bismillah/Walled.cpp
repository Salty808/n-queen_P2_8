#include "nqueens.h"
#include <string>
#include <utility>
#include <vector>


struct Bits {
  unsigned long long a, b, c;
};
static inline bool isEmpty(const Bits& x) { return (x.a | x.b | x.c) == 0; }
static inline Bits andNot(const Bits& x, const Bits& y) {
  Bits r; r.a = x.a & ~y.a; r.b = x.b & ~y.b; r.c = x.c & ~y.c; return r;
}
static inline Bits andB(const Bits& x, const Bits& y) {
  Bits r; r.a = x.a & y.a; r.b = x.b & y.b; r.c = x.c & y.c; return r;
}
static inline int lowest(const Bits& x) {
  if (x.a) return __builtin_ctzll(x.a);
  if (x.b) return 64 + __builtin_ctzll(x.b);
  return 128 + __builtin_ctzll(x.c);
}
static inline void setBit(Bits& x, int i) {
  if (i < 64) x.a |= 1ULL << i; else if (i < 128) x.b |= 1ULL << (i - 64); else x.c |= 1ULL << (i - 128);
}
static inline void clearBit(Bits& x, int i) {
  if (i < 64) x.a &= ~(1ULL << i); else if (i < 128) x.b &= ~(1ULL << (i - 64)); else x.c &= ~(1ULL << (i - 128));
}

// ---------- Jalur cepat: tanpa wall = N-Queens biasa (bitmask) ----------
struct BitSearch {
  int n;
  unsigned full;
  int col[16];
  std::vector<std::vector<std::pair<int, int> > >* out;
  void go(int r, unsigned cols, unsigned d1, unsigned d2) {
    if (r == n) {
      std::vector<std::pair<int, int> > b;
      b.reserve(n);
      for (int i = 0; i < n; ++i) b.push_back(std::make_pair(i, col[i]));
      out->push_back(std::move(b));
      return;
    }
    unsigned avail = full & ~(cols | d1 | d2);
    while (avail) {
      unsigned bit = avail & (~avail + 1u);
      avail ^= bit;
      int c = 0;
      while (!((bit >> c) & 1u)) ++c;
      col[r] = c;
      go(r + 1, cols | bit, ((d1 | bit) << 1) & full, (d2 | bit) >> 1);
    }
  }
};

// ---------- Jalur umum: maximum independent set pada graf konflik ----------
struct WallSearch {
  int n, V, best;
  bool enumerate;
  int cellOf[144];
  Bits adj[144];
  int chosen[160];
  std::vector<std::vector<std::pair<int, int> > > out;

  void build(int n_, const std::vector<std::string>& grid) {
    n = n_;
    int vid[144];
    V = 0;
    for (int i = 0; i < n * n; ++i) vid[i] = -1;
    int walls = 0;
    for (int r = 0; r < n; ++r)
      for (int c = 0; c < n; ++c)
        if (grid[r][c] == 'W') ++walls;
    // Penomoran sel menentukan kualitas batas atas (clique cover greedy).
    // Wall sedikit: urutan baris-demi-baris. Wall banyak: derajat konflik menaik.
    const int dr[4] = {0, 1, 1, 1};
    const int dc[4] = {1, 0, 1, -1};
    int key[144], idx[144], m = 0;
    for (int r = 0; r < n; ++r)
      for (int c = 0; c < n; ++c) {
        if (grid[r][c] == 'W') continue;
        int deg = 0;
        if (walls >= n) {
          for (int f = 0; f < 4; ++f)
            for (int sg = -1; sg <= 1; sg += 2) {
              int rr = r + sg * dr[f], cc = c + sg * dc[f];
              while (rr >= 0 && rr < n && cc >= 0 && cc < n && grid[rr][cc] != 'W') {
                ++deg;
                rr += sg * dr[f];
                cc += sg * dc[f];
              }
            }
        }
        idx[m] = r * n + c;
        key[m] = deg;
        ++m;
      }
    for (int i = 1; i < m; ++i) {  // insertion sort stabil
      int k = key[i], x = idx[i], j = i - 1;
      while (j >= 0 && key[j] > k) { key[j + 1] = key[j]; idx[j + 1] = idx[j]; --j; }
      key[j + 1] = k;
      idx[j + 1] = x;
    }
    for (int i = 0; i < m; ++i) { vid[idx[i]] = V; cellOf[V] = idx[i]; ++V; }
    for (int v = 0; v < V; ++v) { adj[v].a = adj[v].b = adj[v].c = 0; }
    // dua sel konflik jika segaris dan tidak ada wall di antaranya
    for (int v = 0; v < V; ++v) {
      int r = cellOf[v] / n, c = cellOf[v] % n;
      for (int f = 0; f < 4; ++f) {
        for (int s = 1; s <= 2; ++s) {
          int sg = (s == 1) ? 1 : -1;
          int rr = r + sg * dr[f], cc = c + sg * dc[f];
          while (rr >= 0 && rr < n && cc >= 0 && cc < n && grid[rr][cc] != 'W') {
            setBit(adj[v], vid[rr * n + cc]);
            rr += sg * dr[f];
            cc += sg * dc[f];
          }
        }
      }
    }
  }

  void record(int size) {
    std::vector<std::pair<int, int> > b;
    b.reserve(size);
    for (int i = 0; i < size; ++i) b.push_back(std::make_pair(cellOf[chosen[i]] / n, cellOf[chosen[i]] % n));
    out.push_back(std::move(b));
  }

  void expand(Bits P, int size) {
    if (isEmpty(P)) {
      if (!enumerate) {
        if (size > best) best = size;
      } else if (size == best) {
        record(size);
      }
      return;
    }
    unsigned char order[144], col[144];
    int cnt = 0, k = 0;
    Bits U = P;
    while (!isEmpty(U)) {
      ++k;
      Bits Q = U;
      while (!isEmpty(Q)) {
        int v = lowest(Q);
        order[cnt] = (unsigned char)v;
        col[cnt] = (unsigned char)k;
        ++cnt;
        clearBit(U, v);
        Q = andB(Q, adj[v]);
      }
    }
    const int need = enumerate ? best : best + 1;
    for (int i = cnt - 1; i >= 0; --i) {
      if (size + col[i] < need) return;
      int v = order[i];
      chosen[size] = v;
      Bits Q = andNot(P, adj[v]);
      clearBit(Q, v);
      expand(Q, size + 1);
      clearBit(P, v);
    }
  }

  void solve() {
    Bits all; all.a = all.b = all.c = 0;
    for (int v = 0; v < V; ++v) setBit(all, v);
    best = 0;
    enumerate = false;
    expand(all, 0);
    enumerate = true;
    expand(all, 0);
  }
};

std::vector<std::vector<std::pair<int, int> > > solveWalled(
    int n, const std::vector<std::string>& grid) {
  bool anyWall = false;
  for (int r = 0; r < n && !anyWall; ++r)
    for (int c = 0; c < n; ++c)
      if (grid[r][c] == 'W') { anyWall = true; break; }
  if (!anyWall) {
    std::vector<std::vector<std::pair<int, int> > > res;
    BitSearch bs;
    bs.n = n;
    bs.full = (1u << n) - 1u;
    bs.out = &res;
    bs.go(0, 0, 0, 0);
    return res;
  }
  WallSearch ws;
  ws.build(n, grid);
  ws.solve();
  return std::move(ws.out);
}
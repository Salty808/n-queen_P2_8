#include "nqueens.h"
#include <vector>
#include <string>
#include <utility>
#include <algorithm>

namespace {
const int WORDS = 3;
const int MAXV = WORDS * 64;
typedef unsigned long long u64;

struct Set {
    u64 w[WORDS];
};

int N;
int cellOf[MAXV];
Set compat[MAXV];
int best;
bool collect;
int cur[MAXV];
std::vector<std::vector<std::pair<int, int>>>* out;

inline bool empty(const Set& s) {
    for (int i = 0; i < WORDS; i++) if (s.w[i]) return false;
    return true;
}

inline int first(const Set& s) {
    for (int i = 0; i < WORDS; i++) if (s.w[i]) return i * 64 + __builtin_ctzll(s.w[i]);
    return -1;
}

inline void clearBit(Set& s, int v) {
    s.w[v >> 6] &= ~(1ULL << (v & 63));
}

inline void setBit(Set& s, int v) {
    s.w[v >> 6] |= 1ULL << (v & 63);
}

inline bool pruned(int bound) {
    return collect ? bound < best : bound <= best;
}

void leaf(int size) {
    if (!collect) {
        if (size > best) best = size;
        return;
    }
    if (size == best) {
        std::vector<std::pair<int, int>> board;
        board.reserve(size);
        for (int i = 0; i < size; i++) board.emplace_back(cellOf[cur[i]] / N, cellOf[cur[i]] % N);
        out->push_back(std::move(board));
    }
}

void expand(int cnt, Set P) {
    unsigned char order[MAXV];
    unsigned char colr[MAXV];
    int m = 0;
    int kmin = best - cnt + (collect ? 0 : 1);
    Set U = P;
    int k = 0;
    while (!empty(U)) {
        k++;
        Set Q = U;
        while (!empty(Q)) {
            int v = first(Q);
            clearBit(Q, v);
            clearBit(U, v);
            for (int i = 0; i < WORDS; i++) Q.w[i] &= ~compat[v].w[i];
            if (k >= kmin) {
                order[m] = (unsigned char)v;
                colr[m] = (unsigned char)k;
                m++;
            }
        }
    }
    for (int i = m - 1; i >= 0; i--) {
        if (pruned(cnt + colr[i])) return;
        int v = order[i];
        cur[cnt] = v;
        Set np;
        for (int j = 0; j < WORDS; j++) np.w[j] = P.w[j] & compat[v].w[j];
        if (empty(np)) leaf(cnt + 1);
        else expand(cnt + 1, np);
        clearBit(P, v);
    }
}

inline int absv(int x) {
    return x < 0 ? -x : x;
}
}

std::vector<std::vector<std::pair<int, int>>> solveWalled(int n, const std::vector<std::string>& grid) {
    std::vector<std::vector<std::pair<int, int>>> res;
    if (n < 1 || n * n > MAXV) return res;
    N = n;

    int total = n * n;
    std::vector<char> isFree(total);
    int V = 0;
    for (int r = 0; r < n; r++)
        for (int c = 0; c < n; c++) {
            isFree[r * n + c] = grid[r][c] != 'W';
            if (isFree[r * n + c]) V++;
        }
    if (V == 0) {
        res.push_back(std::vector<std::pair<int, int>>());
        return res;
    }

    std::vector<int> rowSeg(total, -1), colSeg(total, -1), d1Seg(total, -1), d2Seg(total, -1);
    int nr = 0, nc = 0, n1 = 0, n2 = 0;
    for (int r = 0; r < n; r++) {
        for (int c = 0; c < n; c++) {
            int i = r * n + c;
            if (!isFree[i]) continue;
            rowSeg[i] = (c > 0 && isFree[i - 1]) ? rowSeg[i - 1] : nr++;
            colSeg[i] = (r > 0 && isFree[i - n]) ? colSeg[i - n] : nc++;
            d1Seg[i] = (r > 0 && c > 0 && isFree[i - n - 1]) ? d1Seg[i - n - 1] : n1++;
            d2Seg[i] = (r > 0 && c < n - 1 && isFree[i - n + 1]) ? d2Seg[i - n + 1] : n2++;
        }
    }

    std::vector<int> cells;
    cells.reserve(V);
    for (int i = 0; i < total; i++)
        if (isFree[i]) cells.push_back(i);

    if (total - V > n / 2) {
        std::stable_sort(cells.begin(), cells.end(), [&](int x, int y) {
            int dx = std::max(absv(2 * (x / n) - (n - 1)), absv(2 * (x % n) - (n - 1)));
            int dy = std::max(absv(2 * (y / n) - (n - 1)), absv(2 * (y % n) - (n - 1)));
            return dx > dy;
        });
    }

    for (int i = 0; i < V; i++) {
        cellOf[i] = cells[i];
        for (int j = 0; j < WORDS; j++) compat[i].w[j] = 0;
    }
    for (int a = 0; a < V; a++) {
        int x = cells[a];
        for (int b = a + 1; b < V; b++) {
            int y = cells[b];
            bool attack = rowSeg[x] == rowSeg[y] || colSeg[x] == colSeg[y] || d1Seg[x] == d1Seg[y] || d2Seg[x] == d2Seg[y];
            if (!attack) {
                setBit(compat[a], b);
                setBit(compat[b], a);
            }
        }
    }

    Set all;
    for (int j = 0; j < WORDS; j++) all.w[j] = 0;
    for (int i = 0; i < V; i++) setBit(all, i);

    out = &res;
    best = 0;
    collect = false;
    expand(0, all);
    collect = true;
    expand(0, all);
    return res;
}
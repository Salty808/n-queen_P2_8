#include "nqueens.h"
#include <vector>
#include <utility>

namespace {
int N;
unsigned FULL;
int fixedCol[16];                 
int path[16];
std::vector<std::vector<int>>* out;

void dfs(int r, unsigned cols, unsigned d1, unsigned d2) {
    if (r == N) {
        out->emplace_back(path, path + N);
        return;
    }
    if (fixedCol[r] >= 0) {       
        path[r] = fixedCol[r];
        dfs(r + 1, cols, d1, d2);
        return;
    }
    unsigned avail = FULL & ~(cols | (d1 >> r) | (d2 >> (N - 1 - r)));
    while (avail) {
        unsigned bit = avail & -avail;
        avail ^= bit;
        int c = __builtin_ctz(bit);
        path[r] = c;
        dfs(r + 1, cols | bit, d1 | (1u << (r + c)), d2 | (1u << (c - r + N - 1)));
    }
}
}

std::vector<std::vector<int>> solvePreplaced(int n, const std::vector<std::pair<int, int>>& preplaced) {
    std::vector<std::vector<int>> res;
    N = n;
    FULL = (1u << n) - 1;
    for (int i = 0; i < n; i++) fixedCol[i] = -1;

    unsigned cols = 0, d1 = 0, d2 = 0;
    for (const auto& p : preplaced) {
        int r = p.first, c = p.second;
        if (r < 0 || r >= n || c < 0 || c >= n) return res;
        if (fixedCol[r] >= 0) return res;                       
        unsigned bc = 1u << c, b1 = 1u << (r + c), b2 = 1u << (c - r + n - 1);
        if ((cols & bc) || (d1 & b1) || (d2 & b2)) return res;  
        cols |= bc; d1 |= b1; d2 |= b2;
        fixedCol[r] = c;
    }
    out = &res;
    dfs(0, cols, d1, d2);
    return res;
}
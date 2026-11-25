#include "nqueens.h"
#include <vector>
#include <utility>

namespace {
int N;
unsigned FULL;
int path[16];
bool mirror;                                  
std::vector<std::vector<int>>* out;

inline void emit() {
    out->emplace_back(path, path + N);
    if (mirror) {
        std::vector<int> m(N);
        for (int i = 0; i < N; i++) m[i] = N - 1 - path[i];
        out->push_back(std::move(m));
    }
}

void dfs(int r, unsigned cols, unsigned d1, unsigned d2) {
    unsigned avail = FULL & ~(cols | d1 | d2);
    if (r == N - 1) {                         
        if (avail) {
            path[r] = __builtin_ctz(avail);
            emit();
        }
        return;
    }
    while (avail) {
        unsigned bit = avail & -avail;
        avail ^= bit;
        path[r] = __builtin_ctz(bit);
        dfs(r + 1, cols | bit, ((d1 | bit) << 1) & FULL, (d2 | bit) >> 1);
    }
}

const int KNOWN[13] = {0, 1, 0, 0, 2, 10, 4, 40, 92, 352, 724, 2680, 14200};  
}

std::vector<std::vector<int>> solveNQueens(int n) {
    std::vector<std::vector<int>> res;
    if (n < 1 || n > 16) return res;
    if (n == 1) { res.push_back(std::vector<int>(1, 0)); return res; }
    if (n <= 3) return res;

    N = n;
    FULL = (1u << n) - 1;
    out = &res;
    if (n <= 12) res.reserve(KNOWN[n]);       

    int half = n / 2;                       
    mirror = true;
    for (int c = 0; c < half; c++) {
        unsigned bit = 1u << c;
        path[0] = c;
        dfs(1, bit, (bit << 1) & FULL, bit >> 1);
    }
    if (n & 1) {                              
        mirror = false;
        unsigned bit = 1u << half;
        path[0] = half;
        dfs(1, bit, (bit << 1) & FULL, bit >> 1);
    }
    return res;
}
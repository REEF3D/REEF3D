// Unit test: the cuts of the weighted partition with the capacities of the ranks (M 20 5; src/dmpart_cuts.h,
// REEF3D M 42, M 43, G 46; DIVEMesh M 36).
//
//   c++ -std=c++20 -O2 -I ../../src partition_cuts_test.cpp -o partition_cuts_test && ./partition_cuts_test
//
// Covers:
//   - cuts_1d with uniform weights and the capacities 1, 1, 1, 0.5: interval lengths 2 : 2 : 2 : 1
//   - cuts_1d with all capacities 1 (the exact test): never above the greedy sweep without capacities, the least
//     bottleneck (brute force); the partition with all capacities equal goes without them (optimise, below)
//   - cuts_1d with capacities: the bottleneck is the least work / capacity (brute force over all cuts, small N)
//   - optimise: 4 x 4 x 1 subdomains, the ranks 12-15 (the x-slab aa = 4) at capacity 0.5: that slab half as wide,
//     the y-cuts equal; all capacities 1 the cuts without capacities, bitwise
//   - optimise: 7 x 2 x 1 (14 ranks, 10 P + 4 E: the x-slabs 6, 7) on a grid with dry cells
//
// Architect: Hans Bihs

#include "dmpart_cuts.h"

#include <cstdio>
#include <cstdlib>
#include <functional>
#include <random>
#include <vector>

using namespace dmpart;

static int fails = 0;
#define CHECK(c, ...)                                        \
    do                                                       \
    {                                                        \
        if(!(c))                                             \
        {                                                    \
            ++fails;                                         \
            std::printf("FAIL %s:%d: ", __FILE__, __LINE__); \
            std::printf(__VA_ARGS__);                        \
            std::printf("\n");                               \
        }                                                    \
    }                                                        \
    while(0)

// the least bottleneck of [0,N) in m intervals, all cuts tried (small N), each at least wmin wide
static double brute(const std::vector<std::vector<double>> &P, int N, int m, int wmin, const std::vector<std::vector<double>> &cs)
{
    double best = 1e300;
    std::vector<int> c(m + 1);
    c[0] = 0;
    c[m] = N;
    // recursive enumeration
    std::function<void(int)> rec = [&](int k)
    {
        if(k == m)
        {
            double mx = 0.0;
            for(int n = 0; n < m; ++n)
                for(size_t s = 0; s < P.size(); ++s)
                {
                    double sm = 0.0;
                    for(int i = c[n]; i < c[n + 1]; ++i)
                        sm += P[s][i];
                    mx = std::max(mx, sm / cs[n][s]);
                }
            best = std::min(best, mx);
            return;
        }
        for(int x = c[k - 1] + wmin; x <= N - (m - k) * wmin; ++x)
        {
            c[k] = x;
            rec(k + 1);
        }
    };
    rec(1);
    return best;
}

int main()
{
    // ---- 2 : 2 : 2 : 1
    {
        const int N = 70, m = 4;
        std::vector<std::vector<double>> P(1, std::vector<double>(N, 1.0));
        std::vector<std::vector<double>> cs = {{1.0}, {1.0}, {1.0}, {0.5}};
        std::vector<int> cut;
        const double h = cuts_1d(P, N, m, 3, cut, &cs);
        CHECK(cut.size() == 5 && cut[0] == 0 && cut[1] == 20 && cut[2] == 40 && cut[3] == 60 && cut[4] == 70,
              "cuts %d %d %d %d %d, expected 0 20 40 60 70", cut.size() > 0 ? cut[0] : -1, cut.size() > 1 ? cut[1] : -1, cut.size() > 2 ? cut[2] : -1,
              cut.size() > 3 ? cut[3] : -1, cut.size() > 4 ? cut[4] : -1);
        CHECK(std::fabs(h - 20.0) < 1e-9, "bottleneck %g, expected 20", h);
        std::printf("2:2:2:1 cuts:");
        for(int c : cut)
            std::printf(" %d", c);
        std::printf("  bottleneck %g\n", h);
    }

    // ---- the exact test with capacities 1 against the greedy sweep without: never a larger bottleneck, the same
    //      bottleneck where the greedy sweep finds the least one (brute force)
    {
        std::mt19937 rng(4711);
        std::uniform_real_distribution<double> U(0.0, 3.0);
        int ok = 0, tried = 0, better = 0;
        for(int t = 0; t < 200; ++t)
        {
            const int N = 9 + int(rng() % 14), m = 2 + int(rng() % 2), ns = 1 + int(rng() % 3);
            if(N < 3 * m)
                continue;
            std::vector<std::vector<double>> P(ns, std::vector<double>(N));
            for(auto &row : P)
                for(double &v : row)
                    v = (rng() % 5 == 0) ? 0.0 : U(rng);
            std::vector<std::vector<double>> ones(m, std::vector<double>(ns, 1.0));
            std::vector<int> c0, c1;
            const double h0 = cuts_1d(P, N, m, 3, c0);
            const double h1 = cuts_1d(P, N, m, 3, c1, &ones);
            const double b = brute(P, N, m, 3, ones);
            ++tried;
            if(h1 <= h0 * (1.0 + 1e-12) && std::fabs(h1 - b) <= 1e-9 * b)
                ++ok;
            if(h1 < h0 * (1.0 - 1e-9))
                ++better;
        }
        CHECK(ok == tried, "exact test: %d of %d not the least bottleneck", tried - ok, tried);
        std::printf("exact test, capacities 1: the least bottleneck in %d of %d (the greedy sweep above it in %d)\n", ok, tried, better);
    }

    // ---- with capacities: the least bottleneck (brute force)
    {
        std::mt19937 rng(17);
        std::uniform_real_distribution<double> U(0.1, 2.0);
        int ok = 0, tried = 0;
        double worst = 0.0;
        for(int t = 0; t < 300; ++t)
        {
            const int N = 9 + int(rng() % 14), m = 2 + int(rng() % 2), ns = 1 + int(rng() % 2);
            if(N < 3 * m)
                continue;
            std::vector<std::vector<double>> P(ns, std::vector<double>(N));
            for(auto &row : P)
                for(double &v : row)
                    v = U(rng);
            std::vector<std::vector<double>> cs(m, std::vector<double>(ns));
            for(int n = 0; n < m; ++n)
                for(int s = 0; s < ns; ++s)
                    cs[n][s] = (rng() % 2) ? 1.0 : 0.3 + 0.6 * double(rng() % 100) / 100.0;
            std::vector<int> cut;
            const double h = cuts_1d(P, N, m, 3, cut, &cs);
            const double b = brute(P, N, m, 3, cs);
            ++tried;
            const double rel = (h - b) / b;
            worst = std::max(worst, rel);
            if((int)cut.size() == m + 1 && rel < 1e-9)
                ++ok;
        }
        CHECK(ok == tried, "least bottleneck missed in %d of %d (worst %.3g above)", tried - ok, tried, worst);
        std::printf("least bottleneck with capacities: %d of %d\n", ok, tried);
    }

    // ---- optimise: 4 x 4 x 1, E-ranks 12-15 (x-slab 4) at 0.5
    {
        const int Nd[3] = {140, 64, 4}, M[3] = {4, 4, 1};
        std::vector<double> w((size_t)Nd[0] * Nd[1] * Nd[2], 1.0);
        std::vector<double> cap(16, 1.0);
        for(int r = 12; r < 16; ++r)
            cap[r] = 0.5;
        std::vector<int> C[3];
        double h0, mean;
        const double h = optimise(w, Nd, M, C, h0, mean, cap);
        CHECK(C[0].size() == 5 && C[0][1] == 40 && C[0][2] == 80 && C[0][3] == 120 && C[0][4] == 140, "x-cuts %d %d %d %d", C[0][1], C[0][2], C[0][3], C[0][4]);
        CHECK(C[1].size() == 5 && C[1][1] == 16 && C[1][2] == 32 && C[1][3] == 48, "y-cuts %d %d %d", C[1][1], C[1][2], C[1][3]);
        CHECK(std::fabs(h / mean - 1.0) < 1e-12, "slowest / mean %g, expected 1", h / mean);
        CHECK(std::fabs(h0 / mean - (140.0 * 64 * 4 / 16 / 0.5) / mean) < 1e-9, "equal widths slowest / mean %g", h0 / mean);
        std::printf("4x4x1, E-slab at 0.5: x-cuts 0 %d %d %d %d, slowest/mean %g (equal widths %g)\n", C[0][1], C[0][2], C[0][3], C[0][4], h / mean, h0 / mean);

        // all capacities 1: bitwise the partition without capacities (random weights)
        std::mt19937 rng(5);
        std::uniform_real_distribution<double> U(0.0, 1.0);
        for(double &v : w)
            v = U(rng) < 0.2 ? 0.0 : U(rng);
        std::vector<int> A[3], B[3];
        double ha0, hb0, ma, mb;
        const double ha = optimise(w, Nd, M, A, ha0, ma);
        const double hb = optimise(w, Nd, M, B, hb0, mb, std::vector<double>(16, 1.0));
        CHECK(A[0] == B[0] && A[1] == B[1] && A[2] == B[2] && ha == hb && ma == mb && ha0 == hb0, "optimise: capacities 1 differ from none");
    }

    // ---- optimise: 14 ranks 7 x 2 x 1, 10 P + 4 E (x-slabs 6, 7 at 0.4), dry cells at the end of the domain
    {
        const int Nd[3] = {300, 40, 1}, M[3] = {7, 2, 1};
        std::vector<double> w((size_t)Nd[0] * Nd[1] * Nd[2], 1.0);
        for(int i = 260; i < 300; ++i)
            for(int j = 0; j < 20; ++j)
                w[(size_t)i * Nd[1] + j] = 0.0;  // a dry corner
        std::vector<double> cap(14, 1.0);
        for(int r = 10; r < 14; ++r)
            cap[r] = 0.4;
        std::vector<int> C[3];
        double h0, mean;
        const double h = optimise(w, Nd, M, C, h0, mean, cap);
        CHECK(h / mean < 1.1, "7x2: slowest / mean %g", h / mean);
        CHECK(h < h0, "7x2: no gain over equal widths");
        std::printf("7x2x1, 4 E at 0.4, dry corner: slowest/mean %g (equal widths %g); x-cuts", h / mean, h0 / mean);
        for(int c : C[0])
            std::printf(" %d", c);
        std::printf("; y-cut %d\n", C[1][1]);
    }

    if(fails)
    {
        std::printf("%d FAILED\n", fails);
        return 1;
    }
    std::printf("all passed\n");
    return 0;
}

// Unit test: the speed of the E-ranks against the P-ranks from the windows of the rebalance (M 42, G 46;
// src/capacity_estimate.h).
//
//   c++ -std=c++20 -O2 -I ../../src capacity_estimate_test.cpp -o capacity_estimate_test && ./capacity_estimate_test
//
// Covers:
//   - 12 P + 4 E, equal work, E-ranks 0.4 times as fast: c = 0.4, spreads 1
//   - unequal work (the capacity-weighted partition: E-subdomains smaller): c from work / busy, not from busy alone
//   - noise within the classes: c of the class means, the spread max/min reported
//   - the smoothing over the windows (alpha 0.5): 0.4, then 0.5 -> 0.45, then 0.5 -> 0.475
//   - not measurable: no E-rank, all ranks E, a rank without work or busy time skipped, a class left empty
//
// Architect: Hans Bihs

#include "capacity_estimate.h"

#include <cmath>
#include <cstdio>
#include <vector>

static int fails = 0;
#define CHECK(c, ...) do { if(!(c)) { ++fails; std::printf("FAIL %s:%d: ", __FILE__, __LINE__); std::printf(__VA_ARGS__); std::printf("\n"); } } while(0)

static bool near(double a, double b, double tol=1e-12) { return std::fabs(a-b)<=tol*std::max(1.0,std::fabs(b)); }

int main()
{
    // ---- equal work, E 0.4 as fast
    {
        capacity_estimate ce;
        std::vector<double> W(16,1000.0), busy(16,1.0);
        for(int r=12; r<16; ++r) busy[r] = 2.5;
        const auto w = ce.update(W,busy,4);
        CHECK(near(w.c,0.4), "c %g, expected 0.4", w.c);
        CHECK(near(ce.cE,0.4) && ce.windows==1, "cE %g windows %d", ce.cE, ce.windows);
        CHECK(near(w.spreadP,1.0) && near(w.spreadE,1.0), "spreads %g %g", w.spreadP, w.spreadE);
    }

    // ---- unequal work: balanced busy times, E-subdomains 0.45 of the P-subdomains
    {
        capacity_estimate ce;
        std::vector<double> W(16,1000.0), busy(16,2.0);
        for(int r=12; r<16; ++r) W[r] = 450.0;
        const auto w = ce.update(W,busy,4);
        CHECK(near(w.c,0.45), "balanced: c %g, expected 0.45", w.c);
    }

    // ---- noise within the classes
    {
        capacity_estimate ce;
        std::vector<double> W(14,1000.0), busy = {1.0,1.1,0.9,1.0,1.05,0.95,1.0,1.0,1.0,1.0, 2.0,2.5,2.2,2.3};
        const auto w = ce.update(W,busy,4);
        double sP = 0.0, sE = 0.0;
        for(int r=0; r<10; ++r) sP += 1000.0/busy[r];
        for(int r=10; r<14; ++r) sE += 1000.0/busy[r];
        CHECK(near(w.c,(sE/4)/(sP/10)), "noise: c %g, expected %g", w.c, (sE/4)/(sP/10));
        CHECK(near(w.spreadP,1.1/0.9) && near(w.spreadE,2.5/2.0), "spreads %g %g", w.spreadP, w.spreadE);
        std::printf("10 P + 4 E with noise: c %.4f, spread P %.3f E %.3f\n",w.c,w.spreadP,w.spreadE);
    }

    // ---- smoothing
    {
        capacity_estimate ce;
        std::vector<double> W(4,100.0), b1 = {1,1,1,2.5}, b2 = {1,1,1,2.0};
        ce.update(W,b1,1);
        ce.update(W,b2,1);
        CHECK(near(ce.cE,0.45) && ce.windows==2, "smoothed %g, expected 0.45", ce.cE);
        ce.update(W,b2,1);
        CHECK(near(ce.cE,0.475) && ce.windows==3, "smoothed %g, expected 0.475", ce.cE);
    }

    // ---- not measurable
    {
        capacity_estimate ce;
        std::vector<double> W(4,100.0), busy(4,1.0);
        CHECK(ce.update(W,busy,0).c==0.0 && ce.windows==0, "no E-rank measured");
        CHECK(ce.update(W,busy,4).c==0.0 && ce.windows==0, "all E-ranks measured");
        std::vector<double> W2 = {100,0,100,100}, b2 = {1,1,1,4};
        const auto w = ce.update(W2,b2,1);
        CHECK(near(w.c,0.25) && ce.windows==1, "rank without work: c %g", w.c);
        std::vector<double> W3 = {100,100,100,0};
        CHECK(ce.update(W3,busy,1).c==0.0 && ce.windows==1 && near(ce.cE,0.25), "empty class: cE kept");
    }

    if(fails)
    {
        std::printf("%d FAILED\n",fails);
        return 1;
    }
    std::printf("all passed\n");
    return 0;
}

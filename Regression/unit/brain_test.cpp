// Standalone check of the brain of a deforming 6DOF body (X 193, src/6DOF_brain.cpp) in closed loop with a
// kinematic swimmer: no MPI, no REEF3D solver.
// Build:  g++ -O2 -std=c++17 -I../../src -I../../ThirdParty/eigen-5.0.0 -DEIGEN_MPL2_ONLY brain_test.cpp -o brain_test
// Run:    ./brain_test
// The swimmer (L = 1 m, f0 = 1 Hz, in the x-z plane): the gait follows the targets of the brain with
// critically damped dynamics, the speed of the gait f a L with 0.8 s, the speed through the water s U_g
// (s = 0.48 L per beat as the 2D swimmer of validation 08; a gliding body slows down with 2.5 s), the turn
// rate 1.3 b U_w/L, a uniform current.
// Checks: the estimate of the stride and of the current (0.2 m/s along x, along x + z, none); station
// holding in a current (the fish faces it and stays at the target; starting close to it, no glide before the
// current is known); a target ahead in still water at
// U_des (the frequency U_des/(s L)); a moving target; the own turn on the lateral line (-0.5 b q at the
// head) is no escape, a stimulus is; at a station in still water the fish glides (amplitude 0) and bursts
// (burst and coast), early enough for the coast to stop it there; the stride prior and the turn rate of a
// body (X 200): a slow swimmer with a weak turn, a body that turns away from its bias; a body still
// accelerating when the brain starts sees no headwind; a body that keeps its speed when its gait slows
// (X 201) sees no tailwind; a restart from the checkpoint of the brain continues it exactly (in a current and
// while gliding).
// Architect: Hans Bihs
#include"../../src/6DOF_brain.cpp"

#include<cstdio>
#include<iostream>
#include<string>
#include<vector>
#include<functional>
#include<array>
#include<algorithm>

static int nfail = 0, npass = 0;
static void check(bool ok, const std::string &what)
{
    if(ok) ++npass;
    else ++nfail;
    std::cout<<(ok ? "  ok    " : "  FAIL  ")<<what<<std::endl;
}

struct swimmer
{
    double x = 1.6, z = 3.0, h = 0.0, Ug = 0.0;
    double g[3] = {1.0, 1.0, 0.0}, gr[3] = {0.0, 0.0, 0.0}, gt[3] = {1.0, 1.0, 0.0}, ts[3] = {0.25, 0.25, 0.25};
    double s_true = 0.48, k_true = 1.3, tau_true = 0.8, tau_down = 0.0, tau_glide = 2.5, C_brake = 0.0, cx = 0.0, cz = 0.0;
    double stim_t = -1.0;      // a pressure pulse across the head at this time (a startle stimulus)
    bool second = false;
    sixdof_brain::decision pend;
    std::vector<sixdof_brain::decision> log;
    std::vector<sixdof_brain::record> rec;      // the record of the brain at each decision
    std::vector<std::array<double,4> > track;   // t, x, z, h
    int nesc = 0;

    void target(const sixdof_brain::decision &d, bool two)
    {
        const double f[3] = {two ? d.f2 : d.f, two ? d.a2 : d.a, two ? d.b2 : d.b};
        const double t = two ? d.ts2 : d.ts;
        for(int k=0; k<3; ++k)
        {
            gt[k] = f[k];
            ts[k] = (t>0.0) ? t : 0.5/std::max(g[0], 0.1);
        }
    }

    // advance to t1 in steps of dt, the brain sensing every step and deciding when due
    void run(sixdof_brain &B, double &t, double t1, double dt = 0.01)
    {
        const double rho = 1000.0;
        while(t<t1 - 1.0e-12)
        {
            for(int k=0; k<3; ++k)
            {
                const double lam = 4.744/ts[k];
                gr[k] += dt*(lam*lam*(gt[k] - g[k]) - 2.0*lam*gr[k]);
                g[k] += dt*gr[k];
            }
            if(second && t>=pend.t2)
            {
                target(pend, true);
                second = false;
            }
            // the speed rises with tau_true; it falls with tau_glide when gliding, with tau_down (0: tau_true)
            // when the gait only slows down
            const double gt_ = g[0]*g[1]*1.0;
            Ug += dt/(g[1]<0.5 ? tau_glide : ((gt_<Ug && tau_down>0.0) ? tau_down : tau_true))*(gt_ - Ug);
            // the bent body brakes the glide like a quadratic drag, dU/dt = -C_brake (|b|/2.5) U^2 (validation 18)
            if(g[1]<0.5)
            Ug -= dt*C_brake*std::min(1.0, std::fabs(g[2])/2.5)*s_true*Ug*std::fabs(Ug);
            const double Uw = s_true*Ug*std::min(1.0, t/1.0);
            h += dt*k_true*g[2]*Uw/1.0;
            const Eigen::Vector3d ex(std::cos(h), 0.0, std::sin(h)), el(-std::sin(h), 0.0, std::cos(h));
            x += dt*(Uw*ex(0) + cx);
            z += dt*(Uw*ex(2) + cz);
            t += dt;
            track.push_back({t, x, z, h});

            // the lateral line: 10 stations on both sides; at the head the signal of the own turn
            // (-0.5 b q_w across) and the stimulus
            const double q = 0.5*rho*std::max(Uw, 0.1)*std::max(Uw, 0.1);
            const double dp = -0.5*g[2]*q + ((stim_t>0.0 && std::fabs(t - stim_t)<0.25) ? 3.0*q : 0.0);
            std::vector<double> llp, lls;
            std::vector<int> lsd;
            for(int s=0; s<10; ++s)
            for(int sd : {1, -1})
            {
                lls.push_back((s + 0.5)/10.0);
                lsd.push_back(sd);
                llp.push_back((s + 0.5)/10.0<0.35 ? 0.5*sd*dp : 0.0);
            }
            const Eigen::Vector3d c(x, 0.0, z), U = Uw*ex + Eigen::Vector3d(cx, 0.0, cz);
            B.sense(dt, ex, el, c, U, &llp, lls, lsd);

            if(B.due(t))
            {
                sixdof_brain::decision d = B.decide(t, g[0], g[1], g[2], nullptr);
                if(d.on)
                {
                    target(d, false);
                    second = d.second;
                    pend = d;
                    log.push_back(d);
                    rec.push_back(B.last());
                    if(B.last().mode==2)
                    ++nesc;
                }
            }
        }
    }
};

static sixdof_brain brain(double Udes, double th, const std::vector<sixdof_brain::target_point> &tp, double s0 = 0.0,
                          double k = 0.0, double bmax = 1.0, double tau_u = 0.0, double tau_d = 0.0)
{
    // U_des, b_max, no vision, th_esc, f0 1 Hz, L 1 m, rho 1000, start 2 s, lateral axis z (2D), X 200 s0 k,
    // X 201 tau_u tau_d
    return sixdof_brain(Udes, bmax, 0.0, th, 1.0, 1.0, 1000.0, 2.0, 2, tp, s0, k, tau_u, tau_d);
}

static double mean_dist(const swimmer &S, double t0, std::function<std::array<double,2>(double)> T, double *mx = nullptr)
{
    double sum = 0.0, m = 0.0;
    int n = 0;
    for(const auto &p : S.track)
    if(p[0]>t0)
    {
        const auto q = T(p[0]);
        const double d = std::hypot(p[1] - q[0], p[2] - q[1]);
        sum += d; m = std::max(m, d); ++n;
    }
    if(mx) *mx = m;
    return n>0 ? sum/n : 0.0;
}

int main()
{
    // 1. station in a current of 0.2 m/s along x: the estimate and the station
    {
        swimmer S; S.cx = 0.2;
        sixdof_brain B = brain(0.5, 1.5, {{0.0, 3.5, 0.0, 3.0}});
        double t = 0.0;
        S.run(B, t, 20.0);
        const auto &r = B.last();
        check(std::fabs(r.s - 0.48)<0.03, "current 0.2: stride "+std::to_string(r.s)+" (0.48)");
        check(std::fabs(r.cx - 0.2)<0.03 && std::fabs(r.cz)<0.03, "current 0.2: current ("+std::to_string(r.cx)+", "+std::to_string(r.cz)+") (0.2, 0)");
        double mx;
        const double md = mean_dist(S, 12.0, [](double){return std::array<double,2>{3.5, 3.0};}, &mx);
        check(md<0.3, "current 0.2: station held after 12 s, mean distance "+std::to_string(md)+" m (max "+std::to_string(mx)+")");
        double hx = 0.0; int n = 0;
        for(const auto &p : S.track) if(p[0]>14.0) {hx += std::cos(p[3]); ++n;}
        check(hx/n< -0.9, "current 0.2: facing the current, mean cos(heading) "+std::to_string(hx/n));
        check(S.nesc==0, "current 0.2: no escape from the own turns ("+std::to_string(S.nesc)+")");
        const double fm = S.log.back().f;
        check(std::fabs(fm*0.48 - 0.2)<0.08, "current 0.2: swimming at the current speed, f = "+std::to_string(fm)+" Hz");
    }

    // 1b. the same current, the fish starting 0.4 m downstream of the station and facing it: it arrives before
    // it has learnt the current and must not glide there (a glide stops the estimate, the fish drifts off)
    {
        swimmer S; S.cx = 0.2; S.x = 3.9; S.h = 3.14159265358979323846;
        sixdof_brain B = brain(0.5, 1.5, {{0.0, 3.5, 0.0, 3.0}});
        double t = 0.0;
        S.run(B, t, 20.0);
        double mx;
        const double md = mean_dist(S, 8.0, [](double){return std::array<double,2>{3.5, 3.0};}, &mx);
        check(md<0.15 && mx<0.3, "current 0.2, start 0.4 m from the station: held after 8 s, mean distance "+std::to_string(md)+" m (max "+std::to_string(mx)+")");
    }

    // 2. a diagonal current of 0.2 m/s
    {
        swimmer S; S.cx = S.cz = 0.2/std::sqrt(2.0);
        sixdof_brain B = brain(0.5, 1.5, {{0.0, 3.5, 0.0, 3.0}});
        double t = 0.0;
        S.run(B, t, 20.0);
        const auto &r = B.last();
        check(std::fabs(r.cx - S.cx)<0.03 && std::fabs(r.cz - S.cz)<0.03, "diagonal current: ("+std::to_string(r.cx)+", "+std::to_string(r.cz)+")");
        const double md = mean_dist(S, 12.0, [](double){return std::array<double,2>{3.5, 3.0};});
        check(md<0.3, "diagonal current: station held, mean distance "+std::to_string(md)+" m");
    }

    // 3. still water, a target 6 m ahead: no current, the frequency U_des/(s L)
    {
        swimmer S;
        sixdof_brain B = brain(0.5, 1.5, {{0.0, 9.0, 0.0, 3.0}});
        double t = 0.0;
        S.run(B, t, 8.0);
        const auto &r = B.last();
        check(std::hypot(r.cx, r.cz)<0.03, "still water: current "+std::to_string(std::hypot(r.cx, r.cz))+" m/s");
        check(std::fabs(S.log.back().f - 0.5/0.48)<0.06, "still water: f = "+std::to_string(S.log.back().f)+" Hz (U_des/(s L) = 1.04)");
        check(std::fabs(S.z - 3.0)<0.05, "still water: straight, z = "+std::to_string(S.z));
    }

    // 3b. still water, a station: the fish glides at it (burst and coast) instead of circling
    {
        swimmer S;
        sixdof_brain B = brain(0.5, 1.5, {{0.0, 3.5, 0.0, 3.0}});
        double t = 0.0;
        int nglide = 0;
        S.run(B, t, 25.0);
        for(size_t k=0; k<S.log.size(); ++k)
        if(S.log[k].a==0.0)
        ++nglide;
        const double md = mean_dist(S, 12.0, [](double){return std::array<double,2>{3.5, 3.0};});
        check(nglide>S.log.size()/3, "still water station: gliding at "+std::to_string(nglide)+" of "+std::to_string(S.log.size())+" decisions");
        double mx;
        const double md2 = mean_dist(S, 12.0, [](double){return std::array<double,2>{3.5, 3.0};}, &mx);
        check(md<0.5 && mx<0.55, "still water station: within the burst distance 0.5 L after 12 s, mean "+std::to_string(md2)+" m, max "+std::to_string(mx));
        const auto &r = B.last();
        check(std::hypot(r.cx, r.cz)<0.05, "still water station: no current from the gliding, |c| = "+std::to_string(std::hypot(r.cx, r.cz))+" m/s");
    }

    // 3c. the early glide: the fish stops beating its tail before the station and coasts to it
    {
        swimmer S;
        sixdof_brain B = brain(0.5, 1.5, {{0.0, 3.5, 0.0, 3.0}});
        double t = 0.0;
        double d0 = -1.0;
        S.run(B, t, 30.0);
        for(size_t k=0; k<S.log.size() && d0<0.0; ++k)
        if(S.log[k].a==0.0)
        {
            const auto &p = *std::lower_bound(S.track.begin(), S.track.end(), std::array<double,4>{S.log[k].t, 0.0, 0.0, 0.0});
            d0 = std::hypot(p[1] - 3.5, p[2] - 3.0);
        }
        check(d0>0.3, "early glide: the first glide starts "+std::to_string(d0)+" m before the station (> 0.25 L)");
        double mx;
        const double md = mean_dist(S, 12.0, [](double){return std::array<double,2>{3.5, 3.0};}, &mx);
        check(md<0.2 && mx<0.25, "early glide: at the station after 12 s, mean distance "+std::to_string(md)+" m, max "+std::to_string(mx));
    }

    // 3d. the stride prior and the turn rate of a body (X 200): a slow swimmer (0.2 L per beat) with a weak
    // turn (0.5 b U/L), the target 2 m ahead and 1 m to the side; with the defaults (0.5, 1.3) the stride
    // is held at 0.25 or more and the fish sees a current that is not there
    for(int body=0; body<2; ++body)
    {
        swimmer S; S.s_true = 0.2; S.k_true = 0.5;
        sixdof_brain B = body ? brain(0.2, 1.5, {{0.0, 3.6, 0.0, 4.0}}, 0.2, 0.5) : brain(0.2, 1.5, {{0.0, 3.6, 0.0, 4.0}});
        double t = 0.0;
        S.run(B, t, 12.0);
        const auto &r = B.last();
        const double c = std::hypot(r.cx, r.cz);
        if(body)
        {
            check(std::fabs(r.s - 0.2)<0.02 && c<0.03, "X 200 0.2 0.5: stride "+std::to_string(r.s)+" (0.2), no current, |c| = "+std::to_string(c)+" m/s");
            S.run(B, t, 30.0);
            const double md = mean_dist(S, 20.0, [](double){return std::array<double,2>{3.6, 4.0};});
            check(md<0.3, "X 200 0.2 0.5: at the target, mean distance "+std::to_string(md)+" m");
        }
        else
        check(r.s>=0.25 && c>0.03, "no X 200: stride "+std::to_string(r.s)+", a false current "+std::to_string(c)+" m/s");
    }

    // 3e. a body that turns away from its bias (k = -0.9, the 3D salmon with the uniform bias, validation
    // 13): with X 200 k -0.9 it reaches a target 2 m ahead and 1.5 m to the side; with the default 1.3 not
    for(int body=0; body<2; ++body)
    {
        swimmer S; S.k_true = -0.9;
        sixdof_brain B = body ? brain(0.4, 1.5, {{0.0, 3.6, 0.0, 4.5}}, 0.0, -0.9) : brain(0.4, 1.5, {{0.0, 3.6, 0.0, 4.5}});
        double t = 0.0;
        S.run(B, t, 25.0);
        const double md = mean_dist(S, 15.0, [](double){return std::array<double,2>{3.6, 4.5};});
        if(body)
        check(md<0.4, "X 200 k -0.9: at the target, mean distance "+std::to_string(md)+" m");
        else
        check(md>1.0, "k -0.9 without X 200: away from the target, mean distance "+std::to_string(md)+" m");
    }

    // 3f. a body that is still accelerating when the brain starts (speed time 5 s instead of 0.8 s, as the
    // 3D salmon): the speed of the gait from the track, no headwind at the start
    {
        swimmer S; S.tau_true = 5.0;
        sixdof_brain B = brain(0.5, 1.5, {{0.0, 9.0, 0.0, 3.0}});
        double t = 0.0, cmax = 0.0;
        while(t<4.0)
        {
            S.run(B, t, t + 0.25);
            if(t>2.4)
            cmax = std::max(cmax, std::hypot(B.last().cx, B.last().cz));
        }
        check(cmax<0.1, "accelerating body: no current at the start (2.5-4 s), |c| at most "+std::to_string(cmax)+" m/s");
    }

    // 3g. a body that keeps its speed for seconds when its gait slows (falling with 4.5 s, as the 3D salmon,
    // validation 17; stride 0.15 L, U_des 0.2 m/s), a target 3 m ahead in still water: the brain lowers the
    // frequency on the approach. With X 201 0.8 4.5 no current; with the default 0.8 s a false tailwind
    for(int body=0; body<2; ++body)
    {
        swimmer S; S.s_true = 0.15; S.tau_down = 4.5;
        const std::vector<sixdof_brain::target_point> tp = {{0.0, 4.6, 0.0, 3.0}};
        sixdof_brain B = body ? brain(0.2, 1.5, tp, 0.15, 0.0, 1.0, 0.8, 4.5) : brain(0.2, 1.5, tp, 0.15);
        double t = 0.0, cmax = 0.0;
        while(t<20.0)
        {
            S.run(B, t, t + 0.25);
            if(t>3.0)
            cmax = std::max(cmax, std::hypot(B.last().cx, B.last().cz));
        }
        double mx;
        const double md = mean_dist(S, 15.0, [](double){return std::array<double,2>{4.6, 3.0};}, &mx);
        if(body)
        {
            check(cmax<0.02, "X 201 0.8 4.5: no current on the approach, |c| at most "+std::to_string(cmax)+" m/s");
            check(md<0.3 && mx<0.5, "X 201 0.8 4.5: at the target after 15 s, mean distance "+std::to_string(md)+" m, max "+std::to_string(mx));
        }
        else
        check(cmax>0.05, "falling with 4.5 s without X 201: a false tailwind "+std::to_string(cmax)+" m/s");
    }

    // 3h. the station in still water does not hinge on the start (validation 20): a body that coasts long (11 s,
    // the 2D swimmer) and brakes when bent (C 3.5 1/m, validation 18); with the brain starting at 1, 1.5 or 2 s
    // and the fish 1.9, 1.7 or 2.1 m before the station it holds the station, and a fish gliding or braking
    // with the station behind does not beat its tail again within twice the burst distance (0.5 L): the bent
    // body turns it back first (a beat with the station behind costs a turn-round of 1.5 L; the brain of
    // 0024: such a beat and a loop out to 1.29 m)
    {
        int nback = 0;
        double mxall = 0.0, mdall = 0.0;
        for(double ts0 : {1.0, 1.5, 2.0})
        for(double x0 : {1.6, 1.8, 1.4})
        {
            swimmer S;
            S.x = x0;
            S.tau_glide = 11.0;
            S.C_brake = 3.5;
            sixdof_brain B(0.5, 1.0, 0.0, 1.5, 1.0, 1.0, 1000.0, ts0, 2, {{0.0, 3.5, 0.0, 3.0}}, 0.0, 0.0, 0.0, 0.0);
            double t = 0.0, mx;
            S.run(B, t, 25.0);
            mdall = std::max(mdall, mean_dist(S, 12.0, [](double){return std::array<double,2>{3.5, 3.0};}, &mx));
            mxall = std::max(mxall, mx);
            for(size_t k=1; k<S.rec.size(); ++k)
            if((S.rec[k-1].mode==5 || S.rec[k-1].mode==6) && S.rec[k].mode==0 && std::fabs(S.rec[k].psi)>2.0944 && S.rec[k].dist<1.0)
            ++nback;
        }
        check(nback==0, "still water station, 3 starts x 3 distances: no beat with the station behind within 2 x 0.5 L, with it more than 120 deg off the heading ("+std::to_string(nback)+")");
        check(mdall<0.6 && mxall<0.7, "still water station, 3 starts x 3 distances: largest mean distance after 12 s "+std::to_string(mdall)+" m, largest "+std::to_string(mxall)+" m");
    }

    // 3i. a current that starts while the fish glides at its station in still water (validation 22): a coasting
    // body does not move backwards through the water, so the drift is a current at least as large; the brain
    // reads it while it still glides (before: the glide stopped the estimate, it read none)
    {
        swimmer S;
        sixdof_brain B = brain(0.5, 1.5, {{0.0, 3.5, 0.0, 3.0}});
        double t = 0.0;
        S.run(B, t, 15.0);
        S.cx = -0.1*std::cos(S.h);          // against the heading of the gliding fish
        S.cz = -0.1*std::sin(S.h);
        double cg = 0.0;
        bool glided = false;
        while(t<20.0 - 1.0e-9)
        {
            S.run(B, t, t + 0.5);
            if(B.last().mode==5 && S.log.back().a==0.0)
            {
                glided = true;
                cg = std::max(cg, -(B.last().cx*std::cos(S.h) + B.last().cz*std::sin(S.h)));
            }
        }
        check(glided && cg>0.05, "current from the drift of a glide: "+std::to_string(cg)+" m/s read against the heading while gliding (true 0.1)");
    }

    // 4. a moving target at 0.25 m/s across, after the fish reached it
    {
        swimmer S;
        std::vector<sixdof_brain::target_point> tp = {{0.0, 3.5, 0.0, 3.0}, {6.0, 3.5, 0.0, 3.0}, {20.0, 3.5, 0.0, 6.5}};
        sixdof_brain B = brain(0.35, 1.5, tp);
        double t = 0.0;
        S.run(B, t, 20.0);
        auto T = [](double t) {return std::array<double,2>{3.5, t<6.0 ? 3.0 : 3.0 + 0.25*(t - 6.0)};};
        double mx;
        const double md = mean_dist(S, 10.0, T, &mx);
        check(md<0.6, "moving target: mean distance "+std::to_string(md)+" m (max "+std::to_string(mx)+")");
    }

    // 5. escape: the own turns (th_esc 0.3, the signal of a turn at b = 1 is 0.5 q) are no escape, a
    // pulse of 3 q across the head is
    {
        swimmer S; S.cx = 0.2; S.stim_t = 15.0;
        sixdof_brain B = brain(0.5, 0.3, {{0.0, 3.5, 0.0, 3.0}});
        double t = 0.0;
        S.run(B, t, 14.5);
        check(S.nesc==0, "lateral line: no escape from the own turns at th_esc 0.3 ("+std::to_string(S.nesc)+")");
        S.run(B, t, 16.0);
        check(S.nesc==1, "lateral line: one escape on the stimulus ("+std::to_string(S.nesc)+")");
    }

    // 6. restart: the brain from its checkpoint continues exactly
    {
        swimmer S; S.cx = 0.2;
        sixdof_brain B = brain(0.5, 1.5, {{0.0, 3.5, 0.0, 3.0}});
        double t = 0.0;
        S.run(B, t, 7.3);
        std::vector<double> v;
        B.ckpt_state(v);
        swimmer S2 = S;
        sixdof_brain B2 = brain(0.5, 1.5, {{0.0, 3.5, 0.0, 3.0}});
        check(B2.ckpt_restore(v), "restart: the checkpoint of the brain is read");
        double t2 = t;
        S.run(B, t, 15.0);
        S2.run(B2, t2, 15.0);
        bool same = S.log.size()==S2.log.size();
        for(size_t k=0; same && k<S.log.size(); ++k)
        same = S.log[k].f==S2.log[k].f && S.log[k].b==S2.log[k].b && S.log[k].t==S2.log[k].t;
        check(same && S.x==S2.x && S.z==S2.z, "restart: the same decisions and track, bitwise ("+std::to_string(S.log.size())+" decisions)");
    }

    // 6b. restart while gliding at a station in still water (the coasting time, the burst)
    {
        swimmer S;
        sixdof_brain B = brain(0.5, 1.5, {{0.0, 3.5, 0.0, 3.0}});
        double t = 0.0;
        S.run(B, t, 9.3);
        std::vector<double> v;
        B.ckpt_state(v);
        swimmer S2 = S;
        sixdof_brain B2 = brain(0.5, 1.5, {{0.0, 3.5, 0.0, 3.0}});
        check(B2.ckpt_restore(v), "restart at a station: the checkpoint of the brain is read");
        double t2 = t;
        S.run(B, t, 20.0);
        S2.run(B2, t2, 20.0);
        bool same = S.log.size()==S2.log.size();
        for(size_t k=0; same && k<S.log.size(); ++k)
        same = S.log[k].f==S2.log[k].f && S.log[k].a==S2.log[k].a && S.log[k].t==S2.log[k].t;
        check(same && S.x==S2.x && S.z==S2.z, "restart at a station: the same decisions and track, bitwise");
    }

    std::cout<<npass<<" passed, "<<nfail<<" failed"<<std::endl;
    return nfail>0 ? 1 : 0;
}

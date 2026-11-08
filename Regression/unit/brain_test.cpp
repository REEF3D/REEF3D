// Standalone check of the brain of a deforming 6DOF body (X 193, src/6DOF_brain.cpp) in closed loop with a
// kinematic swimmer: no MPI, no REEF3D solver.
// Build:  g++ -O2 -std=c++17 -I../../src -I../../ThirdParty/eigen-5.0.0 -DEIGEN_MPL2_ONLY brain_test.cpp -o brain_test
// Run:    ./brain_test
// The swimmer (L = 1 m, f0 = 1 Hz, in the x-z plane): the gait follows the targets of the brain with
// critically damped dynamics, the speed of the gait f a L with 0.8 s, the speed through the water s U_g
// (s = 0.48 L per beat as the 2D swimmer of validation 08), the turn rate 1.3 b U_w/L, a uniform current.
// Checks: the estimate of the stride and of the current (0.2 m/s along x, along x + z, none); station
// holding in a current (the fish faces it and stays at the target); a target ahead in still water at
// U_des (the frequency U_des/(s L)); a moving target; the own turn on the lateral line (-0.5 b q at the
// head) is no escape, a stimulus is; a restart from the checkpoint of the brain continues it exactly.
// Architect: Hans Bihs
#include"../../src/6DOF_brain.cpp"

#include<cstdio>
#include<iostream>
#include<string>
#include<vector>
#include<functional>
#include<array>

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
    double s_true = 0.48, cx = 0.0, cz = 0.0;
    double stim_t = -1.0;      // a pressure pulse across the head at this time (a startle stimulus)
    bool second = false;
    sixdof_brain::decision pend;
    std::vector<sixdof_brain::decision> log;
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
            Ug += dt/0.8*(g[0]*g[1]*1.0 - Ug);
            const double Uw = s_true*Ug*std::min(1.0, t/1.0);
            h += dt*1.3*g[2]*Uw/1.0;
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
                    if(B.last().mode==2)
                    ++nesc;
                }
            }
        }
    }
};

static sixdof_brain brain(double Udes, double th, const std::vector<sixdof_brain::target_point> &tp)
{
    // U_des, b_max 1, no vision, th_esc, f0 1 Hz, L 1 m, rho 1000, start 2 s, lateral axis z (2D)
    return sixdof_brain(Udes, 1.0, 0.0, th, 1.0, 1.0, 1000.0, 2.0, 2, tp);
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

    std::cout<<npass<<" passed, "<<nfail<<" failed"<<std::endl;
    return nfail>0 ? 1 : 0;
}

/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2026 Hans Bihs

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
Architect: Hans Bihs
--------------------------------------------------------------------*/

#ifndef WAVE_LIB_IRREGULAR_2ND_CACHE_H_
#define WAVE_LIB_IRREGULAR_2ND_CACHE_H_

#include<vector>
#include<map>
#include<algorithm>
#include<array>
#include<Eigen/Dense>
#include<cmath>
#include"wave_lib_fastexp.h"
#include"wave_lib_eta_series.h"

/*--------------------------------------------------------------------
Second-order irregular waves (wave_lib_irregular_2nd_a / _b, B 92 32, 33,
42, 43, 52, 53), multidirectional, finite depth h.

First order, M components with amplitude a_n, frequency w_n, wave vector
k_n (direction beta_n), phase T_n = k_n.x - w_n t - e_n:

  eta1 = sum a_n cos T_n
  phi1 = sum b_n cosh(k_n(z+h))/cosh(k_n h) sin T_n,  b_n = w_n a_n/(k_n tanh(k_n h))

Second order (Sharma & Dean 1981; derived here from the free-surface
conditions at z = 0 to second order,
  phi2_tt + g phi2_z = -d/dt |grad phi1|^2 - eta1 d/dz (phi1_tt + g phi1_z)
  eta2 = -(phi2_t + |grad phi1|^2/2 + eta1 phi1_tz)/g ):
for every pair n < m a sum (+) and a difference (-) term, for every n a
self term (2 T_n; its mean-level part is left out),

  phi2 = B cosh(K(z+h))/cosh(K h) sin(T_n +- T_m),  K = |k_n +- k_m|
  eta2 = G cos(T_n +- T_m)
  B = F / (g K tanh(K h) - W^2),  W = w_n +- w_m
  F+ = -b_n b_m (k_n.k_m - |k_n||k_m| t_n t_m) W - (a_n q_m + a_m q_n)/2
  F- = -b_n b_m (k_n.k_m + |k_n||k_m| t_n t_m) W + (a_n q_m - a_m q_n)/2
  F_nn = -b_n^2 k_n^2 (1 - t_n^2) w_n - a_n q_n/2           (K = 2 k_n, W = 2 w_n)
  G = (W B - P)/g
  P+- = b_n b_m (k_n.k_m -+ |k_n||k_m| t_n t_m)/2 - (a_n b_m w_m k_m t_m + a_m b_n w_n k_n t_n)/2
  P_nn = b_n^2 k_n^2 (1 - t_n^2)/4 - a_n b_n w_n k_n t_n/2
  t_n = tanh(k_n h),  q_n = b_n k_n (g k_n - w_n^2 t_n)

Velocities u = grad phi. Checked against Stokes 2nd order (one
component) and against the exact free-surface conditions (residual
O(a^3) instead of O(a^2), validation 16).

Evaluation from per-component quantities: cos / sin of the spatial phase
per registered point (once), of the temporal phase (once per time), angle
addition per call. The components are ordered internally by direction and,
within a direction, by k. Per column (point and time) the pair phases are
multiplied into the term coefficients once (cu = bu cos(T_n +- T_m), ...);
every level of the column then forms the vertical profiles and adds two
products per term:
- pairs of the same direction: K = k_n + k_m or k_m - k_n, the profiles
  are products of P_n = exp(k_n z) and 1/P_n (no exp per pair); for one n
  the partners m are contiguous, the loops SIMD reductions
- pairs of different directions: their coefficients are summed per
  distinct |K| for the column (directional grids share |K|), each level
  one exp and two products per distinct |K|
- directions with at least 96 components: the coefficients of their pairs
  are smooth functions of (k_n, k_m); the sum and difference terms of
  velocities and eta are evaluated from truncated eigen / singular value
  expansions of the coefficient matrices (rank 7-11 at a relative error of
  1e-10), O(R L) per level instead of O(L^2)
  (lr_build, lr_kin, lr_eta); the potential of their difference terms stays
  direct, with the pair phase formed in the summation loop (the coefficient
  array is the only one read per pair)
v is skipped when all wave vectors lie on the x axis. Results agree with
the term-by-term sums to round-off (low-rank blocks: 1e-12 of the field).
Points with the same phases and the levels of a column: wave_lib_zgroup.h.
--------------------------------------------------------------------*/

struct wave_lib_irregular_2nd_terms;

struct wave_lib_irregular_2nd_cache
{
    int M=0;                              // components
    std::vector<double> cS, sS;           // cos / sin of the spatial phase, point-major
    std::vector<double> cT, sT;           // cos / sin of the temporal phase
    double t=-1.0e300;
    std::vector<double> C, S;             // phases of the current point
    int qlast=-1;                         // point of Ci, Si ...
    double tlast=-1.0e300;                // ... at this time
    int cq=-1;                            // point of C, S ...
    double ct=-1.0e300;                   // ... at this time
    std::vector<double> Ci, Si;           // C, S in the internal order of the terms (by direction and k)
    
    // column coefficients for the current Ci, Si (wave_lib_irregular_2nd_terms::column), with the
    // pair phases cp = cos(T_n +- T_m), sp = sin(T_n +- T_m): 1st order u1 C, v1 C, w1 S, b S;
    // 2nd order bu cp, bv cp, bw sp, bf sp
    std::vector<double> c1u, c1v, c1w, c1f, cu, cv, cw, cf;
    std::vector<double> gu, gv, gw, gf;   // other directions: bu cp ... summed per distinct |K|
    std::vector<double> lw, lt;           // low-rank blocks: per-column vectors, per-level sums
    int colmask=0;                        // sets present: 1 u and w, 2 v, 4 f
    
    std::vector<double> P, iP, Q, Eu, Fu; // vertical profiles of the current level (internal order)

    void points(const std::vector<double> &x, const std::vector<double> &y, int m,
                const double *ki, const double *cosb, const double *sinb)
    {
        M = m;
        const int N = int(x.size());
        cS.assign(size_t(N)*M,0.0);
        sS.assign(size_t(N)*M,0.0);

        for(int q=0; q<N; ++q)
        for(int n=0; n<M; ++n)
        {
            const double a = ki[n]*(cosb[n]*x[q] + sinb[n]*y[q]);
            cS[size_t(q)*M+n] = cos(a);
            sS[size_t(q)*M+n] = sin(a);
        }

        cT.assign(M,0.0); sT.assign(M,0.0);
        C.assign(M,0.0); S.assign(M,0.0);
        t = -1.0e300;
        qlast = cq = -1;
        colmask = 0;
    }

    void time(double wt, const double *wi, const double *ei)
    {
        if(wt==t)
        return;

        for(int n=0; n<M; ++n)
        {
            const double b = -wi[n]*wt - ei[n];
            cT[n] = cos(b);
            sT[n] = sin(b);
        }
        t = wt;
    }

    // phases and pair phases of point q at the cached time, reused while q and t stay
    // (the levels of one column)
    inline void phases_pairs(int q, const wave_lib_irregular_2nd_terms &tt);

    // phases of point q at the cached time
    void phases(int q)
    {
        const double *cs=&cS[size_t(q)*M], *ss=&sS[size_t(q)*M];
        cq = q;
        ct = t;
        colmask = 0;

        for(int n=0; n<M; ++n)
        {
            C[n] = cs[n]*cT[n] - ss[n]*sT[n];
            S[n] = ss[n]*cT[n] + cs[n]*sT[n];
        }
    }

    // phases of an arbitrary point (direct evaluation)
    void phases_at(double x, double y, double wt, const double *ki, const double *cosb, const double *sinb,
                   const double *wi, const double *ei)
    {
        C.resize(M); S.resize(M);
        colmask = 0;

        for(int n=0; n<M; ++n)
        {
            const double T = ki[n]*(cosb[n]*x + sinb[n]*y) - wi[n]*wt - ei[n];
            C[n] = cos(T);
            S[n] = sin(T);
        }
    }
};

struct wave_lib_irregular_2nd_terms
{
    // components, in the internal order: by direction, within a direction by k (lib: the
    // index of the wave library); same-direction components are then contiguous
    int M=0;
    double h=0.0, kmax=0.0, uKmax=0.0;
    std::vector<int> lib;
    std::vector<double> a, w, k, kx, ky, b, T, n1;     // T = exp(-2 k h), n1 = 1/(1+T)
    std::vector<double> u1, v1, w1;                    // 1st-order velocity factors
    std::vector<int> be;                               // end of the direction block of component n

    // terms, by class: [0,ns) same-direction sum and self terms, [ns,ns+nd) same-direction
    // difference terms, [ns+nd,N) other directions. Same direction: for each n the terms with
    // m = n .. be[n]-1 (sum, self first) and m = n+1 .. be[n]-1 (difference, k_m >= k_n), consecutive
    int ns=0, nd=0;
    std::vector<int> tn, tm, kind;                     // n, m (internal), kind (0 sum, 1 difference, 2 self)
    std::vector<double> B, G, Kx, Ky, K, nK, eK;       // nK = 1/(1+exp(-2Kh)), eK = exp(-2Kh)
    std::vector<double> bu, bv, bw, bf;                // B Kx nK, B Ky nK, B K nK, B nK
    std::vector<double> sgn;                           // +1 sum and self terms, -1 difference terms
    std::vector<int> iu;                               // other directions: index of its |K| in uK
    std::vector<double> uK, ueK;                       // distinct |K| of the other directions (directional grids share them)
    std::vector<int> us;                               // other directions: terms of uK[j] are [us[j], us[j+1])
    bool has_v=false;                                  // a wave vector off the x axis

    // low-rank evaluation of large same-direction blocks (lr_build): for one direction the pair
    // coefficients are smooth functions of (k_n, k_m) and the profiles and phases products of
    // per-component factors (Z = C + i S, Q = T / P):
    //   sum and self terms   sum_{n,m} C_nm [(P_n Z_n)(P_m Z_m) + (Q_n Z_n)(Q_m Z_m)],  C = pair / 2
    //   difference terms     sum_{n<m} D_nm [(Z_n / P_n)(P_m conj Z_m) + (Z_n P_n / T_n)(Q_m conj Z_m)]
    // with C = sum_r lam_r c_r c_r^T (symmetric eigen) and D = sum_r d_r e_r^T (SVD of D continued
    // smoothly over its diagonal), truncated at a relative Frobenius error of lr_tol: per level
    // O(R L) instead of O(L^2). Velocities and eta; the potential of the difference terms stays
    // direct (its coefficient B grows as 1 / (k_m - k_n) next to the diagonal).
    struct lr_blk
    {
        int s=0, L=0;                                  // components [s, s+L), ascending k
        double ct=1.0, st=0.0;                         // direction
        int Rv=0, Rf=0, Rg=0, Dv=0, Dg=0;              // ranks: sum velocity, potential, eta; difference velocity, eta
        std::vector<double> Uv, lv, Uf, lf, Ug, lg;    // [n*R+r], eigenvalues
        std::vector<double> DUv, DVv, DUg, DVg;        // [n*R+r]
        int ow=0, of=0, od=0, oz=0;                    // offsets of the per-column vectors in the cache (lw)
    };
    std::vector<lr_blk> lrb;
    std::vector<int> lrc;                              // component: its low-rank block or -1
    std::vector<char> lrt;                             // term: 1 low-rank sum term, 2 low-rank difference term, 0 direct
    int lr_work=0, lr_rmax=0;
    static constexpr int lr_min = 96;                  // smallest block (below it the direct sums are faster)
    static constexpr double lr_tol = 1.0e-10;          // relative Frobenius error of the truncated coefficients

    bool on=false;

    void build(int m_, const double *A, const double *wi, const double *ki, const double *cb, const double *sb, double h_, double g)
    {
        M = m_;
        h = h_;

        // directions (equal within 1e-12) in the order of their first component
        std::vector<int> dir(M,-1);
        int nd_=0;
        for(int n=0; n<M; ++n)
        if(dir[n]<0)
        {
            for(int m=n; m<M; ++m)
            if(dir[m]<0 && fabs(cb[n]-cb[m])<1.0e-12 && fabs(sb[n]-sb[m])<1.0e-12)
            dir[m] = nd_;
            ++nd_;
        }

        lib.resize(M);
        for(int n=0; n<M; ++n)
        lib[n] = n;
        std::stable_sort(lib.begin(),lib.end(),[&](int x, int y){ return dir[x]!=dir[y] ? dir[x]<dir[y] : ki[x]<ki[y]; });

        be.resize(M);
        for(int n=M-1; n>=0; --n)
        be[n] = (n==M-1 || dir[lib[n+1]]!=dir[lib[n]]) ? n+1 : be[n+1];

        a.resize(M); w.resize(M); k.resize(M);
        kx.resize(M); ky.resize(M); b.resize(M); T.resize(M); n1.resize(M);
        u1.resize(M); v1.resize(M); w1.resize(M);
        std::vector<double> t(M), q(M);
        kmax = 0.0;
        has_v = false;

        for(int n=0; n<M; ++n)
        {
            const int l = lib[n];
            a[n] = A[l];
            w[n] = wi[l];
            k[n] = ki[l];
            kx[n] = k[n]*cb[l];
            ky[n] = k[n]*sb[l];
            t[n] = tanh(k[n]*h);
            b[n] = w[n]*a[n]/(k[n]*t[n]);
            q[n] = b[n]*k[n]*(g*k[n] - w[n]*w[n]*t[n]);
            T[n] = exp(-2.0*k[n]*h);
            n1[n] = 1.0/(1.0 + T[n]);
            u1[n] = b[n]*kx[n];
            v1[n] = b[n]*ky[n];
            w1[n] = b[n]*k[n];
            kmax = fmax(kmax,k[n]);

            if(ky[n]!=0.0)
            has_v = true;
        }

        tn.clear(); tm.clear(); kind.clear(); sgn.clear();
        B.clear(); G.clear(); Kx.clear(); Ky.clear(); K.clear(); nK.clear(); eK.clear();
        bu.clear(); bv.clear(); bw.clear(); bf.clear();
        std::vector<int> cls;                          // 0 same-direction sum / self, 1 same-direction difference, 2 other

        auto add = [&](int n, int m, int kd)
        {
            const double sg = kd==1 ? -1.0 : 1.0;
            const double kxx = kx[n] + sg*kx[m];
            const double kyy = ky[n] + sg*ky[m];
            const double KK = sqrt(kxx*kxx + kyy*kyy);
            const double W = w[n] + sg*w[m];
            const double kk = kx[n]*kx[m] + ky[n]*ky[m];
            const double kt = k[n]*k[m]*t[n]*t[m];
            double F, P;

            if(kd==2)
            {
                F = -b[n]*b[n]*k[n]*k[n]*(1.0-t[n]*t[n])*w[n] - 0.5*a[n]*q[n];
                P = 0.25*b[n]*b[n]*k[n]*k[n]*(1.0-t[n]*t[n]) - 0.5*a[n]*b[n]*w[n]*k[n]*t[n];
            }
            else
            {
                const double cr = 0.5*(a[n]*b[m]*w[m]*k[m]*t[m] + a[m]*b[n]*w[n]*k[n]*t[n]);

                if(kd==0)
                {
                F = -b[n]*b[m]*(kk - kt)*W - 0.5*(a[n]*q[m] + a[m]*q[n]);
                P = 0.5*b[n]*b[m]*(kk - kt) - cr;
                }
                else
                {
                F = -b[n]*b[m]*(kk + kt)*W + 0.5*(a[n]*q[m] - a[m]*q[n]);
                P = 0.5*b[n]*b[m]*(kk + kt) - cr;
                }
            }

            const double den = g*KK*tanh(KK*h) - W*W;
            double BB = fabs(den)>1.0e-14*(g*KK + W*W + 1.0e-300) ? F/den : 0.0;

            // two components with the same wave vector: the difference term is a constant
            // (mean level), left out like the mean level of the self terms
            const bool steady = kd==1 && KK<1.0e-12 && fabs(W)<1.0e-12;
            if(steady)
            {
            BB = 0.0;
            P = 0.0;
            }

            // same direction: profiles from the per-component exponentials
            const bool same = dir[lib[n]]==dir[lib[m]];

            cls.push_back(same ? (kd==1 ? 1 : 0) : 2);
            tn.push_back(n); tm.push_back(m); kind.push_back(kd);
            sgn.push_back(kd==1 ? -1.0 : 1.0);
            B.push_back(BB);
            G.push_back((W*BB - P)/g);
            Kx.push_back(kxx); Ky.push_back(kyy); K.push_back(KK);
            eK.push_back(exp(-2.0*KK*h));
            nK.push_back(1.0/(1.0 + exp(-2.0*KK*h)));
            bu.push_back(BB*kxx*nK.back());
            bv.push_back(BB*kyy*nK.back());
            bw.push_back(BB*KK*nK.back());
            bf.push_back(BB*nK.back());

            if(kyy!=0.0)
            has_v = true;
        };

        for(int n=0; n<M; ++n)
        {
            add(n,n,2);

            for(int m=n+1; m<M; ++m)
            {
            add(n,m,0);
            add(n,m,1);
            }
        }

        // order by class (stable within a class: by n, then m)
        const int N = int(B.size());
        std::vector<int> ord;
        ord.reserve(N);
        for(int c=0; c<3; ++c)
        for(int i=0; i<N; ++i)
        if(cls[i]==c)
        ord.push_back(i);

        ns = nd = 0;
        for(int i=0; i<N; ++i)
        {
        ns += cls[i]==0;
        nd += cls[i]==1;
        }

        auto perm = [&](auto &v)
        {
            auto o = v;
            for(int i=0; i<N; ++i)
            v[i] = o[ord[i]];
        };
        perm(tn); perm(tm); perm(kind); perm(sgn);
        perm(B); perm(G); perm(Kx); perm(Ky); perm(K); perm(nK); perm(eK);
        perm(bu); perm(bv); perm(bw); perm(bf);

        // distinct |K| of the terms of other directions: with a frequency-direction grid
        // (B 136 1-3) the pairs of one frequency pair and direction difference share |K|
        uK.clear(); ueK.clear(); iu.assign(N,-1);
        uKmax = 0.0;
        std::map<long long,int> km;
        for(int i=ns+nd; i<N; ++i)
        {
            const long long key = llround(K[i]*1.0e12);
            auto it = km.find(key);
            if(it==km.end())
            {
                it = km.emplace(key,int(uK.size())).first;
                uK.push_back(K[i]);
                ueK.push_back(eK[i]);
                uKmax = fmax(uKmax,K[i]);
            }
            iu[i] = it->second;
        }
        
        // other directions ordered by their |K| (stable), so that the terms of one |K| are contiguous
        {
            const int U = int(uK.size());
            for(int i=0; i<ns+nd; ++i)
            ord[i] = i;
            int pos = ns + nd;
            us.assign(U+1,0);
            std::vector<std::vector<int>> grp(U);
            for(int i=ns+nd; i<N; ++i)
            grp[iu[i]].push_back(i);
            for(int j=0; j<U; ++j)
            {
                us[j] = pos;
                for(int i : grp[j])
                ord[pos++] = i;
            }
            us[U] = N;
            
            perm(tn); perm(tm); perm(kind); perm(sgn);
            perm(B); perm(G); perm(Kx); perm(Ky); perm(K); perm(nK); perm(eK);
            perm(bu); perm(bv); perm(bw); perm(bf); perm(iu);
        }
        
        lr_build(g);

        on = true;
    }
    
    // coefficients of a same-direction pair for unit amplitudes: kd 0 sum (K = k1 + k2), 1 difference
    // (K = k2 - k1, signed: smooth over k1 = k2, where the term itself is left out); as in build
    template<typename R>
    static void pair_unit(R k1, R w1, R k2, R w2, R h, R g, int kd, R &cw, R &cf, R &cg)
    {
        using std::tanh; using std::exp;
        const R t1 = tanh(k1*h), t2 = tanh(k2*h);
        const R b1 = w1/(k1*t1), b2 = w2/(k2*t2);
        const R q1 = b1*k1*(g*k1 - w1*w1*t1), q2 = b2*k2*(g*k2 - w2*w2*t2);
        const R kk = k1*k2, kt = kk*t1*t2;
        const R cr = R(0.5)*(b2*w2*k2*t2 + b1*w1*k1*t1);
        R F, P, W, Kp, Kd;
        
        if(kd==0)
        {
        W = w1 + w2; Kp = k1 + k2; Kd = Kp;
        F = -b1*b2*(kk - kt)*W - R(0.5)*(q2 + q1);
        P = R(0.5)*b1*b2*(kk - kt) - cr;
        }
        else
        {
        W = w1 - w2; Kp = k1 - k2; Kd = k2 - k1;
        F = -b1*b2*(kk + kt)*W + R(0.5)*(q2 - q1);
        P = R(0.5)*b1*b2*(kk + kt) - cr;
        }
        
        const R Bc = F/(g*Kp*tanh(Kp*h) - W*W);
        const R nKd = R(1)/(R(1) + exp(-R(2)*Kd*h));
        cw = Bc*Kd*nKd;
        cf = Bc*nKd;
        cg = (W*Bc - P)/g;
    }
    
    // low-rank blocks: same-direction blocks of at least lr_min components with distinct k and
    // k h < 300 (P / T stays finite); kept where every rank is below L / 4
    void lr_build(double g)
    {
        lrb.clear();
        lrc.assign(M,-1);
        lrt.assign(B.size(),0);
        lr_work = lr_rmax = 0;
        
        for(int s=0; s<M; s=be[s])
        {
            const int L = be[s] - s;
            bool ok = L>=lr_min && k[be[s]-1]*h<300.0;
            for(int n=s+1; n<be[s] && ok; ++n)
            ok = k[n]>k[n-1];
            
            if(!ok)
            continue;
            
            Eigen::MatrixXd Cv(L,L), Cf(L,L), Cg(L,L), Dv(L,L), Dg(L,L);
            
            for(int n=0; n<L; ++n)
            for(int m=0; m<L; ++m)
            {
                const int p=s+n, q=s+m;
                const double am = a[p]*a[q];
                double cw, cf, cg;
                pair_unit<double>(k[p],w[p],k[q],w[q],h,g,0,cw,cf,cg);
                Cv(n,m) = 0.5*am*cw;
                Cf(n,m) = 0.5*am*cf;
                Cg(n,m) = 0.5*am*cg;
                
                // difference terms in long double: next to the diagonal B is a ratio of two small
                // differences (cancellation), and the noise would raise the rank
                if(n!=m)
                {
                long double lw, lf, lg;
                pair_unit<long double>(k[p],w[p],k[q],w[q],h,g,1,lw,lf,lg);
                Dv(n,m) = am*double(lw);
                Dg(n,m) = am*double(lg);
                }
            }
            
            // diagonal of the difference coefficients: their limit k2 -> k1 (Richardson, k2 = k1 (1 +- d),
            // omega from the dispersion relation)
            for(int n=0; n<L; ++n)
            {
                const int p=s+n;
                const long double kp=k[p], wp=w[p], hh=h, gg=g;
                long double vw[2], vg[2];
                
                for(int l=0; l<2; ++l)
                {
                    const long double d = l==0 ? 1.0e-3L : 2.0e-3L;
                    vw[l] = vg[l] = 0.0L;
                    
                    for(int sg=-1; sg<=1; sg+=2)
                    {
                        const long double k2 = kp*(1.0L + sg*d);
                        const long double w2 = std::sqrt(gg*k2*std::tanh(k2*hh));
                        long double cw, cf, cg;
                        pair_unit<long double>(kp,wp,k2,w2,hh,gg,1,cw,cf,cg);
                        vw[l] += 0.5L*cw;
                        vg[l] += 0.5L*cg;
                    }
                }
                
                Dv(n,n) = a[p]*a[p]*double((4.0L*vw[0] - vw[1])/3.0L);
                Dg(n,n) = a[p]*a[p]*double((4.0L*vg[0] - vg[1])/3.0L);
            }
            
            lr_blk bl;
            bl.s = s;
            bl.L = L;
            bl.ct = kx[s]/k[s];
            bl.st = ky[s]/k[s];
            
            // symmetric: C = sum lam c c^T, the largest |lam| until the rest is below lr_tol |C|
            auto sym = [&](Eigen::MatrixXd &C, std::vector<double> &U, std::vector<double> &lam) -> int
            {
                C = 0.5*(C + C.transpose());
                Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(C);
                const Eigen::VectorXd &ev = es.eigenvalues();
                std::vector<int> o(L);
                for(int r=0; r<L; ++r)
                o[r] = r;
                std::sort(o.begin(),o.end(),[&](int x, int y){ return fabs(ev(x))>fabs(ev(y)); });
                
                // the error of keeping R: the sum of the dropped squares, summed from the smallest up
                const double ref = C.norm();
                std::vector<double> tail(L+1,0.0);
                for(int r=L-1; r>=0; --r)
                tail[r] = tail[r+1] + ev(o[r])*ev(o[r]);
                
                int R=0;
                while(R<L && sqrt(tail[R])>lr_tol*ref)
                ++R;
                
                U.assign(size_t(L)*R,0.0);
                lam.assign(R,0.0);
                for(int r=0; r<R; ++r)
                {
                lam[r] = ev(o[r]);
                for(int n=0; n<L; ++n)
                U[size_t(n)*R+r] = es.eigenvectors()(n,o[r]);
                }
                return R;
            };
            
            // general: D = sum d e^T on the strict upper triangle (n < m), until the error there is
            // below lr_tol |D_upper|
            auto gen = [&](const Eigen::MatrixXd &D, std::vector<double> &U, std::vector<double> &V) -> int
            {
                Eigen::BDCSVD<Eigen::MatrixXd,Eigen::ComputeThinU|Eigen::ComputeThinV> sv(D);
                Eigen::MatrixXd E = D;
                double ref=0.0;
                for(int n=0; n<L; ++n)
                for(int m=n+1; m<L; ++m)
                ref += D(n,m)*D(n,m);
                ref = sqrt(ref);
                
                int R=0;
                for(;R<L;)
                {
                    double err=0.0;
                    for(int n=0; n<L; ++n)
                    for(int m=n+1; m<L; ++m)
                    err += E(n,m)*E(n,m);
                    
                    if(sqrt(err)<=lr_tol*ref)
                    break;
                    
                    E -= sv.singularValues()(R)*sv.matrixU().col(R)*sv.matrixV().col(R).transpose();
                    ++R;
                }
                
                U.assign(size_t(L)*R,0.0);
                V.assign(size_t(L)*R,0.0);
                for(int r=0; r<R; ++r)
                for(int n=0; n<L; ++n)
                {
                U[size_t(n)*R+r] = sv.singularValues()(r)*sv.matrixU()(n,r);
                V[size_t(n)*R+r] = sv.matrixV()(n,r);
                }
                return R;
            };
            
            bl.Rv = sym(Cv,bl.Uv,bl.lv);
            bl.Rf = sym(Cf,bl.Uf,bl.lf);
            bl.Rg = sym(Cg,bl.Ug,bl.lg);
            bl.Dv = gen(Dv,bl.DUv,bl.DVv);
            bl.Dg = gen(Dg,bl.DUg,bl.DVg);
            
            if(4*std::max(std::max(std::max(bl.Rv,bl.Rf),std::max(bl.Rg,bl.Dv)),bl.Dg)>=L)
            continue;
            
            // per-column vectors (complex, re then im): (U Z) of the sums, (DU Z) and (DV conj Z)
            bl.ow = lr_work;  lr_work += 2*L*bl.Rv;
            bl.of = lr_work;  lr_work += 2*L*bl.Rf;
            bl.od = lr_work;  lr_work += 2*L*bl.Dv;
            bl.oz = lr_work;  lr_work += 2*L*bl.Dv;
            lr_rmax = std::max(lr_rmax,std::max(std::max(bl.Rv,bl.Rf),std::max(std::max(bl.Rg,bl.Dg),bl.Dv)));
            
            for(int n=s; n<s+L; ++n)
            lrc[n] = int(lrb.size());
            lrb.push_back(bl);
        }
        
        // terms of the low-rank blocks
        const int N = int(B.size());
        for(int i=0; i<ns+nd; ++i)
        if(lrc[tn[i]]>=0)
        lrt[i] = i<ns ? 1 : 2;
    }

    // C, S of the wave library order into the internal order
    void map(const double *C, const double *S, std::vector<double> &Ci, std::vector<double> &Si) const
    {
        Ci.resize(M); Si.resize(M);
        
        for(int n=0; n<M; ++n)
        {
        Ci[n] = C[lib[n]];
        Si[n] = S[lib[n]];
        }
    }
    
    // eta from the phases C = cos T, S = sin T in the internal order
    double eta_i(const double *C, const double *S) const
    {
        double e1=0.0, e2=0.0;
        
        #pragma omp simd reduction(+:e1)
        for(int n=0; n<M; ++n)
        e1 += a[n]*C[n];
        
        // same direction: sum and self terms of n, m = n .. be[n]-1; difference terms m = n+1 ..
        int i0=0;
        for(int n=0; n<M; ++n)
        {
            const double Cn = C[n], Sn = S[n];
            const double *g = G.data() + i0, *Cm = C + n, *Sm = S + n;
            const int cnt = be[n] - n;
            
            if(lrc[n]>=0)
            {
            i0 += cnt;
            continue;
            }
            
            #pragma omp simd reduction(+:e2)
            for(int j=0; j<cnt; ++j)
            e2 += g[j]*(Cn*Cm[j] - Sn*Sm[j]);
            
            i0 += cnt;
        }
        
        for(int n=0; n<M; ++n)
        {
            const double Cn = C[n], Sn = S[n];
            const double *g = G.data() + i0, *Cm = C + n + 1, *Sm = S + n + 1;
            const int cnt = be[n] - n - 1;
            
            if(lrc[n]>=0)
            {
            i0 += cnt;
            continue;
            }
            
            #pragma omp simd reduction(+:e2)
            for(int j=0; j<cnt; ++j)
            e2 += g[j]*(Cn*Cm[j] + Sn*Sm[j]);
            
            i0 += cnt;
        }
        
        // other directions
        const int N = int(B.size());
        const int *pn = tn.data(), *pm = tm.data();
        const double *sg = sgn.data(), *g = G.data();
        
        #pragma omp simd reduction(+:e2)
        for(int i=ns+nd; i<N; ++i)
        e2 += g[i]*(C[pn[i]]*C[pm[i]] - sg[i]*S[pn[i]]*S[pm[i]]);
        
        for(const lr_blk &bl : lrb)
        e2 += lr_eta(bl,C,S);
        
        return e1 + e2;
    }
    
    // eta of a low-rank block: Re sum_r lam_r (sum_n U_nr Z_n)^2 + Re sum_m sum_r DV_mr conj Z_m sum_{n<m} DU_nr Z_n
    double lr_eta(const lr_blk &bl, const double *C, const double *S) const
    {
        const int L = bl.L, s = bl.s;
        double e = 0.0;
        
        int R = bl.Rg;
        std::vector<double> tr(R,0.0), ti(R,0.0);
        for(int n=0; n<L; ++n)
        {
            const double cn = C[s+n], sn = S[s+n];
            const double *u = bl.Ug.data() + size_t(n)*R;
            
            #pragma omp simd
            for(int r=0; r<R; ++r)
            {
            tr[r] += u[r]*cn;
            ti[r] += u[r]*sn;
            }
        }
        for(int r=0; r<R; ++r)
        e += bl.lg[r]*(tr[r]*tr[r] - ti[r]*ti[r]);
        
        R = bl.Dg;
        std::vector<double> ar(R,0.0), ai(R,0.0), er(R,0.0);
        for(int m=0; m<L; ++m)
        {
            const double cm = C[s+m], sm = S[s+m];
            const double *dv = bl.DVg.data() + size_t(m)*R, *du = bl.DUg.data() + size_t(m)*R;
            
            #pragma omp simd
            for(int r=0; r<R; ++r)
            {
            er[r] += dv[r]*(cm*ar[r] + sm*ai[r]);
            ar[r] += du[r]*cm;
            ai[r] += du[r]*sm;
            }
        }
        for(int r=0; r<R; ++r)
        e += er[r];
        
        return e;
    }
    
    // eta from the phases C = cos T, S = sin T in the order of the wave library
    double eta(const double *C, const double *S) const
    {
        std::vector<double> Ci, Si;
        map(C,S,Ci,Si);
        return eta_i(Ci.data(),Si.data());
    }
    
    // eta from the phases of the cache
    double eta(wave_lib_irregular_2nd_cache &c) const
    {
        map(c.C.data(),c.S.data(),c.Ci,c.Si);
        return eta_i(c.Ci.data(),c.Si.data());
    }

    // eta at a point for the times tv, with the phases phi_n = k_n.x - e_n (wave library order):
    // the components and the terms as harmonics, G cos(phi_n +- phi_m - (w_n +- w_m) t); the terms
    // of one frequency (within 1e-14 of it) as one harmonic (wave_lib_eta_series.h). The terms of
    // the low-rank blocks per sample from their expansions (lr_series).
    void series(const double *phi, const std::vector<double> &tv, std::vector<double> &ev) const
    {
        const int N = int(B.size());
        std::vector<std::array<double,3>> hm;
        hm.reserve(M+N);
        
        for(int n=0; n<M; ++n)
        {
            const double ph = phi[lib[n]];
            hm.push_back({w[n],a[n]*cos(ph),a[n]*sin(ph)});
        }
        
        for(int i=0; i<N; ++i)
        {
            if(i<ns+nd && lrc[tn[i]]>=0)
            continue;
            
            double ph = phi[lib[tn[i]]] + sgn[i]*phi[lib[tm[i]]];
            double W = w[tn[i]] + sgn[i]*w[tm[i]];
            
            if(W<0.0)
            {
            W = -W;
            ph = -ph;
            }
            hm.push_back({W,G[i]*cos(ph),G[i]*sin(ph)});
        }
        
        std::stable_sort(hm.begin(),hm.end(),[](const std::array<double,3> &x, const std::array<double,3> &y){ return x[0]<y[0]; });
        
        std::vector<double> Wb, ar, ai;
        for(size_t q=0; q<hm.size(); ++q)
        {
            if(Wb.empty() || hm[q][0] - Wb.back() > 1.0e-14*fmax(1.0,Wb.back()))
            {
            Wb.push_back(hm[q][0]);
            ar.push_back(0.0);
            ai.push_back(0.0);
            }
            ar.back() += hm[q][1];
            ai.back() += hm[q][2];
        }
        
        wave_eta_harmonics(Wb,ar,ai,tv,ev);
        
        if(!lrb.empty())
        lr_series(phi,tv,ev);
    }
    
    // adds the terms of the low-rank blocks to the series ev: the phases T_n = phi_n - w_n t rotated
    // from sample to sample (exact every 32 samples and when the step changes, as
    // wave_eta_harmonics), lr_eta per sample: O(R L) instead of O(L^2) harmonics
    void lr_series(const double *phi, const std::vector<double> &tv, std::vector<double> &ev) const
    {
        const size_t NT = tv.size();
        
        if(NT==0)
        return;
        
        std::vector<double> ph(M), C(M), S(M), cd(M), sd(M);
        const double d0 = NT>1 ? tv[1] - tv[0] : 0.0;
        
        for(int n=0; n<M; ++n)
        {
        ph[n] = phi[lib[n]];
        cd[n] = cos(w[n]*d0);
        sd[n] = sin(w[n]*d0);
        }
        
        int since = 32;
        
        for(size_t q=0; q<NT; ++q)
        {
            const bool exact = since>=32 || q==0 || fabs((tv[q] - tv[q-1]) - d0) > 1.0e-9*fabs(d0);
            
            if(exact)
            {
                for(int n=0; n<M; ++n)
                {
                C[n] = cos(ph[n] - w[n]*tv[q]);
                S[n] = sin(ph[n] - w[n]*tv[q]);
                }
                since = 1;
            }
            else
            {
                #pragma omp simd
                for(int n=0; n<M; ++n)
                {
                    const double cn = C[n]*cd[n] + S[n]*sd[n];
                    S[n] = S[n]*cd[n] - C[n]*sd[n];
                    C[n] = cn;
                }
                ++since;
            }
            
            double e = 0.0;
            
            for(const lr_blk &bl : lrb)
            e += lr_eta(bl,C.data(),S.data());
            
            ev[q] += e;
        }
    }
    
    // phases of the cache into the internal order (once per point and time); the column
    // coefficients follow from them (column)
    void pair_phases(wave_lib_irregular_2nd_cache &c) const
    {
        map(c.C.data(),c.S.data(),c.Ci,c.Si);
        c.colmask = 0;
    }
    
    // column coefficients (sets in need: 1 u and w, 2 v, 4 f) from the phases Ci, Si of the cache,
    // with the pair phases cp = C_n C_m -+ S_n S_m, sp = S_n C_m +- C_n S_m (+- the sign of the term).
    // Terms of other directions: summed per distinct |K| (their profiles are the same at every level)
    void column(wave_lib_irregular_2nd_cache &c, int need) const
    {
        need &= ~c.colmask;
        
        if(need==0)
        return;
        
        const int N = int(B.size());
        const int U = int(uK.size());
        
        if(need&1)
        {
        c.c1u.resize(M); c.c1w.resize(M); c.cu.resize(N); c.cw.resize(N); c.gu.resize(U); c.gw.resize(U);
        }
        if(need&2)
        {
        c.c1v.resize(M); c.cv.resize(N); c.gv.resize(U);
        }
        if(need&4)
        {
        c.c1f.resize(M); c.cf.resize(N); c.gf.resize(U);
        }
        
        switch(need)
        {
            case 1: column_t<true,false,false>(c); break;
            case 2: column_t<false,true,false>(c); break;
            case 3: column_t<true,true,false>(c); break;
            case 4: column_t<false,false,true>(c); break;
            case 5: column_t<true,false,true>(c); break;
            case 6: column_t<false,true,true>(c); break;
            default: column_t<true,true,true>(c); break;
        }
        
        c.colmask |= need;
    }
    
    template<bool DW, bool DV, bool DF>
    void column_t(wave_lib_irregular_2nd_cache &c) const
    {
        const double *C = c.Ci.data(), *S = c.Si.data();
        double *c1u = c.c1u.data(), *c1v = c.c1v.data(), *c1w = c.c1w.data(), *c1f = c.c1f.data();
        double *cu = c.cu.data(), *cv = c.cv.data(), *cw = c.cw.data(), *cf = c.cf.data();
        const double *pbu = bu.data(), *pbv = bv.data(), *pbw = bw.data(), *pbf = bf.data();
        
        #pragma omp simd
        for(int n=0; n<M; ++n)
        {
            if constexpr(DW)
            {
            c1u[n] = u1[n]*C[n];
            c1w[n] = w1[n]*S[n];
            }
            if constexpr(DV)
            c1v[n] = v1[n]*C[n];
            if constexpr(DF)
            c1f[n] = b[n]*S[n];
        }
        
        // same direction, sum and self terms of n (m = n .. be[n]-1), then the difference terms
        // (m = n+1 .. be[n]-1)
        int i0=0;
        for(int pass=0; pass<2; ++pass)
        for(int n=0; n<M; ++n)
        {
            const double Cn = C[n], Sn = S[n];
            const double sg = pass==0 ? 1.0 : -1.0;
            const int m0 = n + pass;
            const int cnt = be[n] - m0;
            
            // low-rank block: the potential of its difference terms is summed directly in kin_t
            if(lrc[n]>=0)
            {
                i0 += cnt;
                continue;
            }
            
            #pragma omp simd
            for(int j=0; j<cnt; ++j)
            {
                const int i = i0 + j;
                const double cp = Cn*C[m0+j] - sg*Sn*S[m0+j];
                const double sp = Sn*C[m0+j] + sg*Cn*S[m0+j];
                
                if constexpr(DW)
                {
                cu[i] = pbu[i]*cp;
                cw[i] = pbw[i]*sp;
                }
                if constexpr(DV)
                cv[i] = pbv[i]*cp;
                if constexpr(DF)
                cf[i] = pbf[i]*sp;
            }
            
            i0 += cnt;
        }
        
        // other directions, by distinct |K| (terms [us[j], us[j+1]))
        const int *pn = tn.data(), *pm = tm.data();
        const double *sg = sgn.data();
        const int U = int(uK.size());
        
        for(int j=0; j<U; ++j)
        {
            double au=0.0, av=0.0, aw=0.0, af=0.0;
            
            #pragma omp simd reduction(+:au,av,aw,af)
            for(int i=us[j]; i<us[j+1]; ++i)
            {
                const int n=pn[i], m=pm[i];
                const double cp = C[n]*C[m] - sg[i]*S[n]*S[m];
                const double sp = S[n]*C[m] + sg[i]*C[n]*S[m];
                
                if constexpr(DW)
                {
                au += pbu[i]*cp;
                aw += pbw[i]*sp;
                }
                if constexpr(DV)
                av += pbv[i]*cp;
                if constexpr(DF)
                af += pbf[i]*sp;
            }
            
            if constexpr(DW)
            {
            c.gu[j] = au;
            c.gw[j] = aw;
            }
            if constexpr(DV)
            c.gv[j] = av;
            if constexpr(DF)
            c.gf[j] = af;
        }
        
        // low-rank blocks: (U Z) of the sums, (DU Z) and (DV conj Z) of the difference terms
        if(!lrb.empty())
        {
            c.lw.resize(lr_work);
            c.lt.resize(12*size_t(lr_rmax));
        }
        
        for(const lr_blk &bl : lrb)
        {
            const int L = bl.L, s = bl.s;
            auto mult = [&](const std::vector<double> &U, int R, double *w, double sgn_i)
            {
                double *wr = w, *wi = w + size_t(L)*R;
                for(int n=0; n<L; ++n)
                {
                    const double cn = C[s+n], sn = sgn_i*S[s+n];
                    const double *u = U.data() + size_t(n)*R;
                    
                    #pragma omp simd
                    for(int r=0; r<R; ++r)
                    {
                    wr[size_t(n)*R+r] = u[r]*cn;
                    wi[size_t(n)*R+r] = u[r]*sn;
                    }
                }
            };
            
            if constexpr(DW || DV)
            {
            mult(bl.Uv,bl.Rv,c.lw.data()+bl.ow,1.0);
            mult(bl.DUv,bl.Dv,c.lw.data()+bl.od,1.0);
            mult(bl.DVv,bl.Dv,c.lw.data()+bl.oz,-1.0);
            }
            if constexpr(DF)
            mult(bl.Uf,bl.Rf,c.lw.data()+bl.of,1.0);
        }
    }
    
    // velocities and potential at z (relative to the still water level), with the phases of the
    // cache (phases_pairs / pair_phases).
    // mode: 1 u and w, 2 v, 4 the potential f (sums of these); the others stay 0
    void kin(wave_lib_irregular_2nd_cache &c, double z, int mode, double &u, double &v, double &ww, double &f) const
    {
        const int need = (mode&1) | (has_v ? (mode&2) : 0) | (mode&4);
        
        u=v=ww=f=0.0;
        
        if(need==0)
        return;
        
        column(c,need);
        
        // products of the per-component exponentials only while 1/P stays finite
        const bool prod = kmax*(h + fabs(z)) < 600.0;
        
        c.P.resize(M); c.iP.resize(M); c.Q.resize(M);
        double *P = c.P.data(), *iP = c.iP.data(), *Q = c.Q.data();
        
        if(prod)
        {
            // |k z| < 600: within the range of wave_exp
            #pragma omp simd
            for(int n=0; n<M; ++n)
            {
            P[n] = wave_exp(k[n]*z);
            iP[n] = 1.0/P[n];
            Q[n] = iP[n]*T[n];
            }
            
            // other directions: one exp per distinct |K|
            const int U = int(uK.size());
            c.Eu.resize(U); c.Fu.resize(U);
            double *Eu = c.Eu.data(), *Fu = c.Fu.data();
            
            if(uKmax*(h + fabs(z)) < 600.0)
            {
                #pragma omp simd
                for(int j=0; j<U; ++j)
                {
                Eu[j] = wave_exp(uK[j]*z);
                Fu[j] = ueK[j]/Eu[j];
                }
            }
            else
            for(int j=0; j<U; ++j)
            {
            Eu[j] = exp(uK[j]*z);
            Fu[j] = exp(-uK[j]*(z+2.0*h));
            }
        }
        else
        for(int n=0; n<M; ++n)
        {
            P[n] = exp(k[n]*z);
            iP[n] = 1.0/P[n];
            Q[n] = exp(-k[n]*(z+2.0*h));
        }
        
        switch(need)
        {
            case 1: kin_t<true,false,false>(c,z,prod,u,v,ww,f); break;
            case 2: kin_t<false,true,false>(c,z,prod,u,v,ww,f); break;
            case 3: kin_t<true,true,false>(c,z,prod,u,v,ww,f); break;
            case 4: kin_t<false,false,true>(c,z,prod,u,v,ww,f); break;
            case 5: kin_t<true,false,true>(c,z,prod,u,v,ww,f); break;
            case 6: kin_t<false,true,true>(c,z,prod,u,v,ww,f); break;
            default: kin_t<true,true,true>(c,z,prod,u,v,ww,f); break;
        }
    }
    
    // the sums of kin over the terms; the loops are SIMD reductions (#pragma omp simd, with
    // -fopenmp or -fopenmp-simd; the order of the additions then depends on the vector width)
    template<bool UW, bool DV, bool DF>
    void kin_t(wave_lib_irregular_2nd_cache &c, double z, bool prod, double &u, double &v, double &ww, double &f) const
    {
        const double *P = c.P.data(), *iP = c.iP.data(), *Q = c.Q.data();
        const double *Ci = c.Ci.data(), *Si = c.Si.data();
        const double *c1u = c.c1u.data(), *c1v = c.c1v.data(), *c1w = c.c1w.data(), *c1f = c.c1f.data();
        const double *cu = c.cu.data(), *cv = c.cv.data(), *cw = c.cw.data(), *cf = c.cf.data();
        const int N = int(B.size());
        
        double su=0.0, sv=0.0, sw=0.0, sf=0.0;
        
        // first order: ch = (P + Q) n1, sh = (P - Q) n1
        const double *pn1 = n1.data();
        
        #pragma omp simd reduction(+:su,sv,sw,sf)
        for(int n=0; n<M; ++n)
        {
            const double E = P[n]*pn1[n], F = Q[n]*pn1[n];
            
            if constexpr(UW)
            {
            su += c1u[n]*(E + F);
            sw += c1w[n]*(E - F);
            }
            if constexpr(DV)
            sv += c1v[n]*(E + F);
            if constexpr(DF)
            sf += c1f[n]*(E + F);
        }
        
        if(prod)
        {
            // same direction, sum and self terms of n: K = k_n + k_m, m = n .. be[n]-1,
            // profiles P_n P_m, Q_n Q_m
            int i0 = 0;
            for(int n=0; n<M; ++n)
            {
                const double Pn = P[n], Qn = Q[n];
                const double *Pm = P + n, *Qm = Q + n;
                const int cnt = be[n] - n;
                
                if(lrc[n]>=0)
                {
                i0 += cnt;
                continue;
                }
                
                #pragma omp simd reduction(+:su,sv,sw,sf)
                for(int j=0; j<cnt; ++j)
                {
                    const double E = Pn*Pm[j], F = Qn*Qm[j];
                    const int i = i0 + j;
                    
                    if constexpr(UW)
                    {
                    su += cu[i]*(E + F);
                    sw += cw[i]*(E - F);
                    }
                    if constexpr(DV)
                    sv += cv[i]*(E + F);
                    if constexpr(DF)
                    sf += cf[i]*(E + F);
                }
                
                i0 += cnt;
            }
            
            // same direction, difference terms of n: K = k_m - k_n, m = n+1 .. be[n]-1,
            // profiles P_m / P_n, P_n / P_m exp(-2Kh)
            const double *e = eK.data();
            for(int n=0; n<M; ++n)
            {
                const double Pn = P[n], iPn = iP[n];
                const double *Pm = P + n + 1, *iPm = iP + n + 1;
                const int cnt = be[n] - n - 1;
                
                // low-rank block: only the potential directly, with the pair phase in the loop,
                // B nK sin(T_n - T_m) (P_m / P_n + (P_n / T_n) Q_m), exp(-2Kh) = T_m / T_n: the
                // coefficients B nK are the only array per pair (no column array written and read)
                if(lrc[n]>=0)
                {
                    if constexpr(DF)
                    {
                        const double *Cm = Ci + n + 1, *Sm = Si + n + 1, *Qm = Q + n + 1, *pbf = bf.data() + i0;
                        const double Cn = Ci[n], Sn = Si[n], Bn = Pn/T[n];
                        
                        #pragma omp simd reduction(+:sf)
                        for(int j=0; j<cnt; ++j)
                        sf += pbf[j]*(Sn*Cm[j] - Cn*Sm[j])*(Pm[j]*iPn + Qm[j]*Bn);
                    }
                    
                    i0 += cnt;
                    continue;
                }
                
                #pragma omp simd reduction(+:su,sv,sw,sf)
                for(int j=0; j<cnt; ++j)
                {
                    const int i = i0 + j;
                    const double E = Pm[j]*iPn, F = iPm[j]*Pn*e[i];
                    
                    if constexpr(UW)
                    {
                    su += cu[i]*(E + F);
                    sw += cw[i]*(E - F);
                    }
                    if constexpr(DV)
                    sv += cv[i]*(E + F);
                    if constexpr(DF)
                    sf += cf[i]*(E + F);
                }
                
                i0 += cnt;
            }
            
            // other directions: per distinct |K| (coefficients summed per column)
            const double *Eu = c.Eu.data(), *Fu = c.Fu.data();
            const double *gu = c.gu.data(), *gv = c.gv.data(), *gw = c.gw.data(), *gf = c.gf.data();
            const int U = int(uK.size());
            
            #pragma omp simd reduction(+:su,sv,sw,sf)
            for(int j=0; j<U; ++j)
            {
                const double E = Eu[j], F = Fu[j];
                
                if constexpr(UW)
                {
                su += gu[j]*(E + F);
                sw += gw[j]*(E - F);
                }
                if constexpr(DV)
                sv += gv[j]*(E + F);
                if constexpr(DF)
                sf += gf[j]*(E + F);
            }
        }
        else
        {
            // deep water: exp per term (same direction) and per distinct |K| (other directions)
            for(int i=0; i<ns+nd; ++i)
            {
                if(lrt[i]==1 || (lrt[i]==2 && !DF))
                continue;
                
                const double E = exp(K[i]*z), F = exp(-K[i]*(z+2.0*h));
                
                if(lrt[i]==2)
                {
                if constexpr(DF)
                sf += bf[i]*(Si[tn[i]]*Ci[tm[i]] - Ci[tn[i]]*Si[tm[i]])*(E + F);
                continue;
                }
                
                if constexpr(UW)
                {
                su += cu[i]*(E + F);
                sw += cw[i]*(E - F);
                }
                if constexpr(DV)
                sv += cv[i]*(E + F);
                if constexpr(DF)
                sf += cf[i]*(E + F);
            }
            
            const int U = int(uK.size());
            for(int j=0; j<U; ++j)
            {
                const double E = exp(uK[j]*z), F = exp(-uK[j]*(z+2.0*h));
                
                if constexpr(UW)
                {
                su += c.gu[j]*(E + F);
                sw += c.gw[j]*(E - F);
                }
                if constexpr(DV)
                sv += c.gv[j]*(E + F);
                if constexpr(DF)
                sf += c.gf[j]*(E + F);
            }
        }
        
        for(const lr_blk &bl : lrb)
        lr_kin<UW,DV,DF>(c,bl,su,sv,sw,sf);
        
        u = su;
        v = sv;
        ww = sw;
        f = sf;
    }
    
    // low-rank block at the current level (profiles P, Q, 1/P of the cache)
    template<bool UW, bool DV, bool DF>
    void lr_kin(wave_lib_irregular_2nd_cache &c, const lr_blk &bl, double &su, double &sv, double &sw, double &sf) const
    {
        const int L = bl.L, s = bl.s;
        const double *P = c.P.data() + s, *Q = c.Q.data() + s, *iP = c.iP.data() + s, *Tn = T.data() + s;
        double *t = c.lt.data();
        
        // sums: s1 = sum_n W_nr P_n, s2 = sum_n W_nr Q_n; sum_r lam_r (s1^2 + s2^2)
        auto sums = [&](const double *w, int R, const double *lam, double &X, double &Y, double &Z)
        {
            const double *wr = w, *wi = w + size_t(L)*R;
            double *ar = t, *ai = t + R, *br = t + 2*R, *bi = t + 3*R;
            for(int r=0; r<4*R; ++r)
            t[r] = 0.0;
            
            for(int n=0; n<L; ++n)
            {
                const double p = P[n], q = Q[n];
                const double *xr = wr + size_t(n)*R, *xi = wi + size_t(n)*R;
                
                #pragma omp simd
                for(int r=0; r<R; ++r)
                {
                ar[r] += xr[r]*p;
                ai[r] += xi[r]*p;
                br[r] += xr[r]*q;
                bi[r] += xi[r]*q;
                }
            }
            
            // X = Re(SP + SQ), Y = Im(SP - SQ), Z = Im(SP + SQ)
            X = Y = Z = 0.0;
            for(int r=0; r<R; ++r)
            {
            X += lam[r]*(ar[r]*ar[r] - ai[r]*ai[r] + br[r]*br[r] - bi[r]*bi[r]);
            Y += lam[r]*2.0*(ar[r]*ai[r] - br[r]*bi[r]);
            Z += lam[r]*2.0*(ar[r]*ai[r] + br[r]*bi[r]);
            }
        };
        
        if constexpr(UW || DV)
        {
            double X, Y, Z;
            sums(c.lw.data()+bl.ow,bl.Rv,bl.lv.data(),X,Y,Z);
            
            // difference terms: prefix sums A_r(m) = sum_{n<m} DU_nr Z_n / P_n, B_r(m) = sum_{n<m} DU_nr Z_n P_n / T_n
            const int R = bl.Dv;
            const double *dr = c.lw.data() + bl.od, *di = dr + size_t(L)*R;
            const double *zr = c.lw.data() + bl.oz, *zi = zr + size_t(L)*R;
            double *aAr = t, *aAi = t + R, *aBr = t + 2*R, *aBi = t + 3*R;
            double *ser = t + 4*R, *sei = t + 5*R, *sfr = t + 6*R, *sfi = t + 7*R;
            for(int r=0; r<8*R; ++r)
            t[r] = 0.0;
            
            for(int m=0; m<L; ++m)
            {
                const double p = P[m], q = Q[m], ip = iP[m], pt = P[m]/Tn[m];
                const double *vr = zr + size_t(m)*R, *vi = zi + size_t(m)*R;
                const double *wr = dr + size_t(m)*R, *wi = di + size_t(m)*R;
                
                #pragma omp simd
                for(int r=0; r<R; ++r)
                {
                ser[r] += p*(vr[r]*aAr[r] - vi[r]*aAi[r]);
                sei[r] += p*(vr[r]*aAi[r] + vi[r]*aAr[r]);
                sfr[r] += q*(vr[r]*aBr[r] - vi[r]*aBi[r]);
                sfi[r] += q*(vr[r]*aBi[r] + vi[r]*aBr[r]);
                aAr[r] += wr[r]*ip;
                aAi[r] += wi[r]*ip;
                aBr[r] += wr[r]*pt;
                aBi[r] += wi[r]*pt;
                }
            }
            
            double Xd=0.0, Yd=0.0;
            for(int r=0; r<R; ++r)
            {
            Xd += ser[r] + sfr[r];
            Yd += sei[r] - sfi[r];
            }
            
            // u = cos(beta) (X_sum - X_dif), v = sin(beta) (...), w = Y_sum + Y_dif
            if constexpr(UW)
            {
            su += bl.ct*(X - Xd);
            sw += Y + Yd;
            }
            if constexpr(DV)
            sv += bl.st*(X - Xd);
        }
        
        if constexpr(DF)
        {
            double X, Y, Z;
            sums(c.lw.data()+bl.of,bl.Rf,bl.lf.data(),X,Y,Z);
            sf += Z;
        }
    }
};

inline void wave_lib_irregular_2nd_cache::phases_pairs(int q, const wave_lib_irregular_2nd_terms &tt)
{
    if(q==qlast && t==tlast && q==cq && t==ct)
    return;

    if(q!=cq || t!=ct)
    phases(q);

    tt.pair_phases(*this);
    qlast = q;
    tlast = t;
}

#endif

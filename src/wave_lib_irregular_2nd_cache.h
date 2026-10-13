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
v is skipped when all wave vectors lie on the x axis. Results agree with
the term-by-term sums to round-off.
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

        on = true;
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
        
        return e1 + e2;
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
    // of one frequency (within 1e-14 of it) as one harmonic (wave_lib_eta_series.h)
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
    void kin_t(const wave_lib_irregular_2nd_cache &c, double z, bool prod, double &u, double &v, double &ww, double &f) const
    {
        const double *P = c.P.data(), *iP = c.iP.data(), *Q = c.Q.data();
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
                const double E = exp(K[i]*z), F = exp(-K[i]*(z+2.0*h));
                
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
        
        u = su;
        v = sv;
        ww = sw;
        f = sf;
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

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
#include<cmath>

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
addition per call; vertical profiles from P_n = exp(k_n z): for pairs of
the same direction K = k_n + k_m or |k_n - k_m| and the profiles are
products of P_n, 1/P_n (no exp per pair), other pairs one exp each.
--------------------------------------------------------------------*/

struct wave_lib_irregular_2nd_terms
{
    // components
    int M=0;
    double h=0.0, kmax=0.0;
    std::vector<double> a, w, k, kx, ky, b, T, n1;     // T = exp(-2 k h), n1 = 1/(1+T)
    std::vector<double> u1, v1, w1;                    // 1st-order velocity factors

    // terms: n, m, kind (0 sum, 1 difference, 2 self), coefficients, wave vector, vertical data
    std::vector<int> tn, tm, kind, coll, hi, lo;
    std::vector<double> B, G, Kx, Ky, K, nK, eK;       // nK = 1/(1+exp(-2Kh)), eK = exp(-2Kh)
    std::vector<double> bu, bv, bw, bf;                // B Kx nK, B Ky nK, B K nK, B nK
    std::vector<int> isum, idif, ioth;                 // same-direction sum (and self) / difference terms, others

    bool on=false;

    void build(int m_, const double *A, const double *wi, const double *ki, const double *cb, const double *sb, double h_, double g)
    {
        M = m_;
        h = h_;
        a.assign(A,A+M); w.assign(wi,wi+M); k.assign(ki,ki+M);
        kx.resize(M); ky.resize(M); b.resize(M); T.resize(M); n1.resize(M);
        u1.resize(M); v1.resize(M); w1.resize(M);
        std::vector<double> t(M), q(M);
        kmax = 0.0;

        for(int n=0; n<M; ++n)
        {
            kx[n] = k[n]*cb[n];
            ky[n] = k[n]*sb[n];
            t[n] = tanh(k[n]*h);
            b[n] = w[n]*a[n]/(k[n]*t[n]);
            q[n] = b[n]*k[n]*(g*k[n] - w[n]*w[n]*t[n]);
            T[n] = exp(-2.0*k[n]*h);
            n1[n] = 1.0/(1.0 + T[n]);
            u1[n] = b[n]*kx[n];
            v1[n] = b[n]*ky[n];
            w1[n] = b[n]*k[n];
            kmax = fmax(kmax,k[n]);
        }

        tn.clear(); tm.clear(); kind.clear(); coll.clear(); hi.clear(); lo.clear();
        B.clear(); G.clear(); Kx.clear(); Ky.clear(); K.clear(); nK.clear(); eK.clear();
        bu.clear(); bv.clear(); bw.clear(); bf.clear(); isum.clear(); idif.clear(); ioth.clear();

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
            const bool same = fabs(cb[n]-cb[m])<1.0e-12 && fabs(sb[n]-sb[m])<1.0e-12;

            const int id = int(B.size());
            if(same)
            (kd==1 ? idif : isum).push_back(id);
            else
            ioth.push_back(id);
            
            tn.push_back(n); tm.push_back(m); kind.push_back(kd);
            coll.push_back(same ? 1 : 0);
            hi.push_back(k[n]>=k[m] ? n : m);
            lo.push_back(k[n]>=k[m] ? m : n);
            B.push_back(BB);
            G.push_back((W*BB - P)/g);
            Kx.push_back(kxx); Ky.push_back(kyy); K.push_back(KK);
            eK.push_back(exp(-2.0*KK*h));
            nK.push_back(1.0/(1.0 + exp(-2.0*KK*h)));
            bu.push_back(BB*kxx*nK.back());
            bv.push_back(BB*kyy*nK.back());
            bw.push_back(BB*KK*nK.back());
            bf.push_back(BB*nK.back());
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

        on = true;
    }

    // eta from the phases C = cos T, S = sin T
    double eta(const double *C, const double *S) const
    {
        double e=0.0;

        for(int n=0; n<M; ++n)
        e += a[n]*C[n];

        const int N = int(B.size());
        for(int i=0; i<N; ++i)
        {
            const int n=tn[i], m=tm[i];
            const double cc = C[n]*C[m], ss = S[n]*S[m];
            e += G[i]*(kind[i]==1 ? cc + ss : cc - ss);
        }

        return e;
    }

    // cos / sin of T_n +- T_m of all terms (once per point and time)
    void pair_phases(const double *C, const double *S, std::vector<double> &cp, std::vector<double> &sp) const
    {
        const int N = int(B.size());
        cp.resize(N); sp.resize(N);
        
        for(int i=0; i<N; ++i)
        {
            const int n=tn[i], m=tm[i];
            
            if(kind[i]==1)
            {
            cp[i] = C[n]*C[m] + S[n]*S[m];
            sp[i] = S[n]*C[m] - C[n]*S[m];
            }
            else
            {
            cp[i] = C[n]*C[m] - S[n]*S[m];
            sp[i] = S[n]*C[m] + C[n]*S[m];
            }
        }
    }
    
    // velocities and potential at z (relative to the still water level); cp, sp from pair_phases.
    // mode: 1 u and w, 2 v, 4 the potential f (sums of these); the others stay 0
    void kin(const double *C, const double *S, const double *cp, const double *sp, double z, int mode,
             double &u, double &v, double &ww, double &f,
             std::vector<double> &P, std::vector<double> &iP, std::vector<double> &Q) const
    {
        P.resize(M); iP.resize(M); Q.resize(M);
        
        // products of the per-component exponentials only while 1/P stays finite
        const bool prod = kmax*(h + fabs(z)) < 600.0;
        const bool uw = mode&1, dov = mode&2, dof = mode&4;
        
        u=v=ww=f=0.0;
        
        for(int n=0; n<M; ++n)
        {
            P[n] = exp(k[n]*z);
            iP[n] = 1.0/P[n];
            Q[n] = prod ? iP[n]*T[n] : exp(-k[n]*(z+2.0*h));
        }
        
        for(int n=0; n<M; ++n)
        {
            const double ch = (P[n] + Q[n])*n1[n];
            const double sh = (P[n] - Q[n])*n1[n];
            
            if(uw)
            {
            u += u1[n]*ch*C[n];
            ww += w1[n]*sh*S[n];
            }
            if(dov)
            v += v1[n]*ch*C[n];
            if(dof)
            f += b[n]*ch*S[n];
        }
        
        auto add = [&](int i, double E, double F)
        {
            const double cs = (E + F)*cp[i];
            
            if(uw)
            {
            u += bu[i]*cs;
            ww += bw[i]*(E - F)*sp[i];
            }
            if(dov)
            v += bv[i]*cs;
            if(dof)
            f += bf[i]*(E + F)*sp[i];
        };
        
        if(prod)
        {
            // same direction: products of the per-component exponentials
            for(int i : isum)
            add(i,P[tn[i]]*P[tm[i]],Q[tn[i]]*Q[tm[i]]);
            
            for(int i : idif)
            add(i,P[hi[i]]*iP[lo[i]],iP[hi[i]]*P[lo[i]]*eK[i]);
            
            // other directions: one exp per term
            for(int i : ioth)
            {
                const double E = exp(K[i]*z);
                add(i,E,K[i]*h<600.0 ? eK[i]/E : exp(-K[i]*(z+2.0*h)));
            }
        }
        else
        {
            const int N = int(B.size());
            for(int i=0; i<N; ++i)
            add(i,exp(K[i]*z),exp(-K[i]*(z+2.0*h)));
        }
    }
};

struct wave_lib_irregular_2nd_cache
{
    int M=0;                              // components
    std::vector<double> cS, sS;           // cos / sin of the spatial phase, point-major
    std::vector<double> cT, sT;           // cos / sin of the temporal phase
    double t=-1.0e300;
    std::vector<double> C, S;             // phases of the current point
    std::vector<double> cp, sp;           // pair phases of the current point
    int qlast=-1;                         // point of cp, sp ...
    double tlast=-1.0e300;                // ... at this time
    int cq=-1;                            // point of C, S ...
    double ct=-1.0e300;                   // ... at this time
    std::vector<double> P, iP, Q;         // work arrays of wave_lib_irregular_2nd_terms::kin

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
    void phases_pairs(int q, const wave_lib_irregular_2nd_terms &tt)
    {
        if(q==qlast && t==tlast && q==cq && t==ct)
        return;

        if(q!=cq || t!=ct)
        phases(q);

        tt.pair_phases(C.data(),S.data(),cp,sp);
        qlast = q;
        tlast = t;
    }

    // phases of point q at the cached time
    void phases(int q)
    {
        const double *cs=&cS[size_t(q)*M], *ss=&sS[size_t(q)*M];
        cq = q;
        ct = t;

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

        for(int n=0; n<M; ++n)
        {
            const double T = ki[n]*(cosb[n]*x + sinb[n]*y) - wi[n]*wt - ei[n];
            C[n] = cos(T);
            S[n] = sin(T);
        }
    }
};

#endif

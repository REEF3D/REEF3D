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
Author: Hans Bihs
--------------------------------------------------------------------*/

#ifndef SCOUR_GAUGE_H_
#define SCOUR_GAUGE_H_

#include<fstream>
#include<vector>

class lexer;

// Accumulative scour at a pile from the undisturbed near-bed velocity, no a priori equilibrium depth.
// Shared by FNPF and NHFLOW (fnpf_print_scour, nhflow_print_scour); the solver gives per gauge and time
// step the bed velocity (ub, vb), the velocity at D/2 above the bed (uc, vc, the current) and the depth h.
//
// The wave statistics are running (exponentially weighted, time constant Tw = P 86) means: the mean bed
// velocity, sigma_U^2 = var(ub) + var(vb), Um = sqrt(2) sigma_U, the current Uc = |mean(uc,vc)|,
// Tp (P 86, else the wave input, else 1.28 Tm02 of the bed velocity), KC = Um Tp/D,
// theta_w from Um with f_w = 0.04 (a/k_s)^-1/4, k_s = 2.5 d50, theta_c from Uc with the log law.
//
// P 83, Sumer-Fredsoe incremental (type 0):
//    S_eq/D = 1.3 K_I [1 - exp(-A (KC-B))], A = 0.03 + 0.75 Ucw^2.6, B = 6 exp(-4.7 Ucw)
//             (Sumer & Fredsoe 2001), the current limit blended in for Ucw > 0.75,
//             K_I = min(1, U/Ucr) clear water (Melville & Sutherland 1988), no change below U/Ucr = 0.5
//    T      = T* D^2/sqrt(g(s-1)d50^3), T*_w = 1e-6 (KC/theta_w)^3 (Sumer et al. 1992),
//             T*_c = 1/2000 (h/D) theta_c^-2.2, log-interpolated with Ucw
//    S_n+1  = S_eq + (S_n - S_eq) exp(-dt/T), backfilling (S > S_eq) with P86_cb T
//
// P 84, sediment balance (type 1): inverted cone around the pile, V = pi cot(phi) S^2 (D/2 + S cot(phi)/3)
//    (1-n) dV/dt = k D sqrt(g(s-1)d50^3) [Phi(alpha theta(t)) - Phi(theta(t))],  Phi = 8 (theta-theta_cr)^1.5
//    alpha(S) = alpha_0 exp(-beta S/D),  alpha_0 = amax^F,  F the KC shape of P 83,  beta = ln(amax)/1.3
//    theta(t) instantaneous: wave part 0.5 f_w |u~| u~ (u~ = ub - mean) plus the current stress (vector sum);
//    live bed it tends to S_eq of P 83, clear water to less, the history comes from the transport

class scour_gauge
{
public:
    // type 0: P 83 (Sumer-Fredsoe incremental), 1: P 84 (sediment balance); module "FNPF" or "NHFLOW"
    scour_gauge(lexer*, int type, const char *module);
    ~scour_gauge();

    // smp: per gauge ub, vb, uc, vc, h (scour_gauge::nsmp values, reduced, the same on all ranks);
    // ub < -1e19: the gauge lies outside the domain
    void step(lexer*, const double *smp);

    // checkpoints (P 48, I 45): the state of all gauges (the same on all ranks)
    void ckpt_state(std::vector<double>&) const;
    void ckpt_restore(const std::vector<double>&, size_t &pos);

    static const int nsmp = 5;

    int num;
    const double *x,*y,*D;

private:
    struct stats
    {
        double Um,Uc,Ucw,Tp,KC,fw,theta_w,theta_c,ufc2,F,r,h;
        double ecx,ecy;     // the direction of the current
    };

    void statistics(lexer*, int n, stats&) const;
    double seq_shape(double KC, double Ucw) const;
    void write(lexer*, int n, const stats&, double a, double b);
    double depth_from_volume(double V, double R, double S0) const;

    const int type;
    double g,s,d50,poro,thcr,cotphi,ks,q0;
    double Tw,Tp_in,cb,dtout;
    double amax,beta,kB;

    // the state per gauge: running weights and means, the scour depth, the scour hole volume
    std::vector<double> wsum,mub,mvb,muu,mvv,muc,mvc,maa,mh,ubp,vbp,hasp,tact,S,V;
    double nextprint;
    bool warned;

    std::ofstream *out;
};

#endif

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

// Romero, Melville & Kleiss (2012) breaking model for FNPF.
//
// Physics-based energy dissipation through matching:
//   rho * c->vb * |grad_H(Fi)|^2 = b(S) * rho * c^5 / (g * L_b)
//
// => c->vb = b(S) * c^5 / (g * L_b * |grad_H(Fi)|^2)
//
// with b(S) = a * (S - S0)^(5/2)
//
// Three methods for computing local steepness S:
//   1. Gradient-based:  S = |grad_H(eta)|  (simplest, underestimates)
//   2. Wavenumber-based: S = 0.5 * k * H_local  (via zero-crossing)
//   3. Hilbert-based:   S = k * A  (envelope from Hilbert transform)
//
// Two methods for local phase speed:
//   1. Linear dispersion: c = sqrt(g/k) for deep water,
//                         c = sqrt(g*tanh(kh)/k) general
//   2. Crest tracking:    c = dx_crest / dt  (requires time history)
//

#include"fnpf_breaking.h"
#include"lexer.h"
#include"fdm_fnpf.h"
#include"ghostcell.h"

void fnpf_breaking::breaking_romero(lexer *p, fdm_fnpf *c, ghostcell *pgc, slice &eta, slice &eta_n, slice &Fifsf, double alpha)
{

    // --- Initialize on first call ---
    if (ini_done == 0)
    {

        

        ini_done = 1;
    }

    // --- Step 0: Compute eta_t ---
    SLICELOOP4
    {
        eta_t(i,j) = (c->eta(i,j) - eta_old(i,j)) / p->dt;
    }
    pgc->gcsl_start4(p, eta_t, 1);

    // --- Step 1: Local wavenumber and steepness ---
    compute_local_wavenumber(p, c, pgc);
    compute_steepness(p, c, pgc);

    // --- Step 2: Local phase speed ---
    compute_phase_speed(p, c, pgc);

    // --- Step 3: Detect breaking ---
    detect_breaking(p, c, pgc);

    // --- Step 4: Compute eddy viscosity ---
    compute_eddy_viscosity(p, c, pgc);

    // --- Store eta for next step ---
    SLICELOOP4
    {
        eta_old(i,j) = c->eta(i,j);
    }
    
    
    SLICELOOP4
    c->breaklog(i,j)=0;
    
    // breaklog
    int count=0; 
    
    SLICELOOP4
    if(c->breaking(i,j)>0)
    {
    c->breaklog(i,j)=1;
    ++count;
    }
    
    /*
    LOOP
    {
    if(c->breaking(i,j)>0)
    c->test[IJK] = 1.0;
    
    else
    c->test[IJK] = 0.0;
    }*/
    
    LOOP
    c->test[IJK] = c->vb(i,j);
    
    count=pgc->globalisum(count);
    
    if(p->mpirank==0 && (p->count%p->P12==0))
    cout<<"breaking: "<<count<<endl;
}

// =====================================================================
// Local wavenumber estimation via zero-crossing analysis
// =====================================================================

void fnpf_breaking::compute_local_wavenumber(lexer *p, fdm_fnpf *c,
                                                     ghostcell *pgc)
{
    // Estimate local wavenumber from the free surface.
    // Method: compute second spatial derivative and use the relation
    //   k^2 ~ -nabla^2(eta) / eta  for a monochromatic wave
    // 
    // More robust alternative: zero-crossing distance in x and y.
    // Here we use a gradient + curvature approach that works for
    // multi-directional seas.

    double dx = p->DXM;
    double dy = p->DYM;
    double d2eta_dx2, d2eta_dy2, eta_c;
    double k_est;

    SLICELOOP4
    {
        // Curvature of eta
        d2eta_dx2 = (c->eta(i + 1, j) - 2.0 * c->eta(i,j) + c->eta(i - 1, j))
                    / (dx * dx);

        d2eta_dy2 = (c->eta(i, j + 1) - 2.0 * c->eta(i,j) + c->eta(i, j - 1))
                    / (dy * dy);

        eta_c = c->eta(i,j);

        // k^2 estimate from curvature/amplitude ratio
        // Only valid where |eta| is not near zero
        if (fabs(eta_c) > 1.0e-6)
        {
            k_est = sqrt(fabs(-(d2eta_dx2 + d2eta_dy2) / eta_c));

            // Sanity bounds: k must be positive and physical
            // Deep water: k = omega^2/g, for T=2s: k~1.0, T=20s: k~0.01
            k_est = MAX(k_est, 0.005);
            k_est = MIN(k_est, 5.0);
        }
        else
        {
            k_est = k_local(i,j); // keep previous value
        }

        k_local(i,j) = k_est;
    }

    pgc->gcsl_start4(p, k_local, 1);

    // Smooth the wavenumber field to reduce noise (3-point average)
    // Apply one pass of spatial smoothing
    SLICELOOP4
    {
        k_smooth(i,j) = (1.0 / 6.0) * (k_local(i + 1, j) + k_local(i - 1, j)
                        + k_local(i, j + 1) + k_local(i, j - 1))
                        + (2.0 / 6.0) * k_local(i,j);
    }

    SLICELOOP4
    {
        k_local(i,j) = k_smooth(i,j);
    }

    pgc->gcsl_start4(p, k_local, 1);
}

// =====================================================================
// Steepness computation
// =====================================================================

void fnpf_breaking::compute_steepness(lexer *p, fdm_fnpf *c,
                                              ghostcell *pgc)
{
    double dx = p->DXM;
    double dy = p->DYM;
    double deta_dx, deta_dy;
    double H_local, k;

    SLICELOOP4
    {
        if (steep_method == 1)
        {
            // --- Method 1: Gradient-based ---
            // S = |grad_H(eta)|
            // Simple but tends to underestimate steepness at crests
            // (where gradient = 0) and overestimate on flanks.

            deta_dx = (c->eta(i + 1, j) - c->eta(i - 1, j)) / (2.0 * dx);
            deta_dy = (c->eta(i, j + 1) - c->eta(i, j - 1)) / (2.0 * dy);

            S_local(i,j) = sqrt(deta_dx * deta_dx + deta_dy * deta_dy);
        }
        else if (steep_method == 2)
        {
            // --- Method 2: Wavenumber-based ---
            // S = k * H/2 where H is estimated from local crest-trough
            //
            // For a wave eta = A*cos(kx), the steepness is S = kA.
            // We estimate A from the local amplitude: A ~ |eta| at
            // locations near crests/troughs, and from the envelope
            // otherwise.
            //
            // Practical approach: use the RMS amplitude in a local
            // stencil to estimate A, then S = k * A.

            k = k_local(i,j);

            // Local RMS amplitude from a 5-point stencil
            double eta2_sum = 0.0;
            int count = 0;

            for (int di = -2; di <= 2; di++)
            for (int dj = -2; dj <= 2; dj++)
            {
                eta2_sum += c->eta(i + di, j + dj) * c->eta(i + di, j + dj);
                count++;
            }

            double A_rms = sqrt(2.0 * eta2_sum / (double)count);

            S_local(i,j) = k * A_rms;
        }
        else if (steep_method == 3)
        {
            // --- Method 3: Hilbert-based ---
            // Compute envelope via discrete Hilbert transform along x.
            // This is the most accurate but most expensive method.
            // For now, use the analytic signal approximation:
            //   A(x) = sqrt(eta^2 + H[eta]^2)
            // where H[eta] ~ -d(eta)/dx / k  is a simple approximation
            // of the Hilbert transform for narrowband signals.

            k = k_local(i,j);

            deta_dx = (c->eta(i + 1, j) - c->eta(i - 1, j)) / (2.0 * dx);

            double eta_c = c->eta(i,j);
            double hilbert_approx = 0.0;

            if (k > 1.0e-6)
                hilbert_approx = -deta_dx / k;

            double A_hilbert = sqrt(eta_c * eta_c + hilbert_approx * hilbert_approx);

            S_local(i,j) = k * A_hilbert;
        }

        // Clamp steepness to physical range [0, 0.5]
        // Stokes limiting steepness is ~0.44
        S_local(i,j) = MAX(S_local(i,j), 0.0);
        S_local(i,j) = MIN(S_local(i,j), 0.5);
    }

    pgc->gcsl_start4(p, S_local, 1);
}

// =====================================================================
// Phase speed computation
// =====================================================================

void fnpf_breaking::compute_phase_speed(lexer *p, fdm_fnpf *c,
                                                ghostcell *pgc)
{
    SLICELOOP4
    {
        if (cspeed_method == 1)
        {
            // --- Method 1: Linear dispersion relation ---
            // c = sqrt(g * tanh(k*h) / k)

            double k = k_local(i,j);
            double h_local = p->wd + c->eta(i,j);
            h_local = MAX(h_local, 0.01);

            double kh = k * h_local;

            if (kh > 10.0)
            {
                // Deep water limit: c = sqrt(g/k)
                c_phase(i,j) = sqrt(9.81 / k);
            }
            else
            {
                c_phase(i,j) = sqrt(9.81 * tanh(kh) / k);
            }
        }
        else if (cspeed_method == 2)
        {
            // --- Method 2: Crest tracking ---
            // Estimate c from eta_t / eta_x at points near crests.
            // c ~ -eta_t / eta_x  (from the kinematic relation for
            //                      a propagating surface)
            //
            // This is noisy, so we use it only where |eta_x| is
            // sufficiently large and smooth the result.

            double dx = p->DXM;
            double deta_dx = (c->eta(i + 1, j) - c->eta(i - 1, j)) / (2.0 * dx);

            if (fabs(deta_dx) > 1.0e-4)
            {
                double c_est = -eta_t(i,j) / deta_dx;

                // Only accept physically reasonable values
                // (positive, forward-propagating, bounded)
                if (c_est > 0.1 && c_est < 50.0)
                    c_phase(i,j) = c_est;
                // else keep previous value
            }
        }

        // Safety bounds on phase speed
        c_phase(i,j) = MAX(c_phase(i,j), 0.1);
        c_phase(i,j) = MIN(c_phase(i,j), 50.0);
    }

    pgc->gcsl_start4(p, c_phase, 1);
}

// =====================================================================
// Breaking detection
// =====================================================================

void fnpf_breaking::detect_breaking(lexer *p, fdm_fnpf *c,
                                            ghostcell *pgc)
{
    double h_local;

    SLICELOOP4
    {
        h_local = p->wd + c->eta(i,j);
        h_local = MAX(h_local, 0.01);

        // --- Primary criterion: steepness exceeds threshold ---
        int steep_trigger = (S_local(i,j) > S0) ? 1 : 0;

        // --- Optional secondary criterion: eta_t pre-filter ---
        // If eta_t_thresh > 0, also require eta_t > threshold
        // to avoid false positives on backward-facing slopes
        int etat_trigger = 1; // default: always pass
        if (eta_t_thresh > 0.0)
        {
            double threshold = eta_t_thresh * sqrt(9.81 * h_local);
            etat_trigger = (eta_t(i,j) > threshold) ? 1 : 0;
        }

        // --- Breaking onset ---
        if (c->breaking(i,j)== 0 && steep_trigger == 1 && etat_trigger == 1)
        {
            c->breaking(i,j)= 1;
            t_break(i,j)  = p->simtime;
        }

        // --- Breaking cessation ---
        // Cease when steepness drops well below threshold
        if (c->breaking(i,j)== 1 && S_local(i,j) < 0.5 * S0)
        {
            c->breaking(i,j)= 0;
            B_ramp(i,j)   = 0.0;
            c->vb(i,j)     = 0.0;
        }
    }
}

// =====================================================================
// Eddy viscosity computation via energy dissipation matching
// =====================================================================

void fnpf_breaking::compute_eddy_viscosity(lexer *p, fdm_fnpf *c,
                                                    ghostcell *pgc)
{
    double dx = p->DXM;
    double dy = p->DYM;
    double h_local, k, cp, S, b, Lb;
    double dFidx, dFidy, gradFi2;
    double T_ramp, dt_break;

    SLICELOOP4
    {
        if (c->breaking(i,j)== 1)
        {
            h_local = p->wd + c->eta(i,j);
            h_local = MAX(h_local, 0.01);

            k  = k_local(i,j);
            cp = c_phase(i,j);
            S  = S_local(i,j);

            // ---- Breaking strength: b(S) = a * (S - S0)^(5/2) ----
            if (S > S0)
                b = a_coeff * pow(S - S0, 2.5);
            else
                b = 0.0;

            b_strength(i,j) = b;

            // ---- Breaking zone width ----
            // L_b = alpha_L * wavelength = alpha_L * 2*pi/k
            if (k > 1.0e-6)
                Lb = alpha_L * 2.0 * PI / k;
            else
                Lb = alpha_L * 100.0; // fallback

            Lb = MAX(Lb, dx); // at least one grid cell
            L_break(i,j) = Lb;

            // ---- Surface velocity magnitude from Fi ----
            dFidx = (c->Fifsf(i + 1, j) - c->Fifsf(i - 1, j)) / (2.0 * dx);
            dFidy = (c->Fifsf(i, j + 1) - c->Fifsf(i, j - 1)) / (2.0 * dy);
            gradFi2 = dFidx * dFidx + dFidy * dFidy;

            // ---- Temporal ramp B(t) ----
            // T_ramp = T_star * sqrt(h/g)
            T_ramp = T_star * sqrt(h_local / 9.81);
            dt_break = p->simtime - t_break(i,j);

            if (dt_break < T_ramp && T_ramp > 1.0e-10)
                B_ramp(i,j) = dt_break / T_ramp;
            else
                B_ramp(i,j) = 1.0;

            B_ramp(i,j) = MIN(B_ramp(i,j), 1.0);
            B_ramp(i,j) = MAX(B_ramp(i,j), 0.0);

            // ---- Eddy viscosity via energy matching ----
            // c->vb = B(t) * b(S) * c^5 / (g * L_b * |grad_H(Fi)|^2)
            if (gradFi2 > 1.0e-10 && b > 0.0)
            {
                c->vb(i,j) = B_ramp(i,j) * b * pow(cp, 5.0)
                             / (9.81 * Lb * gradFi2);
            }
            else if (b > 0.0)
            {
                // Fallback when surface velocity is very small:
                // Use the simplified deep-water expression
                // c->vb = B(t) * b * c / (2*pi * alpha_L * S^2)
                if (S > 1.0e-6)
                    c->vb(i,j) = B_ramp(i,j) * b * cp
                                 / (2.0 * PI * alpha_L * S * S);
                else
                    c->vb(i,j) = 0.0;
            }
            else
            {
                c->vb(i,j) = 0.0;
            }

            // ---- Apply limiter ----
            c->vb(i,j) = MIN(c->vb(i,j), vb_max);
            c->vb(i,j) = MAX(c->vb(i,j), 0.0);
        }
        else
        {
            c->vb(i,j)       = 0.0;
            b_strength(i,j) = 0.0;
            B_ramp(i,j)     = 0.0;
        }
    }

    pgc->gcsl_start4(p, c->vb, 1);
}

/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2026 Hans Bihs

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTIBILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
Authors: Hans Bihs
--------------------------------------------------------------------*/

#include"partres.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"sediment_fdm.h"

void partres::stress_schaeffer(lexer *p, ghostcell *pgc, sediment_fdm *s)
{
    /*double Ps = 10.0;
    double beta = 2.0;
    double epsilon = 1.0e-4;
    double Tc = 0.2;
    double Tmax = 1.1*(1.0-p->S24) + 0.1;
    
    //double Tmax = (1.0-p->S24) + 0.0;
    
    double maxTau = 1.0e7;

    ALOOP
    {        
    if(Ts(i,j,k)<=Tc)
    Tau(i,j,k) = 0.0;
        
    // Total stress
    double P_kinetic = (theta_p < theta_p_critical) ? kinetic_stress(theta_p, Theta) : 0.0;
    
    double P_fric = compute_frictional_stress(Ts(i,j,k), shear_rate);

    Tau(i,j,k) = P_kinetic + P_fric;

        
    Tau(i,j,k) = MIN(Tau(i,j,k), maxTau);
        
        //cout<<"Tau: "<<Tau(i,j,k)<<" Ts: "<<Ts(i,j,k)<<" MAXfunc: "<<MAX(Tc-Ts(i,j,k),epsilon*(1.0-Ts(i,j,k)))<<endl;
    }

    pgc->start4a(p,Tau,1);
    pgc->start4a(p,Ts,1);*/
}
/*

// Schaeffer frictional stress
double partres::compute_frictional_stress(double theta_p, double shear_rate) 
{
    double theta_p_critical = 0.50;  // transition to frictional regime
    double theta_p_max = 0.635;      // random close packing
    
    if (theta_p < theta_p_critical) {
        return 0.0;  // no frictional stress in dilute regime
    }
    
    // Frictional pressure
    double phi = 32.0 * PI / 180.0;  // angle of internal friction (degrees->radians)
    double sin_phi = sin(phi);
    double Fr = 0.05;  // dimensionless coefficient
    double n = 2.0;    // exponent
    
    double P_fric = Fr * pow((theta_p - theta_p_critical) / 
                             (theta_p_max - theta_p), n);
    
    // Frictional shear stress
    double tau_fric = sin_phi * P_fric;
    
    return P_fric;
}




// Kinetic stress from granular kinetic theory
double partres::compute_kinetic_stress(double theta_p, double Theta) 
{
    // Theta = granular temperature (m²/s²)
    //       = 1/3 * <v'_i * v'_i> (mean square particle velocity fluctuation)
    
    double e = 0.9;  // coefficient of restitution (0.8-0.95 typical)
    double rho_p = 2650.0;  // particle density
    
    // Radial distribution function at contact
    double g0 = radial_distribution_function(theta_p);
    
    // Kinetic pressure (collisional contribution)
    double P_kinetic = rho_p * theta_p * Theta * 
                       (1.0 + 2.0 * (1.0 + e) * theta_p * g0);
    
    return P_kinetic;
}

// Radial distribution function (Carnahan-Starling)
double partres::radial_distribution_function(double theta_p) {
    double theta_p_max = 0.635;  // maximum packing
    
    // Original form
    return 1.0 / (1.0 - pow(theta_p / theta_p_max, 1.0/3.0));
    
    // Alternative (Ma-Ahmadi form, more stable)
    // double eta = theta_p / theta_p_max;
    // return (2.0 - eta) / (2.0 * pow(1.0 - eta, 3.0));
}



// Compute shear rate from PARTICLE velocity field
double partres::compute_shear_rate_from_particles(Grid* grid, int i, int j, int k) {
    // Use mapped particle velocities (u_p, v_p, w_p)
    // NOT fluid velocities (u_f, v_f, w_f)
    
    double du_p_dx = (grid->u_p[i+1][j][k] - grid->u_p[i-1][j][k]) / (2*dx);
    double du_p_dy = (grid->u_p[i][j+1][k] - grid->u_p[i][j-1][k]) / (2*dy);
    double du_p_dz = (grid->u_p[i][j][k+1] - grid->u_p[i][j][k-1]) / (2*dz);
    
    double dv_p_dx = (grid->v_p[i+1][j][k] - grid->v_p[i-1][j][k]) / (2*dx);
    double dv_p_dy = (grid->v_p[i][j+1][k] - grid->v_p[i][j-1][k]) / (2*dy);
    double dv_p_dz = (grid->v_p[i][j][k+1] - grid->v_p[i][j][k-1]) / (2*dz);
    
    double dw_p_dx = (grid->w_p[i+1][j][k] - grid->w_p[i-1][j][k]) / (2*dx);
    double dw_p_dy = (grid->w_p[i][j+1][k] - grid->w_p[i][j-1][k]) / (2*dy);
    double dw_p_dz = (grid->w_p[i][j][k+1] - grid->w_p[i][j][k-1]) / (2*dz);
    
    // Strain rate tensor for PARTICLE phase
    double S_xx = du_p_dx;
    double S_yy = dv_p_dy;
    double S_zz = dw_p_dz;
    double S_xy = 0.5 * (du_p_dy + dv_p_dx);
    double S_xz = 0.5 * (du_p_dz + dw_p_dx);
    double S_yz = 0.5 * (dv_p_dz + dw_p_dy);
    
    double II = S_xx*S_xx + S_yy*S_yy + S_zz*S_zz +
                2.0 * (S_xy*S_xy + S_xz*S_xz + S_yz*S_yz);
    
    return sqrt(2.0 * II);
}



// Trilinear mapping: particle distributes to 8 surrounding nodes
void partres::map_particles_to_grid_trilinear(Particle* particles, int n_particles, 
                                     Grid* grid) {
    // Initialize
    for (int i = 0; i < grid->nnodes; i++) {
        grid->u_p[i] = 0.0;
        grid->v_p[i] = 0.0;
        grid->w_p[i] = 0.0;
        grid->volume_sum[i] = 0.0;
    }
    
    // Accumulate from particles
    for (int p = 0; p < n_particles; p++) {
        // Find base cell (lower-left-bottom corner)
        int i_base, j_base, k_base;
        find_cell_indices(particles[p].x, grid, &i_base, &j_base, &k_base);
        
        // Fractional position within cell [0, 1]
        double x_local = (particles[p].x[0] - grid->x[i_base]) / grid->dx;
        double y_local = (particles[p].x[1] - grid->y[j_base]) / grid->dy;
        double z_local = (particles[p].x[2] - grid->z[k_base]) / grid->dz;
        
        // Clamp to [0, 1]
        x_local = max(0.0, min(x_local, 1.0));
        y_local = max(0.0, min(y_local, 1.0));
        z_local = max(0.0, min(z_local, 1.0));
        
        // Trilinear weights for 8 corner nodes
        double w[8];
        w[0] = (1-x_local) * (1-y_local) * (1-z_local);  // (i,j,k)
        w[1] = x_local     * (1-y_local) * (1-z_local);  // (i+1,j,k)
        w[2] = (1-x_local) * y_local     * (1-z_local);  // (i,j+1,k)
        w[3] = x_local     * y_local     * (1-z_local);  // (i+1,j+1,k)
        w[4] = (1-x_local) * (1-y_local) * z_local;      // (i,j,k+1)
        w[5] = x_local     * (1-y_local) * z_local;      // (i+1,j,k+1)
        w[6] = (1-x_local) * y_local     * z_local;      // (i,j+1,k+1)
        w[7] = x_local     * y_local     * z_local;      // (i+1,j+1,k+1)
        
        // Node offsets
        int di[8] = {0, 1, 0, 1, 0, 1, 0, 1};
        int dj[8] = {0, 0, 1, 1, 0, 0, 1, 1};
        int dk[8] = {0, 0, 0, 0, 1, 1, 1, 1};
        
        double V_p = particles[p].volume;
        
        // Distribute to 8 nodes
        for (int n = 0; n < 8; n++) {
            int i_node = i_base + di[n];
            int j_node = j_base + dj[n];
            int k_node = k_base + dk[n];
            
            // Check bounds
            if (i_node >= grid->nx || j_node >= grid->ny || k_node >= grid->nz) 
                continue;
            
            int node_idx = i_node + j_node*grid->nx + k_node*grid->nx*grid->ny;
            
            // Add weighted contribution
            grid->u_p[node_idx] += w[n] * particles[p].v[0] * V_p;
            grid->v_p[node_idx] += w[n] * particles[p].v[1] * V_p;
            grid->w_p[node_idx] += w[n] * particles[p].v[2] * V_p;
            grid->volume_sum[node_idx] += w[n] * V_p;
        }
    }
    
    // Normalize
    double V_node = grid->dx * grid->dy * grid->dz;  // volume per node
    for (int i = 0; i < grid->nnodes; i++) {
        if (grid->volume_sum[i] > 1e-12) {
            grid->u_p[i] /= grid->volume_sum[i];
            grid->v_p[i] /= grid->volume_sum[i];
            grid->w_p[i] /= grid->volume_sum[i];
            grid->theta_p[i] = grid->volume_sum[i] / V_node;
        }
    }
}
*/
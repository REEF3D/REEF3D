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

#ifndef SEDIMENT_BEDSHEAR_H_
#define SEDIMENT_BEDSHEAR_H_

#include"norm_vec.h"
#include"slice4.h"

class lexer;
class fdm;
class fdm_nhf;
class fdm2D;
class ghostcell;
class sediment_fdm;
class turbulence;
class sliceint;

using namespace std;

class bedshear :  public norm_vec
{
public:
    bedshear(lexer*,turbulence*);
    virtual ~bedshear();

	void taubed(lexer*, fdm*,ghostcell*,sediment_fdm*);
    void taubed_curv(lexer *p, fdm2D *b, ghostcell *pgc, sediment_fdm *s);     // curvilinear grid (G 50)
	void taucritbed(lexer*, fdm*,ghostcell*,sediment_fdm*);
    
    void taubed(lexer*, fdm_nhf*, ghostcell*, sediment_fdm*);
    void taucritbed(lexer*, fdm_nhf*, ghostcell*, sediment_fdm*);
    
    void taubed(lexer*, fdm2D*,ghostcell*,sediment_fdm*);
    void taucritbed(lexer*, fdm2D*,ghostcell*,sediment_fdm*);

	const double kappa;
    
    // bed shear stress of the bed column (i,j) from the flow above the bed level set, shared by the Eulerian
    // bed shear stress (S 16 10) and the CPM bedload layer (Q 58): log law at the first cell centre at least
    // half a cell above the bed, or href cells above the bed (href > 0), and with kterm max(log law, 0.3 rho k)
    // with k of the first fluid cell; returns false if the column has no fluid cell above the bed
    static bool column_tau(lexer*, fdm*, turbulence*, int i, int j, double zbed, double ks, double href, bool kterm,
                           double rho, double &tx, double &ty, double &hcell, bool &rest);

    // the same for the NHFLOW sigma grid (S 16 10): the first cell centre is half a cell above the bed, href > 0
    // takes the log law href cells of the first cell above the bed (at most half the water depth)
    static bool column_tau(lexer*, fdm_nhf*, int i, int j, double ks, double href, bool kterm,
                           double rho, double &tx, double &ty);

    // the closure of column_tau: rho (u / u+)^2 with u+ = ln(30 z_r/k_s)/kappa at the reference height z_r,
    // with kterm max(log law, 0.3 rho k)
    static double tau_logk(double um, double zr, double ks, double kin, bool kterm, double rho);
    
    // rough-wall log law u/u* = ln(30 z/ks)/kappa, with z0 = ks/30 clamped so that u+ >= 1/kappa
    // (as the NHFLOW wall drag): no singular shear stress when z approaches ks/30 in thin water
    double uplus(double z, double ks) const
    {
    double r = 30.0*z/(ks>1.0e-12?ks:1.0e-12);
    return (1.0/kappa)*log(r>2.718281828459045?r:2.718281828459045);
    }
    
    // S 16 9: wave friction bed shear stress 0.5 rho fw |u|^2 of the near-bed velocity sample (u,v)
    double tau_wave(lexer*, double, double, double, double);
    void tau_wave_step(lexer*);


private:
    turbulence *pturb;
    double tau,tauc;
    double u_abs,u_plus,dist;
    double nx,ny,nz,norm;
    double xip,yip,zip;
    double uvel,vvel,wvel;
	double kinval;
    double scale;
    
    // S 16 9: running mean of the sample (current part), oscillatory part of the previous call and its
    // amplitude, for A = Uhat T/2pi
    slice4 uw_mean,vw_mean,uw_old,vw_old,uw_hat;
    double tw_old,tw_dt;
    int tw_count;
};

#endif

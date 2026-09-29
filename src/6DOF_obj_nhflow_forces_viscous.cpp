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

#include"6DOF_obj.h"
#include"gradient.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"

namespace
{
    // ITTC-1957 model-ship correlation line
    double ittc57_CF(double Re)
    {
        return 0.075/pow(log10(Re) - 2.0, 2.0);
    }
    
    // Local skin friction consistent with ITTC-1957: for a flat plate
    // cf(Re_x) = d(Re_x*CF)/dRe_x = CF*(1 - 2/(ln10*(log10(Re_x) - 2))),
    // so that (1/L) int_0^L cf dx = CF(Re_L) exactly.
    // Close to Schlichting's cf = (2log10(Re_x) - 0.65)^-2.3 for 1e5 < Re_x < 1e9.
    double ittc57_cf_local(double Re)
    {
        return ittc57_CF(Re)*(1.0 - 2.0/(log(10.0)*(log10(Re) - 2.0)));
    }
    
    // Below Re_min (bow region) the local law is replaced by its mean value CF(Re_min),
    // which keeps the plate integral exact. Also the lower Re limit for X 39 1.
    const double Re_min = 1.0e5;
}

void sixdof_obj::hydrodynamic_viscous_forces_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc, slice &WL,
                                            double &Fv_x, double &Fv_y, double &Fv_z, double A_triang,
                                            double xp, double yp, double zp, double nx, double ny, double nz)
{
    double uval, vval, wval;
    double Cf=0.0001;
    
    Fv_x = Fv_y = Fv_z = 0.0;
    
    if(p->X38==1)
    {
    xc = xp + p->X43*nx*p->DXP[IP];
    yc = yp + p->X43*ny*p->DYP[JP];
    zc = zp + p->X43*nz*p->DZP[KP]*WL(i,j);
        
    uval   = p->ccipol4V(d->U, WL, d->bed, xc, yc, zc); 
    vval   = p->ccipol4V(d->V, WL, d->bed, xc, yc, zc); 
    wval   = p->ccipol4V(d->W, WL, d->bed, xc, yc, zc); 
    

	Fv_x = 0.5*p->W1*Cf*uval*fabs(uval);
    Fv_y = 0.5*p->W1*Cf*vval*fabs(vval);
    Fv_z = 0.5*p->W1*Cf*wval*fabs(wval);
    
    if(p->j_dir==0)
    Fv_y = 0.0;
    }
    
    
}

void sixdof_obj::viscous_forces_ittc_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc, slice &WL,
                                            double &Fv_x, double &Fv_y, double &Fv_z,
                                            double &Kv, double &Mv, double &Nv)
{
    // ITTC-1957 friction line with form factor (1+k), applied as a wall shear stress
    // on the wetted sub-triangles collected in force_calc_stl:
    //
    //   tau = 0.5*rho*(1+k)*cf*|u_t|*u_t
    //
    // u_t: fluid velocity relative to the moving hull, sampled X43 mean cell sizes off the
    //      wall and projected onto the hull tangent plane
    // X 39 1: cf = CF(Re_L),  Re_L = |u_t|*Lwl/nu          (global friction line)
    // X 39 2: cf = cf(Re_x),  Re_x = |u_t|*x/nu            (local, ITTC-consistent)
    //         x: distance from the bow along the mean relative flow direction
    //
    // Returns the local (per-rank) force and moment about c_; the caller sums globally.
    
    const int nv = vis_A.size();
    
    vector<double> ut(nv), vt(nv), wt(nv);
    
    Fv_x = Fv_y = Fv_z = 0.0;
    Kv = Mv = Nv = 0.0;
    
    double S=0.0, SU=0.0, SUx=0.0, SUy=0.0, SUz=0.0;
    
    for(int n=0; n<nv; ++n)
    {
        const double nx = vis_nx[n];
        const double ny = vis_ny[n];
        const double nz = vis_nz[n];
        
        // sample point off the wall (DSM: mean horizontal cell size, as for the pressure)
        const double xs = vis_x[n] + p->X43*nx*DSM;
        const double ys = vis_y[n] + p->X43*ny*DSM;
        const double zs = vis_z[n] + p->X43*nz*DSM;
        
        // hull velocity at the wall point
        const double rx = vis_x[n] - c_(0);
        const double ry = vis_y[n] - c_(1);
        const double rz = vis_z[n] - c_(2);
        
        const double ub = u_fb(0) + u_fb(4)*rz - u_fb(5)*ry;
        const double vb = u_fb(1) + u_fb(5)*rx - u_fb(3)*rz;
        const double wb = u_fb(2) + u_fb(3)*ry - u_fb(4)*rx;
        
        double ur = p->ccipol4V(d->U, WL, d->bed, xs, ys, zs) - ub;
        double vr = p->ccipol4V(d->V, WL, d->bed, xs, ys, zs) - vb;
        double wr = p->ccipol4V(d->W, WL, d->bed, xs, ys, zs) - wb;
        
        if(p->j_dir==0)
        vr = 0.0;
        
        // tangential part
        const double un = ur*nx + vr*ny + wr*nz;
        
        ut[n] = ur - un*nx;
        vt[n] = vr - un*ny;
        wt[n] = wr - un*nz;
        
        const double A = vis_A[n];
        
        S   += A;
        SU  += A*sqrt(ut[n]*ut[n] + vt[n]*vt[n] + wt[n]*wt[n]);
        SUx += A*ut[n];
        SUy += A*vt[n];
        SUz += A*wt[n];
    }
    
    S   = pgc->globalsum(S);
    SU  = pgc->globalsum(SU);
    SUx = pgc->globalsum(SUx);
    SUy = pgc->globalsum(SUy);
    SUz = pgc->globalsum(SUz);
    
    // mean relative flow direction over the wetted hull
    const double SUn = sqrt(SUx*SUx + SUy*SUy + SUz*SUz);
    
    if(S<1.0e-20 || SUn<1.0e-12*MAX(S,1.0))
    return;
    
    const double ex = SUx/SUn;
    const double ey = SUy/SUn;
    const double ez = SUz/SUn;
    
    // bow: most upstream wetted point
    double s_bow = 1.0e20;
    
    if(p->X39==2)
    {
        for(int n=0; n<nv; ++n)
        s_bow = MIN(s_bow, vis_x[n]*ex + vis_y[n]*ey + vis_z[n]*ez);
        
        s_bow = pgc->globalmin(s_bow);
    }
    
    const double fac = 0.5*p->W1*(1.0 + p->X39_k);
    
    for(int n=0; n<nv; ++n)
    {
        const double Ut = sqrt(ut[n]*ut[n] + vt[n]*vt[n] + wt[n]*wt[n]);
        
        double cf;
        
        if(p->X39==2)
        {
            const double x  = MAX(vis_x[n]*ex + vis_y[n]*ey + vis_z[n]*ez - s_bow, 0.0);
            const double Re = Ut*x/p->W2;
            
            cf = (Re>Re_min) ? ittc57_cf_local(Re) : ittc57_CF(Re_min);
        }
        else
        {
            const double Re = Ut*p->X39_Lwl/p->W2;
            
            cf = ittc57_CF(MAX(Re, Re_min));
        }
        
        const double tA = fac*cf*Ut*vis_A[n];
        const double fx = tA*ut[n];
        const double fy = tA*vt[n];
        const double fz = tA*wt[n];
        
        const double rx = vis_x[n] - c_(0);
        const double ry = vis_y[n] - c_(1);
        const double rz = vis_z[n] - c_(2);
        
        Fv_x += fx;
        Fv_y += fy;
        Fv_z += fz;
        Kv += ry*fz - rz*fy;
        Mv += rz*fx - rx*fz;
        Nv += rx*fy - ry*fx;
    }
    
    if(p->j_dir==0)
    Fv_y = Kv = Nv = 0.0;
    
    // Diagnostics: effective friction coefficient CF_eff = F.e/((1+k)*0.5*rho*U_m^2*S)
    // vs. the ITTC-1957 line at Re_L = U_m*Lwl/nu (U_m: area-mean |u_t|)
    const double Fe = pgc->globalsum(Fv_x*ex + Fv_y*ey + Fv_z*ez);
    
    if(p->mpirank==0)
    {
        const double Um = SU/S;
        const double CF_eff = Fe/(fac*Um*Um*S);
        
        cout<<"ITTC-57 viscous: S_wet: "<<S<<" U_m: "<<Um<<" F_v: "<<Fe
            <<" CF_eff: "<<CF_eff<<" CF_ITTC(Lwl): "<<ittc57_CF(MAX(Um*p->X39_Lwl/p->W2, Re_min))<<endl;
    }
}

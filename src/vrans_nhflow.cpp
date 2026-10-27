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

#include"vrans_nhflow_f.h"
#include"nhflow_por_ghost.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"

vrans_nhflow_f::vrans_nhflow_f(lexer *p, fdm_nhf *d, ghostcell *pgc) : nhflow_geometry(p,d,pgc), Cval(p->B264)
{
    p->Darray(UN,p->imax*p->jmax*(p->kmax+2));
    p->Darray(VN,p->imax*p->jmax*(p->kmax+2));
    p->Darray(WN,p->imax*p->jmax*(p->kmax+2));
    
    p->Darray(P,p->imax*p->jmax*(p->kmax+2));
    
    // porous layers (B 202)
    LDEP = LN = LD50 = LALPHA = LBETA = nullptr;
    
    if(p->B202>0)
    {
    p->Darray(LDEP,p->imax*p->jmax*(p->kmax+2));
    p->Darray(LN,p->imax*p->jmax*(p->kmax+2));
    p->Darray(LD50,p->imax*p->jmax*(p->kmax+2));
    p->Darray(LALPHA,p->imax*p->jmax*(p->kmax+2));
    p->Darray(LBETA,p->imax*p->jmax*(p->kmax+2));
    
    // depth of the inner boundary of each layer
    double T=0.0;
    
    for(int q=0; q<p->B202; ++q)
    {
        if(p->B202_t[q]<=0.0 && p->mpirank==0)
        cout<<"VRANS porous layers: !!! B 202 layer "<<q+1<<" has thickness "<<p->B202_t[q]<<" and is ignored !!!"<<endl;
        
        T += MAX(p->B202_t[q],0.0);
        layer_T.push_back(T);
    }
    }
    
    print_force_ini(p,d,pgc);
    
    // UN is set at the first force_calc, dU/dt starts with the second call
    cmfac=0.0;
}

vrans_nhflow_f::~vrans_nhflow_f()
{
    delete [] LDEP;
    delete [] LN;
    delete [] LD50;
    delete [] LALPHA;
    delete [] LBETA;
}

void vrans_nhflow_f::update(lexer *p, fdm_nhf *d, ghostcell *pgc, double alpha, int val)
{
    ray_cast(p, d, pgc, d->PORSTRUC);
    reini_RK2(p, d, pgc, d->PORSTRUC);
    
    // porous layers (B 202): depth below the exposed surface -> layered n, d50, alpha, beta
    if(p->B202>0)
    {
    layer_dist(p, d, layer_T.back(), d->PORSTRUC, LDEP);
    layer_prop(p, d);
    }

    LOOP
    {
    H = Hporface(p,d,0,0,0);
    
    porval  = (p->B202>0) ? LN[IJK]   : p->B201_n;
    partval = (p->B202>0) ? LD50[IJK] : p->B201_d50;
    
    d->POR[IJK]     = H*porval + (1.0-H)*1.0;
	d->PORPART[IJK] = H*partval;
    }
    
    // porous floating body (X 16): POR is reset above, so re-apply the moving body porosity
    // n = 1 - H_fb(1 - n_fb), same Heaviside as sixdof_obj::Hsolidface_nhflow()
    if(p->X10>0 && p->X16==1)
    LOOP
    {
        double psi,Hfb;
        
        if(p->j_dir==0)
        psi = p->A526*p->DXN[IP];
        
        if(p->j_dir==1)
        psi = p->A526*0.5*(p->DXN[IP] + p->DYN[JP]);
        
        Hfb = 0.5*(1.0 + (-d->FB[IJK])/psi + (1.0/PI)*sin((PI*(-d->FB[IJK]))/psi));
        
        if(-d->FB[IJK] > psi)
        Hfb = 1.0;
        
        if(-d->FB[IJK] < -psi)
        Hfb = 0.0;
        
        d->POR[IJK] = MIN(d->POR[IJK], 1.0 - Hfb*(1.0 - p->X16_n));
    }
    
    // porosity of the DEM particles (E 28), 1 otherwise
    LOOP
    d->POR[IJK] *= d->PORDEM[IJK];

    pgc->start5Vfull(p,d->POR,1);
    pgc->start5Vfull(p,d->PORPART,1);
    
    // below the bed / above the surface of all columns incl. halo and ghost columns (nhflow_por_ghost.h)
    nhflow_por_ghost(p,d);
    
    // porosity at the free surface: storage coefficient of the continuity
    por_surface(p,d,pgc);
    
    // print force
    if(p->B208==1)
    force_calc(p, d, pgc, alpha, val);
}

// Porous layers (B 202 t n d50 alpha beta, counted from the exposed surface inwards; B 201 is the core).
// Layer q covers the depth T_(q-1) < dep < T_q. The properties are blended from the core outwards
// with the sine Heaviside of Hporface at every layer boundary, width psi = B209 dx:
//   f = f_core,   f = w_q f_q + (1 - w_q) f,  q = last..0,   w_q = H(T_q - dep)
void vrans_nhflow_f::layer_prop(lexer *p, fdm_nhf *d)
{
    double psi,dep,x,w;
    double nval,dval,aval,bval;
    
    LOOP
    {
        psi = porpsi(p,d);
        
        dep = LDEP[IJK];
        
        nval = p->B201_n;
        dval = p->B201_d50;
        aval = p->B201_alpha;
        bval = p->B201_beta;
        
        for(int q=p->B202-1; q>=0; --q)
        {
            if(p->B202_t[q]<=0.0)
            continue;
            
            x = layer_T[q] - dep;
            
            if(x > psi)
            w = 1.0;
            
            else if(x < -psi)
            w = 0.0;
            
            else
            w = 0.5*(1.0 + x/psi + (1.0/PI)*sin((PI*x)/psi));
            
            nval = w*p->B202_n[q]     + (1.0-w)*nval;
            dval = w*p->B202_d50[q]   + (1.0-w)*dval;
            aval = w*p->B202_alpha[q] + (1.0-w)*aval;
            bval = w*p->B202_beta[q]  + (1.0-w)*bval;
        }
        
        LN[IJK]     = nval;
        LD50[IJK]   = dval;
        LALPHA[IJK] = aval;
        LBETA[IJK]  = bval;
    }
}

// Porosity at the free surface n(eta), the storage coefficient of the continuity
//   d/dt int_{-h}^{eta} n dz + div(Q) = 0  ->  n(eta) deta/dt + div(Q) = 0   (n static)
// The surface cell centre lies half a sigma cell below eta. With its porosity, a free surface within
// dz/2 above the structure surface was stored with n_armour instead of 1 (deta/dt 1/n too large), and
// the storage jumped when the top cell centre crossed the Heaviside band during run-up and run-down.
// The static-structure level set is extrapolated from the two top cells to z = eta (vertical slope
// limited to the distance-function bound |d phi/dz| <= 1), and the factors of porous floating bodies
// and DEM particles in POR of the top cell are kept.
void vrans_nhflow_f::por_surface(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    double phit,phib,zt,zb,zs,g,phis,psi,Hs,Ht,nl,nst,nss;
    
    SLICELOOP4
    {
        k = p->knoz-1;
        
        d->PORS(i,j) = d->POR[IJK];
        
        if(p->wet[IJ]==0 || p->knoz<2)
        continue;
        
        phit = d->PORSTRUC[IJK];
        phib = d->PORSTRUC[IJKm1];
        zt   = p->ZSP[IJK];
        zb   = p->ZSP[IJKm1];
        zs   = d->WL(i,j) + d->bed(i,j);
        
        g = (phit-phib)/MAX(zt-zb,1.0e-20);
        g = MAX(-1.0, MIN(1.0, g));
        
        phis = phit + g*(zs-zt);
        
        psi = porpsi(p,d);
        
        Hs = 0.5*(1.0 - phis/psi + (1.0/PI)*sin((PI*(-phis))/psi));
        
        if(-phis > psi)
        Hs = 1.0;
        
        if(-phis < -psi)
        Hs = 0.0;
        
        Ht = Hporface(p,d,0,0,0);
        
        nl  = (p->B202>0) ? LN[IJK] : p->B201_n;
        nst = 1.0 - Ht*(1.0-nl);
        nss = 1.0 - Hs*(1.0-nl);
        
        d->PORS(i,j) = d->POR[IJK]*nss/MAX(nst,1.0e-6);
    }
    
    pgc->gcsl_start4(p,d->PORS,50);
}

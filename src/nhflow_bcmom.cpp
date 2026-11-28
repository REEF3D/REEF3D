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

#include"nhflow_bcmom.h"
#include"nhflow_wall.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"slide2D.h"
#include"ghostcell.h"
#include"turbulence.h"
#include"increment.h"

nhflow_bcmom::nhflow_bcmom(lexer* p):roughness(p),kappa(0.4)
{
}

nhflow_bcmom::~nhflow_bcmom()
{
}

// wall shear stress of one wall face on one velocity component, per unit volume, for the WL-weighted
// momentum equation: |u_t| u / (u+^2 dn) WL, with the log law at the cell centre (z0 = dn/2 >= ks/30)
static inline double nhflow_wall_drag(double vel, double ut, double dn, double ks, double WLval)
{
    const double kappa=0.4;
    double z0 = 0.5*dn;
    
    if(ks<=0.0)
    ks=0.0001;
    
    if(30.0*z0<ks)
    z0=ks/30.0;
    
    const double uplus = (1.0/kappa)*MAX(1.0,log(30.0*(z0/ks)));
    
    return ut*vel*WLval/(uplus*uplus*dn);
}

// The bed shear stress of the wall-function bed (A 518 1, A 519 1/2), implicit in time. The friction rate of the
// bottom cell is gamma = |u| / (u+^2 dn), d(UH)/dt = -gamma UH, linearised in the velocity of the stage input.
//
//  - with the implicit diffusion (A 512 2, or A 512 0 with a turbulence model): on the diagonal of the bed cell
//    of the diffusion system (nhflow_idiff), for the velocity of the stage X + dt F. The bed shear and the
//    vertical diffusion balance each other in the bed cell and belong into one system: as a source in F next to
//    a diffusion solved on its own, the steady state had the splitting error dt D(friction), of the order
//    nu_t dt/dz^2 of the bed shear (uniform open-channel flow, k-epsilon, h 1 m, dt 1 s: the velocity of the bed
//    cell 12 % above the log law with 10 layers and 29 % with 20, the normal depth 8 % and 17 % too small).
//  - without it (A 512 0 and no turbulence model, A 512 1): after the stage update (bed_friction_update),
//    UH_new (1 + fac dt gamma) = UH_stage, as the Darcy-Forchheimer drag (vrans_nhflow_f::implicit_drag). The
//    former drag in F, c |u| u WL / (1 + dt c |u|), was the explicit drag divided by 1 + dt gamma: stable in a
//    thin layer at a wet/dry front, but it left the steady bed shear too small by that factor.
bool nhflow_bcmom::bed_implicit(lexer *p)
{
    return p->A519>=1 && p->A518!=2 && (p->A512==2 || (p->A512==0 && p->A560>0));
}

double nhflow_bcmom::bed_gamma(lexer *p, fdm_nhf *d, int i, int j, int k, double WLval)
{
    const double kappa=0.4;
    const int marge = increment::marge;   // for KP
    const double dn = p->DZN[KP]*WLval;
    double ks = d->ks(i,j);
    
    if(dn<1.0e-12)
    return 0.0;
    
    double z0 = 0.5*dn;
    
    if(ks<=0.0)
    ks=0.0001;
    
    if(30.0*z0<ks)
    z0=ks/30.0;
    
    const double uplus = (1.0/kappa)*MAX(1.0,log(30.0*(z0/ks)));
    const double U = d->U[IJK];
    const double V = d->V[IJK];
    
    return sqrt(U*U + V*V)/(uplus*uplus*dn);
}

// bed friction after the stage update, where the diffusion system does not carry it; fac: stage weight of the
// explicit terms
void nhflow_bcmom::bed_friction_update(lexer *p, fdm_nhf *d, double *UH, slice &WL, double fac, int comp)
{
    if(p->A519<1 || p->A518==2 || comp>1 || bed_implicit(p))
    return;
    
    k=0;
    
    SLICELOOP4
    if(p->DF[IJK]>0 && p->wet[IJ]==1)
    UH[IJK] /= (1.0 + fac*MAX(p->dt,0.0)*bed_gamma(p,d,i,j,k,WL(i,j)));
}

// comp 0: U into F, 1: V into G, 2: W into H
// A519 1: bed friction in the bottom cell layer on u and v: implicit, not in F (bed_gamma: in the diffusion
//         system or after the stage update).
// A519 2: friction at every wall face (nhflow_wall.h, the same walls as the turbulence wall functions),
//         on the velocity components tangential to the face; ks = B57 at side and top walls, the bed
//         roughness d->ks at the bed. Faces add up (corner cells feel both walls).
// The bed friction belongs to the wall-function bed (A 518 1). The resolved no-slip bed (A 518 2) carries the
// bed shear in the diffusion (u = 0 on the bed face): no bed friction on top of it.
void nhflow_bcmom::wall_friction(lexer *p, fdm_nhf *d, int comp, double *VEL, double *F, slice &WL)
{
    double U,V,W,dz;
    int w[6],n;
    
    const bool bedfric = p->A518!=2;
    
    if(p->A519==2)
    LOOP
    if(p->DF[IJK]>0)
    {
    nhflow_wall_faces(p,i,j,k,w);
    
    U=d->U[IJK];
    V=d->V[IJK];
    W=d->W[IJK];
    dz = p->DZN[KP]*WL(i,j);
    
        // x walls: v, w
        n = w[0]+w[1];
        if(n>0 && comp!=0)
        F[IJK] -= n*nhflow_wall_drag(VEL[IJK], sqrt(V*V + W*W), p->DXN[IP], p->B57, WL(i,j));
        
        // y walls: u, w
        n = w[2]+w[3];
        if(n>0 && comp!=1)
        F[IJK] -= n*nhflow_wall_drag(VEL[IJK], sqrt(U*U + W*W), p->DYN[JP], p->B57, WL(i,j));
        
        // bed: u, v. The bottom layer (k = 0) is implicit (bed_gamma); the bed of a cell above a solid here
        if(w[4] && comp!=2 && bedfric && k>0)
        F[IJK] -= nhflow_wall_drag(VEL[IJK], sqrt(U*U + V*V), dz, d->ks(i,j), WL(i,j));
        
        // solid above: u, v
        if(w[5] && comp!=2)
        F[IJK] -= nhflow_wall_drag(VEL[IJK], sqrt(U*U + V*V), dz, p->B57, WL(i,j));
    }
    
    // sharp thin bodies (X 330 'mobility sharp') with a turbulence model: wall friction at the faces at the body
    // (the turbulence model sees them as walls with a wall function, nhflow_wall.h)
    if(d->thinbody!=nullptr && p->A560>0)
    LOOP
    if(p->DF[IJK]>0 && p->wet[IJ]==1)
    {
    d->thinbody->wall_faces(p,i,j,k,w);
    
    U=d->U[IJK];
    V=d->V[IJK];
    W=d->W[IJK];
    dz = p->DZN[KP]*WL(i,j);
    
        n = w[0]+w[1];
        if(n>0 && comp!=0)
        F[IJK] -= n*nhflow_wall_drag(VEL[IJK], sqrt(V*V + W*W), p->DXN[IP], p->B57, WL(i,j));
        
        n = w[2]+w[3];
        if(n>0 && comp!=1)
        F[IJK] -= n*nhflow_wall_drag(VEL[IJK], sqrt(U*U + W*W), p->DYN[JP], p->B57, WL(i,j));
        
        n = w[4]+w[5];
        if(n>0 && comp!=2)
        F[IJK] -= n*nhflow_wall_drag(VEL[IJK], sqrt(U*U + V*V), dz, p->B57, WL(i,j));
    }
}

void nhflow_bcmom::roughness_u(lexer* p, fdm_nhf *d, double *U, double *F, slice &WL)
{
    wall_friction(p,d,0,U,F,WL);
    
    // the slide layer (W 150, W 155): the drag of the slide on the bottom water layer
    if(d->slide!=nullptr)
    d->slide->drag_nhflow(p,d,0,U,F,WL);
}

void nhflow_bcmom::roughness_v(lexer* p, fdm_nhf *d, double *V, double *G, slice &WL)
{
    wall_friction(p,d,1,V,G,WL);
    
    if(d->slide!=nullptr)
    d->slide->drag_nhflow(p,d,1,V,G,WL);
}

void nhflow_bcmom::roughness_w(lexer* p, fdm_nhf *d, double *W, double *H, slice &WL)
{
    wall_friction(p,d,2,W,H,WL);
}

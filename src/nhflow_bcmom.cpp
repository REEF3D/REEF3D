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
#include"ghostcell.h"
#include"turbulence.h"

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

// comp 0: U into F, 1: V into G, 2: W into H
// A519 1: bed friction in the bottom cell layer on u and v.
// A519 2: friction at every wall face (nhflow_wall.h, the same walls as the turbulence wall functions),
//         on the velocity components tangential to the face; ks = B57 at side and top walls, the bed
//         roughness d->ks at the bed. Faces add up (corner cells feel both walls).
void nhflow_bcmom::wall_friction(lexer *p, fdm_nhf *d, int comp, double *VEL, double *F, slice &WL)
{
    double U,V,W,dz;
    int w[6],n;
    
    if(p->A519==1 && comp<2)
    {
    k=0;
    
    SLICELOOP4
    if(p->DF[IJK]>0)
    {
    U=d->U[IJK];
    V=d->V[IJK];
    
    F[IJK] -= nhflow_wall_drag(VEL[IJK], sqrt(U*U + V*V), p->DZN[KP]*WL(i,j), d->ks(i,j), WL(i,j));
    }
    }
    
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
        
        // bed: u, v
        if(w[4] && comp!=2)
        F[IJK] -= nhflow_wall_drag(VEL[IJK], sqrt(U*U + V*V), dz, d->ks(i,j), WL(i,j));
        
        // solid above: u, v
        if(w[5] && comp!=2)
        F[IJK] -= nhflow_wall_drag(VEL[IJK], sqrt(U*U + V*V), dz, p->B57, WL(i,j));
    }
}

void nhflow_bcmom::roughness_u(lexer* p, fdm_nhf *d, double *U, double *F, slice &WL)
{
    wall_friction(p,d,0,U,F,WL);
}

void nhflow_bcmom::roughness_v(lexer* p, fdm_nhf *d, double *V, double *G, slice &WL)
{
    wall_friction(p,d,1,V,G,WL);
}

void nhflow_bcmom::roughness_w(lexer* p, fdm_nhf *d, double *W, double *H, slice &WL)
{
    wall_friction(p,d,2,W,H,WL);
}

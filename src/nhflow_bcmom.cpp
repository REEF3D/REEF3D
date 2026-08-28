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

void nhflow_bcmom::roughness_u(lexer* p, fdm_nhf *d, double *U, double *F, slice &WL)
{
    if(p->A519==1)
    {
	k=0;
    
    SLICELOOP4
    if(p->DF[IJK]>0)
    {
    deltaZ = p->DZN[KP]*WL(i,j);
    z0 = 0.5*deltaZ;
	
	ks=d->ks(i,j);


		if(30.0*z0<ks)
		z0=ks/30.0;

		uplus = (1.0/kappa)*log(30.0*(z0/ks));

	if(fabs(z0*uplus)>1.0e-10)
	F[IJK] -= (fabs(U[IJK])*U[IJK]*WL(i,j))/(uplus*uplus*deltaZ);
    }
    }
    
    int check;
    
    if(p->A519==2)
    LOOP
    {
        if(p->DF[IJK]>0)
        {
        deltaZ = p->DZN[KP]*WL(i,j);
        
        
        check=0;
            
            if(p->DF[IJK]>0)
            {
                
            if((p->flag4[Im1JK]<0 || p->DF[Im1JK]<0) && i+p->origin_i != 0)
            {
            deltaZ = p->DXN[IP];
            ks=p->B57;
            check=1;
            }

            if((p->flag4[Ip1JK]<0 || p->DF[Ip1JK]<0) && i+p->origin_i != p->gknox-1)
            {
            deltaZ = p->DXN[IP];
            ks=p->B57;
            check=1;
            }

            if((p->flag4[IJm1K]<0 || p->DF[IJm1K]<0) && p->j_dir==1)
            {
            deltaZ = p->DYN[JP];
            ks=p->B57;
            check=1;
            }
                
            if((p->flag4[IJp1K]<0 || p->DF[IJp1K]<0) && p->j_dir==1)
            {
            deltaZ = p->DYN[JP];
            ks=p->B57;
            check=1;
            }
                
            if(p->flag4[IJKm1]<0 || p->DF[IJKm1]<0 || k==0)
            {
            deltaZ = p->DZN[KP]*d->WL(i,j);
            ks=p->B50;
            check=1;
            }

            if((p->flag4[IJKp1]<0 || p->DF[IJKp1]<0) && k!=p->knoz-1)
            {
            deltaZ = p->DZN[KP]*d->WL(i,j);
            ks=p->B57;
            check=1;
            }
        
                if(check==1)
                {
                z0 = 0.5*deltaZ;
                ks=d->ks(i,j);
                
                if(k==0 && p->S10>0)
                ks = p->S20*p->S21;


                if(30.0*z0<ks)
                z0=ks/30.0;

                uplus = (1.0/kappa)*log(30.0*(z0/ks));

                if(fabs(z0*uplus)>1.0e-10)
                F[IJK] -= (fabs(U[IJK])*U[IJK]*WL(i,j))/(uplus*uplus*deltaZ);
                }
            }
        }
    }
}

void nhflow_bcmom::roughness_v(lexer* p, fdm_nhf *d, double *V, double *G, slice &WL)
{
    if(p->A519==1)
    {
    k=0;
    
    SLICELOOP4
    if(p->DF[IJK]>0)
    {
    deltaZ = p->DZN[KP]*WL(i,j);
    z0 = 0.5*deltaZ;
	
	ks=d->ks(i,j);


		if(30.0*z0<ks)
		z0=ks/30.0;

		uplus = (1.0/kappa)*log(30.0*(z0/ks));

	if(fabs(z0*uplus)>1.0e-10)
	G[IJK] -= (fabs(V[IJK])*V[IJK]*WL(i,j))/(uplus*uplus*deltaZ);
    }
    }
	
    
    int check;
    
    if(p->A519==2)
    LOOP
    {
        if(p->DF[IJK]>0)
        {
        deltaZ = p->DZN[KP]*WL(i,j);
        
        
        check=0;
            
            if(p->DF[IJK]>0)
            {
                
            if((p->flag4[Im1JK]<0 || p->DF[Im1JK]<0) && i+p->origin_i != 0)
            {
            deltaZ = p->DXN[IP];
            ks=p->B57;
            check=1;
            }

            if((p->flag4[Ip1JK]<0 || p->DF[Ip1JK]<0) && i+p->origin_i != p->gknox-1)
            {
            deltaZ = p->DXN[IP];
            ks=p->B57;
            check=1;
            }

            if((p->flag4[IJm1K]<0 || p->DF[IJm1K]<0) && p->j_dir==1)
            {
            deltaZ = p->DYN[JP];
            ks=p->B57;
            check=1;
            }
                
            if((p->flag4[IJp1K]<0 || p->DF[IJp1K]<0) && p->j_dir==1)
            {
            deltaZ = p->DYN[JP];
            ks=p->B57;
            check=1;
            }
                
            if(p->flag4[IJKm1]<0 || p->DF[IJKm1]<0 || k==0)
            {
            deltaZ = p->DZN[KP]*d->WL(i,j);
            ks=p->B50;
            check=1;
            }

            if((p->flag4[IJKp1]<0 || p->DF[IJKp1]<0) && k!=p->knoz-1)
            {
            deltaZ = p->DZN[KP]*d->WL(i,j);
            ks=p->B57;
            check=1;
            }
        
                if(check==1)
                {
                z0 = 0.5*deltaZ;
                ks=d->ks(i,j);
                
                if(k==0 && p->S10>0)
                ks = p->S20*p->S21;


                if(30.0*z0<ks)
                z0=ks/30.0;

                uplus = (1.0/kappa)*log(30.0*(z0/ks));

                if(fabs(z0*uplus)>1.0e-10)
                G[IJK] -= (fabs(V[IJK])*V[IJK]*WL(i,j))/(uplus*uplus*deltaZ);
                }
            }
        }
    }
}

void nhflow_bcmom::roughness_w(lexer* p, fdm_nhf *d, double *W, double *H, slice &WL)
{
    int check;
    
    if(p->A519==2)
    LOOP
    {
        if(p->DF[IJK]>0)
        {
        deltaZ = p->DZN[KP]*WL(i,j);
        
        
        check=0;
            
            if(p->DF[IJK]>0)
            {
                
            if((p->flag4[Im1JK]<0 || p->DF[Im1JK]<0) && i+p->origin_i != 0)
            {
            deltaZ = p->DXN[IP];
            ks=p->B57;
            check=1;
            }

            if((p->flag4[Ip1JK]<0 || p->DF[Ip1JK]<0) && i+p->origin_i != p->gknox-1)
            {
            deltaZ = p->DXN[IP];
            ks=p->B57;
            check=1;
            }

            if((p->flag4[IJm1K]<0 || p->DF[IJm1K]<0) && p->j_dir==1)
            {
            deltaZ = p->DYN[JP];
            ks=p->B57;
            check=1;
            }
                
            if((p->flag4[IJp1K]<0 || p->DF[IJp1K]<0) && p->j_dir==1)
            {
            deltaZ = p->DYN[JP];
            ks=p->B57;
            check=1;
            }
                
            if(p->flag4[IJKm1]<0 || p->DF[IJKm1]<0 || k==0)
            {
            deltaZ = p->DZN[KP]*d->WL(i,j);
            ks=p->B50;
            check=1;
            }

            if((p->flag4[IJKp1]<0 || p->DF[IJKp1]<0) && k!=p->knoz-1)
            {
            deltaZ = p->DZN[KP]*d->WL(i,j);
            ks=p->B57;
            check=1;
            }
        
                if(check==1)
                {
                z0 = 0.5*deltaZ;
                ks=d->ks(i,j);
                
                if(k==0 && p->S10>0)
                ks = p->S20*p->S21;


                if(30.0*z0<ks)
                z0=ks/30.0;

                uplus = (1.0/kappa)*log(30.0*(z0/ks));

                if(fabs(z0*uplus)>1.0e-10)
                H[IJK] -= (fabs(W[IJK])*W[IJK]*WL(i,j))/(uplus*uplus*deltaZ);
                }
            }
        }
    }
}





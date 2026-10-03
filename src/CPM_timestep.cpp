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
Author: Alexander Hanke
--------------------------------------------------------------------*/

#include"CPM.h"
#include"lexer.h"
#include"ghostcell.h"

void CPM::timestep(lexer *p, ghostcell *pgc)
{
    double maxVelU=0.0, maxVelV=0.0, maxVelW=0.0;
    double maxvz=0.0;

    for(size_t n=0;n<P.index;n++)
    if(P.Flag[n]==ACTIVE)
    {
        maxVelU = MAX(maxVelU,fabs(P.U[n]));
        maxVelV = MAX(maxVelV,fabs(P.V[n]));
        maxVelW = MAX(maxVelW,fabs(P.W[n]));
    }

    maxVelU = pgc->globalmax(maxVelU);
    maxVelV = pgc->globalmax(maxVelV);
    maxVelW = pgc->globalmax(maxVelW);
    
    maxvz = MAX(maxVelU,maxVelV);
    maxvz = MAX(maxvz,maxVelW);
    
    // MP-PIC: synchronised with the fluid, adaptive sub-steps in CPM::substep_size
    if(p->Q11==2)
    p->dtsed = p->dt;
    
    // grid-limited step: rejected moves in the previous step and largest occupancy
    int nrej = pgc->globalisum(nrej_step);
    int nclip = pgc->globalisum(nclip_step);
    int nit = pgc->globalimax(nit_step);
    nrej_step = 0;
    nclip_step = 0;
    nit_step = 0;
    double omax = (p->Q11==2 && p->Q19==1) ? occupancy_max(p,pgc) : 0.0;
    
    if(p->Q11==1)
    {
        if(timestep_ini<5)
        {
            maxvz = 1000.0;
            ++timestep_ini;
        }

        if(p->S15==0)
            p->dtsed=MIN(p->S13, (p->S14*p->DXM)/(fabs(maxvz)>1.0e-15?maxvz:1.0e-15));
        else if(p->S15==1)
            p->dtsed=MIN(p->dt, (p->S14*p->DXM)/(fabs(maxvz)>1.0e-15?maxvz:1.0e-15));
        else if(p->S15==2)
            p->dtsed=p->S13;

        p->dtsed=pgc->timesync(p->dtsed);
    }

    p->sedtime+=p->dtsed;

    if(p->mpirank==0)
    {
        cout<<"Sediment Iter: "<<p->sediter<<"  Sediment Time: "<<setprecision(4)<<p->sedtime<<"  Sediment Timestep: "<<setprecision(4)<<p->dtsed<<endl;
        
        if(p->Q11==2)
        cout<<"CPM sub-steps (previous step): "<<nsub<<"  mean dt_sub: "<<setprecision(4)<<dtsub<<"  stress wave speed: "<<setprecision(4)<<cmax<<endl;
        
        if(p->Q11==2 && p->Q19==1)
        cout<<"CPM grid-limited step: rejected moves "<<nrej<<"  shortened moves "<<nclip<<"  max occupancy "<<setprecision(4)<<omax<<" (capacity "<<MAX(theta_max, theta_0 + P.ParcelFactor*Vp/(p->DXM*p->DXM*(p->j_dir==1?p->DXM:p->DYN[marge])))<<")  fix-up passes "<<nit<<endl;
        
        cout<<"Up_max: "<<setprecision(4)<<maxVelU<<"  Vp_max: "<<maxVelV<<"  Wp_max: "<<maxVelW<<endl;
        cout<<defaultfloat;
    }
}

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
Architect: Hans Bihs
--------------------------------------------------------------------*/

#include"6DOF_obj_nhflow.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"

void sixdof_fluid_nhflow::velocity(int n, const double *xyz, double *uvw)
{
    // each point is interpolated by the rank that owns it (x, y of the subdomain; NHFLOW is
    // decomposed horizontally), the others contribute 0
    for(int q=0; q<n; ++q)
    {
        const double x = xyz[3*q], y = xyz[3*q+1], z = xyz[3*q+2];
        
        uvw[3*q] = uvw[3*q+1] = uvw[3*q+2] = 0.0;
        
        const bool own = (x>=p->originx && x<p->endx) && (p->j_dir==0 || (y>=p->originy && y<p->endy));
        
        if(own)
        {
            uvw[3*q]   = p->ccipol4V(d->U, d->WL, d->bed, x, y, z);
            uvw[3*q+1] = p->ccipol4V(d->V, d->WL, d->bed, x, y, z);
            uvw[3*q+2] = p->ccipol4V(d->W, d->WL, d->bed, x, y, z);
        }
    }
    
    for(int q=0; q<3*n; ++q)
    uvw[q] = pgc->globalsum(uvw[q]);
}

void sixdof_obj_nhflow::actuator_source(lexer *p, fdm_nhf *d, ghostcell *pgc, slice &WL, int comp, double *F)
{
    if(pload.empty())
    return;
    
    vector<sixdof_actuator_disk> disks;
    
    for(size_t ql=0; ql<pload.size(); ++ql)
    pload[ql]->actuator_disks(disks);
    
    for(size_t qd=0; qd<disks.size(); ++qd)
    {
        const sixdof_actuator_disk &ad = disks[qd];
        double wa, wt, r;
        Eigen::Vector3d et;
        
        // discrete integrals of the Hough-Ordway weights over the cells of the disk: the
        // force and the torque given to the fluid are exactly T and Q on this grid
        // only the fluid part of a cell gets the source (1 - FHB, FHB: Heaviside of the bodies of the
        // direct forcing), so no thrust is put into the forcing zone of the hull, where the forcing
        // would take it out again and hand it to the hull as a resistance
        double SA=0.0, ST=0.0, SA0=0.0;
        
        LOOP
        {
            const Eigen::Vector3d x(p->XP[IP], p->YP[JP], p->ZSP[IJK]);
            
            if(ad.weights(x,wa,wt,et,r))
            {
                const double V = p->DXN[IP]*p->DYN[JP]*p->DZN[KP]*WL(i,j);
                const double fl = 1.0 - MIN(MAX(d->FHB[IJK],0.0),1.0);
                SA0 += wa*V;
                SA += fl*wa*V;
                ST += fl*wt*r*V;
            }
        }
        
        SA0 = pgc->globalsum(SA0);
        SA = pgc->globalsum(SA);
        ST = pgc->globalsum(ST);
        
        if(SA0>0.0 && SA<0.95*SA0 && !actuator_warned)
        {
            actuator_warned = true;
            
            if(p->mpirank==0)
            cout<<"6DOF: "<<int(100.0*(1.0-SA/SA0)+0.5)<<" % of the actuator disk lie in the direct forcing zone of a body;"
                <<" the source acts on the fluid part only, place the propeller further from the hull"<<endl;
        }
        
        if(SA<=0.0)
        {
            if(p->mpirank==0 && comp==0)
            cout<<"6DOF: actuator disk without fluid cells (centre "<<ad.centre.transpose()<<"), no momentum source"<<endl;
            
            continue;
        }
        
        LOOP
        {
            const Eigen::Vector3d x(p->XP[IP], p->YP[JP], p->ZSP[IJK]);
            
            if(ad.weights(x,wa,wt,et,r))
            {
                // force density on the fluid [N/m^3]: thrust reaction along -axis, swirl along et
                const double fl = 1.0 - MIN(MAX(d->FHB[IJK],0.0),1.0);
                Eigen::Vector3d f = -fl*ad.T*wa/SA*ad.axis;
                
                if(ST>0.0)
                f += fl*ad.Q*wt/ST*et;
                
                F[IJK] += WL(i,j)*f(comp)/p->W1;
            }
        }
    }
}

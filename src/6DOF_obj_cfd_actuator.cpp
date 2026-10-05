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

#include"6DOF_obj_cfd.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

void sixdof_fluid_cfd::velocity(int n, const double *xyz, double *uvw)
{
    // each point is interpolated by the rank that owns it (CFD is decomposed in x, y and z),
    // the others contribute 0
    for(int q=0; q<n; ++q)
    {
        const double x = xyz[3*q], y = xyz[3*q+1], z = xyz[3*q+2];
        
        uvw[3*q] = uvw[3*q+1] = uvw[3*q+2] = 0.0;
        
        const bool own = (x>=p->originx && x<p->endx) && (p->j_dir==0 || (y>=p->originy && y<p->endy))
                      && (z>=p->originz && z<p->endz);
        
        if(own)
        {
            uvw[3*q]   = p->ccipol1(a->u, x, y, z);
            uvw[3*q+1] = p->ccipol2(a->v, x, y, z);
            uvw[3*q+2] = p->ccipol3(a->w, x, y, z);
        }
    }
    
    for(int q=0; q<3*n; ++q)
    uvw[q] = pgc->globalsum(uvw[q]);
}

// share of the momentum source a velocity point (i,j,k) of component comp gets: the part outside
// the bodies of the direct forcing (fbh: Heaviside of the bodies), times the water part from the
// face density; so no thrust is put into the forcing zone of the hull, where the forcing would take
// it out again and hand it to the hull as a resistance, and none into air (a disk that breaks the
// surface pushes only the water part)
double sixdof_obj_cfd::actuator_fluid_fraction(lexer *p, fdm *a, int comp)
{
    double h, rof;
    
    if(comp==0)
    {
        h = a->fbh1(i,j,k);
        rof = 0.5*(a->ro(i,j,k) + a->ro(i+1,j,k));
    }
    else if(comp==1)
    {
        h = a->fbh2(i,j,k);
        rof = 0.5*(a->ro(i,j,k) + a->ro(i,j+1,k));
    }
    else
    {
        h = a->fbh3(i,j,k);
        rof = 0.5*(a->ro(i,j,k) + a->ro(i,j,k+1));
    }
    
    double water = 1.0;
    
    if(fabs(p->W1-p->W3)>1.0e-10)
    water = MIN(MAX((rof - p->W3)/(p->W1 - p->W3),0.0),1.0);
    
    return (1.0 - MIN(MAX(h,0.0),1.0))*water;
}

// source point and volume of velocity point (i,j,k) of component comp (staggered grid)
void sixdof_obj_cfd::actuator_point(lexer *p, int comp, Eigen::Vector3d &x, double &V)
{
    if(comp==0)
    {
        x << p->pos1_x(), p->pos1_y(), p->pos1_z();
        V = p->DXP[IP]*p->DYN[JP]*p->DZN[KP];
    }
    else if(comp==1)
    {
        x << p->pos2_x(), p->pos2_y(), p->pos2_z();
        V = p->DXN[IP]*p->DYP[JP]*p->DZN[KP];
    }
    else
    {
        x << p->pos3_x(), p->pos3_y(), p->pos3_z();
        V = p->DXN[IP]*p->DYN[JP]*p->DZP[KP];
    }
}

void sixdof_obj_cfd::actuator_forcing(lexer *p, fdm *a, ghostcell *pgc, field &fx, field &fy, field &fz)
{
    if(pload.empty())
    return;
    
    vector<sixdof_actuator_disk> disks;
    
    for(size_t ql=0; ql<pload.size(); ++ql)
    pload[ql]->actuator_disks(disks);
    
    Eigen::Vector3d x;
    double V;
    
    for(size_t qd=0; qd<disks.size(); ++qd)
    {
        const sixdof_actuator_disk &ad = disks[qd];
        
        // discrete sums of the distribution, each component on its own grid (staggered velocity
        // points), so that the force given to the fluid is exactly T along the axis and the torque
        // exactly Q on this grid
        sixdof_actuator_disk::sums S[3];
        
        ULOOP
        {
            actuator_point(p,0,x,V);
            ad.accumulate(x,V,actuator_fluid_fraction(p,a,0),0,S[0]);
        }
        
        VLOOP
        {
            actuator_point(p,1,x,V);
            ad.accumulate(x,V,actuator_fluid_fraction(p,a,1),1,S[1]);
        }
        
        WLOOP
        {
            actuator_point(p,2,x,V);
            ad.accumulate(x,V,actuator_fluid_fraction(p,a,2),2,S[2]);
        }
        
        for(int c=0; c<3; ++c)
        {
            S[c].SA0 = pgc->globalsum(S[c].SA0);
            S[c].SA  = pgc->globalsum(S[c].SA);
            S[c].SW  = pgc->globalsum(S[c].SW);
            S[c].SE  = pgc->globalsum(S[c].SE);
            S[c].G1  = pgc->globalsum(S[c].G1);
            S[c].G0  = pgc->globalsum(S[c].G0);
            S[c].GA  = pgc->globalsum(S[c].GA);
        }
        
        const double SA0 = S[0].SA0, SA = S[0].SA;
        
        if(SA0>0.0 && SA<0.95*SA0 && !actuator_warned)
        {
            actuator_warned = true;
            
            if(p->mpirank==0)
            cout<<"6DOF: "<<int(100.0*(1.0-SA/SA0)+0.5)<<" % of the actuator disk lie in the direct forcing zone of a body or in air;"
                <<" the source acts on the water part only, place the propeller further from the hull and the free surface"<<endl;
        }
        
        if(SA<=0.0)
        {
            if(p->mpirank==0)
            cout<<"6DOF: actuator disk without water cells (centre "<<ad.centre.transpose()<<"), no momentum source"<<endl;
            
            continue;
        }
        
        const double kappa = ad.swirl_factor(S);
        
        // force density [N/m^3] divided by the face density: the forcing terms are accelerations
        ULOOP
        {
            actuator_point(p,0,x,V);
            fx(i,j,k) += ad.force(x,actuator_fluid_fraction(p,a,0),0,S[0],kappa)/MAX(0.5*(a->ro(i,j,k) + a->ro(i+1,j,k)),1.0e-10);
        }
        
        VLOOP
        {
            actuator_point(p,1,x,V);
            fy(i,j,k) += ad.force(x,actuator_fluid_fraction(p,a,1),1,S[1],kappa)/MAX(0.5*(a->ro(i,j,k) + a->ro(i,j+1,k)),1.0e-10);
        }
        
        WLOOP
        {
            actuator_point(p,2,x,V);
            fz(i,j,k) += ad.force(x,actuator_fluid_fraction(p,a,2),2,S[2],kappa)/MAX(0.5*(a->ro(i,j,k) + a->ro(i,j,k+1)),1.0e-10);
        }
    }
}

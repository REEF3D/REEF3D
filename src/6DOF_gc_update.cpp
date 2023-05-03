/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2022 Hans Bihs

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
Author: Hans Bihs, Tobias Martin
--------------------------------------------------------------------*/

#include"6DOF_gc.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"ioflow.h"
#include"freesurface.h"


void sixdof_gc::solidUpdate
(
	lexer *p,
	fdm* a, 
	ghostcell *pgc, 
	const std::vector<double>& ek
)
{
	// Print quaternion
	
	if(p->mpirank==0)
    {
		cout<<"Quaternion "<<ek[0]<<" "<<ek[1]<<" "<<ek[2]<<" "<<ek[3]<<endl;
    }
	
	
    // Update transformation matrix from Shivarama thesis p. 19
	
	R_[0][0] = ek[0]*ek[0] + ek[1]*ek[1] - ek[2]*ek[2] - ek[3]*ek[3]; 
	R_[0][1] = 2.0*ek[1]*ek[2] - 2.0*ek[0]*ek[3]; 
	R_[0][2] = 2.0*ek[0]*ek[2] + 2.0*ek[1]*ek[3];
	R_[1][0] = 2.0*ek[0]*ek[3] + 2.0*ek[1]*ek[2]; 
	R_[1][1] = ek[0]*ek[0] - ek[1]*ek[1] + ek[2]*ek[2] - ek[3]*ek[3];
	R_[1][2] = 2.0*ek[2]*ek[3] - 2.0*ek[0]*ek[1];
	R_[2][0] = 2.0*ek[1]*ek[3] - 2.0*ek[0]*ek[2]; 
	R_[2][1] = 2.0*ek[0]*ek[1] + 2.0*ek[2]*ek[3]; 
	R_[2][2] = ek[0]*ek[0] - ek[1]*ek[1] - ek[2]*ek[2] + ek[3]*ek[3];

    quatRotMat(0,0) = R_[0][0]; quatRotMat(0,1) = R_[0][1]; quatRotMat(0,2) = R_[0][2];
    quatRotMat(1,0) = R_[1][0]; quatRotMat(1,1) = R_[1][1]; quatRotMat(1,2) = R_[1][2];
    quatRotMat(2,0) = R_[2][0]; quatRotMat(2,1) = R_[2][1]; quatRotMat(2,2) = R_[2][2];


    // Update angular velocities in body coordinates
	
	Ps = 
		(ek[6]*(I_[0][1]*I_[1][2] - I_[0][2]*I_[1][1]))/(I_[0][0]*I_[1][1]*I_[2][2] - I_[0][0]*I_[1][2]*I_[2][1] - I_[0][1]*I_[1][0]*I_[2][2] + I_[0][1]*I_[1][2]*I_[2][0] + I_[0][2]*I_[1][0]*I_[2][1] - I_[0][2]*I_[1][1]*I_[2][0]) 
		- (ek[5]*(I_[0][1]*I_[2][2] - I_[0][2]*I_[2][1]))/(I_[0][0]*I_[1][1]*I_[2][2] - I_[0][0]*I_[1][2]*I_[2][1] - I_[0][1]*I_[1][0]*I_[2][2] + I_[0][1]*I_[1][2]*I_[2][0] + I_[0][2]*I_[1][0]*I_[2][1] - I_[0][2]*I_[1][1]*I_[2][0]) 
		+ (ek[4]*(I_[1][1]*I_[2][2] - I_[1][2]*I_[2][1]))/(I_[0][0]*I_[1][1]*I_[2][2] - I_[0][0]*I_[1][2]*I_[2][1] - I_[0][1]*I_[1][0]*I_[2][2] + I_[0][1]*I_[1][2]*I_[2][0] + I_[0][2]*I_[1][0]*I_[2][1] - I_[0][2]*I_[1][1]*I_[2][0]);
	Qs = 
		(ek[5]*(I_[0][0]*I_[2][2] - I_[0][2]*I_[2][0]))/(I_[0][0]*I_[1][1]*I_[2][2] - I_[0][0]*I_[1][2]*I_[2][1] - I_[0][1]*I_[1][0]*I_[2][2] + I_[0][1]*I_[1][2]*I_[2][0] + I_[0][2]*I_[1][0]*I_[2][1] - I_[0][2]*I_[1][1]*I_[2][0]) 
		- (ek[6]*(I_[0][0]*I_[1][2] - I_[0][2]*I_[1][0]))/(I_[0][0]*I_[1][1]*I_[2][2] - I_[0][0]*I_[1][2]*I_[2][1] - I_[0][1]*I_[1][0]*I_[2][2] + I_[0][1]*I_[1][2]*I_[2][0] + I_[0][2]*I_[1][0]*I_[2][1] - I_[0][2]*I_[1][1]*I_[2][0]) 
		- (ek[4]*(I_[1][0]*I_[2][2] - I_[1][2]*I_[2][0]))/(I_[0][0]*I_[1][1]*I_[2][2] - I_[0][0]*I_[1][2]*I_[2][1] - I_[0][1]*I_[1][0]*I_[2][2] + I_[0][1]*I_[1][2]*I_[2][0] + I_[0][2]*I_[1][0]*I_[2][1] - I_[0][2]*I_[1][1]*I_[2][0]);
		
	Rs = 
		(ek[6]*(I_[0][0]*I_[1][1] - I_[0][1]*I_[1][0]))/(I_[0][0]*I_[1][1]*I_[2][2] - I_[0][0]*I_[1][2]*I_[2][1] - I_[0][1]*I_[1][0]*I_[2][2] + I_[0][1]*I_[1][2]*I_[2][0] + I_[0][2]*I_[1][0]*I_[2][1] - I_[0][2]*I_[1][1]*I_[2][0]) 
		- (ek[5]*(I_[0][0]*I_[2][1] - I_[0][1]*I_[2][0]))/(I_[0][0]*I_[1][1]*I_[2][2] - I_[0][0]*I_[1][2]*I_[2][1] - I_[0][1]*I_[1][0]*I_[2][2] + I_[0][1]*I_[1][2]*I_[2][0] + I_[0][2]*I_[1][0]*I_[2][1] - I_[0][2]*I_[1][1]*I_[2][0]) 
		+ (ek[4]*(I_[1][0]*I_[2][1] - I_[1][1]*I_[2][0]))/(I_[0][0]*I_[1][1]*I_[2][2] - I_[0][0]*I_[1][2]*I_[2][1] - I_[0][1]*I_[1][0]*I_[2][2] + I_[0][1]*I_[1][2]*I_[2][0] + I_[0][2]*I_[1][0]*I_[2][1] - I_[0][2]*I_[1][1]*I_[2][0]);
	
	
	// Update linear velocities and position
	
	xg = ek[7];
	yg = ek[8];
	zg = ek[9];
    Ue = ek[10]/Mfb;
    Ve = ek[11]/Mfb;
    We = ek[12]/Mfb;
	
	
	// Shrink motions

	motion_ext_quaternion(p,a,pgc);
	
	
	// Update Euler angles and position of triangles
	
	fb_position_quaternion(p,a,pgc,ek);	

	// Update floating field
	
	ray_cast(p,a,pgc);
	reini_AB2(p,a,pgc,a->fb);
    pgc->start4a(p,a->fb,50);
	
	
	// Save body data for fluid solver
	
	transform_angle_SE(Ps,Qs,Rs,Pe,Qe,Re);		
	interface(p,true);
	maxvel(p,a,pgc);


	// Update ghostcells
	
	double starttime1=pgc->timer();
	pgc->gcfb_update(p,a);
	double endtime1 = pgc->timer()-starttime1;

	if (p->mpirank == 0)
	{
		cout<<"Ue: "<<p->ufbi<<" Ve: "<<p->vfbi<<" We: "<<p->wfbi<<" Pe: "<<p->pfbi<<" Qe: "<<p->qfbi<<" Re: "<<p->rfbi<<endl;
		cout<<"6DOF update time: "<<endtime1<<endl;	
	}
}

void sixdof_gc::forceUpdate(lexer *p, fdm* a, ghostcell *pgc, vrans *pvrans, vector<net*>& pnet)
{
    if (p->X12 == 1)
	{
		fluidForces(p,a,pgc);
	}

	if (p->X310 > 0)
	{
		mooringForces(p,a,pgc);
		
        for (int i=0; i<p->mooring_count; i++)
		{
			if(p->X11_u==1)
			Xe += Xme[i];
			if(p->X11_v==1)
			Ye += Yme[i];
			if(p->X11_w==1)
			Ze += Zme[i];
			
			if(p->X11_p==1)
			Ke += Kme[i];
			if(p->X11_q==1)
			Me += Mme[i];
			if(p->X11_r==1)
			Ne += Nme[i];
		}
	} 
	
	if (p->X320 > 0)
	{
		netForces(p,a,pgc,1.0,pvrans,pnet);
        
        for (int i=0; i<p->net_count; i++)
		{
			if(p->X11_u==1)
			Xe += Xne[i];
			if(p->X11_v==1)
			Ye += Yne[i];
			if(p->X11_w==1)
			Ze += Zne[i];
			
			if(p->X11_p==1)
			Ke += Kne[i];
			if(p->X11_q==1)
			Me += Mne[i];
			if(p->X11_r==1)
			Ne += Nne[i];
		}
	}
	
    preventMotion(p);

	// Transform to body coordinate system
	transform_vec_ES(Xe,Ye,Ze,Xs,Ys,Zs);
	transform_vec_ES(Ke,Me,Ne,Ks,Ms,Ns);	
}

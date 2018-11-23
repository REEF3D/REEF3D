/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2018 Hans Bihs

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
--------------------------------------------------------------------*/

#include"6DOF_f.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"momentum_FSI.h"
#include"ioflow.h"
#include"freesurface.h"


void sixdof_f::solidUpdate
(
	lexer *p,
	fdm* a, 
	ghostcell *pgc, 
	const std::vector<double>& ek,
	bool final
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
	interface(p,final);
	maxvel(p,a,pgc);


	// Update ghostcells
	
	double starttime1=pgc->timer();
	pgc->gcfb_update(p,a);
	double endtime1 = pgc->timer()-starttime1;


	if (final == true && p->mpirank == 0)
	{
		cout<<"Ue: "<<p->ufbi<<" Ve: "<<p->vfbi<<" We: "<<p->wfbi<<" Pe: "<<p->pfbi<<" Qe: "<<p->qfbi<<" Re: "<<p->rfbi<<endl;
		cout<<"6DOF update time: "<<endtime1<<endl;	
	}
}


void sixdof_f::fluidUpdate
(
	lexer *p,
	fdm *a, 
	ghostcell *pgc, 
	momentum *pmom,
	ioflow *pflow,
	freesurface *pfsf,
	convection *pfsfdisc,
	solver *psolv,
	reini *preini,
	particlecorr *ppart,
	bool final,
	int pcfInd_
)
{
	// Has to be used with N 40 0 and D 30 4
	
	
	pflow->u_relax(p,a,pgc,a->u);
	pflow->v_relax(p,a,pgc,a->v);
	pflow->w_relax(p,a,pgc,a->w);
	pfsf->update(p,a,pgc,a->phi);
	
	momentum_FSI::pcfInd = pcfInd_;
	pmom->start(p,a,pgc,pmom);
	
	// Update phi if not final loop												// ghostcells update necessary?
	if (final == false)
	{
		LOOP
		{
			a->phi(i,j,k) = phin(i,j,k);
		}
		
		pfsf->start(a,p, pfsfdisc,psolv,pgc,pflow,preini,ppart,a->phi);
	}
	else
	{
		LOOP
		{
			phin(i,j,k) = a->phi(i,j,k);
		}
	}
}


void sixdof_f::forceUpdate(lexer *p, fdm* a, ghostcell *pgc)
{
    if (p->X12 == 1)
	{
		fluidForces(p,a,pgc);
	}

	if (p->X310 > 0)
	{
		mooringForces(p,a,pgc);
	} 
	
	if (p->X320 > 0)
	{
		netForces(p,a,pgc);
	}
		
	// Add mooring forces
	if (p->X310 > 0)
	{
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
	
	preventMotion(p);
}
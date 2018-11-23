/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2018 Hans Bihs

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTIBILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
--------------------------------------------------------------------*/

#include"6DOF_f.h"
#include"net_QuasiStatic.h"
#include"net_void.h"
#include"mooring_void.h"
#include"mooring_DGSEM.h"
#include"mooring_QuasiStatic.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"reinidisc.h"
#include"reinidisc_f.h"
#include"reinidisc_fsf.h"
#include"momentum_FSI.h"

sixdof_f::sixdof_f
(
	lexer *p, 
	fdm *a, 
	ghostcell *pgc, 
	momentum *pmom,
	ioflow* pflow,
	freesurface* pfsf,
	convection* pfsfdisc,
	solver* psolv,
	reini* preini,
	particlecorr* ppart
) : gradient(p), cutl(p), cutr(p), epsifb(1.6*p->dx),f(p),dt(p),frk1(p),frk2(p),L(p),eta(p),phin(p),vertice(p),nodeflag(p)
{
	p->printcount_sixdof=0;
	
	prdisc = new reinidisc_fsf(p);
    
    Xe = 0.0;
	Ye = 0.0;
	Ze = 0.0;
	Ke = 0.0;
	Me = 0.0;
	Ne = 0.0;
	
	epsi=1.6*p->dx;
    zero=0.0;
	
	LOOP
	{
		phin(i,j,k) = a->phi(i,j,k);
	}
}

sixdof_f::~sixdof_f()
{
}

void sixdof_f::start
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
	particlecorr *ppart
)
{
	starttime=pgc->timer();
		
	// Main loop
	if (p->X13 == 0)
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
		
        solve(p,a,pgc);
        motion_ext(p,a,pgc);
        fb_position(p,a,pgc);
		
		ray_cast(p,a,pgc);
		reini_AB2(p,a,pgc,a->fb);
		pgc->start4a(p,a->fb,50);
		
		interface(p,true);
		maxvel(p,a,pgc);
		
		if(p->mpirank==0)
		cout<<"Ue: "<<p->ufbi<<" Ve: "<<p->vfbi<<" We: "<<p->wfbi<<" Pe: "<<p->pfbi<<" Qe: "<<p->qfbi<<" Re: "<<p->rfbi<<endl;
		
		double starttime1=pgc->timer();
		pgc->gcfb_update(p,a);
		double endtime1 = pgc->timer()-starttime1;
		
		print_stl(p,a,pgc);
		print_E_position(p,a,pgc);
		print_E_velocity(p,a,pgc);
		print_E_force(p,a,pgc);
		print_S_force(p,a,pgc);
		
		if(p->mpirank==0)
		cout<<"6DOF time: "<<setprecision(3)<<pgc->timer()-starttime<<"  update time: "<<endtime1<<endl;	
	}
    else if (p->X13 == 1)
    {
		// Initialise fluid
		//pgc->start4(p,a->press,40);
		forceUpdate(p,a,pgc);
		//fluidUpdate(p,a,pgc,pmom,pflow,pfsf,pfsfdisc,psolv,preini,ppart,false,0);
		
		// Initialise forces

		
		// FSI loop
		solve_quaternion(p,a,pgc,pmom,pflow,pfsf,pfsfdisc,psolv,preini,ppart);

		// Finalise solid 
		solidUpdate(p,a,pgc,e_,true);
		//forceUpdate(p,a,pgc);
		
		// Finalise fluid 
		//fluidUpdate(p,a,pgc,pmom,pflow,pfsf,pfsfdisc,psolv,preini,ppart,true,2);		
		
		// Print
		print_stl(p,a,pgc);
		print_E_position(p,a,pgc);
		print_E_velocity(p,a,pgc);
		print_E_force(p,a,pgc);
		
		if(p->mpirank==0)
		cout<<"6DOF time: "<<setprecision(3)<<pgc->timer()-starttime<<endl;
    }   
}
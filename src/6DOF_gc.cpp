/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2023 Hans Bihs

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
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"6DOF_gc.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"reinidisc.h"
#include"reinidisc_f.h"
#include"reinidisc_fsf.h"

sixdof_gc::sixdof_gc(
	lexer *p,
	fdm *a,
	ghostcell *pgc
) : gradient(p), cutl(p), cutr(p), fbio(p), epsifb(1.6*p->DXM), epsi(1.6),f(p),dt(p),frk1(p),frk2(p),L(p),eta(p),phin(p),vertice(p),nodeflag(p)
{
    if(p->mpirank==0)
    cout<<"6DOF_gc startup ..."<<endl;
	p->printcount_sixdof=0;

	prdisc = new reinidisc_fsf(p);

    Xe = 0.0;
	Ye = 0.0;
	Ze = 0.0;
	Ke = 0.0;
	Me = 0.0;
	Ne = 0.0;
    
    zero=0.0;
}

sixdof_gc::~sixdof_gc()
{
}


void sixdof_gc::start
(
	lexer *p, 
	fdm *a, 
	ghostcell *pgc,
    double alpha,
    vrans *pvrans,
    vector<net*>& pnet
)
{
	// Main loop
	if (p->X13 == 0)
	{
        start_Euler(p,a,pgc,pvrans,pnet);
    }
    else if (p->X13 == 1)
    {
        start_Quaternion(p,a,pgc,pvrans,pnet);
    }
}

void sixdof_gc::start_Euler
(
	lexer *p, 
	fdm *a, 
	ghostcell *pgc, 
    vrans *pvrans,
    vector<net*>& pnet
)
{
	starttime=pgc->timer();

    forceUpdate(p,a,pgc,pvrans,pnet);

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

    if(p->X50==1)
    print_vtp(p,a,pgc);
    
    if(p->X50==2)
    print_stl(p,a,pgc);
    
    print_E_position(p,a,pgc);
    print_E_velocity(p,a,pgc);
    print_E_force(p,a,pgc);
    print_S_force(p,a,pgc);

    if(p->mpirank==0)
    cout<<"6DOF time: "<<setprecision(3)<<pgc->timer()-starttime<<"  update time: "<<endtime1<<endl;	
}

void sixdof_gc::start_Quaternion
(
	lexer *p, 
	fdm *a, 
	ghostcell *pgc, 
    vrans *pvrans,
    vector<net*>& pnet
)
{
	starttime=pgc->timer();
    
    forceUpdate(p,a,pgc,pvrans,pnet);
    
    solve_quaternion();

    solidUpdate(p,a,pgc,e_);
    
    if(p->X50==1)
    print_vtp(p,a,pgc);
    
    if(p->X50==2)
    print_stl(p,a,pgc);
    
    print_E_position(p,a,pgc);
    print_E_velocity(p,a,pgc);
    print_E_force(p,a,pgc);
    
    if(p->mpirank==0)
    cout<<"6DOF time: "<<setprecision(3)<<pgc->timer()-starttime<<endl;
}

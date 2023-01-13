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

void sixdof_gc::motion_ext(lexer *p, fdm *a, ghostcell *pgc)
{
    if(p->X210==1 || p->X211==1)
	motion_fixed(p,a,pgc);
    
    if(p->X221==1)
	motion_vec(p,a,pgc);
    
    if(p->X210==1 || p->X211==1 || p->X221==1)
    {
	if(p->X11_u==2)
	dxg = p->dt*Uext;
	
	if(p->X11_v==2)
	dyg = p->dt*Vext;

	if(p->X11_w==2)
	dzg = p->dt*Wext;
	
	if(p->X11_p==2)
	dphi = p->dt*Pext;
	
	if(p->X11_q==2)
	dtheta = p->dt*Qext;
	
	if(p->X11_r==2)
	dpsi = p->dt*Rext;
    }
	
}


void sixdof_gc::motion_ext_quaternion(lexer *p, fdm *a, ghostcell *pgc)
{
    if(p->X210==1 || p->X211==1)
	motion_fixed(p,a,pgc);
    
    if(p->X221==1)
	motion_vec(p,a,pgc);
    
    if(p->X210==1 || p->X211==1 || p->X221==1)
    {
		if(p->X11_u==2)
		dxg = p->dt*Uext;
		
		if(p->X11_v==2)
		dyg = p->dt*Vext;

		if(p->X11_w==2)
		dzg = p->dt*Wext;
		
		if(p->X11_p==2)
		cout<<"not implemented yet"<<endl;
		
		if(p->X11_q==2)
		cout<<"not implemented yet"<<endl;
		
		if(p->X11_r==2)
		cout<<"not implemented yet"<<endl;
    }
	
}


void sixdof_gc::preventMotion(lexer *p)
{
	if(p->X11_u != 1) Xe = 0.0;
	if(p->X11_v != 1) Ye = 0.0;
	if(p->X11_w != 1) Ze = 0.0;
	if(p->X11_p != 1) Ke = 0.0;
	if(p->X11_q != 1) Me = 0.0;
	if(p->X11_r != 1) Ne = 0.0;
}


void sixdof_gc::motion_fixed(lexer *p, fdm *a, ghostcell *pgc)
{
	if(p->X210==1)
	{
	Uext = p->X210_u;
	Vext = p->X210_v;
	Wext = p->X210_w;
	}
	
	if(p->X211==1)
	{
	Pext = p->X211_p;
	Qext = p->X211_q;
	Rext = p->X211_r;
	}
}


void sixdof_gc::motion_vec(lexer *p, fdm* a, ghostcell* pgc)
{
    Uext=Vext=Wext=0.0;
    
    if(p->simtime>=ts_motion && p->simtime<=te_motion)
    {
	
    double vel = (motion[tcount_motion+1][1]-motion[tcount_motion][1])/(motion[tcount_motion+1][0]-motion[tcount_motion][0]);
	
    Uext = vel*nvecx;
	Vext = vel*nvecy;
	Wext = vel*nvecz;
    
	if(p->simtime>motion[tcount_motion+1][0])
	++tcount_motion;
    }
	
}

void sixdof_gc::read_motionvec(lexer *p, fdm* a, ghostcell* pgc)
{
	char name[100];
	double val,val0,val1;
    double sign,beta,s;
	int count,ptnum;
	
	sprintf(name,"motionfile_vec.dat");

// open file------------
	ifstream file(name, ios_base::in);
	
	if(!file)
	cout<<endl<<("no 'motionfile_vec.dat' file found")<<endl<<endl;
	
	count=0;
	while(!file.eof())
	{
	file>>val0>>val1;
	++count;
	}
	
	file.close();

    
    ptnum=count;
	
	p->Darray(motion,ptnum,3);
	
	file.open ("motionfile_vec.dat", ios_base::in);
	
	count=0;
	while(!file.eof())
	{
        file>>val0>>val1;
        
        motion[count][0] = val0;
        motion[count][1] = val1;
        ++count;
	}
	
	ts_motion = motion[0][0];
	te_motion = motion[ptnum-1][0];
    
    
    // vec
    vecx = p->X221_xe - p->X221_xs;
    vecy = p->X221_ye - p->X221_ys;
    vecz = p->X221_ze - p->X221_zs;
    
    double length = sqrt(vecx*vecx + vecy*vecy + vecz*vecz);
    
    length = length>1.0e-20?length:1.0e20;
    
    nvecx = vecx/length;
    nvecy = vecy/length;
    nvecz = vecz/length;
    
    tcount_motion=0;
    
}

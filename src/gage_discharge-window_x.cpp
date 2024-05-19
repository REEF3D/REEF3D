/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

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

#include"gage_discharge_window_x.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include<sys/stat.h>
#include<sys/types.h>

gage_discharge_window_x::gage_discharge_window_x(lexer *p, fdm* a, ghostcell *pgc)
{
	p->Iarray(iloc,p->P168);
	p->Iarray(flag,p->P168);
	p->Darray(q,p->P168);
	
	// Create Folder
	if(p->mpirank==0)
	mkdir("./REEF3D_Log",0777);
	
    if(p->mpirank==0 && p->P168>0)
    {
    // open file
	qout.open("./REEF3D_Log/REEF3D-discharge_window_x.dat");

    qout<<"number of x-discharge window gages:  "<<p->P168<<endl<<endl;
    qout<<"x_coord   zs   ze"<<endl;
    for(n=0;n<p->P168;++n)
    qout<<n+1<<"\t "<<p->P168_x[n]<<" "<<p->P168_zs[n]<<" "<<p->P168_ze[n]<<endl;

    qout<<endl<<endl;

    qout<<"iter";
    for(n=0;n<p->P168;++n)
    qout<<"\t P"<<n+1;

    qout<<endl<<endl;
    }
	
	ini_location(p,a,pgc);
}

gage_discharge_window_x::~gage_discharge_window_x()
{
    qout.close();
}

void gage_discharge_window_x::start(lexer *p, fdm *a, ghostcell *pgc)
{
    double epsi,H;

    for(n=0;n<p->P168;++n)
    q[n]=0.0;

	
    for(n=0;n<p->P168;++n)
    {
	area=0.0;
	Ai=0.0;

    i=iloc[n];
		
        /*
        if(flag[n]==1)
        JLOOP
        KLOOP
        PCHECK
        {
			area=0.0;
            if(a->phi(i,j,k)>-0.5*p->DZN[KP]-1.0e-20 && a->topo(i,j,k)>0.0)
			{
            if(a->phi(i,j,k)>=0.5*p->DZN[KP])
            area=p->DYN[JP]*p->DZN[KP];

            if(a->phi(i,j,k)<0.5*p->DZN[KP] && a->phi(i,j,k)>0.0)
            area=p->DYN[JP]*(p->DZN[KP]*0.5 + a->phi(i,j,k));
			
			if(a->phi(i,j,k)>=-0.5*p->DZN[KP] -1.0e-20 && a->phi(i,j,k)<=0.0*p->DZN[KP])
            area=p->DYN[JP]*(p->DZN[KP]*0.5 - fabs(a->phi(i,j,k)));

            q[n]+=area*0.5*(a->u(i,j,k) + a->u(i-1,j,k));
			}
        }*/
        
        if(flag[n]==1)
        JLOOP
        KLOOP
        if(p->ZP[KP]>p->P168_zs[n] && p->ZP[KP]<=p->P168_ze[n])
        PCHECK
        {
            epsi = 1.6*p->DXN[KP];
            
            if(a->phi(i,j,k)>epsi)
            H=1.0;

            if(a->phi(i,j,k)<-epsi)
            H=0.0;

            if(fabs(a->phi(i,j,k))<=epsi)
            H=0.5*(1.0 + a->phi(i,j,k)/epsi + (1.0/PI)*sin((PI*a->phi(i,j,k))/epsi));

            area=H*p->DYN[JP]*p->DZN[KP];

            q[n]+=area*0.5*(a->u(i,j,k) + a->u(i-1,j,k));
			
        }
	
    }
	
	
    for(n=0;n<p->P168;++n)
    q[n]=pgc->globalsum(q[n]);
	
	if(p->mpirank==0 && p->P166==1)
	for(n=0;n<p->P168;++n)
	cout<<n+1<<setprecision(6)<<" Qi: "<<q[n]<<endl;  

    // write to file
    if(p->mpirank==0)
    {
    qout<<setprecision(9)<<p->count<<"\t";
    for(n=0;n<p->P168;++n)
    qout<<setprecision(6)<<q[n]<<"  \t  ";
    qout<<endl;
    }
}

void gage_discharge_window_x::ini_location(lexer *p, fdm *a, ghostcell *pgc)
{
	
	for(n=0;n<p->P168;++n)
	flag[n]=0;

    for(n=0;n<p->P168;++n)
    {
    iloc[n] = p->posc_i(p->P168_x[n]);
	
	if(iloc[n]>=0 && iloc[n]<p->knox)
	flag[n]=1;
	}
}


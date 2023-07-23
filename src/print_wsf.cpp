/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2023 Hans Bihs

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

#include"print_wsf.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include<sys/stat.h>
#include<sys/types.h>

print_wsf::print_wsf(lexer *p, fdm* a, ghostcell *pgc, int num)
{
	gauge_num = p->P51;
	x = p->P51_x;
	y = p->P51_y;
	
    if(p->P51>0 && num==0)
	{
	gauge_num = p->P51;
	x = p->P51_x;
	y = p->P51_y;
	}
	
	if(p->P351>0 && num==1)
	{
	gauge_num = p->P351;
	x = p->P351_x;
	y = p->P351_y;
	}
	
	if(p->P352>0 && num==2)
	{
	gauge_num = p->P352;
	x = p->P352_x;
	y = p->P352_y;
	}
	
	// Create Folder
	if(p->mpirank==0 && p->P14==1)
	mkdir("./REEF3D_CFD_WSF",0777);
	
    if(p->mpirank==0 && p->P51>0 && num==0)
    {
    // open file
	if(p->P14==0)
    wsfout.open("REEF3D-CFD-WSF-HG.dat");
	
	if(p->P14==1)
	wsfout.open("./REEF3D_CFD_WSF/REEF3D-CFD-WSF-HG.dat");

    wsfout<<"number of gauges:  "<<gauge_num<<endl<<endl;
    wsfout<<"x_coord     y_coord"<<endl;
    for(n=0;n<gauge_num;++n)
    wsfout<<n+1<<"\t "<<x[n]<<"\t "<<y[n]<<endl;

    wsfout<<endl<<endl;

    wsfout<<"time";
    for(n=0;n<gauge_num;++n)
    wsfout<<"\t P"<<n+1;

    wsfout<<endl<<endl;
    }
	
	//-------------------
	
	if(p->mpirank==0 && p->P351>0 && num==1)
    {
    // open file
	if(p->P14==0)
    wsfout.open("REEF3D-CFD-WSF-HG-1.dat");
	
	if(p->P14==1)
	wsfout.open("./REEF3D_CFD_WSF/REEF3D-CFD-WSF-HG-1.dat");

    wsfout<<"number of gauges:  "<<gauge_num<<endl<<endl;
    wsfout<<"x_coord     y_coord"<<endl;
    for(n=0;n<gauge_num;++n)
    wsfout<<n+1<<"\t "<<x[n]<<"\t "<<y[n]<<endl;

    wsfout<<endl<<endl;

    wsfout<<"time";
    for(n=0;n<gauge_num;++n)
    wsfout<<"\t P"<<n+1;

    wsfout<<endl<<endl;
    }
	
	//-------------------
	
	if(p->mpirank==0 && p->P352>0 && num==2)
    {
    // open file
	if(p->P14==0)
    wsfout.open("REEF3D-CFD-WSF-HG-2.dat");
	
	if(p->P14==1)
	wsfout.open("./REEF3D_CFD_WSF/REEF3D-CFD-WSF-HG-2.dat");

    wsfout<<"number of gauges:  "<<gauge_num<<endl<<endl;
    wsfout<<"x_coord     y_coord"<<endl;
    for(n=0;n<gauge_num;++n)
    wsfout<<n+1<<"\t "<<x[n]<<"\t "<<y[n]<<endl;

    wsfout<<endl<<endl;

    wsfout<<"time";
    for(n=0;n<gauge_num;++n)
    wsfout<<"\t P"<<n+1;

    wsfout<<endl<<endl;
    }
	
	p->Iarray(iloc,gauge_num);
	p->Iarray(jloc,gauge_num);
	p->Iarray(flag,gauge_num);
	p->Darray(wsf,gauge_num);

    ini_location(p,a,pgc);
}

print_wsf::~print_wsf()
{
    wsfout.close();
}

void print_wsf::height_gauge(lexer *p, fdm *a, ghostcell *pgc, field &f)
{
    double zval=0.0;

    for(n=0;n<gauge_num;++n)
    wsf[n]=-1.0e20;

	if(p->A10==6)
    for(n=0;n<gauge_num;++n)
    if(flag[n]>0)
    {
    zval=0.0;

    i=iloc[n];
    j=jloc[n];
	
        KLOOP
        PCHECK
        {
            if(f(i,j,k)>=0.0 && f(i,j,k+1)<0.0)
            wsf[n]=MAX(wsf[n],-(f(i,j,k)*p->DZP[KP])/(f(i,j,k+1)-f(i,j,k)) + p->pos_z());
        }
    }
    
    if(p->A10==5 || p->A10==4)
    for(n=0;n<gauge_num;++n)
    if(flag[n]>0)
    {
    zval=0.0;

    i=iloc[n];
    j=jloc[n];
	
			wsf[n] = a->eta(i,j);

    }
	
    for(n=0;n<gauge_num;++n)
    wsf[n]=pgc->globalmax(wsf[n]);

    // write to file
    if(p->mpirank==0)
    {
    wsfout<<setprecision(9)<<p->simtime<<"\t";
    for(n=0;n<gauge_num;++n)
    wsfout<<setprecision(9)<<wsf[n]<<"  \t  ";
    wsfout<<endl;
    }
}

void print_wsf::write(lexer *p, fdm *a, ghostcell *pgc)
{
}

void print_wsf::ini_location(lexer *p, fdm *a, ghostcell *pgc)
{
    int check;

    for(n=0;n<gauge_num;++n)
    {
    iloc[n] = p->posc_i(x[n]); 
    
    if(p->j_dir==0)
    jloc[n]=0;
    
    if(p->j_dir==1)
    jloc[n] = p->posc_j(y[n]); 

    check=ij_boundcheck(p,a,iloc[n],jloc[n],0);

    if(check==1)
    flag[n]=1;
    
    //cout<<p->mpirank<<" n: "<<n<<" flag: "<<flag[n]<<" x: "<<x[n]<<" y: "<<y[n]<<" iloc: "<<iloc[n]<<" jloc: "<<jloc[n]<<endl;
    }
}



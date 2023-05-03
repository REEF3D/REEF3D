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

#include"fnpf_print_wsf.h"
#include"lexer.h"
#include"fdm2D.h"
#include"ghostcell.h"
#include<sys/stat.h>
#include<sys/types.h>

fnpf_print_wsf::fnpf_print_wsf(lexer *p, fdm_fnpf *c)
{

	gauge_num = p->P51;
	x = p->P51_x;
	y = p->P51_y;

	
	// Create Folder
	if(p->mpirank==0 && p->P14==1)
	mkdir("./REEF3D_FNPF_WSF",0777);
	
    if(p->mpirank==0 && p->P51>0)
    {
    // open WSF file
	if(p->P14==0)
    wsfout.open("REEF3D-FNPF-WSF-HG.dat");
	
	if(p->P14==1)
	wsfout.open("./REEF3D_FNPF_WSF/REEF3D-FNPF-WSF-HG.dat");

    wsfout<<"number of gauges:  "<<gauge_num<<endl<<endl;
    wsfout<<"x_coord     y_coord"<<endl;
    for(n=0;n<gauge_num;++n)
    wsfout<<n+1<<"\t "<<x[n]<<"\t "<<y[n]<<endl;

    wsfout<<endl<<endl;

    wsfout<<"time";
    for(n=0;n<gauge_num;++n)
    wsfout<<"\t P"<<n+1;

    wsfout<<endl<<endl;
    
    
        if(p->P57==1 || p->P57==3)
        {
        // open rise velovity file
        if(p->P14==0)
        detaout.open("REEF3D-FNPF-Rise-Velocity.dat");
        
        if(p->P14==1)
        detaout.open("./REEF3D_FNPF_WSF/REEF3D-FNPF-Rise-Velocity.dat");

        detaout<<"number of rise velocity gauges:  "<<gauge_num<<endl<<endl;
        detaout<<"x_coord     y_coord"<<endl;
        for(n=0;n<gauge_num;++n)
        detaout<<n+1<<"\t "<<x[n]<<"\t "<<y[n]<<endl;

        detaout<<endl<<endl;

        detaout<<"time";
        for(n=0;n<gauge_num;++n)
        detaout<<"\t P"<<n+1;

        detaout<<endl<<endl;
        }
        
        if(p->P57==2 || p->P57==3)
        {
        // open horizontal velocity file
        if(p->P14==0)
        Uhorzout.open("REEF3D-FNPF-Horiontal-FSF-Velocity.dat");
        
        if(p->P14==1)
        Uhorzout.open("./REEF3D_FNPF_WSF/REEF3D-FNPF-Horiontal-FSF-Velocity.dat");

        Uhorzout<<"number of horizontal fsf velocity gauges:  "<<gauge_num<<endl<<endl;
        Uhorzout<<"x_coord     y_coord"<<endl;
        for(n=0;n<gauge_num;++n)
        Uhorzout<<n+1<<"\t "<<x[n]<<"\t "<<y[n]<<endl;

        Uhorzout<<endl<<endl;

        Uhorzout<<"time";
        for(n=0;n<gauge_num;++n)
        Uhorzout<<"\t P"<<n+1;

        Uhorzout<<endl<<endl;
        }
    }
	
	//-------------------
	
	
	p->Iarray(iloc,gauge_num);
	p->Iarray(jloc,gauge_num);
	p->Iarray(flag,gauge_num);
	p->Darray(wsf,gauge_num);
    p->Darray(deta,gauge_num);
    p->Darray(Uhorz,gauge_num);

    ini_location(p,c);
}

fnpf_print_wsf::~fnpf_print_wsf()
{
    wsfout.close();
}

void fnpf_print_wsf::height_gauge(lexer *p, fdm_fnpf *c, ghostcell *pgc, slice &f)
{
    
    fill_eta(p,c,pgc,f);
    
    // write to file
    if(p->mpirank==0)
    {
    wsfout<<setprecision(9)<<p->simtime<<"\t";
    for(n=0;n<gauge_num;++n)
    wsfout<<setprecision(9)<<wsf[n]<<"  \t  ";
    wsfout<<endl;
    }
    
    // Rise Velocity
    
    
    
    // Horionztal Velocity
}

void fnpf_print_wsf::fill_eta(lexer *p, fdm_fnpf *c, ghostcell *pgc, slice &f)
{
    double zval=0.0;

    for(n=0;n<gauge_num;++n)
    wsf[n]=-1.0e20;

	
    for(n=0;n<gauge_num;++n)
    if(flag[n]>0)
    {
    zval=0.0;

    i=iloc[n];
    j=jloc[n];
	
			wsf[n] = f(i,j);

    }
	
    for(n=0;n<gauge_num;++n)
    wsf[n]=pgc->globalmax(wsf[n]);
}

void fnpf_print_wsf::fill_deta(lexer *p, fdm_fnpf *c, ghostcell *pgc, slice &f)
{
    
}
    
void fnpf_print_wsf::fill_Uhorz(lexer *p, fdm_fnpf *c, ghostcell *pgc, slice &f)
{
    
}


void fnpf_print_wsf::ini_location(lexer *p, fdm_fnpf *c)
{
    for(n=0;n<gauge_num;++n)
    {
    iloc[n] = p->posc_i(x[n]); 
    
    if(p->j_dir==0)
    jloc[n] = 0; 
    
    if(p->j_dir==1)
    jloc[n] = p->posc_j(y[n]); 

    if(iloc[n]>=0 && iloc[n]<p->knox)
    if(jloc[n]>=0 && jloc[n]<p->knoy)
    flag[n]=1;
    }
}

int fnpf_print_wsf::conv(double a)
{

    int b,c;
    double d,diff;

    c= int( a);
    d=double(c);
    diff=a-d;

    b=c;

    if(diff>0.5)
    b=c+1;

    if(diff<=-0.5)
    b=c-1;

    return b;

}

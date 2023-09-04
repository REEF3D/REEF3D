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

#include"wave_lib_hdc.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

wave_lib_hdc::wave_lib_hdc(lexer *p, ghostcell *pgc) : wave_lib_parameters(p,pgc) 
{ 
    // read header
    read_header(p,pgc);
    allocate(p,pgc);
    
    // time_interpol
    if(p->mpirank==0)
    {
    cout<<"Wave Tank: wave coupling FNPF->CFD "<<endl;
    cout<<" HDC Nx: "<<Nx<<" Ny: "<<Ny<<" Nz: "<<Nz<<" . jdir: "<<jdir<<endl;
    cout<<" HDC Xs: "<<Xstart<<" Xe: "<<Xend<<" Ys: "<<Ystart<<" Ye: "<<Yend<<endl;
    cout<<" HDC numiter: "<<numiter<<" t_start: "<<t_start<<" t_end: "<<t_end<<endl;
    cout<<" HDC simtime[0]: "<<simtime[0]<<" simtime[numiter-1]: "<<simtime[numiter-1]<<endl;
    }
    
    startup=0;
    endseries=0;
}

wave_lib_hdc::~wave_lib_hdc()
{
}

double wave_lib_hdc::wave_u(lexer *p, double x, double y, double z)
{
    double vel=0.0;
    
    if(p->B125==1)
    y=p->B125_y;
    
    if(endseries==0)
    vel = space_interpol(p,U,x,y,z);
    
    return vel;
}

double wave_lib_hdc::wave_v(lexer *p, double x, double y, double z)
{
    double vel=0.0;
    
    if(p->B125==1)
    y=p->B125_y;
    
    if(endseries==0 && p->B125==0 && p->B127==0)
    vel = space_interpol(p,V,x,y,z);

    return vel;
}

double wave_lib_hdc::wave_w(lexer *p, double x, double y, double z)
{
    double vel=0.0;
    
    if(p->B125==1)
    y=p->B125_y;
    
    if(endseries==0)
    vel = space_interpol(p,W,x,y,z);

    return vel;
}

double wave_lib_hdc::wave_eta(lexer *p, double x, double y)
{
    double eta=0.0;
    
    if(p->B125==1)
    y=p->B125_y;
    
    if(endseries==0)
    eta = plane_interpol(p,E,x,y);
    
    return eta;
}

double wave_lib_hdc::wave_fi(lexer *p, double x, double y, double z)
{
    double fi=0.0;
    
    return fi;
}

void wave_lib_hdc::parameters(lexer *p, ghostcell *pgc)
{
}

void wave_lib_hdc::wave_prestep(lexer *p, ghostcell *pgc)
{
    // only at startup
    if(startup==0)
    {
        deltaT = simtime[1]-simtime[0];
        
        t1 = (simtime[1]-(p->simtime+p->I241))/deltaT;
        t2 = ((p->simtime+p->I241)-simtime[0])/deltaT;
        
        q1 = diter;
        q2 = diter+1;
        
        if(file_type==2)
        {
        filename_continuous(p,pgc);
        result.open(name, ios::binary);
        }
    
        read_result(p,pgc,E1,U1,V1,W1,q1);
        read_result(p,pgc,E2,U2,V2,W2,q2);
        startup=1;
        
        
        // find q1
        while(simtime[q1+1-diter]<=p->simtime+p->I241)
        {
        ++q1;
        
        //cout<<"HDC ++q1: "<<q1<<endl;
        if(file_type==2)
        read_result_continuous(p,pgc,E1,U1,V1,W1,q1);
        }
            
        // find q2
        while(simtime[q2-diter]<p->simtime+p->I241)
        {
        ++q2;
        
        //cout<<"HDC ++q2: "<<q2<<endl;
        if(file_type==2 )
        read_result_continuous(p,pgc,E2,U2,V2,W2,q2);
        }
    }
    
    q1n=q1;
    q2n=q2;
    
    // check: open next timestep
           
    // find q1
    while(simtime[q1+1-diter]<=p->simtime+p->I241)
    ++q1;
        
    // find q2
    while(simtime[q2-diter]<p->simtime+p->I241)
    ++q2;
    
    if(q2>=numiter+diter)
    endseries=1;
        
    q1=MIN(q1,numiter+diter);
    q2=MIN(q2,numiter+diter);
    
    if(q1==q2)
    ++q2;
    
    
    // single file read 
    if(file_type==1)
    {
        if(q1!=q1n)
        {
        // Open File 1
        filename_single(p,pgc,q1);
        read_result(p,pgc,E1,U1,V1,W1,q1);
        }
        
        if(q2!=q2n)
        {
        // Open File 2
        filename_single(p,pgc,q2);
        read_result(p,pgc,E2,U2,V2,W2,q2);
        }
    }
        
    // contiuous file read
    if(file_type==2)
    {
        if(q1!=q1n)
        fill_result_continuous(p,pgc);
        
        if(q2!=q2n)
        read_result_continuous(p,pgc,E2,U2,V2,W2,q2);
    }
        

        deltaT = simtime[q2-diter]-simtime[q1-diter];
        
        if(p->mpirank==0)
        cout<<"HDC  q1: "<<q1<<" q2: "<<q2<<" t1: "<<t1<<" t2: "<<t2<<" deltaT: "<<deltaT<<" simtime[q1]: "<<simtime[q1-diter]<<" simtime[q2]: "<<simtime[q2-diter]<<endl;

        t1 = (simtime[q2-diter]-(p->simtime+p->I241))/deltaT;
        t2 = ((p->simtime+p->I241)-simtime[q1-diter])/deltaT;
        
        
    // time interpolation
    //if(q1!=q1n || q2!=q2n)
    time_interpol(p);
    
}

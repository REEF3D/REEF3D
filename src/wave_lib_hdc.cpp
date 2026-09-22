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
    
    if(p->A10==3)
    allocate_fnpf(p,pgc);
    
    if(p->A10==6)
    allocate_cfd(p,pgc);
    
    // time_interpol
    if(p->mpirank==0)
    {
    cout<<"Wave_Lib: hydrodynamic coupling "<<endl;
    cout<<" HDC Nx: "<<Nx<<" Ny: "<<Ny<<" Nz: "<<Nz<<" . jdir: "<<jdir<<endl;
    cout<<" HDC Xs: "<<Xstart<<" Xe: "<<Xend<<" Ys: "<<Ystart<<" Ye: "<<Yend<<endl;
    cout<<" HDC numiter: "<<numiter<<" t_start: "<<t_start<<" t_end: "<<t_end<<endl;
    cout<<" HDC simtime[0]: "<<simtime[0]<<" simtime[numiter-1]: "<<simtime[numiter-1]<<endl;
    }
    
    startup=0;
    endseries=0;
    
    /*
    // testfile
    if(p->mpirank==0)
    {
    ifstream testfile;
    
    testfile.open("testfile.dat", ios::binary);
    
        for(int qn=0;qn<100;++qn)
        {
        testfile.read((char*)&ffn, sizeof (float));
        cout<<"ffn: "<<ffn<<endl;
        }
    }*/
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
    
    if(p->B125==1)
    y=p->B125_y;
    
    if(endseries==0)
    fi = plane_interpol(p,F,x,y);
    
    return fi;
}

void wave_lib_hdc::parameters(lexer *p, ghostcell *pgc)
{
}

void wave_lib_hdc::wave_prestep(lexer *p, ghostcell *pgc)
{
    if(p->A10==3)
    wave_prestep_fnpf(p,pgc);
    
    if(p->A10==6)
    wave_prestep_cfd(p,pgc);
}

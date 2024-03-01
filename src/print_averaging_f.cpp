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

#include"print_averaging_f.h"
#include"lexer.h"
#include"fdm.h"
#include"heat.h"
#include"ghostcell.h"
#include<iostream>
#include<fstream>
#include<sys/stat.h>
#include<sys/types.h>

print_averaging_f::print_averaging_f(lexer *p, fdm* a, ghostcell *pgc) : um(p),vm(p),wm(p),pm(p),Tm(p)
{
    ULOOP
    um(i,j,k) = 0.0;
    
    VLOOP
    vm(i,j,k) = 0.0;
    
    WLOOP
    wm(i,j,k) = 0.0;
    
    LOOP
    {
    pm(i,j,k) = 0.0;
    Tm(i,j,k) = 0.0;
    }
    
    stime = p->P22;
}

print_averaging_f::~print_averaging_f()
{

}

void print_averaging_f::averaging(lexer *p, fdm *a, ghostcell *pgc, heat *pheat)
{
    // u,v,w,p,T
    
    if(p->simtime>stime)
    {
        ULOOP
        um(i,j,k) = um(i,j,k) + p->dt*a->u(i,j,k);
        
        VLOOP
        vm(i,j,k) = vm(i,j,k) + p->dt*a->v(i,j,k);
        
        WLOOP
        wm(i,j,k) = wm(i,j,k) + p->dt*a->w(i,j,k);
        
        LOOP
        pm(i,j,k) = pm(i,j,k) + p->dt*a->press(i,j,k);
        
        LOOP
        Tm(i,j,k) = Tm(i,j,k) + p->dt*pheat->val(i,j,k);
    }
}

void print_averaging_f::offset_vtu(lexer *p, fdm *a, ghostcell *pgc, ofstream &result, int *offset, int &n)
{
    // velocity
	offset[n]=offset[n-1]+4*(p->pointnum)*3+4;
	++n;
    // pressure
    offset[n]=offset[n-1]+4*(p->pointnum)+4;
	++n;
    
    // heat
    if(p->H10>0)
    {
    offset[n]=offset[n-1]+4*(p->pointnum)+4;
	++n;
    }
}

void print_averaging_f::name_vtu(lexer *p, fdm *a, ghostcell *pgc, ofstream &result, int *offset, int &n)
{
    result<<"<DataArray type=\"Float32\" Name=\"velocity_mean\" NumberOfComponents=\"3\" format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
    
	result<<"<DataArray type=\"Float32\" Name=\"pressure_mean\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
    
	if(p->H10>0)
    {
    result<<"<DataArray type=\"Float32\" Name=\"T_mean\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
    }
}

void print_averaging_f::name_pvtu(lexer *p, fdm *a, ghostcell *pgc, ofstream &result)
{
    result<<"<PDataArray type=\"Float32\" Name=\"velocity_mean\" NumberOfComponents=\"3\"/>"<<endl;
    result<<"<PDataArray type=\"Float32\" Name=\"pressure_mean\"/>"<<endl;
    
    if(p->H10>0)
    result<<"<PDataArray type=\"Float32\" Name=\"T_mean\"/>"<<endl;
}

void print_averaging_f::print_3D(lexer* p, fdm *a, ghostcell *pgc, ofstream &result)
{
    pgc->start1(p,um,110);
    pgc->start2(p,vm,111);
	pgc->start3(p,wm,112);
    pgc->start4(p,pm,40);
    pgc->start4(p,Tm,1);
    
    
    //  Velocities
    if(p->simtime<=stime+1.0e-8)
    {
    iin=3*4*(p->pointnum);
	result.write((char*)&iin, sizeof (int));
    TPLOOP
	{
	ffn=0.0;
	result.write((char*)&ffn, sizeof (float));

	ffn=0.0;
	result.write((char*)&ffn, sizeof (float));

	ffn=0.0;
	result.write((char*)&ffn, sizeof (float));
	}
    
    //  Pressure
	iin=4*(p->pointnum);
	result.write((char*)&iin, sizeof (int));
	TPLOOP
	{
	ffn=0.0;
	result.write((char*)&ffn, sizeof (float));
	}
    
    //  Temperature
    if(p->H10>0)
    {
	iin=4*(p->pointnum);
	result.write((char*)&iin, sizeof (int));
	TPLOOP
	{
	ffn=0.0;
	result.write((char*)&ffn, sizeof (float));
	}
    }
    
    }
    
    if(p->simtime>stime+1.0e-8)
    {
    iin=3*4*(p->pointnum);
	result.write((char*)&iin, sizeof (int));
    TPLOOP
	{
	ffn=float(p->ipol1(um)/(p->simtime-stime));
	result.write((char*)&ffn, sizeof (float));

	ffn=float(p->ipol2(vm)/(p->simtime-stime));
	result.write((char*)&ffn, sizeof (float));

	ffn=float(p->ipol3(wm)/(p->simtime-stime));
	result.write((char*)&ffn, sizeof (float));
	}
    
    //  Pressure
	iin=4*(p->pointnum);
	result.write((char*)&iin, sizeof (int));
	TPLOOP
	{
	ffn=float(p->ipol4press(pm)/(p->simtime-stime));
	result.write((char*)&ffn, sizeof (float));
	}
    
    //  Temperature
    if(p->H10>0)
    {
	iin=4*(p->pointnum);
	result.write((char*)&iin, sizeof (int));
	TPLOOP
	{
	ffn=float(p->ipol4(Tm)/(p->simtime-stime));
	result.write((char*)&ffn, sizeof (float));
	}
    }
    
    }

}


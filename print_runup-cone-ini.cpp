/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2020 Hans Bihs

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

#include"print_runup.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

void print_runup::cone_ini(lexer *p, fdm* a, ghostcell *pgc)
{	
	double U,ds,phi;
	double xm,ym,z1,z2,r1,r2;
	double rmax;
	int snum,n;
	int vertice_mem, center1_num,center2_num;
	int rank=0;
	
	xm=p->P101_xm;
    ym=p->P101_ym;
	
	z1=p->P101_zs;
	z2=p->P101_ze;
	
    r1=p->P101_r1;
	r2=p->P101_r2;    
	
	rmax = MAX(r1,r2);
	
	U = 2.0 * PI * rmax;
	
	ds = 0.1*(U*p->dx);
	
	snum = int(U/ds);
	
	//ds = U/double(snum);

// Vertices	
	ds = (2.0*PI)/double(snum);
	
	++snum;
	
	//cout<<"snum: "<<snum<<" U: "<<U<<endl;
	if(p->mpirank==0)
	cout<<"RUNUP: xm,ym: "<<xm<<"   "<<ym<<endl;
	p->Darray(line,snum*3,6);
	
	//side lines
	phi=0.0;
	line_num=0;
	for(n=0;n<snum;++n)
	{
	line[line_num][0] = xm + r1*cos(phi);
	line[line_num][1] = ym + r1*sin(phi);
	line[line_num][2] = z1;
	line[line_num][3] = xm + r2*cos(phi);
	line[line_num][4] = ym + r2*sin(phi);
	line[line_num][5] = z2;
	phi+=ds;
	++line_num;
	}
	
	vertice_mem = line_num;
/*
	// top lines
	phi=0.0;
	for(n=0;n<snum;++n)
	{
	line[line_num][0] = xm + r2*cos(phi);
	line[line_num][1] = ym + r2*sin(phi);
	line[line_num][2] = z2;
	line[line_num][3] = xm;
	line[line_num][4] = ym;
	line[line_num][5] = z2;
	phi+=ds;
	++line_num;
	}
*/
	line_num_all = line_num*p->M10;
	
	
	p->Darray(cut_x,line_num);
	p->Darray(cut_y,line_num);
	p->Darray(cut_z,line_num);
	p->Iarray(cut_ID,line_num);
	p->Iarray(cut_active,line_num);
	
	p->Darray(cutall_x,line_num_all);
	p->Darray(cutall_y,line_num_all);
	p->Darray(cutall_z,line_num_all);
	p->Iarray(cutall_ID,line_num_all);
	p->Iarray(cutall_count,p->M10);
	
	p->Darray(runup_x,line_num);
	p->Darray(runup_y,line_num);
	p->Darray(runup_z,line_num);
	p->Iarray(runup_active,line_num);
	
	p->Iarray(displ,p->M10);
	
	
	for(n=0;n<line_num;++n)
	{
	runup_x[n]=0.0;
	runup_y[n]=0.0;
	runup_z[n]=-10000000.0;
	runup_active[n]=-1;
	}
	
}
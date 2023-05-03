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

#include"print_wsf_theory.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"ioflow.h"
#include<sys/stat.h>
#include<sys/types.h>

print_wsf_theory::print_wsf_theory(lexer *p, fdm* a, ghostcell *pgc, int num)
{
	gauge_num = p->P50;
	x = p->P50_x;
	y = p->P50_y;
	
    if(p->P50>0 && num==0)
	{
	gauge_num = p->P50;
	x = p->P50_x;
	y = p->P50_y;
	}
	
	// Create Folder
	if(p->mpirank==0 && p->P14==1)
	mkdir("./REEF3D_CFD_WSF",0777);
	
    if(p->mpirank==0 && p->P50>0 && num==0)
    {
    // open file
	if(p->P14==0)
    wsfout.open("REEF3D-CFD-WSF-HG-THEORY.dat");
	
	if(p->P14==1)
	wsfout.open("./REEF3D_CFD_WSF/REEF3D-CFD-WSF-HG-THEORY.dat");

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
	
	
}

print_wsf_theory::~print_wsf_theory()
{
    wsfout.close();
}

void print_wsf_theory::height_gauge(lexer *p, fdm *a, ghostcell *pgc, ioflow *pflow, field &f)
{
	
    // write to file
    if(p->mpirank==0)
    {
    wsfout<<setprecision(9)<<p->simtime<<"\t";
    for(n=0;n<gauge_num;++n)
    wsfout<<setprecision(9)<<pflow->wave_fsf(p,pgc,x[n])<<"  \t  ";
    wsfout<<endl;
    }
}





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

#include<iomanip>
#include"print_wsfline_y.h"
#include"wsfline_core.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"ioflow.h"
#include"wave_interface.h"
#include<sys/stat.h>
#include<sys/types.h>

print_wsfline_y::print_wsfline_y(lexer *p, fdm* a, ghostcell *pgc)
{
	p->Iarray(iloc,p->P56);

    pcore = new wsfline_core(p,pgc,p->P56,p->knoy);

    ini_location(p,a,pgc);

	// Create Folder
	if(p->mpirank==0)
	mkdir("./REEF3D_CFD_WSFLINE_Y",0777);
}

print_wsfline_y::~print_wsfline_y()
{
    wsfout.close();
    delete pcore;
}

void print_wsfline_y::wsfline(lexer *p, fdm *a, ghostcell *pgc, ioflow *pflow)
{
    char name[250];
    int num;

    num = p->count;

    if(p->mpirank==0)
    {
		// open file
		snprintf(name,sizeof(name),"./REEF3D_CFD_WSFLINE_Y/REEF3D-CFD-wsfline_y-%08i.dat",num);

		wsfout.open(name);

		wsfout<<"simtime:  "<<p->simtime<<endl;
		wsfout<<"number of wsf-lines_y:  "<<p->P56<<endl<<endl;
		wsfout<<"line_No     x_coord"<<endl;
		for(q=0;q<p->P56;++q)
		wsfout<<q+1<<"\t "<<p->P56_x[q]<<endl;

		if(p->P53==1)
		wsfout<<q+1<<"\t "<<" Wave Theory "<<endl;

		wsfout<<endl<<endl;

		wsfout<<"y_coord";
		for(q=0;q<p->P56+p->P53;++q)
		wsfout<<"\t \t P "<<q+1;

		wsfout<<endl<<endl;
    }

    //-------------------

    pcore->reset();

    for(q=0;q<p->P56;++q)
    {
        JLOOP
        if(pcore->flag[q][j]>0)
        {
		
        i=iloc[q];

            KLOOP
            PCHECK
            {
                if(a->phi(i,j,k)>=0.0 && a->phi(i,j,k+1)<0.0)
                {
                pcore->wsf[q][j]=MAX(pcore->wsf[q][j],-(a->phi(i,j,k)*p->DZP[KP])/(a->phi(i,j,k+1)-a->phi(i,j,k)) + p->pos_z());
                pcore->loc[q][j]=p->pos_y();
                }
            }
        }
    }

    // gather, sort and write the rows
    pcore->write_rows(p,pgc,wsfout,5,
                      [&](double x){return x;},
                      p->P53==1 ? function<double(double)>([&](double x){return pflow->wave_fsf(p,pgc,x);}) : function<double(double)>(),
                      " \t \t  ");

    if(p->mpirank==0)
    {
    wsfout.close();
    }
}

void print_wsfline_y::ini_location(lexer *p, fdm *a, ghostcell *pgc)
{
    int check,count;

    for(q=0;q<p->P56;++q)
    {
        count=0;
        JLOOP
        {
        iloc[q]=p->posc_i(p->P56_x[q]);

        check=ij_boundcheck(p,iloc[q],j,0);

        if(check==1)
        pcore->flag[q][count]=1;

        ++count;
        }
    }
}

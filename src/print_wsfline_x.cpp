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
#include"print_wsfline_x.h"
#include"wsfline_core.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"ioflow.h"
#include"wave_interface.h"
#include<sys/stat.h>
#include<sys/types.h>
#include"runlog.h"

print_wsfline_x::print_wsfline_x(lexer *p, fdm* a, ghostcell *pgc)
{
	p->Iarray(jloc,p->P52);

    pcore = new wsfline_core(p,pgc,p->P52,p->knox);

    ini_location(p,a,pgc);

	// Create Folder
	if(p->mpirank==0)
	mkdir("./REEF3D_CFD_WSFLINE",0777);
}

print_wsfline_x::~print_wsfline_x()
{
    wsfout.close();
    delete pcore;
}

void print_wsfline_x::wsfline(lexer *p, fdm *a, ghostcell *pgc, ioflow *pflow)
{
    char name[250];
    int num;

    num = p->count;

    if(p->mpirank==0)
    {
		// open file
		snprintf(name,sizeof(name),"./REEF3D_CFD_WSFLINE/REEF3D-CFD-wsfline-%08i.dat",num);
        
		wsfout.open(name);

		wsfout<<"simtime:  "<<p->simtime<<endl;
		wsfout<<"number of wsf-lines:  "<<p->P52<<endl<<endl;
		wsfout<<"line_No     y_coord"<<endl;
		for(q=0;q<p->P52;++q)
		wsfout<<q+1<<"\t "<<p->P52_y[q]<<endl;

		if(p->P53==1)
		wsfout<<q+1<<"\t "<<" Wave Theory "<<endl;

		wsfout<<endl<<endl;

		
		for(q=0;q<p->P52;++q)
		{
		wsfout<<"X "<<q+1;
		wsfout<<"\t P "<<q+1<<" \t \t ";
		if(p->P53==1)
		wsfout<<"\t \t W "<<q+1;
		}

		wsfout<<endl<<endl;
    }

    //-------------------

    pcore->reset();

    for(q=0;q<p->P52;++q)
    {
        ILOOP
        if(pcore->flag[q][i]>0)
        {
        j=jloc[q];

            KLOOP
            PCHECK
            {
                if(a->phi(i,j,k)>=0.0 && a->phi(i,j,k+1)<0.0)
                {
                pcore->wsf[q][i]=MAX(pcore->wsf[q][i],-(a->phi(i,j,k)*p->DZP[KP])/(a->phi(i,j,k+1)-a->phi(i,j,k)) + p->pos_z());
                pcore->loc[q][i]=p->pos_x();
				
				
                }
            }
        }
    }
	
	
    // gather, sort and write the rows
    pcore->write_rows(p,pgc,wsfout,5,
                      [&](double x){return x;},
                      p->P53==1 ? function<double(double)>([&](double x){return pflow->wave_fsf(p,pgc,x);}) : function<double(double)>(),
                      " \t ");

    if(p->mpirank==0)
    {
    wsfout.close();
    if(p->plog)
    p->plog->written(p,num,"wsfline","profiles",name,0);
    }
}

void print_wsfline_x::ini_location(lexer *p, fdm *a, ghostcell *pgc)
{
    int check,count;

    for(q=0;q<p->P52;++q)
    {
        count=0;
        ILOOP
        {
        if(p->j_dir==0)
        jloc[q]=0;
        
        if(p->j_dir==1)
        jloc[q]=p->posc_j(p->P52_y[q]);

        check=ij_boundcheck(p,i,jloc[q],0);

        if(check==1)
        pcore->flag[q][count]=1;

        ++count;
        }
    }
}

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
#include"nhflow_print_wsfline.h"
#include"wsfline_core.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include"ioflow.h"
#include"wave_interface.h"
#include<sys/stat.h>
#include<sys/types.h>
#include"runlog.h"

nhflow_print_wsfline::nhflow_print_wsfline(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
	p->Iarray(jloc,p->P52);

    pcore = new wsfline_core(p,pgc,p->P52,p->knox);

    ini_location(p,d,pgc);

	// Create Folder
	if(p->mpirank==0)
	mkdir("./REEF3D_NHFLOW_WSFLINE",0777);
}

nhflow_print_wsfline::~nhflow_print_wsfline()
{
    wsfout.close();
    delete pcore;
}

void nhflow_print_wsfline::start(lexer *p, fdm_nhf *d, ghostcell *pgc, ioflow *pflow, slice &f)
{
    char name[250];
    int num;

    num = p->count;

    if(p->mpirank==0)
    {
		// open file
		sprintf(name,"./REEF3D_NHFLOW_WSFLINE/REEF3D-NHFLOW-wsfline-%08i.dat",num);

		
		wsfout.open(name);

		wsfout<<"simtime:  "<<p->simtime<<endl;
		wsfout<<"number of wsf-lines:  "<<p->P52<<endl<<endl;
		wsfout<<"line_No     y_coord"<<endl;
		for(q=0;q<p->P52;++q)
		wsfout<<q+1<<"\t "<<p->Yout(0.0,p->P52_y[q])<<endl;

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

		wsfout<<endl;
    }

    //-------------------

    pcore->reset();

    for(q=0;q<p->P52;++q)
    {
        ILOOP
        if(pcore->flag[q][i]>0)
        {
        j=jloc[q];


            pcore->wsf[q][i]=f(i,j)+p->phimean;
            pcore->loc[q][i]=p->pos_x();
            
            //cout<<p->mpirank<<" pcore->wsf[q][i]: "<<pcore->wsf[q][i]<<" "<<f(i,j)<<" "<<p->phimean<<" "<<p->wd<<endl;
        }
    }
	
	
    // gather, sort and write the rows
    pcore->write_rows(p,pgc,wsfout,12,
                      [&](double x){return p->Xout(x,0.0);},
                      p->P53==1 ? function<double(double)>([&](double x){return pflow->wave_fsf(p,pgc,x);}) : function<double(double)>(),
                      " \t ");

    if(p->mpirank==0)
    {
    wsfout.close();
    if(p->plog)
    p->plog->written(p,num,"wsfline","profiles",name,0);
    }
}

void nhflow_print_wsfline::ini_location(lexer *p, fdm_nhf *d, ghostcell *pgc)
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


        if(jloc[q]>=0 && jloc[q]<p->knoy)
        pcore->flag[q][count]=1;

        ++count;
        }
    }
}

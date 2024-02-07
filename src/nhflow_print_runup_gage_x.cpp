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

#include<iomanip>
#include"nhflow_print_runup_gage_x.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include"ioflow.h"
#include"wave_interface.h"
#include<sys/stat.h>
#include<sys/types.h>

nhflow_print_runup_gage_x::nhflow_print_runup_gage_x(lexer *p, fdm_nhf *d, ghostcell *pgc)
{	
	p->Iarray(jloc,p->P133);

    maxknox=pgc->globalimax(p->knox);
    sumknox=pgc->globalisum(maxknox);
	
    p->Darray(xloc,p->P133+1);
    p->Darray(yloc,p->P133+1);
    p->Darray(zloc,p->P133+1);
    
    p->Darray(xloc_all,p->P133+1,p->M10+1);
    p->Darray(zloc_all,p->P133+1,p->M10+1);

    p->Iarray(flag,p->P133+1);
    p->Iarray(flag_all,p->P133+1,p->M10+1);


    ini_location(p,d,pgc);
	
	// Create Folder
	if(p->mpirank==0 && p->P14==1)
	mkdir("./REEF3D_NHFLOW_RUNUP",0777);
    
    
    
    if(p->mpirank==0)
    {
		// open file
		if(p->P14==0)
		sprintf(name,"REEF3D-NHFLOW-runup-x.dat");


		if(p->P14==1)
		sprintf(name,"./REEF3D_NHFLOW_RUNUP/REEF3D-NHFLOW-runup-x.dat");

		
		wsfout.open(name);

		wsfout<<"number of runup-probes:  "<<p->P133<<endl<<endl;
		wsfout<<"line_No     y_coord"<<endl;
		for(q=0;q<p->P133;++q)
		wsfout<<q+1<<"\t "<<p->P133_y[q]<<endl;


		wsfout<<endl<<endl;

		
		for(q=0;q<p->P133;++q)
		{
		wsfout<<"X "<<q+1;
		wsfout<<"\t P "<<q+1<<" \t \t ";
		}

		wsfout<<endl<<endl;
    }
}

nhflow_print_runup_gage_x::~nhflow_print_runup_gage_x()
{
    wsfout.close();
}

void nhflow_print_runup_gage_x::start(lexer *p, fdm_nhf *d, ghostcell *pgc, ioflow *pflow, slice &f)
{
    double zval=0.0;
    int check;

    //-------------------

    for(q=0;q<p->P133;++q)
    {
    xloc[q] = -1.0e20;
    zloc[q] = -1.0e20;
    flag[q] = 0;
    }
    
    for(q=0;q<p->P133;++q)
    for(n=0;n<p->M10;++n)
    {
    xloc_all[q][n]=0.0;
    zloc_all[q][n]=0.0;
    }

    for(q=0;q<p->P133;++q)
    {
        j=jloc[q];
        
        ILOOP
        if(p->wet[IJ]==1 && p->wet[Ip1J]==0)
        {
            if(p->XN[IP1]>xloc[q])
            {
            xloc[q] = p->XN[IP1];
            zloc[q] = 0.5*(f(i,j)+f(i+1,j))+p->phimean;
            flag[q] = 1;
            //cout<<p->mpirank<<" xloc[q]: "<<xloc[q]<<" zloc[q]: "<<zloc[q]<<endl;
            }
        }
    }
	
	
    // gather
    for(q=0;q<p->P133;++q)
    {
    pgc->gather_double(&xloc[q],1,xloc_all[q],1);
    pgc->gather_double(&zloc[q],1,zloc_all[q],1);
	pgc->gather_int(&flag[q],1,flag_all[q],1);
    }
    
    if(p->mpirank==0)
    {
        for(q=0;q<p->P133;++q)
        for(n=0;n<p->M10;++n)
        {
        //cout<<p->mpirank<<" xloc_all[q]: "<<xloc_all[q][n]<<endl;

        if(xloc_all[q][n]>xloc[q])
        {
        xloc[q] = xloc_all[q][n];
        zloc[q] = zloc_all[q][n];
        }
        }
        
    }
    
    
	
    // write to file
    if(p->mpirank==0)
    {
    wsfout<<setprecision(9)<<p->simtime<<"\t";
    for(q=0;q<p->P133;++q)
    wsfout<<setprecision(9)<<xloc[q]<<"  \t  "<<zloc[q];
    wsfout<<endl;
    }
}

void nhflow_print_runup_gage_x::ini_location(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    int check,count;
    
    
    for(q=0;q<p->P133;++q)
    {
        count=0;
        ILOOP
        {
        
        if(p->j_dir==0)
        jloc[q]=0;
        
        if(p->j_dir==1)
        jloc[q]=p->posc_j(p->P133_y[q]);

        ++count;
        }
    }
}


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
#include"bedprobe_line_x.h"
#include"wsfline_core.h"
#include"lexer.h"
#include"sediment_fdm.h"
#include"ghostcell.h"
#include"ioflow.h"
#include"wave_interface.h"
#include<sys/stat.h>
#include<sys/types.h>

bedprobe_line_x::bedprobe_line_x(lexer *p, ghostcell *pgc, sediment_fdm *s)
{	
	p->Iarray(jloc,p->P123);

    pcore = new wsfline_core(p,pgc,p->P123,p->knox);

    ini_location(p,pgc,s);
	
	// Create Folder     
    if(p->mpirank==0 && p->A10==2)
    {
	mkdir("./REEF3D_SFLOW_Sediment",0777);
    mkdir("./REEF3D_SFLOW_Sediment/Line",0777);
    }
    
    if(p->mpirank==0 && p->A10==5)
    {
	mkdir("./REEF3D_NHFLOW_Sediment",0777);
    mkdir("./REEF3D_NHFLOW_Sediment/Line",0777);
    }
    
    if(p->mpirank==0 && p->A10==6)
    {
	mkdir("./REEF3D_CFD_Sediment",0777);
    mkdir("./REEF3D_CFD_Sediment/Line",0777);
    }

}

bedprobe_line_x::~bedprobe_line_x()
{
    wsfout.close();
    delete pcore;
}

void bedprobe_line_x::start(lexer *p, ghostcell *pgc, sediment_fdm *s, ioflow *pflow)
{
	
    char name[250];
    int num;
	
    num = p->count;

    if(p->mpirank==0)
    {
		// open file
        if(p->A10==2)
        snprintf(name,sizeof(name),"./REEF3D_SFLOW_Sediment/Line/REEF3D-SFLOW-bedprobe_line_x-%06i.dat",num);
        
        if(p->A10==5)
        snprintf(name,sizeof(name),"./REEF3D_NHFLOW_Sediment/Line/REEF3D-NHFLOW-bedprobe_line_x-%06i.dat",num);
        
        if(p->A10==6)
        snprintf(name,sizeof(name),"./REEF3D_CFD_Sediment/Line/REEF3D-CFD-bedprobe_line_x-%06i.dat",num);

		wsfout.open(name);

		wsfout<<"sedtime:  "<<p->sedtime<<endl;
		wsfout<<"simtime:  "<<p->simtime<<endl;
		wsfout<<"number of topo-x-lines:  "<<p->P123<<endl<<endl;
		wsfout<<"line_No     y_coord"<<endl;
		for(q=0;q<p->P123;++q)
		wsfout<<q+1<<"\t "<<p->P123_y[q]<<endl;


		wsfout<<endl<<endl;

		
		for(q=0;q<p->P123;++q)
		{
		wsfout<<"X "<<q+1;
		wsfout<<"\t P "<<q+1<<" \t \t ";
		}

		wsfout<<endl<<endl;
    }

    //-------------------

    pcore->reset();

    for(q=0;q<p->P123;++q)
    {
        ILOOP
        if(pcore->flag[q][i]>0)
        {
        j=jloc[q];

        pcore->wsf[q][i] = MAX(pcore->wsf[q][i], s->bedzh(i,j));
        pcore->loc[q][i]=p->XP[IP];
        }
    }
	
	
    // gather, sort and write the rows
    pcore->write_rows(p,pgc,wsfout,5,
                      [&](double x){return x;},
                      function<double(double)>(),
                      " \t ");

    if(p->mpirank==0)
    {
    wsfout.close();
    }
}

void bedprobe_line_x::ini_location(lexer *p, ghostcell *pgc, sediment_fdm *s)
{
    int check,count;

    for(q=0;q<p->P123;++q)
    {
        count=0;
        ILOOP
        {
        jloc[q]=p->posc_j(p->P123_y[q]);

        check=ij_boundcheck_topo(p,i,jloc[q],0);

        if(check==1)
        pcore->flag[q][count]=1;

        ++count;
        }
    }
}



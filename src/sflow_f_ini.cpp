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

#include"sflow_f.h"
#include"lexer.h"
#include"fdm2D.h"
#include"ghostcell.h"
#include"iowave.h"
#include"sediment.h"
#include"hypre_struct2D.h"
#include"sflow_etimestep.h"
#include"sflow_weno_flux.h"
#include"sflow_eta.h"
#include"sflow_hydrostatic.h"
#include"sflow_potential.h"
#include"sflow_vtp_fsf.h"
#include"sflow_vtp_bed.h"
#include"6DOF_sflow.h"

void sflow_f::ini(lexer *p, fdm2D* b, ghostcell* pgc)
{
    p->count=0;
	p->printcount=0;

    int count=0;
	p->pointnum2D=0;
	p->cellnum2D=0;
	p->polygon_sum=0;


    TPSLICELOOP
	{
	++count;
	++p->pointnum2D;
	b->nodeval(i,j)=count;
    }

	SLICEBASELOOP
	++p->polygon_sum;

	p->polygon_sum*=2;

	SLICELOOP4
	++p->cellnum2D;

    p->cellnumtot2D=pgc->globalisum(p->cellnum2D);


    if(p->mpirank==0)
    cout<<"number of cells: "<<p->cellnumtot2D<<endl;


    if(p->mpirank==0)
    cout<<"starting driver_ini"<<endl;

    ptime->ini(p,b ,pgc);
    
    
    // bed ini
	ILOOP
    JLOOP
	b->bed(i,j) = p->bed[IJ];
    
    ILOOP
    JLOOP
	b->solidbed(i,j) = p->solidbed[IJ];

    ILOOP
    JLOOP
	b->topobed(i,j) = p->topobed[IJ];


    pgc->gcsl_start4(p,b->bed,50);
    pgc->gcsl_start4(p,b->solidbed,50);
    pgc->gcsl_start4(p,b->topobed,50);

    
	for(int qn=0; qn<p->A209;++qn)
    {
	SLICELOOP4
	b->bed(i,j) = b->bed0(i,j) = 0.5*b->bed(i,j) + 0.125*(b->bed(i-1,j) +b->bed(i+1,j) +b->bed(i,j-1) +b->bed(i,j+1) );
    
    pgc->gcsl_start4(p,b->bed,50);
    }
    

    // FSF ini
    ini_fsf(p,b,pgc);

    SLICELOOP4
    p->wet[IJ]=1;


    SLICELOOP4
    b->breaking(i,j)=0;

	pgc->gcsl_start4(p,b->depth,50);

	SLICELOOP4
	b->ws(i,j) = 0.0;


    pgc->gcsl_start4(p,b->eta_n,50);

    SLICELOOP4
	b->hp(i,j) = MAX(b->eta(i,j) + p->wd - b->bed(i,j),0.0);

     pflow->ini2D(p,b,pgc);
     
     // P,Q ini
	pflow->um_relax(p,pgc,b->P,b->bed,b->eta);
	pflow->vm_relax(p,pgc,b->Q,b->bed,b->eta);

	pgc->gcsl_start1(p,b->P,10);
	pgc->gcsl_start2(p,b->Q,11);

	pgc->gcsl_start1(p,b->P,10);
	pgc->gcsl_start2(p,b->Q,11);
	pgc->gcsl_start4(p,b->eta,50);
    pgc->gcsl_start4(p,b->hp,50);
    pgc->gcsl_start4(p,b->bed,50);
    
    pfsf->depth_update(p,b,pgc,b->P,b->Q,b->ws,b->eta);
    
    // potential flow ini
    potflow->start(p,b,ppoissonsolv,pgc);
    
    // FSF ini
    ini_fsf_2(p,b,pgc);

    pfsf->depth_update(p,b,pgc,b->P,b->Q,b->ws,b->eta);

    //roughness ini
    SLICELOOP4
    b->ks(i,j) = p->B50;

    pgc->gcsl_start4(p,b->ks,50);

    //sediment ini
    psed->ini_sflow(p,b,pgc);

    //6DOF ini
    p6dof_sflow->ini(p,pgc);

    // print
    log_ini(p);
	print_debug(p,b,pgc);
    pprint->start(p,b,pgc,pflow,pturb,psed);

	pprintbed->start(p,b,pgc,psed);
}

void sflow_f::ini_fsf(lexer *p, fdm2D* b, ghostcell* pgc)
{
    int istart,iend,jstart,jend;
    p->phimean = 0.0;
    p->phiout = 0.0; 
    
    // depth
    if(p->F60>-1.0e20)
    {
    p->phimean=p->F60;
    p->phiout=p->F60;
    p->wd=p->F60;
    }

    // eta plain
    SLICELOOP4
    b->eta(i,j)=0.0;

    // eta slope
    if(p->A251==1)
    SLICELOOP4
    {
    b->eta(i,j)= -p->A251_val*(p->XP[IP]-p->global_xmin);
    }

    // eta box area
    for(int qn=0;qn<p->F72;++qn)
    {
		istart = p->posc_i(p->F72_xs[qn]);
        iend = p->posc_i(p->F72_xe[qn]);

        jstart = p->posc_j(p->F72_ys[qn]);
        jend = p->posc_j(p->F72_ye[qn]);

        SLICELOOP4
        if(i>=istart && i<iend && j>=jstart && j<jend)
        b->eta(i,j) = p->F72_h[qn]-p->wd;
	}

        // fix inflow fsf
        for(n=0;n<p->gcslin_count;n++)
        {
        i=p->gcslin[n][0];
        j=p->gcslin[n][1];

        b->eta(i-1,j) = 0.0;
        b->eta(i-2,j) = 0.0;
        b->eta(i-3,j) = 0.0;

        b->hp(i-1,j) = MAX(b->eta(i-1,j) + p->wd - b->bed(i,j),0.0);
        b->hp(i-2,j) = MAX(b->eta(i-2,j) + p->wd - b->bed(i,j),0.0);
        b->hp(i-3,j) = MAX(b->eta(i-3,j) + p->wd - b->bed(i,j),0.0);
        }
        
        // fix outflow fsf
        for(n=0;n<p->gcslout_count;n++)
        {
        i=p->gcslout[n][0];
        j=p->gcslout[n][1];

        b->eta(i+1,j) = 0.0;
        b->eta(i+2,j) = 0.0;
        b->eta(i+3,j) = 0.0;

        b->hp(i+1,j) = MAX(b->eta(i+1,j) + p->wd - b->bed(i,j),0.0);
        b->hp(i+2,j) = MAX(b->eta(i+2,j) + p->wd - b->bed(i,j),0.0);
        b->hp(i+3,j) = MAX(b->eta(i+3,j) + p->wd - b->bed(i,j),0.0);
        }
    

        GCSL1LOOP
        {
        i = p->gcbsl1[n][0];
        j = p->gcbsl1[n][1];
                
        
            if(p->gcbsl1[n][4]==1)
            {
            b->hx(i-1,j) = MAX(p->wd - b->bed(i,j),0.0);
            b->hx(i-2,j) = MAX(p->wd - b->bed(i,j),0.0);
            b->hx(i-3,j) = MAX(p->wd - b->bed(i,j),0.0);
            }
            
            if(p->gcbsl1[n][4]==2)
            {
            b->hx(i+1,j) = MAX(p->wd - b->bed(i,j),0.0);
            b->hx(i+2,j) = MAX(p->wd - b->bed(i,j),0.0);
            b->hx(i+3,j) = MAX(p->wd - b->bed(i,j),0.0);
            }
        }

        GCSL2LOOP
        {
        i = p->gcbsl2[n][0];
        j = p->gcbsl2[n][1];

            if(p->gcbsl2[n][4]==1)
            {
            b->hy(i-1,j) = MAX(p->wd - b->bed(i,j),0.0);
            b->hy(i-2,j) = MAX(p->wd - b->bed(i,j),0.0);
            b->hy(i-3,j) = MAX(p->wd - b->bed(i,j),0.0);
            }
            
            if(p->gcbsl2[n][4]==2)
            {
            b->hy(i+1,j) = MAX(p->wd - b->bed(i,j),0.0);
            b->hy(i+2,j) = MAX(p->wd - b->bed(i,j),0.0);
            b->hy(i+3,j) = MAX(p->wd - b->bed(i,j),0.0);
            }
        }
        
        
    if(p->F61>-1.0e20)
    {
        GCSL1LOOP
        {
        i = p->gcbsl1[n][0];
        j = p->gcbsl1[n][1];
                
        
            if(p->gcbsl1[n][4]==1)
            {
            b->hx(i-1,j) = MAX(p->F61 - b->bed(i,j),0.0);
            b->hx(i-2,j) = MAX(p->F61 - b->bed(i,j),0.0);
            b->hx(i-3,j) = MAX(p->F61 - b->bed(i,j),0.0);
            }
        }

        GCSL2LOOP
        {
        i = p->gcbsl2[n][0];
        j = p->gcbsl2[n][1];

            if(p->gcbsl2[n][4]==1)
            {
            b->hy(i-1,j) = MAX(p->F61 - b->bed(i,j),0.0);
            b->hy(i-2,j) = MAX(p->F61 - b->bed(i,j),0.0);
            b->hy(i-3,j) = MAX(p->F61 - b->bed(i,j),0.0);
            }
        }
        
    }
    
    
    
    if(p->F62>-1.0e20)
    {
        GCSL1LOOP
        {
        i = p->gcbsl1[n][0];
        j = p->gcbsl1[n][1];
                
            if(p->gcbsl1[n][4]==2)
            {
            b->hx(i+1,j) = MAX(p->F62 - b->bed(i,j),0.0);
            b->hx(i+2,j) = MAX(p->F62 - b->bed(i,j),0.0);
            b->hx(i+3,j) = MAX(p->F62 - b->bed(i,j),0.0);
            }
        }

        GCSL2LOOP
        {
        i = p->gcbsl2[n][0];
        j = p->gcbsl2[n][1];
            
            if(p->gcbsl2[n][4]==2)
            {
            b->hy(i+1,j) = MAX(p->F62 - b->bed(i,j),0.0);
            b->hy(i+2,j) = MAX(p->F62 - b->bed(i,j),0.0);
            b->hy(i+3,j) = MAX(p->F62 - b->bed(i,j),0.0);
            }
        }
        
        for(n=0;n<p->gcslout_count;n++)
        {
        i=p->gcslout[n][0];
        j=p->gcslout[n][1];
        
        b->eta(i,j) = p->F62-p->wd;
        b->eta(i+1,j) = p->F62-p->wd;
        b->eta(i+2,j) = p->F62-p->wd;
        b->eta(i+3,j) = p->F62-p->wd;

        b->hp(i,j) = MAX(b->eta(i+1,j) + p->wd - b->bed(i,j),0.0);
        b->hp(i+1,j) = MAX(b->eta(i+1,j) + p->wd - b->bed(i,j),0.0);
        b->hp(i+2,j) = MAX(b->eta(i+2,j) + p->wd - b->bed(i,j),0.0);
        b->hp(i+3,j) = MAX(b->eta(i+3,j) + p->wd - b->bed(i,j),0.0);
        }
    }
      
	pfsf->depth_update(p,b,pgc,b->P,b->Q,b->ws,b->eta);

}

void sflow_f::ini_fsf_2(lexer *p, fdm2D* b, ghostcell* pgc)
{
    // eta ini
	pflow->eta_relax(p,pgc,b->eta);
    pgc->gcsl_start4(p,b->eta,50);

    // eta_n ini
    SLICELOOP4
    b->eta_n(i,j) = b->eta(i,j);
    
    pgc->gcsl_start4(p,b->eta_n,50);
}

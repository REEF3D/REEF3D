/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2022 Hans Bihs

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

#include"driver.h"
#include"ghostcell.h"
#include"lexer.h"
#include"fdm.h"

void driver::driver_ini_nsewave()
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
    a->eta(i,j)=0.0;

    // eta slope
    if(p->A251==1)
    SLICELOOP4
    {
    a->eta(i,j)= -p->A251_val*(p->XP[IP]-p->global_xmin);
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
        a->eta(i,j) = p->F72_h[qn]-p->wd;
	}
/*
        // fix inflow fsf
        for(n=0;n<p->gcslin_count;n++)
        {
        i=p->gcslin[n][0];
        j=p->gcslin[n][1];

        a->eta(i-1,j) = 0.0;
        a->eta(i-2,j) = 0.0;
        a->eta(i-3,j) = 0.0;

        a->hp(i-1,j) = MAX(a->eta(i-1,j) + p->wd - a->bed(i,j),0.0);
        a->hp(i-2,j) = MAX(a->eta(i-2,j) + p->wd - a->bed(i,j),0.0);
        a->hp(i-3,j) = MAX(a->eta(i-3,j) + p->wd - a->bed(i,j),0.0);
        }
        
        // fix outflow fsf
        for(n=0;n<p->gcslout_count;n++)
        {
        i=p->gcslout[n][0];
        j=p->gcslout[n][1];

        a->eta(i+1,j) = 0.0;
        a->eta(i+2,j) = 0.0;
        a->eta(i+3,j) = 0.0;

        a->hp(i+1,j) = MAX(a->eta(i+1,j) + p->wd - a->bed(i,j),0.0);
        a->hp(i+2,j) = MAX(a->eta(i+2,j) + p->wd - a->bed(i,j),0.0);
        a->hp(i+3,j) = MAX(a->eta(i+3,j) + p->wd - a->bed(i,j),0.0);
        }
    

        GCSL1LOOP
        {
        i = p->gcbsl1[n][0];
        j = p->gcbsl1[n][1];
                
        
            if(p->gcbsl1[n][4]==1)
            {
            a->hx(i-1,j) = MAX(p->wd - a->bed(i,j),0.0);
            a->hx(i-2,j) = MAX(p->wd - a->bed(i,j),0.0);
            a->hx(i-3,j) = MAX(p->wd - a->bed(i,j),0.0);
            }
            
            if(p->gcbsl1[n][4]==2)
            {
            a->hx(i+1,j) = MAX(p->wd - a->bed(i,j),0.0);
            a->hx(i+2,j) = MAX(p->wd - a->bed(i,j),0.0);
            a->hx(i+3,j) = MAX(p->wd - a->bed(i,j),0.0);
            }
        }

        GCSL2LOOP
        {
        i = p->gcbsl2[n][0];
        j = p->gcbsl2[n][1];

            if(p->gcbsl2[n][4]==1)
            {
            a->hy(i-1,j) = MAX(p->wd - a->bed(i,j),0.0);
            a->hy(i-2,j) = MAX(p->wd - a->bed(i,j),0.0);
            a->hy(i-3,j) = MAX(p->wd - a->bed(i,j),0.0);
            }
            
            if(p->gcbsl2[n][4]==2)
            {
            a->hy(i+1,j) = MAX(p->wd - a->bed(i,j),0.0);
            a->hy(i+2,j) = MAX(p->wd - a->bed(i,j),0.0);
            a->hy(i+3,j) = MAX(p->wd - a->bed(i,j),0.0);
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
            a->hx(i-1,j) = MAX(p->F61 - a->bed(i,j),0.0);
            a->hx(i-2,j) = MAX(p->F61 - a->bed(i,j),0.0);
            a->hx(i-3,j) = MAX(p->F61 - a->bed(i,j),0.0);
            }
        }

        GCSL2LOOP
        {
        i = p->gcbsl2[n][0];
        j = p->gcbsl2[n][1];

            if(p->gcbsl2[n][4]==1)
            {
            a->hy(i-1,j) = MAX(p->F61 - a->bed(i,j),0.0);
            a->hy(i-2,j) = MAX(p->F61 - a->bed(i,j),0.0);
            a->hy(i-3,j) = MAX(p->F61 - a->bed(i,j),0.0);
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
            a->hx(i+1,j) = MAX(p->F62 - a->bed(i,j),0.0);
            a->hx(i+2,j) = MAX(p->F62 - a->bed(i,j),0.0);
            a->hx(i+3,j) = MAX(p->F62 - a->bed(i,j),0.0);
            }
        }

        GCSL2LOOP
        {
        i = p->gcbsl2[n][0];
        j = p->gcbsl2[n][1];
            
            if(p->gcbsl2[n][4]==2)
            {
            a->hy(i+1,j) = MAX(p->F62 - a->bed(i,j),0.0);
            a->hy(i+2,j) = MAX(p->F62 - a->bed(i,j),0.0);
            a->hy(i+3,j) = MAX(p->F62 - a->bed(i,j),0.0);
            }
        }
    }
      
	//pfsf->depth_update(p,b,pgc,a->P,a->Q,a->ws,a->eta);*/
    
    int gcval_phi;
    
    if(p->F50==1)
	gcval_phi=51;

	if(p->F50==2)
	gcval_phi=52;

	if(p->F50==3)
	gcval_phi=53;

	if(p->F50==4)
	gcval_phi=54;
    
    LOOP
    a->phi(i,j,k) = a->eta(i,j) + p->phimean - p->pos_z();
    
    pgc->start4(p,a->phi,gcval_phi);

}
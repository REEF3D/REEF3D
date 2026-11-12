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

#include"lexer.h"
#include"ghostcell.h"
#include<vector>
#include<cmath>

static void periodic_ghost_nodes(lexer*, ghostcell*);

void lexer::gridini(ghostcell *pgc)
{
    if(G502==1)
    grid::sigma_coord_ini();
    
    periodic_ghost_nodes(this,pgc);

    grid::gridspacing(pgc);

    // horizontal geometry layer: 2D nodes and metrics (not read by the solvers yet)
    geometry_ini(pgc);

    gcd_ini(pgc);
    
    // solids and topography of the grid file: fields, bed levels, ghost cell estimates
    grid_solids(pgc);

    // curvilinear grid (G 50): bed levels of its nodes
    geometry_bed();
}

void lexer::flagini()
{
    control_calc();
	gridsize();

	Iarray(flag1,imax*jmax*kmax);
	Iarray(flag2,imax*jmax*kmax);
	Iarray(flag3,imax*jmax*kmax);
    Iarray(flag5,imax*jmax*kmax);

    Iarray(flagsf1,imax*jmax*kmax);
	Iarray(flagsf2,imax*jmax*kmax);
	Iarray(flagsf3,imax*jmax*kmax);
	Iarray(flagsf4,imax*jmax*kmax);

    for(i=0; i<knox; ++i)
    for(j=0; j<knoy; ++j)
    for(k=0; k<knoz; ++k)
    {
    flagsf1[(i-imin)*jmax*kmax + (j-jmin)*kmax + k-kmin]=1;
    flagsf2[(i-imin)*jmax*kmax + (j-jmin)*kmax + k-kmin]=1;
    flagsf3[(i-imin)*jmax*kmax + (j-jmin)*kmax + k-kmin]=1;
    flagsf4[(i-imin)*jmax*kmax + (j-jmin)*kmax + k-kmin]=1;
    }

    // boundary conditions
    Iarray(IO,imax*jmax*kmax);
    Iarray(IOSL,imax*jmax);
    Iarray(DF,imax*jmax*kmax);
    Iarray(DF1,imax*jmax*kmax);
    Iarray(DF2,imax*jmax*kmax);
    Iarray(DF3,imax*jmax*kmax);
    Iarray(DFF,imax*jmax*(kmax+1));

    Iarray(DFBED,imax*jmax);

    for(i=-margin; i<knox+margin; ++i)
    for(j=-margin; j<knoy+margin; ++j)
    DFBED[(i-imin)*jmax + j-jmin] = 1;

    // flag
	makeflag(flag1);
	makeflag(flag2);
	makeflag(flag3);

    for(i=-margin; i<knox+margin; ++i)
    for(j=-margin; j<knoy+margin; ++j)
    for(k=-margin; k<knoz+margin; ++k)
    IO[(i-imin)*jmax*kmax + (j-jmin)*kmax + k-kmin] = 0;

    for(i=-margin; i<knox+margin; ++i)
    for(j=-margin; j<knoy+margin; ++j)
    IOSL[(i-imin)*jmax + j-jmin] = 0;

    for(i=-margin; i<knox+margin; ++i)
    for(j=-margin; j<knoy+margin; ++j)
    for(k=-margin; k<knoz+margin; ++k)
    DF[(i-imin)*jmax*kmax + (j-jmin)*kmax + k-kmin] = 1;

	x_dir=y_dir=z_dir=1.0;

	if(i_dir==0)
	x_dir=0.0;

	if(j_dir==0)
	y_dir=0.0;

	if(k_dir==0)
	z_dir=0.0;


	if(B98>=3)
	for(n=0;n<gcb4_count;++n)
	if(gcb4[n][4]==6)
	gcb4[n][4]=1;

    // gcdf
    gcdf1_count=gcdf2_count=gcdf3_count=gcdf4_count=1;

    Iarray(gcdf1,gcdf1_count,6);
    Iarray(gcdf2,gcdf2_count,6);
    Iarray(gcdf3,gcdf3_count,6);
    Iarray(gcdf4,gcdf4_count,6);

    // gcsldf
    gcsldfeta4_count=1;
    gcsldfbed4_count=1;

    Iarray(gcsldfeta4,gcsldfeta4_count,6);
    Iarray(gcsldfbed4,gcsldfbed4_count,6);
}

int lexer::conv(double a)
{
	int b,c;
	double d,diff;

	c= int( a);
	d=double(c);
	diff=a-d;

	b=c;

	if(diff>0.5)
	b=c+1;

	if(diff<=-0.5)
	b=c-1;

	return b;
}

void lexer::gcd_ini(ghostcell *pgc)
{
    for(int q=0;q<gcb4_count;++q)
    {
        if(gcb4[q][3]==1 || gcb4[q][3]==4)
        {
            i=gcb4[q][0];
            gcd4[q] = 0.5*DXP[IP];
        }
        else if(gcb4[q][3]==2 || gcb4[q][3]==3)
        {
            j=gcb4[q][1];
            gcd4[q] = 0.5*DYP[JP];
        }
        else if(gcb4[q][3]==5 || gcb4[q][3]==6)
        {
            k=gcb4[q][2];
            gcd4[q] = 0.5*DZP[KP];
        }
    }
}

// periodic sides (DIVEMesh C 21, C 22): the ghost nodes beyond a periodic side are the nodes of the cells at the
// opposite end. The grid file extrapolates them with the edge cell of the own side, so where the two ends have
// different cell sizes the ghost cells across the periodic side did not have the size of the cells they mirror
// (a channel with 2.5 mm cells at one end and 5 mm at the other, validation 14). Ends of equal cell size keep the
// nodes of the grid file.
static void periodic_ghost_nodes(lexer *p, ghostcell *pgc)
{
    auto side = [&](double *N, int kno, int origin, int gkno)
    {
        const int m = p->marge;
        const bool lo = origin==0;
        const bool hi = origin+kno==gkno;
        
        // the cells at the two ends of the domain: the first m (low end) and the last m (high end)
        std::vector<double> hl(m+1,0.0), hh(m+1,0.0);
        
        for(int q=1;q<=m && q<=kno;++q)
        {
            if(lo)
            hl[q] = N[q+p->marge] - N[q-1+p->marge];
            
            if(hi)
            hh[q] = N[kno-q+1+p->marge] - N[kno-q+p->marge];
        }
        
        for(int q=1;q<=m;++q)
        {
            hl[q] = pgc->globalmax(hl[q]);
            hh[q] = pgc->globalmax(hh[q]);
        }
        
        // ghost nodes from the cells of the opposite end, applied only where the two ends differ
        double dev=0.0, hmin=1.0e20;
        std::vector<double> gl(m+1), gh(m+1);
        
        gl[0] = N[0+p->marge];
        gh[0] = N[kno+p->marge];
        
        for(int q=1;q<=m;++q)
        {
            gl[q] = gl[q-1] - hh[q];
            gh[q] = gh[q-1] + hl[q];
            hmin = MIN(hmin, MIN(hh[q], hl[q]));
            
            if(lo)
            dev = MAX(dev, fabs(gl[q] - N[-q+p->marge]));
            
            if(hi)
            dev = MAX(dev, fabs(gh[q] - N[kno+q+p->marge]));
        }
        
        dev = pgc->globalmax(dev);
        
        if(hmin<=0.0 || dev <= 1.0e-6*hmin)
        return;
        
        if(lo)
        for(int q=1;q<=m;++q)
        N[-q+p->marge] = gl[q];
        
        if(hi)
        for(int q=1;q<=m;++q)
        N[kno+q+p->marge] = gh[q];
        
        if(p->mpirank==0)
        cout<<"periodic side with different end cells: ghost nodes from the opposite end"<<endl;
    };
    
    if(p->periodic1>0)
    side(p->XN,p->knox,p->origin_i,p->gknox);
    
    if(p->periodic2>0 && p->j_dir==1)
    side(p->YN,p->knoy,p->origin_j,p->gknoy);
}

// FNPF reads the 3D cell flags flag4, flag5 and flag7 in its time steps (and the printers). The flags
// of the CFD velocity components (flag1-3), of the solid forcing (flagsf1-4, DF, DF1-3, DFF) and IO
// are read by the grid set-up only. For FNPF they are given back after the initialisation: 13 int
// arrays of the 3D grid (with the ghost cells), about 0.65 GB each for 1e8 cells. Any later access
// fails at once (nullptr).
void lexer::flagfree_fnpf()
{
    const int n3 = imax*jmax*kmax;
    const int nF = imax*jmax*(kmax+1);
    
    del_Iarray(flag1,n3);
    del_Iarray(flag2,n3);
    del_Iarray(flag3,n3);
    
    del_Iarray(flagsf1,n3);
    del_Iarray(flagsf2,n3);
    del_Iarray(flagsf3,n3);
    del_Iarray(flagsf4,n3);
    
    del_Iarray(IO,n3);
    del_Iarray(DF,n3);
    del_Iarray(DF1,n3);
    del_Iarray(DF2,n3);
    del_Iarray(DF3,n3);
    del_Iarray(DFF,nF);
}

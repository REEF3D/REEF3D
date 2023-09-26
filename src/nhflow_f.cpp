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
#include"nhflow_f.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"

nhflow_f::nhflow_f(lexer *p, fdm_nhf *d, ghostcell *pgc) 
{
    margin=3;
}

nhflow_f::~nhflow_f()
{
}

void nhflow_f::ini(lexer *p, fdm_nhf *d, ghostcell *pgc, ioflow *pflow)
{
    // count cells 3D
    int count=0;
    p->count=0;
	p->pointnum=0;
	p->cellnum=0;
	p->tpcellnum=0;

	TPLOOP
	{
	++count;
	++p->pointnum;
    d->NODEVAL[IJK]=count;
	}

	LOOP
	++p->cellnum;
    
    LOOP
    ++p->tpcellnum;
    
    p->count=0;
    
    // 2D mesh
    count=0;
	p->pointnum2D=0;
	p->cellnum2D=0;
	p->polygon_sum=0;
    
   
    TPSLICELOOP  
	{
	++count;
	++p->pointnum2D;
	d->nodeval2D(i,j)=count;
    }
	
	SLICEBASELOOP
	++p->polygon_sum;
	
	p->polygon_sum *=2;

	SLICELOOP4
	++p->cellnum2D;
    
    SLICELOOP4
	++p->cellnum2D;

    p->cellnumtot2D=pgc->globalisum(p->cellnum2D);
    
// --
    p->cellnumtot=pgc->globalisum(p->cellnum);
    p->pointnumtot=pgc->globalisum(p->pointnum);


    if(p->mpirank==0)
    cout<<"number of cells: "<<p->cellnumtot<<endl;
    
    // SIGMA grid
    // bed ini
    SLICELOOP4
	d->bed(i,j) = p->bed[IJ];
    
    pgc->gcsl_start4(p,d->bed,50);
    
    SLICELOOP4
	d->depth(i,j) = p->wd - d->bed(i,j);
    
    pgc->gcsl_start4(p,d->depth,50);
    
    // eta ini
	SLICELOOP4
    {
	d->eta(i,j) = 0.0;
    p->deep[IJ] = 1;
    p->wet[IJ]=1;
    d->breaking(i,j)=0;
    }

    pgc->gcsl_start4Vint(p,p->wet,50);
    pgc->gcsl_start4Vint(p,p->deep,50);
    pgc->gcsl_start4(p,d->eta,50);

    
    ALOOP
    d->porosity[IJK]=1.0;
    
    pgc->start4V(p,d->porosity,1);
    
    SLICELOOP4
    d->WL(i,j) = d->eta(i,j) + d->depth(i,j);
    
    SLICELOOP4
    if(d->WL(i,j)<p->A544)
    p->wet[IJ]=0;
    
    SLICELOOP4
    d->WL(i,j) = MAX(p->A544,d->eta(i,j) + d->depth(i,j));

    pgc->gcsl_start4(p,d->WL,50);
}


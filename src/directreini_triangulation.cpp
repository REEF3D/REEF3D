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

#include"directreini.h"
#include"lexer.h"
#include"fdm.h"

void directreini::triangulation(lexer *p,fdm* a, field& b, fieldint& nodeflag, fieldint& vertice)
{

    LOOP
    {
    vertice(i,j,k)=-1;

    vertice(i+1,j,k)=-1;
    vertice(i+1,j+1,k)=-1;
    vertice(i,j+1,k)=-1;
    vertice(i,j,k+1)=-1;
    vertice(i+1,j,k+1)=-1;
    vertice(i+1,j+1,k+1)=-1;
    vertice(i,j+1,k+1)=-1;

    vertice(i-1,j,k)=-1;
    vertice(i-1,j-1,k)=-1;
    vertice(i,j-1,k)=-1;
    vertice(i,j,k-1)=-1;
    vertice(i-1,j,k-1)=-1;
    vertice(i-1,j-1,k-1)=-1;
    vertice(i,j-1,k-1)=-1;
    }

    LOOP
    {
    nodeflag(i,j,k)=0;

    nodeflag(i+1,j,k)=0;
    nodeflag(i+1,j+1,k)=0;
    nodeflag(i,j+1,k)=0;
    nodeflag(i,j,k+1)=0;
    nodeflag(i+1,j,k+1)=0;
    nodeflag(i+1,j+1,k+1)=0;
    nodeflag(i,j+1,k+1)=0;

    nodeflag(i-1,j,k)=0;
    nodeflag(i-1,j-1,k)=0;
    nodeflag(i,j-1,k)=0;
    nodeflag(i,j,k-1)=0;
    nodeflag(i-1,j,k-1)=0;
    nodeflag(i-1,j-1,k-1)=0;
    nodeflag(i,j-1,k-1)=0;
    }

    LOOP
    {
        check=1;

        if(b(i,j,k)<zero && b(i+1,j,k)<zero && b(i+1,j+1,k)<zero && b(i,j+1,k)<zero &&
           b(i,j,k+1)<zero && b(i+1,j,k+1)<zero && b(i+1,j+1,k+1)<zero && b(i,j+1,k+1)<zero)
        check=0;
		
		if(b(i,j,k)>zero && b(i+1,j,k)>zero && b(i+1,j+1,k)>zero && b(i,j+1,k)>zero &&
           b(i,j,k+1)>zero && b(i+1,j,k+1)>zero && b(i+1,j+1,k+1)>zero && b(i,j+1,k+1)>zero)
        check=0;

        if(check==1)
        {
        nodeflag(i,j,k)=1;
        nodeflag(i+1,j,k)=1;
        nodeflag(i+1,j+1,k)=1;
        nodeflag(i,j+1,k)=1;
        nodeflag(i,j,k+1)=1;
        nodeflag(i+1,j,k+1)=1;
        nodeflag(i+1,j+1,k+1)=1;
        nodeflag(i,j+1,k+1)=1;
        }
    }

    LOOP
    {
        check=1;

        if(b(i,j,k)<zero && b(i-1,j,k)<zero && b(i-1,j-1,k)<zero && b(i,j-1,k)<zero &&
           b(i,j,k-1)<zero && b(i-1,j,k-1)<zero && b(i-1,j-1,k-1)<zero && b(i,j-1,k-1)<zero)
        check=0;
		
		if(b(i,j,k)>zero && b(i-1,j,k)>zero && b(i-1,j-1,k)>zero && b(i,j-1,k)>zero &&
           b(i,j,k-1)>zero && b(i-1,j,k-1)>zero && b(i-1,j-1,k-1)>zero && b(i,j-1,k-1)>zero)
        check=0;

        if(check==1)
        {
        nodeflag(i,j,k)=1;
        nodeflag(i-1,j,k)=1;
        nodeflag(i-1,j-1,k)=1;
        nodeflag(i,j-1,k)=1;
        nodeflag(i,j,k-1)=1;
        nodeflag(i-1,j,k-1)=1;
        nodeflag(i-1,j-1,k-1)=1;
        nodeflag(i,j-1,k-1)=1;
        }
    }

    countM=0;
    BLOOP
    if(nodeflag(i,j,k)==1)
    ++countM;

    numtri = 6*(countM);
    numvert = countM;

    numtri_mem = numtri;
    numvert_mem = numvert;

    Iarray(tri,numtri,4);

    Darray(pt,numvert,3);
    Iarray(ijk,numvert,3);
    Darray(ls,numvert);
	Darray(ls0,numvert);
	Darray(ls1,numvert);
	Darray(lsvert,numvert);
	Darray(lsfac,numvert);
	Iarray(reiniflag,numvert);

    Iarray(facet,numtri,4);
    Iarray(confac,numtri);
    Iarray(numfac,numtri);
    Darray(ccpt,numtri*4,3);


    countM=0;
    BLOOP
    if(nodeflag(i,j,k)==1)
    {
    pt[countM][0] = p->pos_x();
    pt[countM][1] = p->pos_y();
    pt[countM][2] = p->pos_z();

    ijk[countM][0] = i;
    ijk[countM][1] = j;
    ijk[countM][2] = k;

    ls[countM] = b(i,j,k);
	ls0[countM] = b(i,j,k);
	lsvert[countM] = b(i,j,k);
	lsfac[countM] = b(i,j,k);

    vertice(i,j,k) = countM;

    ++countM;
    }

	// p. 725, 956
    count=0;
    countM=0;
    TPLOOP
    if(nodeflag(i,j,k)==1)
    if(nodeflag(i+1,j,k)==1)
    if(nodeflag(i+1,j+1,k)==1)
    if(nodeflag(i,j+1,k)==1)
    if(nodeflag(i,j,k+1)==1)
    if(nodeflag(i+1,j,k+1)==1)
    if(nodeflag(i+1,j+1,k+1)==1)
    if(nodeflag(i,j+1,k+1)==1)
    {
		
    // 1
    tri[count][0] = vertice(i,j,k);
    tri[count][1] = vertice(i,j+1,k);
    tri[count][2] = vertice(i,j,k+1);
    tri[count][3] = vertice(i+1,j,k+1);
    ++count;

    // 2
    tri[count][0] = vertice(i,j,k);
    tri[count][1] = vertice(i+1,j,k);
    tri[count][2] = vertice(i,j+1,k);
    tri[count][3] = vertice(i+1,j,k+1);
    ++count;

    // 3
    tri[count][0] = vertice(i,j+1,k);
    tri[count][1] = vertice(i+1,j+1,k);
    tri[count][2] = vertice(i+1,j,k);
    tri[count][3] = vertice(i+1,j,k+1);
    ++count;

    // 4
    tri[count][0] = vertice(i+1,j+1,k);
    tri[count][1] = vertice(i,j+1,k);
    tri[count][2] = vertice(i+1,j,k+1);
    tri[count][3] = vertice(i+1,j+1,k+1);
    ++count;

    // 5
	tri[count][0] = vertice(i,j+1,k);
    tri[count][1] = vertice(i,j+1,k+1);
    tri[count][2] = vertice(i+1,j+1,k+1);
    tri[count][3] = vertice(i+1,j,k+1);
    ++count;

    // 6
    tri[count][0] = vertice(i,j+1,k);
    tri[count][1] = vertice(i,j,k+1);
    tri[count][2] = vertice(i+1,j,k+1);
    tri[count][3] = vertice(i,j+1,k+1);
    ++count;
    }
    numtri=count;
    //cout<<p->mpirank<<"  COUNTnumtri:  "<<count<<endl;
}


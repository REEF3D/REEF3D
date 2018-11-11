/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2018 Hans Bihs

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

#include"fsf_vtp.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

void fsf_vtp::triangulation(lexer *p,fdm* a, ghostcell *pgc, field& f) 
{
	int negcount, poscount;
	nodefill4(p,a,pgc,f,eta);

	
    NDBASELOOP
    vertice(i,j,k)=-1;

    NDBASELOOP
    nodeflag(i,j,k)=0;
	
	for(n=0;n<p->facetnum;n++)
	{
		i=p->facet[n][0];
		j=p->facet[n][1];
		k=p->facet[n][2];
		
		vertice(i,j,k)=1;
	}

    BASELOOP
    {
        epsi=interfac*(1.0/3.0)*(p->DXN[IP] + p->DYN[JP] + p->DZN[KP]);
        
        if(vertice(i,j,k)<0 && fabs(a->phi(i,j,k))<epsi)
        {
            check=1;

            if(eta(i,j,k)<zero && eta(i-1,j,k)<zero && eta(i-1,j-1,k)<zero && eta(i,j-1,k)<zero &&
               eta(i,j,k-1)<zero && eta(i-1,j,k-1)<zero && eta(i-1,j-1,k-1)<zero && eta(i,j-1,k-1)<zero)
            check=0;
            
            if(eta(i,j,k)>zero && eta(i-1,j,k)>zero && eta(i-1,j-1,k)>zero && eta(i,j-1,k)>zero &&
               eta(i,j,k-1)>zero && eta(i-1,j,k-1)>zero && eta(i-1,j-1,k-1)>zero && eta(i,j-1,k-1)>zero)
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
    }
	
	NDBASELOOP
    vertice(i,j,k)=-1;
	
	// CC
	
	for(n=0;n<p->ccellnum;++n)
	for(q=0;q<ccnode[n];++q)
	{
		if(ccid[n][q]==0)	
		lscc[n][q] = p->ccipol4(a->phi,ccell[n][q][0],ccell[n][q][1],ccell[n][q][2]);
		
		if(ccid[n][q]==1)
		{
		i = ccijk[n][q][0];
		j = ccijk[n][q][1];
		k = ccijk[n][q][2];
		lscc[n][q] = eta(i,j,k);
		}
	}
	
	for(n=0;n<p->ccellnum;++n)
	ccflag[n]=0;
	
	countCC=0;
	for(n=0;n<p->ccellnum;++n)
	{
		check=0;
		negcount=0;
		poscount=0;
		for(q=0;q<ccnode[n];++q)
		{
			if(lscc[n][q]<zero)
			++negcount;
			
			if(lscc[n][q]>zero)
			++poscount;
		}
		if(negcount!=ccnode[n] && poscount!=ccnode[n])
		{
		++countCC;
		ccflag[n]=1;
		check=1;
		}

	}
	
	
	//------
    countM=0;
    NDBASELOOP
    if(nodeflag(i,j,k)==1)
    ++countM;

    numtri = 6*(countM+countCC);
    numvert = countM+countCC*8;

    numtri_mem = numtri;
    numvert_mem = numvert;

    p->Iarray(tri,numtri,4);
    p->Darray(pt,numvert,3);
    p->Darray(ls,numvert);
    p->Iarray(facet,numtri,4);
    p->Iarray(confac,numtri);
    p->Iarray(numfac,numtri);
	p->Iarray(numpt,numtri);
    p->Darray(ccpt,numtri*4,3);


    countM=0;
    NDBASELOOP
    if(nodeflag(i,j,k)==1)
    {
    pt[countM][0] = p->posnode_x();
    pt[countM][1] = p->posnode_y();
    pt[countM][2] = p->posnode_z();

    ls[countM] = eta(i,j,k);

    vertice(i,j,k) = countM;

    ++countM;
    }

	// p. 725, 956
    count=0;
    BASELOOP
    if(nodeflag(i,j,k)==1)
    if(nodeflag(i-1,j,k)==1)
    if(nodeflag(i-1,j-1,k)==1)
    if(nodeflag(i,j-1,k)==1)
    if(nodeflag(i,j,k-1)==1)
    if(nodeflag(i-1,j,k-1)==1)
    if(nodeflag(i-1,j-1,k-1)==1)
    if(nodeflag(i,j-1,k-1)==1)
    {
    // 1
    tri[count][0] = vertice(i-1,j-1,k-1);
    tri[count][1] = vertice(i-1,j,k-1);
    tri[count][2] = vertice(i-1,j-1,k);
    tri[count][3] = vertice(i,j-1,k);
    ++count;

    // 2
    tri[count][0] = vertice(i-1,j-1,k-1);
    tri[count][1] = vertice(i,j-1,k-1);
    tri[count][2] = vertice(i-1,j,k-1);
    tri[count][3] = vertice(i,j-1,k);
    ++count;

    // 3
    tri[count][0] = vertice(i-1,j,k-1);
    tri[count][1] = vertice(i,j,k-1);
    tri[count][2] = vertice(i,j-1,k-1);
    tri[count][3] = vertice(i,j-1,k);
    ++count;

    // 4
    tri[count][0] = vertice(i,j,k-1);
    tri[count][1] = vertice(i-1,j,k-1);
    tri[count][2] = vertice(i,j-1,k);
    tri[count][3] = vertice(i,j,k);
    ++count;

    // 5
	tri[count][0] = vertice(i-1,j,k-1);
    tri[count][1] = vertice(i-1,j,k);
    tri[count][2] = vertice(i,j,k);
    tri[count][3] = vertice(i,j-1,k);
    ++count;

    // 6
    tri[count][0] = vertice(i-1,j,k-1);
    tri[count][1] = vertice(i-1,j-1,k);
    tri[count][2] = vertice(i,j-1,k);
    tri[count][3] = vertice(i-1,j,k);
    ++count;
    }

	
	// CC tri
	for(n=0;n<p->ccellnum;++n)
	if(ccflag[n]==1)
    {	
		for(q=0;q<ccnode[n];++q)
		{
		pt[countM][0] = ccell[n][q][0];
		pt[countM][1] = ccell[n][q][1];
		pt[countM][2] = ccell[n][q][2];
		ls[countM] = lscc[n][q];
		vertice_cc[n][q]=countM;
		++countM;
		}
    }
	
	for(n=0;n<p->ccellnum;++n)
	if(ccflag[n]==1)
	{	
		if(ccnode[n]==3)
		{
		// 1
		tri[count][0] = vertice_cc[n][0];
		tri[count][1] = vertice_cc[n][1];
		tri[count][2] = vertice_cc[n][2];
		tri[count][3] = vertice_cc[n][3];
		++count;
		}
		
		if(ccnode[n]==5)
		{
		// 1
		tri[count][0] = vertice_cc[n][0];
		tri[count][1] = vertice_cc[n][1];
		tri[count][2] = vertice_cc[n][2];
		tri[count][3] = vertice_cc[n][4];
		++count;

		// 2		
		tri[count][0] = vertice_cc[n][0];
		tri[count][1] = vertice_cc[n][2];
		tri[count][2] = vertice_cc[n][3];
		tri[count][3] = vertice_cc[n][4];
		++count;
		}
		
		if(ccnode[n]==6)
		{
		// 1
		tri[count][0] = vertice_cc[n][0];
		tri[count][1] = vertice_cc[n][1];
		tri[count][2] = vertice_cc[n][2];
		tri[count][3] = vertice_cc[n][4];
		++count;

		// 2		
		tri[count][0] = vertice_cc[n][0];
		tri[count][1] = vertice_cc[n][2];
		tri[count][2] = vertice_cc[n][3];
		tri[count][3] = vertice_cc[n][4];
		++count;

		// 3
		tri[count][0] = vertice_cc[n][3];
		tri[count][1] = vertice_cc[n][4];
		tri[count][2] = vertice_cc[n][5];
		tri[count][3] = vertice_cc[n][2];
		++count;
		}
		
		if(ccnode[n]==8)
		{
		// 1
		tri[count][0] = vertice_cc[n][0];
		tri[count][1] = vertice_cc[n][3];
		tri[count][2] = vertice_cc[n][4];
		tri[count][3] = vertice_cc[n][5];
		++count;

		// 2		
		tri[count][0] = vertice_cc[n][0];
		tri[count][1] = vertice_cc[n][1];
		tri[count][2] = vertice_cc[n][3];
		tri[count][3] = vertice_cc[n][5];
		++count;

		// 3
		tri[count][0] = vertice_cc[n][3];
		tri[count][1] = vertice_cc[n][2];
		tri[count][2] = vertice_cc[n][1];
		tri[count][3] = vertice_cc[n][5];
		++count;

		// 4
		tri[count][0] = vertice_cc[n][2];
		tri[count][1] = vertice_cc[n][3];
		tri[count][2] = vertice_cc[n][5];
		tri[count][3] = vertice_cc[n][6];
		++count;

		// 5
		tri[count][0] = vertice_cc[n][3];
		tri[count][1] = vertice_cc[n][7];
		tri[count][2] = vertice_cc[n][6];
		tri[count][3] = vertice_cc[n][5];
		++count;

		// 6
		tri[count][0] = vertice_cc[n][3];
		tri[count][1] = vertice_cc[n][4];
		tri[count][2] = vertice_cc[n][5];
		tri[count][3] = vertice_cc[n][7];
		++count;		
		}
		
		
	}
	//cout<<" "<<p->mpirank<<"  numtri_cc: "<<count<<endl;
    numtri=count;
}


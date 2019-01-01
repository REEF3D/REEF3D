/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2019 Hans Bihs

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

#include"pvccparse.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

pvccparse::pvccparse(lexer* p, fdm* a, ghostcell *pgc) : dx(p->dx), eps(1.0e-5)
{
}

pvccparse::~pvccparse()
{
}

void pvccparse::start(lexer* p, fdm* a, ghostcell *pgc)
{
    count=0;
	int countold;
    for(n=0;n<p->facetnum;++n)
    {
		countold = count;
        ini(p,a,pgc);
        cellnodes(p,a,pgc);
        pointcheck(p,a,pgc);
        collectpoints(p,a,pgc);

        pointcount=0;
        for(int qn=0;qn<8;++qn)
        if(pt[qn]>0)
        ++pointcount;
		
		
		
        //Paraview
        if(pcount==4)
		{
		face_tetra(p,a,pgc);
        cell_tetra(p,a,pgc);
		}

        if(pcount==10 && pointcount==7 && clcount==3) 
        cell_reverse_tetra(p,a,pgc);

        if(pcount==9 && pointcount==7 && clcount==2) 
        cell_reverse_tetra_a(p,a,pgc);

        if(pcount==8 && pointcount==7 && clcount==1)
        cell_reverse_tetra_b(p,a,pgc);

        if(pcount==6) 
        {
		face_wedge(p,a,pgc);
		cell_wedge(p,a,pgc);
		}
		
		if(pcount>8 && pointcount==6 && clcount==4) 
        cell_divide4(p,a,pgc);
		
		if(pcount==9 && pointcount==6 && clcount==3) 
        cell_divide4a(p,a,pgc);
		
        if(pcount==7)
        cell_divide4b(p,a,pgc);

		if(pcount==8 && pointcount==6 && clcount==2)
        cell_divide4c(p,a,pgc);
		
        if(pcount==8 && ((pointcount==4 && clcount==4) || (pointcount==5 && clcount==3) || (pointcount==6 && clcount==2 ))) // f
		{
		face_hexahedron(p,a,pgc);
        cell_hexahedron(p,a,pgc);
		}

        if(pcount==8 && ((pointcount==3 && clcount==5) || (pointcount==4 && clcount==4))) 
        cell_divide5a(p,a,pgc);

		if(pcount==10 && pointcount==5 && clcount==5) 
        cell_divide5b(p,a,pgc);

        if(pcount==9 && pointcount==5 && clcount==4)
        cell_divide5d(p,a,pgc);

        if(pcount>8 && pointcount==4 && clcount==6)
        cell_divide6(p,a,pgc);
		
		int diff;
		
		diff = count-countold;

		if(diff<=0 && p->P17==1)
		{
		cout<<"Points: "<<pcount<<" Pts: "<<pointcount<<" Clcount: "<<clcount;
		cout<<" | CSx: "<<p->facet[n][3]<<" CSy: "<<p->facet[n][4]<<" CSz: "<<p->facet[n][5]<<" |  ";
		
		for(int qnn=0; qnn<8; ++qnn)
		if(pt[qnn]==1)
		cout<<qnn<<"  ";
		
		cout<<"  :   ";
		
		for(int qnn=0; qnn<12; ++qnn)
		if(cl[qnn]<=0)
		cout<<qnn<<"  ";
		
		cout<<endl;
		}
    }

    p->ccellnum=count;

    p->ccedgenum=0;
    for(n=0;n<p->ccellnum;++n)
    p->ccedgenum+=a->ccedge[n];

}


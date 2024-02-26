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

#include"ghostcell.h"
#include"lexer.h"
#include"fdm.h"


void ghostcell::gcb_generic(lexer* p,field& f,int *gcb_count, int ***gcb)
{ 
    int aa,bb,cc;
	int r;
    double gravity;
    
    for(n=0;n<6;++n)
    for(q=0;q<gcb_count[n];++q)
    {
    i=gcb[n][q][0];
    j=gcb[n][q][1];
    k=gcb[n][q][2];
    
    aa=bb=cc=0;
    
		for(r=1;r<=2;++r)
		{
		if(n==0)
		aa=r;

		if(n==1)
		bb=r;
		
		if(n==2)
		bb=-r;
		
		if(n==3)
		aa=-r;
		
		if(n==4)
		cc=-r;
		
		if(n==5)
		cc=r;
        
        
		f(i+aa,j+bb,k+cc)=f(i,j,k);
        
        //if(k==5)
        //cout<<p->mpirank<<" i: "<<i<<" n: "<<n<<" f(i,j,k): "<<f(i,j,k)<<" f(i+aa,j+bb,k+cc): "<<f(i+aa,j+bb,k+cc)<<" aa: "<<aa<<endl;
		}
    }

}

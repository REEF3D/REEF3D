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

#include"lexer.h"
#include"ghostcell.h"

void lexer::gridini(ghostcell *pgc)
{        
    if(G2==1)
    sigma_coord_ini(this);
    
    lexer_gridspacing(pgc);
	parse();	
}

void lexer::flagini()
{
    control_calc();

	gridsize();
	vellast();
	
	Iarray(flag1,imax*jmax*kmax);
	Iarray(flag2,imax*jmax*kmax);
	Iarray(flag3,imax*jmax*kmax);
    Iarray(flag5,imax*jmax*kmax);
    Iarray(flag,imax*jmax*kmax);
	
	Iarray(tpflag,imax*jmax*kmax);
    Iarray(ndbaseflag,imax*jmax*kmax);


	makeflag(flag1);
	makeflag(flag2);
	makeflag(flag3);
	makeflag(tpflag);

	
	if(N5==0)
	i_dir=j_dir=k_dir=1;
	
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
}

void lexer::gridini_outflow()
{
	int istart,iend,jstart,jend,kstart,kend,qn;
	int count=0;
	
    for(qn=0;qn<G95;++qn)
    {

        istart = posc_i(G95_xs[qn]);
        iend = posc_i(G95_xe[qn]);
        
        jstart = posc_j(G95_ys[qn]);
        jend = posc_j(G95_ye[qn]);
        
        kstart = posc_k(G95_zs[qn]);
        kend = posc_k(G95_ze[qn]);
        
        
        for(n=0;n<gcb4_count;++n)
		{
		i=gcb4[n][0];
		j=gcb4[n][1];
		k=gcb4[n][2];
		
			if(i>=istart && i<iend && j>=jstart && j<jend && k>=kstart && k<kend && gcb4[n][3]==5 && (gcb4[n][4]==21||gcb4[n][4]==22))
			{
			++count;
			gcb4[n][4]=2;
			}
		}
    }
	count+=gcout_count;
	
	Iresize(gcout, gcout_count,count,6,6);
    Iresize(gcout6, gcout6_count,count,6,6);
    /*
    gcout_count = count;
    cout<<mpirank<<" "<<gcout_count<<" "<<count<<endl;
    gcout6_count = count;*/
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

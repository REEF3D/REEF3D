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

#include"lexer.h"
#include"ghostcell.h"

void lexer::gridini(ghostcell *pgc)
{        
    if(G2==1)
    sigma_coord_ini();
    
    lexer_gridspacing(pgc);
	parse();	
    gcd_ini(pgc);
}

void lexer::flagini()
{
    control_calc();
	gridsize();
	
	Iarray(flag1,imax*jmax*kmax);
	Iarray(flag2,imax*jmax*kmax);
	Iarray(flag3,imax*jmax*kmax);
    Iarray(flag5,imax*jmax*kmax);
    Iarray(flag,imax*jmax*kmax);
    
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
    Iarray(DF,imax*jmax*kmax);
    
    // flag
	Iarray(tpflag,imax*jmax*kmax);
    Iarray(ndbaseflag,imax*jmax*kmax);


	makeflag(flag1);
	makeflag(flag2);
	makeflag(flag3);
	makeflag(tpflag);
    
    for(i=-margin; i<knox+margin; ++i)
    for(j=-margin; j<knoy+margin; ++j)
    for(k=-margin; k<knoz+margin; ++k)
    IO[(i-imin)*jmax*kmax + (j-jmin)*kmax + k-kmin] = 0;
    
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
}

void lexer::gridini_patchBC()
{
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
        i=gcb4[q][0];
		j=gcb4[q][1];
		k=gcb4[q][2];
	
    if(gcb4[q][3]==1 || gcb4[q][3]==4)
    gcd4[q] = 0.5*DXP[IP];
    
    if(gcb4[q][3]==2 || gcb4[q][3]==3)
    gcd4[q] = 0.5*DYP[JP];
    
    if(gcb4[q][3]==5 || gcb4[q][3]==6)
    gcd4[q] = 0.5*DZP[KP];
	}
}

void lexer::sigma_coord_ini()
{
    double L, ZN0temp;
    
    L = ZN[knoz+marge] - ZN[0+marge];
    
    ZN0temp = ZN[0+marge];
    
    for(k=-marge;k<knoz+marge;++k)
    ZN[KP] = (ZN[KP]-ZN0temp)/L;
}

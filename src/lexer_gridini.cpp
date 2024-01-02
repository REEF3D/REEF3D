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
	
    
    //cout<<mpirank<<" imax: "<<imax<<" jmax: "<<jmax<<" kmax: "<<kmax<<endl;
	
	Iarray(flag1,imax*jmax*kmax);
	Iarray(flag2,imax*jmax*kmax);
	Iarray(flag3,imax*jmax*kmax);
    Iarray(flag5,imax*jmax*kmax);
    Iarray(flag,imax*jmax*kmax);
    
    //cout<<mpirank<<" flagini: "<<imax*jmax*kmax<<endl;
	
    Iarray(flagsf1,imax*jmax*kmax);
	Iarray(flagsf2,imax*jmax*kmax);
	Iarray(flagsf3,imax*jmax*kmax);
	Iarray(flagsf4,imax*jmax*kmax);
    
    Iarray(BC,imax*jmax*kmax);
    
	Iarray(tpflag,imax*jmax*kmax);
    Iarray(ndbaseflag,imax*jmax*kmax);


	makeflag(flag1);
	makeflag(flag2);
	makeflag(flag3);
	makeflag(tpflag);
	
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

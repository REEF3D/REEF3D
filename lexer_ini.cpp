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

void lexer::lexer_ini()
{
    umax=0.0;
    vmax=0.0;
    wmax=0.0;
    kinmax=0.0;
    epsmax=0.0;
    pressmax=0.0;

    originx+=global_xmin;
    originy+=global_ymin;
    originz+=global_zmin;

    utime=vtime=wtime=0.0;
    kintime=epstime=poissontime=lsmtime=susptime=topotime=printouttime=0.0;

    uiter=viter=witer=0;
    kiniter=epsiter=poissoniter=lsmiter=suspiter=topoiter=0;
    count_statestart=-1;

    phimean=0.0;
    phiout=0.0;
    phiin=0.0;

    gcextra1=gcextra2=gcextra3=gcextra4=gcextra4a=0;

    dtsed=0.0;
    sedtime=0.0;
    sediter=0;
    slidecells=0;
    bedmin=bedmax=0.0;
	
	maxdt=mindt=0.0;

    G1=0;
    if(S10>0||(G50>0 && G51>0)||G60>0||G61>0||G39==1)
    G1=1;

    wT=0.0;
    wV=0.0;
    wH=0.0;
    wL=0.0;
    wd=0.0;
	
	velcorr=1;
	
	ufbmax=0.0;
	vfbmax=0.0;
	wfbmax=0.0;
	
	ufbi=vfbi=wfbi=0.0;
	pfbi=qfbi=rfbi=0.0;
		
}

void lexer::makeflag( int *field)
{
    int n;
	for(n=0;n<imax*jmax*kmax;++n)
	field[n]=OBJ;
}

void lexer::parse()
{
	

	if(F80>0)
	{
	D32=7;
	}
    
    if(F80>0 && F35>0)
    F35=0;
	
	if(I10==1)
    {
    I11=1;
    I12=2;
    I13=1;
    }
	
	if(I40>0)
	{
	I10=0;
    I11=0;
    I12=0;
    I13=0;
    }
	
	
	if(N15==-1)
	{
		if(I12==0)
		{
		N15=1;
		N16=1;
		}
	}
	
	else
	if(N15==1)
	{
		if(N16==0)
		N16=N46;
	}

    if(T10==0)
    I13=0;
	
	if(S10==1 || G50==1)
	P27=1;
	
	if(F350>0)
	F50=F350;
	
	
}

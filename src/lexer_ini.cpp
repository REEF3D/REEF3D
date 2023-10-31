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
#include"lexer.h"

void lexer::lexer_ini()
{
    umax=0.0;
    vmax=0.0;
    wmax=0.0;
    kinmax=0.0;
    epsmax=0.0;
    pressmax=0.0;
    omegamax=0.0;

    utime=vtime=wtime=0.0;
    kintime=epstime=poissontime=lsmtime=susptime=printouttime=0.0;
    recontime=fsftime=0.0;

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
    solver_status=0;
	
	maxdt=mindt=0.0;

    G1=0;
    if(S10>0 || toporead>0 || solidread==1)
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
	fbmax=0.0;
    sfmax=0.0;
    pressgage=0.0;
    
	ufbi=vfbi=wfbi=0.0;
	pfbi=qfbi=rfbi=0.0;
    
    if(B98==1)
    B98=2;
    
    if(A10==3 || A10==5)
    G2=1;
		
}

void lexer::makeflag( int *field)
{
    int n;
	for(n=0;n<imax*jmax*kmax;++n)
	field[n]=OBJ;
}

void lexer::parse()
{
    if(F80>0 && F35>0)
    F35=0;
	
	if(I10==1)
    {
    I11=1;
    I12=2;
    I13=1;
    }
    
    if(I10==2)
    {
    I11=2;
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


    if(T10==0)
    I13=0;
	
	if(S10>=1 || toporead==1)
	P27=1;
	
	
}

/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2026 Hans Bihs

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
    kintime=epstime=lsmtime=susptime=printouttime=dftime=0.0;
    poissontime=laplacetime=matrixtime=ptime=0.0;
    recontime=fsftime=0.0;

    uiter=viter=witer=0;
    kiniter=epsiter=poissoniter=lsmiter=suspiter=topoiter=0;
    count_statestart=-1;

    phimean=0.0;
    phiout=0.0;
    phiin=0.0;

    dtsed=0.0;
    sedtime=0.0;
    sediter=0;
    slidecells=0;
    bedmin=bedmax=0.0;
    solver_status=0;
    solver_error=0;
	
	maxdt=mindt=0.0;
    RK_alpha=0.0;
    wavetime=0.0;

    wT=0.0;
    wV=0.0;
    wH=0.0;
    wL=0.0;
    wd=0.0;
    wC=0.0;
	
	velcorr=1;
	
	ufbmax=0.0;
	vfbmax=0.0;
	wfbmax=0.0;
	fbmax=0.0;
    sfmax=0.0;
    pressgage=0.0;
    
	ufbi=vfbi=wfbi=0.0;
	pfbi=qfbi=rfbi=0.0;
}

void lexer::makeflag( int *field)
{
    int n;
	for(n=0;n<imax*jmax*kmax;++n)
	field[n]=OBJ_FLAG;
}

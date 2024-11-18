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
#include"fdm_nhf.h"
#include"ghostcell.h"
#include"nhflow_reinidisc_fsf.h"
#include"cpt.h"

nhflow_reinidisc_fsf::nhflow_reinidisc_fsf(lexer *p) :  ddweno_nug_sig(p)
{
}

nhflow_reinidisc_fsf::~nhflow_reinidisc_fsf()
{
}

void nhflow_reinidisc_fsf::start(lexer *p, fdm_nhf *d, ghostcell *pgc, double *F, double *L)
{
    LOOP
    disc(p,d,pgc,F,L);
}

void nhflow_reinidisc_fsf::disc(lexer *p, fdm_nhf *d, ghostcell *pgc, double *F, double *L)
{	
    
    L[IJK]=0.0;
    
	if((F[IJK]>=0.0 && F[Ip1JK]>=0.0 && F[Im1JK]>=0.0 && F[IJp1K]>=0.0 && F[IJm1K]>=0.0 && F[IJKp1]>=0.0 && F[IJKm1]>=0.0) 
	|| (F[IJK]<0.0  && F[Ip1JK]<0.0  && F[Im1JK]<0.0  && F[IJp1K]<0.0  && F[IJm1K]<0.0   && F[IJKp1]<0.0  && F[IJKm1]<0.0)) 
	{
	dx=0.0;
	dy=0.0;
	dz=0.0;
	lsv=F[IJK];
    lsSig=lsv/sqrt(lsv*lsv);

    if(fabs(lsv)<1.0e-8)
    lsSig=1.0;

// x
	xmin=(lsv-F[Im1JK])/p->DXP[IM1];
	xplus=(F[Ip1JK]-lsv)/p->DXP[IP];
	
	if(xmin*lsSig>0.0 && xplus*lsSig>-xmin*lsSig)
	dx=ddwenox(F,1.0);

	if(xplus*lsSig<0.0 && xmin*lsSig<-xplus*lsSig)
	dx=ddwenox(F,-1.0);

	if(xplus*lsSig>0.0 && xmin*lsSig<0.0)
	dx=0.0;

// y
    if(p->j_dir==1)
    {
	ymin=(lsv-F[IJm1K])/p->DYP[JM1];
	yplus=(F[IJp1K]-lsv)/p->DYP[JP];
	
	if(ymin*lsSig>0.0 && yplus*lsSig>-ymin*lsSig)
	dy=ddwenoy(F,1.0);

	if(yplus*lsSig<0.0 && ymin*lsSig<-yplus*lsSig)
	dy=ddwenoy(F,-1.0);

	if(yplus*lsSig>0.0 && ymin*lsSig<0.0)
	dy=0.0;
    }

// z
	zmin=(lsv-F[IJKm1])/p->DZP[KM1];
	zplus=(F[IJKp1]-lsv)/p->DZP[KP];
	
	if(zmin*lsSig>0.0 && zplus*lsSig>-zmin*lsSig)
	dz=ddwenoz(F,1.0);

	if(zplus*lsSig<0.0 && zmin*lsSig<-zplus*lsSig)
	dz=ddwenoz(F,-1.0);

	if(zplus*lsSig>0.0 && zmin*lsSig<0.0)
	dz=0.0;	
					

	dnorm=sqrt(dx*dx + dy*dy + dz*dz);
    

    if(p->j_dir==0)
    deltax = (1.0/2.0)*(p->DXN[IP] + p->DZN[KP]*d->WL(i,j));
	
    if(p->j_dir==1)
    deltax = (1.0/3.0)*(p->DXN[IP] + p->DYN[JP] + p->DZN[KP]*d->WL(i,j));
	
	sign=lsv/sqrt(lsv*lsv + dnorm*dnorm*deltax*deltax);
    
    if(sign!=sign)
    sign=1.0;
    
	L[IJK] = -(sign*dnorm - sign);
    }
}

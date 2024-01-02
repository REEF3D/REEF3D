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

#include"particle_pls.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include<math.h>

void particle_pls::reseed(lexer* p, fdm* a, ghostcell* pgc, double pnum_coeff)
{
	double lsc, maxpos, maxneg;
	double distfac=1.0;
	int qn;
	
	if(p->count>0)
		distfac=2.0;
	
	reseeded=0;
	
	LOOP
	if(fabs(a->phi(i,j,k))<epsi)
	{
		lsc = a->phi(i,j,k);
		
		// POS
		if(lsc<0.5*p->DXM && lsc>-0.5*p->DXM)
		maxpos = (0.5 + lsc/p->DXM)*double(pnum)*pnum_coeff;
		
		if(lsc>=0.5*p->DXM)
		maxpos = double(pnum)*pnum_coeff;
		
		if(lsc<=-0.5*p->DXM)
		maxpos = 0.0;
		
		qn=0;
		while(posnum(i,j,k)<maxpos && qn<pnum)
		{		
		check=posseed(p,a,pgc,distfac);
		
		if(check==1)
		posnum(i,j,k)+=1.0;
		
		++qn;
		}    
	
		// NEG		
		if(lsc>-0.5*p->DXM && lsc<0.5*p->DXM)
		maxneg = (0.5 - lsc/p->DXM)*double(pnum)*pnum_coeff;
		
		if(lsc<=-0.5*p->DXM)
		maxneg = double(pnum)*pnum_coeff;
		
		if(lsc>=0.5*p->DXM)
		maxneg = 0.0;
		
		//cout<<"maxpos: "<<maxpos<<"  maxneg: "<<maxneg<<"      lsc: "<<lsc/p->DXM<<endl;
		
		qn=0;
		while(negnum(i,j,k)<maxneg && qn<pnum)		
		{
		check=negseed(p,a,pgc,distfac);
		
		if(check==1)
		negnum(i,j,k)+=1.0;
		
		++qn;
		}    
    }	
}

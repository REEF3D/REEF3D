/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTIBILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"vrans_f.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

void vrans_f::u_source(lexer *p, fdm *a)
{
	// VRANS porosity
    count=0;
    if(p->B270>0 || p->B274>0 || p->B281>1 || p->B282>1 || p->B291>1)
		ULOOP
		{
			porval = 0.5*(a->porosity(i,j,k) + a->porosity(i+1,j,k));
			partval = 0.5*(porpart(i,j,k) + porpart(i+1,j,k));
			alphaval = 0.5*(alpha(i,j,k) + alpha(i+1,j,k));
			betaval = 0.5*(beta(i,j,k) + beta(i+1,j,k));
			viscval = 0.5*(a->visc(i,j,k)+a->visc(i+1,j,k));
			
			
			Aporval = Apor(porval,partval,alphaval,viscval);
			Bporval = Bpor(porval,partval,betaval);
				

			porousterm = Aporval*a->u(i,j,k) + Bporval*a->u(i,j,k)*fabs(a->u(i,j,k)); 
			
			a->rhsvec.V[count] -= porousterm;
			++count;
		}
	
}

void vrans_f::v_source(lexer *p, fdm *a)
{
	
	// VRANS porosity
    count=0;
    if(p->B270>0 || p->B274>0 || p->B281>1 || p->B282>1 || p->B291>1)
		VLOOP
		{
			porval = 0.5*(a->porosity(i,j,k) + a->porosity(i,j+1,k));
			partval = 0.5*(porpart(i,j,k) + porpart(i,j+1,k));
			alphaval = 0.5*(alpha(i,j,k) + alpha(i,j+1,k));
			betaval = 0.5*(beta(i,j,k) + beta(i,j+1,k));
			viscval = 0.5*(a->visc(i,j,k)+a->visc(i,j+1,k));
		
			Aporval = Apor(porval,partval,alphaval,viscval);
			Bporval = Bpor(porval,partval,betaval);

			porousterm = Aporval*a->v(i,j,k) + Bporval*a->v(i,j,k)*fabs(a->v(i,j,k));

		
			a->rhsvec.V[count] -= porousterm;
			++count;
		}
	
}

void vrans_f::w_source(lexer *p, fdm *a)
{
	// VRANS porosity
    count=0;
    if(p->B270>0 || p->B274>0 || p->B281>1 || p->B282>1 || p->B291>1)
		WLOOP
		{
			porval = 0.5*(a->porosity(i,j,k) + a->porosity(i,j,k+1));
			partval = 0.5*(porpart(i,j,k) + porpart(i,j,k+1));
			alphaval = 0.5*(alpha(i,j,k) + alpha(i,j,k+1));
			betaval = 0.5*(beta(i,j,k) + beta(i,j,k+1));
			viscval = 0.5*(a->visc(i,j,k)+a->visc(i,j,k+1));
			
			Aporval = Apor(porval,partval,alphaval,viscval);
			Bporval = Bpor(porval,partval,betaval);
			

			porousterm = Aporval*a->w(i,j,k) + Bporval*a->w(i,j,k)*fabs(a->w(i,j,k));

		
			a->rhsvec.V[count] -= porousterm;
			++count;
		}
		
}

double vrans_f::Apor(double por, double part, double alpha, double visc)
{
	val = alpha*(pow(1.0-por,2.0)/pow(por,3.0))*(viscval/pow(part,2.0));
	
	return val;
}

double vrans_f::Bpor(double por, double part, double beta)
{
	val = beta*(1.0 + 7.5/Cval)*((1.0-por)/pow(por,3.0))*(1.0/part);
	
	return val;
}

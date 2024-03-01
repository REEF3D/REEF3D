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

#include"surftens.h"
#include"lexer.h"
#include"fdm.h"

surftens::surftens(lexer* p):gradient(p),epsi(p->F45*p->DXM)
{
	tension=p->W5;
}

surftens::~surftens()
{
}

void surftens::surface_tension(fdm* a,lexer*p,field& surf,int gcval)
{
	n=0;
	
	if(gcval==10 && tension>1.0e-10)
	ULOOP
	{
	dirac=0.0;
		if(fabs(0.5*(a->phi(i,j,k)+a->phi(i+1,j,k)))<epsi && (p->F30!=0 || p->F40!=0))
		{
		dirac = (0.5/epsi)*(1.0 + cos((PI*0.5*(a->phi(i,j,k)+a->phi(i+1,j,k)))/epsi));

		curv = (pow(xdx(a,a->phi),2.0)*(xdyy(a,a->phi) + xdzz(a,a->phi))
				+pow(xdy(a,a->phi),2.0)*(xdxx(a,a->phi) + xdzz(a,a->phi))
				+pow(xdz(a,a->phi),2.0)*(xdxx(a,a->phi) + xdyy(a,a->phi))
				-2.0*(xdx(a,a->phi)*xdy(a,a->phi)*xdxy(a,a->phi)
					 +xdx(a,a->phi)*xdz(a,a->phi)*xdxz(a,a->phi)
					 +xdy(a,a->phi)*xdz(a,a->phi)*xdyz(a,a->phi)))
			  /(pow(xdx(a,a->phi)*xdx(a,a->phi)+xdy(a,a->phi)*xdy(a,a->phi)+xdz(a,a->phi)*xdz(a,a->phi),1.5)+1.0e-20);


		a->rhsvec.V[n]-=(curv*tension*dirac*xdx(a,a->phi))/(0.5*(a->ro(i,j,k)+a->ro(i+1,j,k)));
		}
	++n;
	}

	if(gcval==11 && tension>1.0e-10)
	VLOOP
	{
	dirac=0.0;
		if( fabs(0.5*(a->phi(i,j,k)+a->phi(i,j+1,k)))<epsi && (p->F30!=0 || p->F40!=0))
		{
		dirac = (0.5/epsi)*(1.0 + cos((PI*0.5*(a->phi(i,j,k)+a->phi(i,j+1,k)))/epsi));

		curv = (pow(ydx(a,a->phi),2.0)*(ydyy(a,a->phi) + ydzz(a,a->phi))
				+pow(ydy(a,a->phi),2.0)*(ydxx(a,a->phi) + ydzz(a,a->phi))
				+pow(ydz(a,a->phi),2.0)*(ydxx(a,a->phi) + ydyy(a,a->phi))
				-2.0*(ydx(a,a->phi)*ydy(a,a->phi)*ydxy(a,a->phi)
					 +ydx(a,a->phi)*ydz(a,a->phi)*ydxz(a,a->phi)
					 +ydy(a,a->phi)*ydz(a,a->phi)*ydyz(a,a->phi)))
			  /(pow(ydx(a,a->phi)*ydx(a,a->phi)+ydy(a,a->phi)*ydy(a,a->phi)+ydz(a,a->phi)*ydz(a,a->phi),1.5)+1.0e-20);


		a->rhsvec.V[n]-=(curv*tension*dirac*ydy(a,a->phi))/(0.5*(a->ro(i,j,k)+a->ro(i,j+1,k)));
		}
	++n;
	}

	if(gcval==12 && tension>1.0e-10)
	WLOOP
	{
	dirac=0.0;
		if(fabs(0.5*(a->phi(i,j,k)+a->phi(i,j,k+1)))<epsi && (p->F30!=0 || p->F40!=0))
		{
		dirac = (0.5/epsi)*(1.0 + cos((PI*0.5*(a->phi(i,j,k)+a->phi(i,j,k+1)))/epsi));

		curv = (pow(zdx(a,a->phi),2.0)*(zdyy(a,a->phi) + zdzz(a,a->phi))
				+pow(zdy(a,a->phi),2.0)*(zdxx(a,a->phi) + zdzz(a,a->phi))
				+pow(zdz(a,a->phi),2.0)*(zdxx(a,a->phi) + zdyy(a,a->phi))
				-2.0*(zdx(a,a->phi)*zdy(a,a->phi)*zdxy(a,a->phi)
					 +zdx(a,a->phi)*zdz(a,a->phi)*zdxz(a,a->phi)
					 +zdy(a,a->phi)*zdz(a,a->phi)*zdyz(a,a->phi)))
			  /(pow(zdx(a,a->phi)*zdx(a,a->phi)+zdy(a,a->phi)*zdy(a,a->phi)+zdz(a,a->phi)*zdz(a,a->phi),1.5)+1.0e-20);


		a->rhsvec.V[n]-=(curv*tension*dirac*zdz(a,a->phi))/(0.5*(a->ro(i,j,k)+a->ro(i,j,k+1)));
		}
	++n;
	}
}



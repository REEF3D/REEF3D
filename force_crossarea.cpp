/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2020 Hans Bihs

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

#include"force.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

void force::crossarea(lexer* p, fdm *a, ghostcell *pgc)
{
    int *iloc,*jloc,*flag;
    double *wsf,*xpos,*ypos;
    double xl,yl;
    int n;

    xl = xe-xs;
    yl = ye-ys;

    Ai=Aj=0.0;

	p->Iarray(iloc,4);
	p->Iarray(jloc,4);
	p->Iarray(flag,4);
	p->Darray(wsf,4);
	p->Darray(xpos,4);
	p->Darray(ypos,4);


    for(n=0;n<4;++n)
    {
    iloc[n]=0;
    jloc[n]=0;
    flag[n]=0;
    wsf[n]=0.0;
    }

    xpos[0] = xs + 0.5*xl;
    ypos[0] = ys + 0.125*yl;

    xpos[1] = xs + 0.5*xl;
    ypos[1] = ys + 0.875*yl;

    xpos[2] = xs + 0.125*xl;
    ypos[2] = ys + 0.5*yl;

    xpos[3] = xs + 0.875*xl;
    ypos[3] = ys + 0.5*yl;

    // -

    int check;

    for(n=0;n<4;++n)
    {
    iloc[n]=p->posc_i(xpos[n]);
    jloc[n]=p->posc_j(ypos[n]);

    check=ij_boundcheck(p,a,iloc[n],jloc[n],0);

    if(check==1)
    flag[n]=1;
    }


    // -

    double zval=0.0;

    for(n=0;n<4;++n)
    wsf[n]=-1.0e20;


    for(n=0;n<4;++n)
    if(flag[n]>0)
    {
    zval=0.0;

    i=iloc[n];
    j=jloc[n];

        KLOOP
        PCHECK
        {
            if(a->phi(i,j,k)>=0.0 && a->phi(i,j,k+1)<0.0)
            wsf[n]=MAX(wsf[n],-(a->phi(i,j,k)*p->dx)/(a->phi(i,j,k+1)-a->phi(i,j,k)) + p->pos_z());
        }
    }

    for(n=0;n<4;++n)
    wsf[n]=pgc->globalmax(wsf[n]);


    Ai = 0.5*(wsf[0] + wsf[1])*p->P82_y;
    Aj = 0.5*(wsf[2] + wsf[3])*p->P82_x;

    rcyl = 0.5*p->P84;

    hcyl = 0.25*(wsf[0] + wsf[1] + wsf[2] + wsf[3]);
    Vcyl = PI*rcyl*rcyl*hcyl;
	
	Vrect = p->P86_x*p->P86_y*(hcyl-p->P90);

}


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

#include"interface_width.h"
#include"picard_f.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"heaviside.h"

picard_f::picard_f(lexer *p) : gradient(p), vol1(0.0), vol2(0.0), dvol2(0.0)
{
}

picard_f::~picard_f()
{
}

void picard_f::volcalc(lexer *p, fdm *a, ghostcell *pgc, field& b)
{
    double H = 0.0;
    vol1=0.0;

    LOOP
	{
		H = heaviside(b(i,j,k),interface_width(p,b,i,j,k));

		vol1+=p->DXN[IP]*p->DYN[JP]*p->DZN[KP]*H;
	}

	vol1 = pgc->globalsum(vol1);
}

void picard_f::volcalc2(lexer *p, fdm *a, ghostcell *pgc, field& b)
{
    // volume and its derivative with respect to a uniform shift of the level set,
    // dV/dw = sum delta(phi) dV (the interface area for a signed distance)
    double V = 0.0;
    vol2=0.0;
    dvol2=0.0;

    LOOP
	{
        V = p->DXN[IP]*p->DYN[JP]*p->DZN[KP];

        double epsi = interface_width(p,b,i,j,k);

		vol2  += V*heaviside(b(i,j,k),epsi);
		dvol2 += V*heaviside_delta(b(i,j,k),epsi);
	}

	vol2  = pgc->globalsum(vol2);
	dvol2 = pgc->globalsum(dvol2);

}


void picard_f::correct_ls(lexer *p, fdm *a, ghostcell *pgc, field& b)
{
    double w;

    // no reference volume yet (volcalc not called)
    if(vol1<=0.0)
    return;

    // Newton iterations for a uniform shift w of the level set that restores the
    // reference volume: w = (V_ref - V) / (dV/dw)
    for(int n=0;n<p->F47;++n)
    {
    volcalc2(p,a,pgc,b);

    if(dvol2<=1.0e-20)
    break;

    w = (vol1-vol2)/dvol2;

    LOOP
    b(i,j,k)+=w;
    }
}

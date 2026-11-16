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

#include"diff_dynamic.h"
#include"lexer.h"
#include"fdm.h"
#include"field4.h"

diff_dynamic::diff_dynamic(lexer *p) : mode(p->D24>0 ? p->D24 : 0), mu(nullptr), zero(nullptr), size(0)
{
    if(mode>0)
    {
    mu = new field4(p);
    zero = new field4(p);
    size = p->imax*p->jmax*p->kmax;
    }
}

diff_dynamic::~diff_dynamic()
{
    delete mu;
    delete zero;
}

void diff_dynamic::update(lexer *p, fdm *a)
{
    if(mode==0)
    return;

    // every cell of the arrays, ghost and solid cells included: rho and nu are set there by the
    // fluid update and its ghost-cell update
    for(int n=0; n<size; ++n)
    mu->V[n] = a->ro.V[n]*(a->visc.V[n] + a->eddyv.V[n]);
}

field& diff_dynamic::visc(fdm *a)
{
    if(mode==0)
    return a->visc;

    return *mu;
}

field& diff_dynamic::eddyv(fdm *a)
{
    if(mode==0)
    return a->eddyv;

    return *zero;
}

double diff_dynamic::irho(fdm *a, int c, int i, int j, int k) const
{
    if(mode==0)
    return 1.0;

    const double r1 = a->ro(i,j,k);
    const double r2 = a->ro(i+(c==0),j+(c==1),k+(c==2));

    // a solid or unset neighbour: the density of the cell itself
    if(r2<=0.0)
    return r1>0.0 ? 1.0/r1 : 1.0;

    if(r1<=0.0)
    return 1.0/r2;

    return 2.0/(r1 + r2);
}

double diff_dynamic::corner(double m1, double m2, double m3, double m4) const
{
    if(mode!=2)
    return 0.25*(m1 + m2 + m3 + m4);

    // harmonic mean: the shear stress through layers of different viscosity (a film under air);
    // cells without a viscosity (solid cells without ghost values) do not count
    const double m[4] = {m1,m2,m3,m4};
    double s = 0.0;
    int n = 0;

    for(int q=0; q<4; ++q)
    if(m[q]>0.0)
    {
    s += 1.0/m[q];
    ++n;
    }

    return n>0 ? double(n)/s : 0.0;
}

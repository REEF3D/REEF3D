/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

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

#include"strain.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"fieldint.h"

double strain::qij(lexer *p, fdm *a, int ii, int jj)
{
	double q=0.0;

    if((ii==1 && jj==1) || (ii==2 && jj==2) || (ii==3 && jj==3))
        q = 0.0;

    if(ii==1 && jj==2)
        q = pudy(p,a) - pvdx(p,a);

    if(ii==2 && jj==1)
        q = -pudy(p,a) + pvdx(p,a);

    if(ii==1 && jj==3)
        q = pudz(p,a) - pwdx(p,a);

    if(ii==3 && jj==1)
        q = -pudz(p,a) + pwdx(p,a);

    if(ii==2 && jj==3)
        q = pvdz(p,a) - pwdy(p,a);

    if(ii==3 && jj==2)
        q = -pvdz(p,a) + pwdy(p,a);

	return q;
}

double strain::Qij2(lexer *p, fdm *a)
{    
    skewSymmetricStrainRateTensor(p,a->u,a->v,a->w);

    double r = sqrt(r12*r12 + r13*r13 + r23*r23);
    
    r=r*r;

	return r;
}

double strain::rotationterm(lexer *p, fdm *a)
{
    return rotationterm(p,a->u,a->v,a->w);
}

double strain::rotationterm(lexer *p, field &u, field &v, field &w)
{    
    skewSymmetricStrainRateTensor(p,u,v,w);

    double r = sqrt(r12*r12 + r13*r13 + r23*r23);

	return r;
}

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

#include"fluxlim.h"
#include"increment.h"

#ifndef SUPERBEE_H_
#define SUPERBEE_H_

using namespace std;

class superbee : public fluxlim, public increment
{
public:
	superbee (lexer *);
	virtual ~superbee();

	virtual double iphi(field&,int,int,int,int);
	virtual double jphi(field&,int,int,int,int);
	virtual double kphi(field&,int,int,int,int);

private:
    double r, phi,denom;
	double dx,dy,dz;
	double L;
};

#endif


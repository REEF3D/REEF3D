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

#ifndef LIMO3_H_
#define LIMO3_H_

#include"fluxlim.h"
#include"increment.h"

using namespace std;

class limo3 : public fluxlim, public increment
{

public:

	limo3 (lexer *);
	virtual ~limo3();

	virtual double iphi(field&,int,int,int,int);
	virtual double jphi(field&,int,int,int,int);
	virtual double kphi(field&,int,int,int,int);

private:
	
	double max(double,double,double);
	double max(double,double);
	double min(double,double,double);
	double min(double,double);

    double ul,ur,vl,vr,wl,wr;
    double r, phi;
	double eta,phihat,d1,d2;
	double dx,dy,dz;
	const double delta,radius,eps;
	double L;
};

#endif


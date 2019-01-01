/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2019 Hans Bihs

This file is part of REEF3D.

REEF3D is fra->eps software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Fra->eps Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. Sa->eps the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, sa->eps <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
--------------------------------------------------------------------*/

#include"boundarycheck.h"

class fdm;
class lexer;
class field;

#ifndef INTERPOLATION_H_
#define INTERPOLATION_H_

using namespace std;

class interpolation : virtual public boundarycheck
{
public:
    interpolation(lexer*);
	virtual ~interpolation();
    
    double ccipol1(field&,double,double,double);
    double ccipol2(field&,double,double,double);
    double ccipol3(field&,double,double,double);
    double ccipol4(field&,double,double,double);
    double ccipol4phi(fdm*,field&,double,double,double);
    double ccipol4press(fdm*,field&,double,double,double);
	double ccipol1_a(field&,double,double,double);
    double ccipol2_a(field&,double,double,double);
    double ccipol3_a(field&,double,double,double);
    double ccipol4_a(field&,double,double,double);
    
    double cctripol4_a(fdm*,field&,double,double,double);
	double cint4a(double,double,double,double,double);
	double tricubic4a(lexer*,fdm*,field&,int&,int&,int&,double,double,double);
    
    double ipol1(field&);
    double ipol2(field&);
    double ipol3(field&);
    double ipol4(field&);
    
	double ipol4ro(fdm*,field&);
    double ipol4phi(fdm*,field&);
    double ipol4press(field&);
    double ipol4_a(field&);
    
    double lint(field&,int&,int&,int&,double,double,double);
    double lint1(field&,int&,int&,int&,double,double,double);
    double lint2(field&,int&,int&,int&,double,double,double);
    double lint3(field&,int&,int&,int&,double,double,double);
    double lint4(field&,int&,int&,int&,double,double,double);
    double lint4phi(fdm*,field&,int&,int&,int&,double,double,double);
    double lint_a(field&,int&,int&,int&,double,double,double);
    
    double tricubic4a(field&,int&,int&,int&,double,double,double);

private:
    double v1,v2,v3,v4,v5,v6,v7,v8;
    double x1,x2,x3,x4,y1,y2;
    double wa,wb,wc,wx,wy,wz;
    double value;
    double epphi;
    
    int ii,jj,kk;
    
    lexer *p;
};

#endif
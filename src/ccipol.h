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
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"increment.h"

class lexer;
class fdm;
class field;


#ifndef CCIPOL_H_
#define CCIPOL_H_

using namespace std;

class ccipol : public increment
{

public:
	ccipol(lexer*);
	virtual ~ccipol();

    double ccipol1(fdm*,field&,double&,double&,double&);
    double ccipol2(fdm*,field&,double&,double&,double&);
    double ccipol3(fdm*,field&,double&,double&,double&);
    double ccipol4(fdm*,field&,double&,double&,double&);
    double ccipol4_a(fdm*,field&,double&,double&,double&);
    double ipol1(fdm*,lexer*,field&);
    double ipol2(fdm*,lexer*,field&);
    double ipol3(fdm*,lexer*,field&);
    double ipol4(fdm*,lexer*,field&);
    double ipol4_a(fdm*,lexer*,field&);
    double lint(fdm*,field&,int&,int&,int&,double,double,double);
    double lint_a(field&,int&,int&,int&,double,double,double);

    int ii,jj,kk;
    double wa,wb,wc,wx,wy,wz;
    double value;
    double x1,x2,x3,x4,y1,y2;
    const double dx;
    double val;
    double v1,v2,v3,v4,v5,v6,v7,v8;

};

#endif

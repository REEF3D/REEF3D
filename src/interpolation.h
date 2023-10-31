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

#include"boundarycheck.h"

class fdm;
class lexer;
class field;
class slice;
class sliceint;

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
    double ccipol4_b(field&,double,double,double);
    double ccipol4_kin(field&,double,double,double);
    
    double cctripol4_a(fdm*,field&,double,double,double);
	double cint4a(double,double,double,double,double);
	double tricubic4a(lexer*,fdm*,field&,int&,int&,int&,double,double,double);
    
    double ipol1(field&);
    double ipol2(field&);
    double ipol3(field&);
    double ipol4(field&);
    
	double ipol4ro(fdm*,field&);
    double ipol4phi(fdm*,field&);
    double ipol4topo(fdm*,field&);
    double ipol4press(field&);
    double ipol4_a(field&);

    
    double lint(field&,int&,int&,int&,double,double,double);
    double lint1(field&,int&,int&,int&,double,double,double);
    double lint2(field&,int&,int&,int&,double,double,double);
    double lint3(field&,int&,int&,int&,double,double,double);
    double lint4(field&,int&,int&,int&,double,double,double);
    double lint4phi(fdm*,field&,int&,int&,int&,double,double,double);
    double lint_a(field&,int&,int&,int&,double,double,double);
    double lint4b(field&,int&,int&,int&,double,double,double);
    double lint4kin(field&,int&,int&,int&,double,double,double);
    
    double lint1_2D(field&,int&,int&,int&,double,double,double);
    double lint2_2D(field&,int&,int&,int&,double,double,double);
    double lint3_2D(field&,int&,int&,int&,double,double,double);
    double lint4_2D(field&,int&,int&,int&,double,double,double);
    double lint_a_2D(field&,int&,int&,int&,double,double,double);
    double lint4phi_2D(fdm*,field&,int&,int&,int&,double,double,double);
    
    double tricubic4a(field&,int&,int&,int&,double,double,double);
    
    
    // slice
    double sl_ipol1(slice&);
    double sl_ipol2(slice&);
	double sl_ipol1a(slice&);
    double sl_ipol2a(slice&);
    double sl_ipol4(slice&);
    double sl_ipol1eta(int*,slice&,slice&);
    double sl_ipol2eta(int*,slice&,slice&);
    double sl_ipol4eta(int*,slice&,slice&);
    double sl_ipol4eta_wd(int*,slice&,slice&);
    double sl_ipolint(sliceint&);
    
    double ccslipol1(slice&,double,double);
    double ccslipol2(slice&,double,double);
    double ccslipol4(slice&,double,double);
    
    double lintsl1(slice&,int&,int&,double,double);
    double lintsl2(slice&,int&,int&,double,double);
    double lintsl4(slice&,int&,int&,double,double);
    
    

private:
    double v1,v2,v3,v4,v5,v6,v7,v8;
    double x1,x2,x3,x4,y1,y2;
    int c1,c2,c3,c4;
    double wa,wb,wc,wx,wy,wz;
    double value;

    
    int ii,jj,kk;
    
    lexer *p;
};

#endif

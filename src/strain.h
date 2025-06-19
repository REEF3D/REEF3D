/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTIBILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
Author: Hans Bihs
--------------------------------------------------------------------*/

#ifndef STRAIN_H_
#define STRAIN_H_

#include"gradient.h"

class lexer;
class fdm;
class ghostcell;
class field;
class fieldint;

class strain : public gradient
{

public:
    strain(lexer*);
    virtual ~strain()=default;

	double sij(lexer*,fdm*,int,int);
	double qij(lexer*,fdm*,int,int);
	double pk(lexer*,fdm*,field&);
    double pk_b(lexer*,fdm*,field&);

    void wallf_update(lexer*,fdm*,ghostcell*,fieldint&);
    double strainterm(lexer*,fdm*);
    double strainterm(lexer*,field&,field&,field&);
    double rotationterm(lexer*,fdm*);
    double rotationterm(lexer*,field&,field&,field&);
    double magSqrSd(lexer*,fdm*);
    double magSqrSd(lexer*,field&,field&,field&);
    double strainplain(lexer*,fdm*);

private:
    void symmetricStrainRateTensor(lexer*,field&,field&,field&);
    void skewSymmetricStrainRateTensor(lexer*,field&,field&,field&);
    
    double s11,s22,s33,s12,s13,s23;
    double r11,r22,r33,r12,r13,r23;
    double ss11,ss22,ss33,ss12,ss13,ss23;
    double rr11,rr22,rr33,rr12,rr13,rr23;
    double q11,q22,q33,q12,q13,q23;
	double pkterm,s,q,val;
	const double epsi;

};

#endif

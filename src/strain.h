/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2023 Hans Bihs

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

#include"gradient.h"
#include"field1.h"
#include"field2.h"
#include"field3.h"
#include"field4.h"

class ghostcell;
class fieldint;

#ifndef STRAIN_H_
#define STRAIN_H_

using namespace std;

class strain : public gradient
{

public:
	strain (lexer*,fdm*);
	virtual ~strain();

	double sij(lexer*,fdm*,int,int);
	double qij(lexer*,fdm*,int,int);
	double pk(lexer*,fdm*);
	double pk_k(lexer*,fdm*);
	double pk_w(lexer*,fdm*);
	void Pk_update(lexer*,fdm*,ghostcell*);
	void wallf_update(lexer*,fdm*,ghostcell*,fieldint&);
	virtual double strainterm(lexer*,fdm*);
    virtual double strainterm(lexer*,field&,field&,field&);
	virtual double rotationterm(lexer*,fdm*);
    virtual double rotationterm(lexer*,field&,field&,field&);
	virtual double magSqrSd(lexer*,fdm*);
    virtual double magSqrSd(lexer*,field&,field&,field&);
	double strainplain(lexer*,fdm*);
	field4 Pk;

private:
    double s11,s22,s33,s12,s13,s23;
    double r11,r22,r33,r12,r13,r23;
    double ss11,ss22,ss33,ss12,ss13,ss23;
    double rr11,rr22,rr33,rr12,rr13,rr23;
    double q11,q22,q33,q12,q13,q23;
	double pkterm,s,q,val;
	const double epsi;

    

};

#endif

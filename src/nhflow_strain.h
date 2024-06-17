/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

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

#ifndef NHFLOW_STRAIN_H_
#define NHFLOW_STRAIN_H_

#include"nhflow_gradient.h"
#include"field1.h"
#include"field2.h"
#include"field3.h"
#include"field4.h"

class ghostcell;
class fieldint;

using namespace std;

class nhflow_strain : public nhflow_gradient
{

public:
	nhflow_strain (lexer*,fdm_nhf*);
	virtual ~nhflow_strain();

	double sij(lexer*,fdm_nhf*,int,int);
	double qij(lexer*,fdm_nhf*,int,int);
	double pk(lexer*,fdm_nhf*);
	double pk_k(lexer*,fdm_nhf*);
	double pk_w(lexer*,fdm_nhf*);
	void Pk_update(lexer*,fdm_nhf*,ghostcell*);
	void wallf_update(lexer*,fdm_nhf*,ghostcell*,int*);
	virtual double strainterm(lexer*,fdm_nhf*);
    virtual double strainterm(lexer*,double*,double*,double*);
	virtual double rotationterm(lexer*,fdm_nhf*);
    virtual double rotationterm(lexer*,double*,double*,double*);
	virtual double magSqrSd(lexer*,fdm_nhf*);
    virtual double magSqrSd(lexer*,double*,double*,double*);
	double strainplain(lexer*,fdm_nhf*);
    
	double *PK;

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

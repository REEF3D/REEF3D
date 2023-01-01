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
--------------------------------------------------------------------*/

#include"strain.h"
#include"field4.h"

class fdm;

using namespace std;

#ifndef EARSM_H_
#define EARSM_H_

class EARSM : private strain
{
public:
	EARSM(lexer *,fdm*);
	virtual ~EARSM();
	virtual void sq(lexer*,fdm*);
	virtual void invar();
	virtual void beta(fdm*);
	virtual void terms();
	virtual void isource(lexer*, fdm*);
	virtual void jsource(lexer*, fdm*);
	virtual void ksource(lexer*, fdm*);
	double tk;

    field4 rs11, rs22, rs33, rs12, rs13, rs23;
	double IIs,IIq,IV,V;
	double rep;
	double b1,b3,b4,b6,b9;
	double P1,P2,N,Q,cmu;
	double phi1,phi2,D;
	double T1_11, T1_12, T1_13, T1_22, T1_23, T1_33;
	double T3_11, T3_12, T3_13, T3_22, T3_23, T3_33;
	double T4_11, T4_12, T4_13, T4_22, T4_23, T4_33;
	double T6_11, T6_12, T6_13, T6_22, T6_23, T6_33;
	double T9_11, T9_12, T9_13, T9_22, T9_23, T9_33;
	double s11,s22,s33,s12,s13,s23;
	double ss11,ss22,ss33,ss12,ss13,ss23;
	double q12, q13, q23;
	double qq12,qq13,qq23;
	const double c_1,c1,tt,ot,ctau;
	double signP;
	double repacos,t;
	int gcval_earsm;
};

#endif



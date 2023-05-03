/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2019 Hans Bihs

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

#include"EARSM.h"
#include"fdm.h"
#include"lexer.h"

EARSM::EARSM(lexer* p, fdm* a) : strain(p,a),rs11(p),rs22(p),rs33(p),rs12(p),rs13(p),rs23(p),
c_1(1.8), c1(1.8),tt(2.0/3.0),ot(1.0/3.0),ctau(6.0)
{
}

EARSM::~EARSM()
{
}

void EARSM::sq(lexer *p, fdm* a)
{
    s11 = tk*pudx(p,a);
	s22 = tk*pvdy(p,a);
	s33 = tk*pwdz(p,a);

	s12 = tk*0.5*(pudy(p,a) + pvdx(p,a));
	s13 = tk*0.5*(pudz(p,a) + pwdx(p,a));
	s23 = tk*0.5*(pvdz(p,a) + pwdy(p,a));

	q12 = tk*0.5*(pudy(p,a) - pvdx(p,a));
	q13 = tk*0.5*(pudz(p,a) - pwdx(p,a));
	q23 = tk*0.5*(pvdz(p,a) - pwdy(p,a));

	ss11=s11*s11;
	ss22=s22*s22;
	ss33=s33*s33;
	ss12=s12*s12;
	ss13=s13*s13;
	ss23=s23*s23;

	qq12=q12*q12;
	qq13=q13*q13;
	qq23=q23*q23;
}

void EARSM::invar()
{
		IIs = ss11+ss22+ss33 + 2.0*(ss12+ss13+ss23);

		IIq = -2.0*(qq12+qq13+qq23);

		IV  = -s11*(qq12+qq13)-s22*(qq12+qq23)
			  -s33*(qq13+qq23)
			  +2.0*(-s12*q13*q23 + s13*q12*q23 - s23*q12*q13);

		V   = 2.0*(-(s12*s13 + s22*s23 + s23*s33)*q12*q13
				   +(s11*s13 + s12*s23 + s13*s33)*q12*q23
				   -(s11*s12 + s12*s22 + s13*s23)*q13*q23)
			   -(ss11+ss12+ss13)*(qq12+qq13)
			   -(ss12+ss22+ss23)*(qq12+qq23)
			   -(ss13+ss23+ss33)*(qq13+qq23);

}

void EARSM::beta(fdm* a)
{

        P1 = ((1.0/27.0)*c_1*c_1 + (9.0/20.0)*IIs - (2.0/3.0)*IIq)*c_1;
		P2 = P1*P1 - pow(((c_1*c_1)/9.0 + 0.9*IIs + (2.0/3.0)*IIq),3.0);
		rep=0.0;

		if(P2>=0.0)
		{
			signP = (P1 - sqrt(P2))/fabs((P1 - sqrt(P2)));
			N = (c_1/3.0) + pow(fabs(P1 + sqrt(P2)),(1.0/3.0)) + signP*pow(fabs(P1 - sqrt(P2)),(1.0/3.0));
		}
		if(P2<0.0)
		{
			rep = (P1/sqrt(P1*P1-P2));

				if(rep>=1.0)
				rep = 1.0;

				if(rep<=-1.0)
				rep = -1.0;

				repacos = acos(rep);

			N = (c_1/3.0) + 2.0*pow(P1*P1 - P2,(1.0/6.0)) * cos((1.0/3.0)*repacos);
		}

		phi1=IV*IV;
		phi2=V-IIs-IIq*0.5;

        D = 20*pow(N,4.0)*(N-0.5*c_1)-IIq*(10.0*pow(N,3.0)+15.0*c_1*pow(N,2.0))
          + 10.0*c_1*pow(IIq,2.0);

        N += (162.0*(phi1+phi2*pow(N,2.0)))/D;


		Q = (5.0/6.0)*(N*N - 2.0*IIq)*(2.0*N*N-IIq);

		b1 = -(N*(2.0*N*N - 7.0*IIq))/Q;
		b3 = -(12.0*IV)/(Q*N);
		b4 = -(2.0*(N*N - 2.0*IIq))/Q;
		b6 = -(6.0*N)/Q;
		b9 = 6.0/Q;
}

void EARSM::terms()
{
        T1_11 = b1*s11;
        T1_22 = b1*s22;
        T1_33 = b1*s33;
        T1_12 = b1*s12;
        T1_13 = b1*s13;
        T1_23 = b1*s23;

		T3_11 = b3*(-qq12-qq13-ot*IIq);
		T3_22 = b3*(-qq12-qq23-ot*IIq);
		T3_12 = b3*(-q13*q23);
		T3_13 = b3*(q12*q23);
		T3_23 = b3*(-q12*q13);

		T4_11 = b4*(-2.0*(s12*q12+s13*q13));
		T4_22 = b4*(2.0*(s12*q12-s23*q23));
		T4_12 = b4*((s11-s22)*q12 - s23*q13 - s13*q23);
		T4_13 = b4*(-s23*q12 + (s11-s33)*q13 +s12*q23);
		T4_23 = b4*(s13*q12 + s12*q13 + (s22-s33)*q23);

		T6_11 = b6*(-2.0*((s12*q13 - s13*q12)*q23 + s11*(qq12+qq13)) - tt*IV - IIq*s11);
		T6_22 = b6*(-2.0*((s23*q12 - s12*q23)*q13 + s22*(qq12+qq23)) - tt*IV - IIq*s22);
		T6_12 = b6*(-s12*(2.0*qq12 + qq13 + qq23) - (s13*q13 - s23*q23)*q12 - (s11+s22)*q13*q23 - IIq*s12);
		T6_13 = b6*(-s13*(qq12 + 2.0*qq13 + qq23) - (s12*q12 + s23*q23)*q13 + (s11+s33)*q12*q23 - IIq*s13);
		T6_23 = b6*(-s23*(qq12 + qq13 + 2.0*qq23) + (s12*q12 - s13*q13)*q23 - (s22+s33)*q12*q13 - IIq*s23);

		T9_11 = b9*(-2.0*((s12*q12 + s13*q13 - s23*q23)*qq12
						 +(s12*q12 + s13*q13 + s23*q23)*qq13
						 +(s22-s33)*q12*q13*q23));
		T9_22 = b9*(-2.0*((-s12*q12 - s13*q13 + s23*q23)*qq12
						 +(-s12*q12 + s13*q13 + s23*q23)*qq23
						 +(-s11+s33)*q12*q13*q23));
		T9_12 = b9*(((s11-s22)*q12 - 2.0*(s13*q23 + s23*q13))*qq12
				   +((s11-s33)*q12 - 2.0*s13*q23)*qq13
				   +((s33-s22)*q12 - 2.0*s23*q13)*qq23);
		T9_13 = b9*(((s11-s22)*q13 + 2.0*s13*q23)*qq12
				   +((s11-s33)*q13 + 2.0*(s12*q23 - s23*q12))*qq13
				   +((s22-s33)*q13 - 2.0*s23*q12)*qq23);
		T9_23 = b9*(((s22-s11)*q23 + 2.0*s12*q13)*qq12
				   +((s11-s33)*q23 + 2.0*s13*q12)*qq13
				   +((s22-s33)*q23 + 2.0*(s12*q13 + s13*q12))*qq23);
}

void EARSM::isource(lexer *p, fdm* a)
{
	ULOOP
	{
	a->F(i,j,k) = -(xdx(a,rs11)+xdy(a,rs12)+xdz(a,rs13));
    a->maxF=MAX(fabs(a->F(i,j,k)),a->maxF);
	}
}
void EARSM::jsource(lexer *p, fdm* a)
{
	VLOOP
	{
	a->G(i,j,k) = -(ydx(a,rs12)+ydy(a,rs22)+ydz(a,rs23));
    a->maxG=MAX(fabs(a->G(i,j,k)),a->maxG);
	}
}

void EARSM::ksource(lexer *p, fdm* a)
{
	WLOOP
	{
	a->H(i,j,k) = -(zdx(a,rs13)+zdy(a,rs23)+zdz(a,rs33));
    a->maxH=MAX(fabs(a->H(i,j,k)),a->maxH);
	}
}






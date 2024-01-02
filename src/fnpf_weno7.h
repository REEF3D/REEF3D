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

#include"fnpf_convection.h"
#include"increment.h"
#include"ddweno_f_nug.h"

#ifndef FNPF_WENO7_H_
#define FNPF_WENO7_H_

using namespace std;

class fnpf_weno7 : public fnpf_convection, public increment, public ddweno_f_nug
{
public:
	fnpf_weno7(lexer*);
	virtual ~fnpf_weno7();

    virtual double fx(lexer*, field&, double, double);
	virtual double fy(lexer*, field&, double, double);
	virtual double fz(lexer*, field&, double, double);
    
    virtual double sx(lexer*, slice&, double);
	virtual double sy(lexer*, slice&, double);
    virtual double sz(lexer*, double*);

private:
    
    double is1,is2,is3,is4;
	double alpha1,alpha2,alpha3,alpha4;
	double w1,w2,w3,w4;
	double q0,q1,q2,q3,q4,q5,q6,q7;
	double gradx, grady, gradz;
    
    double ivel1,ivel2,jvel1,jvel2;


	void is();
	void alpha();
	void weight();
    
    const double epsilon;
    
    void iqmin(lexer*, slice&);
    void jqmin(lexer*, slice&);
    void iqmax(lexer*, slice&);
    void jqmax(lexer*, slice&);
    

};

#endif

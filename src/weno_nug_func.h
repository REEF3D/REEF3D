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

#include"increment.h"

class lexer;
class fdm;
class field;

#ifndef WENO_NUG_FUNC_H_
#define WENO_NUG_FUNC_H_

using namespace std;

class weno_nug_func : public increment
{
public:
	weno_nug_func(lexer*);
	virtual ~weno_nug_func();

	void precalc_qf(lexer*);
    void precalc_cf(lexer*);
    void precalc_isf(lexer*);
    
    void ini(lexer*);
    


	void is_min_x();
    void is_min_y();
    void is_min_z();
    
    void is_max_x();
    void is_max_y();
    void is_max_z();
    
	void weight_min_x();
    void weight_min_y();
    void weight_min_z();
    
    void weight_max_x();
    void weight_max_y();
    void weight_max_z();
    
    
    static int* ggcmem;
    
    static double ****qfx,****qfy,****qfz;
    static double ***cfx,***cfy,***cfz;
    static double ****isfx,****isfy,****isfz;
    
	static int iniflag;
    
    
    
    double q1,q2,q3,q4,q5;
    
    const double epsilon,psi;
	double is1x,is2x,is3x;
    double is1y,is2y,is3y;
    double is1z,is2z,is3z;
	double w1x,w2x,w3x;
    double w1y,w2y,w3y;
    double w1z,w2z,w3z;
    
    int uf,vf,wf;

private:
    lexer *pp;

    
    
};

#endif

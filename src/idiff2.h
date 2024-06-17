/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

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

#ifndef IDIFF2_H_
#define IDIFF2_H_

#include"gradient.h"
#include"diffusion.h"

using namespace std;

class rheology;


class idiff2 : public diffusion, public gradient
{

public:

	idiff2(lexer*);
	virtual ~idiff2();

	virtual void diff_scalar(lexer*, fdm*, ghostcell*, solver*, field&, field&, field&, double, double);
	virtual void diff_scalar(lexer*, fdm*, ghostcell*, solver*, field&, field&, field&, field&, double, double){};
    virtual void idiff_scalar(lexer*, fdm*, ghostcell*, solver*, field&, field&, double, double);
    
	virtual void diff_u(lexer*, fdm*, ghostcell*, solver*, field&, field&, field&, double);
	virtual void diff_v(lexer*, fdm*, ghostcell*, solver*, field&, field&, field&, double);
	virtual void diff_w(lexer*, fdm*, ghostcell*, solver*, field&, field&, field&, double);
    
    virtual void diff_u(lexer*, fdm*, ghostcell*, solver*, field&, field&, field&, field&, field&, double){};
	virtual void diff_v(lexer*, fdm*, ghostcell*, solver*, field&, field&, field&, field&, field&, double){};
	virtual void diff_w(lexer*, fdm*, ghostcell*, solver*, field&, field&, field&, field&, field&, double){};

private:
    rheology *prheo;
    
	double D;
	int count,n;
	double b_ijk,ev_ijk,visc_ijk;
	double b_im_j_k, b_ip_j_k, b_i_jm_k, b_i_jp_k, b_i_j_km, b_i_j_kp;
	double ev_im_j_k, ev_ip_j_k, ev_i_jm_k, ev_i_jp_k, ev_i_j_km, ev_i_j_kp;
	double visc_im_j_k, visc_ip_j_k, visc_i_jm_k, visc_i_jp_k, visc_i_j_km, visc_i_j_kp;
};
#endif

/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2023 Hans Bihs

This file is part of REEF3D.

REEF3D is fra->eps software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Fra->eps Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. Sa->eps the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, sa->eps <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
Author: Tobias Martin
--------------------------------------------------------------------*/

#include"increment.h"
#include"diffusion.h"

class rheology;
class heat;
class density;
class concentration;

using namespace std;

#ifndef IDIFF_IMEX_2D_H_
#define IDIFF_IMEX_2D_H_


class idiff_IMEX_2D : public diffusion, public increment
{

public:

	idiff_IMEX_2D(lexer*);
	idiff_IMEX_2D(lexer*,heat*,concentration*);
	virtual ~idiff_IMEX_2D();

	virtual void diff_scalar(lexer*, fdm*, ghostcell*, solver*, field&, field&, field&, double, double);
	virtual void diff_scalar(lexer*, fdm*, ghostcell*, solver*, field&, field&, field&, field&, double, double){};
    virtual void idiff_scalar(lexer*, fdm*, ghostcell*, solver*, field&, field&, double, double);
    
	virtual void diff_u(lexer*, fdm*, ghostcell*, solver*, field&, field&, field&, double);
	virtual void diff_v(lexer*, fdm*, ghostcell*, solver*, field&, field&, field&, double);
	virtual void diff_w(lexer*, fdm*, ghostcell*, solver*, field&, field&, field&, double);
    
    virtual void diff_u(lexer*, fdm*, ghostcell*, solver*, field&, field&, field&, field&, double){};
	virtual void diff_v(lexer*, fdm*, ghostcell*, solver*, field&, field&, field&, field&, double){};
	virtual void diff_w(lexer*, fdm*, ghostcell*, solver*, field&, field&, field&, field&, double){};
	
private:
    rheology *prheo;
    density *pd;
    
	double D;
	double time,starttime,endtime;
	int count,q;
	int gcval_u,gcval_v,gcval_w;
	double b_ijk,ev_ijk,visc_ijk;
	double b_im_j_k, b_ip_j_k, b_i_jm_k, b_i_jp_k, b_i_j_km, b_i_j_kp;
	double ev_im_j_k, ev_ip_j_k, ev_i_jm_k, ev_i_jp_k, ev_i_j_km, ev_i_j_kp;
	double visc_im_j_k, visc_ip_j_k, visc_i_jm_k, visc_i_jp_k, visc_i_j_km, visc_i_j_kp;
};
#endif


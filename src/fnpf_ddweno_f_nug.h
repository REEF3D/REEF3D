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

#ifndef FNPF_DDWENO_F_NUG_H_
#define FNPF_DDWENO_F_NUG_H_

#include"increment.h"
#include"weno_nug_func.h"

class fdm_fnpf;
class field;
class slice;
class lexer;
class ghostcell;
class vec;
class cpt;

using namespace std;

class fnpf_ddweno_f_nug : public weno_nug_func
{
public:

	 fnpf_ddweno_f_nug(lexer*,fdm_fnpf*);
	 ~fnpf_ddweno_f_nug();
    
    // field
    double ddwenox(field&, double);
    double ddwenoy(field&, double);
    double ddwenoz(field&, double);
     
    void iqmin(field&);
    void jqmin(field&);
    void kqmin(field&);
	void iqmax(field&);
	void jqmax(field&);
	void kqmax(field&);
     
    // slice
    double dswenox(slice&, double);
    double dswenoy(slice&, double);

    void isqmin(slice&);
	void jsqmin(slice&);
	void isqmax(slice&);
	void jsqmax(slice&);
    
    void is_wd_x_min();
    void is_wd_x_max();
    void is_wd_y_min();
    void is_wd_y_max();

    
    double grad;
    double *DX,*DY,*DZ;
    
private:
    lexer *p;
    fdm_fnpf *c;
};

#endif

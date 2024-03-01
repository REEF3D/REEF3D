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

#include"increment.h"
#include"weno_nug_func.h"

class fdm;
class field;
class lexer;
class ghostcell;
class vec;
class cpt;

#ifndef DDWENO_NUG_SF_H_
#define DDWENO_NUG_SF_H_

using namespace std;

class ddweno_nug_sf : public weno_nug_func
{
public:

	 ddweno_nug_sf(lexer*);
	 ~ddweno_nug_sf();

	 double ddwenox(fdm*, vec&, double, int, cpt&);
	 double ddwenoy(fdm*, vec&, double, int, cpt&);
	 double ddwenoz(fdm*, vec&, double, int, cpt&);

	void iqmin0(fdm*, vec&, cpt&);
	void jqmin0(fdm*, vec&, cpt&);
	void kqmin0(fdm*, vec&, cpt&);
	void iqmax0(fdm*, vec&, cpt&);
	void jqmax0(fdm*, vec&, cpt&);
	void kqmax0(fdm*, vec&, cpt&);
    
    void iqmin1(fdm*, vec&, cpt&);
	void jqmin1(fdm*, vec&, cpt&);
	void kqmin1(fdm*, vec&, cpt&);
	void iqmax1(fdm*, vec&, cpt&);
	void jqmax1(fdm*, vec&, cpt&);
	void kqmax1(fdm*, vec&, cpt&);
    
    void iqmin2(fdm*, vec&, cpt&);
	void jqmin2(fdm*, vec&, cpt&);
	void kqmin2(fdm*, vec&, cpt&);
	void iqmax2(fdm*, vec&, cpt&);
	void jqmax2(fdm*, vec&, cpt&);
	void kqmax2(fdm*, vec&, cpt&);
    
    void iqmin3(fdm*, vec&, cpt&);
	void jqmin3(fdm*, vec&, cpt&);
	void kqmin3(fdm*, vec&, cpt&);
	void iqmax3(fdm*, vec&, cpt&);
	void jqmax3(fdm*, vec&, cpt&);
	void kqmax3(fdm*, vec&, cpt&);
    
    void weight_min_sfcheck_x(fdm*);
    void weight_max_sfcheck_x(fdm*);

    double grad;
    double *DX,*DY,*DZ;
    
    int check1,check2,check3;
    
    int modus;
    
private:
    lexer *p;
};

#endif

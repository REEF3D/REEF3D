/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

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

	 double ddwenox(fdm*, field&, double);
	 double ddwenoy(fdm*, field&, double);
	 double ddwenoz(fdm*, field&, double);

	void iqmin0(fdm*, field&);
	void jqmin0(fdm*, field&);
	void kqmin0(fdm*, field&);
	void iqmax0(fdm*, field&);
	void jqmax0(fdm*, field&);
	void kqmax0(fdm*, field&);
    
    void iqmin1(fdm*, field&);
	void jqmin1(fdm*, field&);
	void kqmin1(fdm*, field&);
	void iqmax1(fdm*, field&);
	void jqmax1(fdm*, field&);
	void kqmax1(fdm*, field&);
    
    void iqmin2(fdm*, field&);
	void jqmin2(fdm*, field&);
	void kqmin2(fdm*, field&);
	void iqmax2(fdm*, field&);
	void jqmax2(fdm*, field&);
	void kqmax2(fdm*, field&);
    
    void iqmin3(fdm*, field&);
	void jqmin3(fdm*, field&);
	void kqmin3(fdm*, field&);
	void iqmax3(fdm*, field&);
	void jqmax3(fdm*, field&);
	void kqmax3(fdm*, field&);

    double grad;
    double *DX,*DY,*DZ;
    
    int check1,check2,check3;
    
    int modus;
    
private:
    lexer *p;
};

#endif

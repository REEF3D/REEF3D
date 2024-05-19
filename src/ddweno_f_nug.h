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
class slice;
class lexer;
class ghostcell;
class vec;
class cpt;

#ifndef DDWENO_F_NUG_H_
#define DDWENO_F_NUG_H_

using namespace std;

class ddweno_f_nug : public weno_nug_func
{
public:

	 ddweno_f_nug(lexer*);
	 ~ddweno_f_nug();

	 double ddwenox(field&, double);
	 double ddwenoy(field&, double);
	 double ddwenoz(field&, double);
     
     double dswenox(slice&, double);
	 double dswenoy(slice&, double);


	void iqmin(lexer*, field&);
	void jqmin(lexer*, field&);
	void kqmin(lexer*, field&);
	void iqmax(lexer*, field&);
	void jqmax(lexer*, field&);
	void kqmax(lexer*, field&);
    
    void isqmin(lexer*, slice&);
	void jsqmin(lexer*, slice&);
	void isqmax(lexer*, slice&);
	void jsqmax(lexer*, slice&);

    
    double grad;
    double *DX,*DY,*DZ;
    
private:
    lexer *p;
};

#endif

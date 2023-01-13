/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2023 Hans Bihs

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

#ifndef DDWENO_NUG_H_
#define DDWENO_NUG_H_

using namespace std;

class ddweno_nug : public weno_nug_func
{
public:

	 ddweno_nug(lexer*);
	 ~ddweno_nug();

	 double ddwenox(fdm*, vec&, double, int, cpt&);
	 double ddwenoy(fdm*, vec&, double, int, cpt&);
	 double ddwenoz(fdm*, vec&, double, int, cpt&);


	void iqmin(vec&, cpt&);
	void jqmin(vec&, cpt&);
	void kqmin(vec&, cpt&);
	void iqmax(vec&, cpt&);
	void jqmax(vec&, cpt&);
	void kqmax(vec&, cpt&);

    double grad;
    double *DX,*DY,*DZ;
    
private:
    lexer *p;
};

#endif

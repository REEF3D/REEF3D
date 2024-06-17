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

#ifndef DDWENO_NUG_SIG_H_
#define DDWENO_NUG_SIG_H_

#include"increment.h"
#include"weno_nug_func.h"

class fdm;
class field;
class lexer;
class ghostcell;
class vec;
class cpt;

using namespace std;

class ddweno_nug_sig : public weno_nug_func
{
public:

	 ddweno_nug_sig(lexer*);
	 ~ddweno_nug_sig();

	 double ddwenox(double*, double);
	 double ddwenoy(double*, double);
	 double ddwenoz(double*, double);


	void iqmin(double*);
	void jqmin(double*);
	void kqmin(double*);
	void iqmax(double*);
	void jqmax(double*);
	void kqmax(double*);

    double grad;
    double *DX,*DY,*DZ;
    
private:
    lexer *p;
};

#endif

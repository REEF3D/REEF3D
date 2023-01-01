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

class fdm;
class field;
class lexer;
class ghostcell;
class vec;
class cpt;

#ifndef DDWENO_H_
#define DDWENO_H_

using namespace std;

class ddweno : virtual public increment
{
public:

	 ddweno(lexer*);
	 ~ddweno();

	 double ddwenox(fdm*, vec&, double, cpt&);
	 double ddwenoy(fdm*, vec&, double, cpt&);
	 double ddwenoz(fdm*, vec&, double, cpt&);
	 

	 const double dx;

	double grad;


	void iqmin(vec&, double, cpt&);
	void jqmin(vec&, double, cpt&);
	void kqmin(vec&, double, cpt&);
	void iqmax(vec&, double, cpt&);
	void jqmax(vec&, double, cpt&);
	void kqmax(vec&, double, cpt&);

	const double tttw,fourth,third,sevsix,elvsix,sixth,fivsix,tenth;
	const double sixten,treten;
	const double epsilon,smallnum;
	double is1,is2,is3;
	double alpha1,alpha2,alpha3;
	double w1,w2,w3;
	double q1,q2,q3,q4,q5;

	void is();
	void alpha();
	void weight();
};

#endif

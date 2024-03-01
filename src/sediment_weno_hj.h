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

#include"sediment_exnerdisc.h"
#include"increment.h"
#include"ddweno_f_nug.h"

#ifndef SEDIMENT_WENO_HJ_H_
#define SEDIMENT_WENO_HJ_H_

using namespace std;

class sediment_weno_hj : public sediment_exnerdisc, public increment, public ddweno_f_nug
{
public:
	sediment_weno_hj(lexer*);
	virtual ~sediment_weno_hj();

    virtual double sx(lexer*, slice&, double,double);
	virtual double sy(lexer*, slice&, double,double);

private:
    double **ckz;
    
    double ivel,jvel;
};

#endif

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

#ifndef SEDIMENT_FOU_H_
#define SEDIMENT_FOU_H_

#include"sediment_exnerdisc.h"
#include"increment.h"
#include"weno_nug_func.h"

using namespace std;

class sediment_fou : public sediment_exnerdisc, public increment
{
public:
	sediment_fou(lexer*);
	virtual ~sediment_fou();

    virtual double sx(lexer*, slice&, double, double);
	virtual double sy(lexer*, slice&, double, double);


private:

    double ivel1,ivel2,jvel1,jvel2;
    double grad;
    
    double fu1,fu2,fv1,fv2;

};

#endif

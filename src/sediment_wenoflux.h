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
#include"weno_nug_func.h"

#ifndef SEDIMENT_WENOFLUX_H_
#define SEDIMENT_WENOFLUX_H_

using namespace std;

class sediment_wenoflux : public sediment_exnerdisc, public increment, public weno_nug_func
{
public:
	sediment_wenoflux(lexer*);
	virtual ~sediment_wenoflux();

    virtual double sx(lexer*, slice&, double, double);
	virtual double sy(lexer*, slice&, double, double);


private:
    double ffx(lexer *p, slice &f, double advec);
    double ffy(lexer *p, slice &f, double advec);

    
    void iqmin(lexer*, slice&);
	void jqmin(lexer*, slice&);
	void iqmax(lexer*, slice&);
	void jqmax(lexer*, slice&);
    
    
    double **ckz;
    double grad;
    
    double fu1,fu2,fv1,fv2;

};

#endif

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

#include"fnpf_etadisc.h"
#include"increment.h"
#include"weno_nug_func.h"

#ifndef FNPF_WENOFLUX_H_
#define FNPF_WENOFLUX_H_

using namespace std;

class fnpf_wenoflux : public fnpf_etadisc, public increment, public weno_nug_func
{
public:
	fnpf_wenoflux(lexer*);
	virtual ~fnpf_wenoflux();

    virtual double sx(lexer*, slice&, slice&);
	virtual double sy(lexer*, slice&, slice&);


private:
    double ffx(lexer *p, slice &f, double advec);
    double ffy(lexer *p, slice &f, double advec);

    
    void iqmin(lexer*, slice&);
	void jqmin(lexer*, slice&);
	void iqmax(lexer*, slice&);
	void jqmax(lexer*, slice&);
    
    
    double **ckz;
    double ivel1,ivel2,jvel1,jvel2;
    double grad;
    
    double fu1,fu2,fv1,fv2;

};

#endif

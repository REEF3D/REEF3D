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

#include"fnpf_ddx.h"
#include"increment.h"

using namespace std;

#ifndef FNPF_DDX_CDS4_H_
#define FNPF_DDX_CDS4_H_

class fnpf_ddx_cds4 : public fnpf_ddx, public increment
{
public:
    fnpf_ddx_cds4(lexer*);
	virtual ~fnpf_ddx_cds4();

    virtual double sxx(lexer*, slice&);
	virtual double syy(lexer*, slice&);
    
private:
    double X1,X2,X3,X4,X0;
    double Y1,Y2,Y3,Y4,Y0;
    double grad;

};

#endif

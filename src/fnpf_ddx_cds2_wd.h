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

#ifndef FNPF_DDX_CDS2_WD_H_
#define FNPF_DDX_CDS2_WD_H_

#include"fnpf_ddx.h"
#include"increment.h"

class fdm_fnpf;

using namespace std;

class fnpf_ddx_cds2_wd : public fnpf_ddx, public increment
{
public:
    fnpf_ddx_cds2_wd(lexer*,fdm_fnpf*);
	virtual ~fnpf_ddx_cds2_wd();

    virtual double sxx(lexer*, slice&);
	virtual double syy(lexer*, slice&);
    
private:
    fdm_fnpf *c;

};

#endif

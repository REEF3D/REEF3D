/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2023 Hans Bihs

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

#include"fnpf_laplace.h"
#include"increment.h"

class fnpf_bed_update;

#ifndef LAPLACE_FNPF_CDS2_H_
#define LAPLACE_FNPF_CDS2_H_

using namespace std;

class fnpf_laplace_cds2 : public fnpf_laplace, public increment
{
public:
    fnpf_laplace_cds2 (lexer*);
	virtual ~fnpf_laplace_cds2();

    virtual void start(lexer *,fdm_fnpf*,ghostcell*,solver*,fnpf_fsf*,double*,slice&);
    
private:
    
    fnpf_bed_update *pbed;
    double denom;

};

#endif

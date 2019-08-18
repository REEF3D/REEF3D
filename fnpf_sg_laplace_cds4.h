/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2019 Hans Bihs

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
--------------------------------------------------------------------*/

#include"fnpf_sg_laplace.h"
#include"increment.h"

#ifndef FNPF_SG_LAPLACE_CDS4_H_
#define FNPF_SG_LAPLACE_CDS4_H_

using namespace std;

class fnpf_sg_laplace_cds4 : public fnpf_sg_laplace, public increment
{
public:
    fnpf_sg_laplace_cds4 (lexer*);
	virtual ~fnpf_sg_laplace_cds4();

    virtual void start(lexer *,fdm_fnpf*,ghostcell*,solver*,fnpf_sg_fsf*,double*);
    
private:
    
    double **ckx,**cky,**ckz;
};

#endif
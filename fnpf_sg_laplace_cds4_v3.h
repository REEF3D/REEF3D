/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2018 Hans Bihs

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

#ifndef FNPF_SG_LAPLACE_CDS4_V3_H_
#define FNPF_SG_LAPLACE_CDS4_V3_H_

using namespace std;

class fnpf_sg_laplace_cds4_v3 : public fnpf_sg_laplace, public increment
{
public:
    fnpf_sg_laplace_cds4_v3 ();
	virtual ~fnpf_sg_laplace_cds4_v3();

    virtual void start(lexer *,fdm_fnpf*,ghostcell*,solver*,fnpf_sg_fsfbc*,double*);
    
private:
    
    double X1,X2,X3,X4,X0;
    double Y1,Y2,Y3,Y4,Y0;
    double Z1,Z2,Z3,Z4,Z0;
};

#endif
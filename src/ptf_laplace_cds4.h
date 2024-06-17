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

#ifndef PTF_LAPLACE_CDS4_H_
#define PTF_LAPLACE_CDS4_H_

#include"ptf_laplace.h"
#include"increment.h"

using namespace std;

class ptf_laplace_cds4 : public ptf_laplace, public increment
{
public:
    ptf_laplace_cds4 ();
	virtual ~ptf_laplace_cds4();

    virtual void start(lexer *,fdm*,ghostcell*,solver*,field&,slice&);
    
private:
    
    double X1,X2,X3,X4,X0;
    double Y1,Y2,Y3,Y4,Y0;
    double Z1,Z2,Z3,Z4,Z0;

};

#endif

/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

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

#ifndef SOLVER_VOID_H_
#define SOLVER_VOID_H_

#include"solver.h"

class solver_void : public solver
{
public:
    solver_void();
    virtual ~solver_void() = default;
    void start(lexer*,fdm*,ghostcell*,field&,vec&,int) override;
    void startf(lexer*,ghostcell*,field&,vec&,matrix_diag&,int) override;
    void startF(lexer*,ghostcell*,double*,vec&,matrix_diag&,int) override;
    void startV(lexer*,ghostcell*,double*,vec&,matrix_diag&,int) override;
    void startM(lexer*,ghostcell*,double*,double*,double*,int) override;
};

#endif

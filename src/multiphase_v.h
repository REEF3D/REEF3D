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

#ifndef MULTIPHASE_V_H_
#define MULTIPHASE_V_H_

#include"multiphase.h"
#include<fstream>

class multiphase_v : public multiphase
{
public:
    multiphase_v();
    virtual ~multiphase_v() = default;
    void start(lexer*,fdm*,ghostcell*,convection*,solver*,ioflow*,reini*,particle_corr*,printer*) override;
    void ini(lexer*,fdm*,ghostcell*,ioflow*,printer*,convection*,solver*) override;
    virtual void update(lexer*,fdm*,ghostcell*) override;
    
    void print_3D(lexer*,fdm*,ghostcell*,ofstream&) override;
    void print_file(lexer*,fdm*,ghostcell*) override;
    double ls1val(int,int,int) override;
    double ls2val(int,int,int) override;
    double ccipol_ls1val(lexer*,ghostcell*,double,double,double) override;
    double ccipol_ls2val(lexer*,ghostcell*,double,double,double) override;
    void ls1get(int,int,int,double) override;
    void ls2get(int,int,int,double) override;

    void name_pvtu(lexer*,fdm*,ghostcell*,ofstream&) override;
    void name_vtu(lexer*,fdm*,ghostcell*,ofstream&,int*,int&) override;
    void offset_vtu(lexer*,fdm*,ghostcell*,ofstream&,int*,int&) override;
};

#endif

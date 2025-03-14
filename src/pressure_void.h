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

#ifndef PRESSURE_VOID_H_
#define PRESSURE_VOID_H_

#include"pressure.h"
#include"increment.h"

class pressure_void : public pressure, public increment
{

public:
    pressure_void() = default;
    virtual ~pressure_void() = default;

    void start(lexer*,fdm*,ghostcell*,ioflow*,solver*,field&,field&,field&,double) override;
    void ini(lexer*,fdm*,ghostcell*) override;
    void upgrad(lexer*,fdm*,slice&,slice&) override;
    void vpgrad(lexer*,fdm*,slice&,slice&) override;
    void wpgrad(lexer*,fdm*,slice&,slice&) override;
    void ucorr(lexer*,fdm*,field&,double) override;
    void vcorr(lexer*,fdm*,field&,double) override;
    void wcorr(lexer*,fdm*,field&,double) override;
};

#endif

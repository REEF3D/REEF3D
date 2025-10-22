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

#ifndef NHFLOW_VTP_BED_H_
#define NHFLOW_VTP_BED_H_

#include"increment.h"

class lexer;
class fdm_nhf;
class ghostcell;
class sediment;

using namespace std;

class nhflow_vtp_bed : public increment
{
public:
    nhflow_vtp_bed(lexer*);
    virtual ~nhflow_vtp_bed() = default;

    void start(lexer*,fdm_nhf*,ghostcell*,sediment*);

private:
    void print2D(lexer*,fdm_nhf*,ghostcell*,sediment*);
    void pvtp(lexer*,sediment*,int);

    char name[200];
    int n,iin,offset[200];
    float ffn;

    int printcount;
};

#endif

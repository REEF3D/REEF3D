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

#ifndef SFLOW_VTP_BED_H_
#define SFLOW_VTP_BED_H_

#include"increment.h"

class lexer;
class fdm2D;
class ghostcell;
class sediment;

using namespace std;

class sflow_vtp_bed : public increment
{
public:
    sflow_vtp_bed(lexer*);
    virtual ~sflow_vtp_bed() = default;
    void start(lexer*,fdm2D*,ghostcell*,sediment*);

private:
    void print2D(lexer*,fdm2D*,ghostcell*,sediment*);
    void pvtp(lexer*,sediment*,int);

    char name[200];
    int n,iin,offset[200];
    float ffn;
    int printbedcount;
    double printbedtime;
};
#endif

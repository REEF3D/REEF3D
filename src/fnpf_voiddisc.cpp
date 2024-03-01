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

#include"fnpf_voiddisc.h"
#include"lexer.h"
#include"vec.h"

fnpf_voiddisc::fnpf_voiddisc(lexer* p)
{
}

fnpf_voiddisc::~fnpf_voiddisc()
{
}

double fnpf_voiddisc::fx(lexer *p, field &f, double ivel1, double ivel2)
{
    return 0.0;
}

double fnpf_voiddisc::fy(lexer *p, field &f, double jvel1, double jvel2)
{
    return 0.0;
}

double fnpf_voiddisc::fz(lexer *p, field &f, double kvel1, double kvel2)
{
    return 0.0;
}

double fnpf_voiddisc::sx(lexer *p, slice &f, double ivel)
{
    return 0.0; 
}

double fnpf_voiddisc::sy(lexer *p, slice &f, double jvel)
{
    return 0.0;  
}

double fnpf_voiddisc::sz(lexer *p, double *f)
{
    return 0.0;
}

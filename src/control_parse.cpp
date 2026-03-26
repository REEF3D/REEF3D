/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2026 Hans Bihs

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
Author: Alexander Hanke
--------------------------------------------------------------------*/

#include "control.h"
#include "lexer.h"

void control::parse(lexer* p)
{
    if(B98==1)
    B98=2;
    else if(B98==3 && A10==5)
    B98=4;

    if(X10>0)
    D22=2;

    if(F80>0 && F35>0)
    F35=0;

    G1=0;
    if(S10>0 || p->toporead>0 || p->solidread==1)
    G1=1;

    if(A10==3 || A10==5)
    G2=1;

    if(I10==1)
    {
        I11=1;
        I12=2;
        I13=1;
    }
    else if(I10==2)
    {
        I11=2;
        I12=2;
        I13=1;
    }

    if(I40>0)
    {
        I10=0;
        I11=0;
        I12=0;
        I13=0;
    }

    if(T10==0)
    I13=0;

    if((N40==3 || N40==23 ) && X10>0)
    N40=4;
    else if(N40==13 && X10>0)
    N40=14;

    if(S10>=1 || p->toporead==1)
    P27=1;
}

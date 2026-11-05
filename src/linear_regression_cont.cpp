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
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"linear_regression_cont.h"
#include"lexer.h"
#include"ghostcell.h"

linear_regression_cont::linear_regression_cont (lexer *p)
{
    num=1.0;
    xm=ym=0.0;
    SS_xy=SS_xx=0.0; 
    xsum=ysum=0.0;
    xx_sum=xy_sum=0.0;
    
}

linear_regression_cont::~linear_regression_cont()
{
}


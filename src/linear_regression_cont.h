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

#ifndef LINEAR_REGRESSION_CONT_H_
#define LINEAR_REGRESSION_CONT_H_

#include"increment.h"

class lexer;
class ghostcell;
class flux;

using namespace std;

class linear_regression_cont : public increment
{

public:

	linear_regression_cont (lexer *);
	virtual ~linear_regression_cont();

	void linreg_cont_func(lexer*,ghostcell*,double, double, double &b0, double &b1);

private:
    double num;
    double xm,ym;
    double xsum,ysum;
    double xy_sum,xx_sum;
    double SS_xy,SS_xx;
  
};

#endif


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
Author: Filip Hahs
--------------------------------------------------------------------*/

#ifndef SIXDOF_FORCESEXT_FILE_H_
#define SIXDOF_FORCESEXT_FILE_H_

#include "6DOF_forcesext.h"
#include <Eigen/Dense>
#include <fstream>

class lexer;
class ghostcell;

using namespace std;

class sixdof_forcesext_file final : public sixdof_forcesext
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

    void forcesext_trans(lexer *, ghostcell *) override final;
    void forcesext_rot(lexer *, ghostcell *) override final;
    void ini(lexer *, ghostcell *) override final;
    sixdof_forcesext_file(lexer *, ghostcell *);
    // virtual ~sixdof_forcesext_file();

private:
    void read_format(lexer *, ghostcell *);
    ofstream file;
    char name[200];
    int qn, count, ptnum;
    int rowcount, colcount;
    int colnum;
    double val;
    double **data;
    double ts, te;
    int timecount, timecount_old;
};
#endif

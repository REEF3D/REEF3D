/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

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

#include"6DOF_motionext.h"
#include<fstream>
#include <Eigen/Dense>

class lexer;
class fdm;
class fdm_nhf;
class fdm2D;
class ghostcell;
class vrans;
class net;
class field;

using namespace std;

#ifndef SIXDOF_MOTIONEXT_FILE_H_
#define SIXDOF_MOTIONEXT_FILE_H_

class sixdof_motionext_file : public sixdof_motionext
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
    
    virtual void motionext_trans(lexer*, ghostcell*, Eigen::Vector3d&, Eigen::Vector3d&);
    virtual void motionext_rot(lexer*, Eigen::Vector3d&, Eigen::Vector3d&, Eigen::Vector4d&, Eigen::Matrix<double, 3, 4>&,  Eigen::Matrix3d&);

    virtual void ini(lexer*,ghostcell*);
    
    
    
    sixdof_motionext_file(lexer*, ghostcell*);
	virtual ~sixdof_motionext_file();
    
private:
    double ramp_vel(lexer*);
    double ramp_draft(lexer*);
    
    void read_format_1(lexer*,ghostcell*);
    void read_format_2(lexer*,ghostcell*);
    
    ofstream file;
    char name[200];
    int qn,count,ptnum;
    int rowcount,colcount;
    int colnum;
    double val;
    double **data;
    double ts,te;
    int timecount,timecount_old;
    
    double Uext, Vext, Wext, Pext, Qext, Rext;
};

#endif

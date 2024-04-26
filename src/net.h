/*--------------------------------------------------------------------
REEF3D
Copyright 2018-2024 Tobias Martin

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
Author: Tobias Martin
--------------------------------------------------------------------*/

#include<vector>
#include <Eigen/Dense>

class lexer;
class fdm;
class ghostcell;
class sixdof;

using namespace std;

#ifndef NET_H_
#define NET_H_

class net
{    
public:
    
    typedef vector<Eigen::Vector3d> EigenMat;
        
	virtual void start(lexer*, fdm*, ghostcell*, double,Eigen::Matrix3d)=0;
	virtual void initialize(lexer*, fdm*, ghostcell*)=0;	
	virtual void netForces(lexer*, double&, double&, double&, double&, double&, double&)=0;
    
    virtual const EigenMat& getLagrangePoints()=0;
    virtual const EigenMat& getLagrangeForces()=0;
    virtual const EigenMat& getCollarVel()=0;
    virtual const EigenMat& getCollarPoints()=0;
};

#endif

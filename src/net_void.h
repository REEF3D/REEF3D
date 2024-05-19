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

#include"net.h"
#include"field1.h"
#include"field2.h"
#include"field3.h"
#include"field4.h"
#include"field5.h"
#include"fieldint5.h"
#include"vec.h"
#include<fstream>
#include<iostream>
#include<vector>

using namespace std;

#ifndef NET_VOID_H_
#define NET_VOID_H_

class net_void : public net
{
public:

	virtual void start(lexer*, fdm*, ghostcell*, double,Eigen::Matrix3d);
	virtual void initialize(lexer*, fdm*, ghostcell*);
	virtual void netForces(lexer*, double&, double&, double&, double&, double&, double&);
    virtual const EigenMat& getLagrangePoints(){return lagrangePoints;} 
    virtual const EigenMat& getLagrangeForces(){return lagrangeForces;} 
    virtual const EigenMat& getCollarVel(){return collarVel;} 
    virtual const EigenMat& getCollarPoints(){return collarPoints;} 
 

private:
  
    vector<Eigen::Vector3d> lagrangePoints;    
    vector<Eigen::Vector3d> lagrangeForces;    
    vector<Eigen::Vector3d> collarVel;
    vector<Eigen::Vector3d> collarPoints;
};

#endif

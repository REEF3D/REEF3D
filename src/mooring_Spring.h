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

#include"mooring.h"
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

#ifndef MOORING_SPRING_H_
#define MOORING_SPRING_H_

class mooring_Spring : public mooring
{
public:
    mooring_Spring(int);
    virtual ~mooring_Spring();
    
    virtual void start(lexer*, ghostcell*);
    virtual void initialize(lexer*, ghostcell*);
    virtual void mooringForces(double&, double&, double&);
    
private:    

    // Print
    void print(lexer*);
    
    // --------------------
    
    // Line number
    int line;
    
    // Material constants
    double L0, k;
    
    // Mesh
    double dx, dy, dz, L;

    // Forces
    double T0, T;    
    double Xme_, Yme_, Zme_;
    
    // Print
    char name[100];
    ofstream eTout;
    
    // Breaking
    bool broken;
    double breakTension, breakTime, curr_time;
    double printtime;
};

#endif

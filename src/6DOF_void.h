/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2023 Hans Bihs

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

#include"6DOF.h"
#include<vector>
#include<fstream>
#include<iostream>
#include <Eigen/Dense>

class lexer;
class fdm;
class ghostcell;
class mooring;
class net;

using namespace std;

#ifndef SIXDOF_VOID_H_
#define SIXDOF_VOID_H_

class sixdof_void : public sixdof
{
public:
	sixdof_void();
	virtual ~sixdof_void();
    
    virtual void start(lexer*,ghostcell*);
	virtual void ini(lexer*,ghostcell*);

    
    virtual void isource(lexer*,fdm*,ghostcell*);
    virtual void jsource(lexer*,fdm*,ghostcell*);
    virtual void ksource(lexer*,fdm*,ghostcell*);
    
    virtual void isource(lexer*,fdm_nhf*,ghostcell*);
    virtual void jsource(lexer*,fdm_nhf*,ghostcell*);
    virtual void ksource(lexer*,fdm_nhf*,ghostcell*);
    
    virtual void isource2D(lexer*,fdm2D*,ghostcell*);
    virtual void jsource2D(lexer*,fdm2D*,ghostcell*);
    
private:
    Eigen::Matrix3d quatRotMat;

    vector<mooring*> pmooring;

	vector<double> Xme, Yme, Zme, Kme, Mme, Nme;    
	vector<double> Xne, Yne, Zne, Kne, Mne, Nne;    
};

#endif

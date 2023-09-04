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

#include"FSI.h"
#include<vector>
#include <Eigen/Dense>

class lexer;
class fdm;
class ghostcell;
class fsi_strip;

using namespace std;

#ifndef FSI_STRIPS_H_
#define FSI_STRIPS_H_

class fsi_strips : public fsi
{
public:
	
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
	
    fsi_strips(lexer*,ghostcell*);
	virtual ~fsi_strips();
	virtual void start(lexer*,fdm*,ghostcell*);
	virtual void initialize(lexer*,fdm*,ghostcell*);
    virtual void forcing(lexer*,fdm*,ghostcell*,double,field&,field&,field&,field&,field&,field&,bool);
    
    
private:
    int numberStrips;
    double starttime, endtime;
    double starttime0;
	
    vector<fsi_strip*> pstrip;
};

#endif

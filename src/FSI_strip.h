/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2021 Hans Bihs

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
#include"beam.h"

class lexer;
class fdm;
class ghostcell;

using namespace std;

#ifndef FSI_STRIP_H_
#define FSI_STRIP_H_

class fsi_strip : public beam
{
public:
	
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
	
    fsi_strip(int);
	virtual ~fsi_strip();
	virtual void start(lexer*,fdm*,ghostcell*);
	virtual void initialize(lexer*,fdm*,ghostcell*);
    
    
private:
    int nstrip;
};

#endif

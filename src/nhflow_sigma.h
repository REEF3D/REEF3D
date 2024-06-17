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

#ifndef NHFLOW_SIGMA_H_
#define NHFLOW_SIGMA_H_

#include"fnpf.h"
#include"nhflow_gradient.h"

class lexer;
class fdm_nhf;
class ghostcell;
class field;
class nhflow_sigma_data;
class slice;

using namespace std;

class nhflow_sigma : public nhflow_gradient
{
public:
	nhflow_sigma(lexer*);
	virtual ~nhflow_sigma();
    
    virtual void sigma_coord_ini(lexer*);
    virtual void sigma_ini(lexer*, fdm_nhf*, ghostcell*, slice&);
    virtual void sigma_update(lexer*, fdm_nhf*, ghostcell*, slice&);
    
    void omega_update(lexer*,fdm_nhf*,ghostcell*,slice&,double*,double*,double*);

        
private:
    
    void disc_bed(lexer*, fdm_nhf*, ghostcell*);
    void disc_eta(lexer*, fdm_nhf*, ghostcell*);
    
    double sig;
};

#endif

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

#include"nhflow_diffusion.h"
#include"increment.h"

#ifndef NHFLOW_IDIFF_H_
#define NHFLOW_IDIFF_H_

using namespace std;

class nhflow_idiff : public nhflow_diffusion, public increment
{
public:
    nhflow_idiff(lexer*);
	virtual ~nhflow_idiff();

	virtual void diff_u(lexer*, fdm_nhf*, ghostcell*, solver*, double*, double*, double*, double*, double*, slice&, double);
	virtual void diff_v(lexer*, fdm_nhf*, ghostcell*, solver*, double*, double*, double*, double*, double*, slice&, double);
    virtual void diff_w(lexer*, fdm_nhf*, ghostcell*, solver*, double*, double*, double*, double*, double*, slice&, double);
    virtual void diff_scalar(lexer*, fdm_nhf*, ghostcell*, solver*, double*, double);
    
private:
    int gcval_u,gcval_v,gcval_w;
    int gcval_uh,gcval_vh,gcval_wh;
    
    double time,starttime,endtime;
    
    double visc;
    double sigxyz2;

};

#endif

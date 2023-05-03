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
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"nhflow_hxy_disc.h"
#include"nhflow_hxy_flux.h"
#include"increment.h"

#ifndef NHFLOW_HXY_FOU_H_
#define NHFLOW_HXY_FOU_H_

using namespace std;

class nhflow_hxy_fou : public nhflow_hxy_disc, nhflow_hxy_flux
{
public:
	nhflow_hxy_fou(lexer*,patchBC_interface*);
	virtual ~nhflow_hxy_fou();

	virtual void start(lexer*,slice&,slice&,slice&,int*,slice&,slice&,slice&);

private:

    double ivel1,ivel2,jvel1,jvel2;
    int qq;
    
    patchBC_interface *pBC;
};

#endif

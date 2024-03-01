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

#include"increment.h"
#include<vector>

class lexer;
class fdm;
class ghostcell;
class field;
class turbulence;
class sixdof;
class vrans;
class mooring;
class net;
class fsi;

#ifndef MOMENTUM_FORCING_H_
#define MOMENTUM_FORCING_H_

using namespace std;

class momentum_forcing : public increment
{
public:
	momentum_forcing(lexer*);
	virtual ~momentum_forcing();
	void momentum_forcing_start(fdm*,lexer*,ghostcell*, sixdof*, vrans*, vector<net*>&, fsi*,
                                field&,field&,field&,field&,field&,field&,int,double,bool);

private:
	double uplus,ks_plus,dist,ks,ustar;
	int ii,jj,kk;
	double value;
	int gcval_u,gcval_v,gcval_w;
    double starttime, endtime;
};
#endif

/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2020 Hans Bihs

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

class lexer;
class fdm;
class ghostcell;
class momentum;
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
	virtual void start(lexer*,fdm*,ghostcell*,vrans*,vector<net*>&);
	virtual void initialize(lexer*,fdm*,ghostcell*,vector<net*>&);
    
private:

    vector<mooring*> pmooring;

	vector<double> Xme, Yme, Zme, Kme, Mme, Nme;    
	vector<double> Xne, Yne, Zne, Kne, Mne, Nne;    
};

#endif

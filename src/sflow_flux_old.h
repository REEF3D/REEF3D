/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2023 Hans Bihs

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
/*
#include"increment.h"
#include"fdm2D.h"

class lexer;
class slice;

#ifndef SFLOW_FLUX_H_
#define SFLOW_FLUX_H_

using namespace std;

class sflow_flux : public increment
{
public:

	sflow_flux (lexer*);
	virtual ~sflow_flux();
    
    void iflux(int,slice&);
	void jflux(int,slice&);

	void ifluxC(fdm2D*,int,slice&);
	void jfluxC(fdm2D*,int,slice&);
	
	void ifluxHJ(int,slice&);
	void jfluxHJ(int,slice&);
    

	double iadvec,jadvec;
	double ivel1,jvel1;
	double ivel2,jvel2;

};

#endif
*/

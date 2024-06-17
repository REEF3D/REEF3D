/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

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

#ifndef DENSITY_HEAT_H_
#define DENSITY_HEAT_H_

#include"density.h"
#include"increment.h"

class fdm;
class lexer;
class heat;


using namespace std;

class density_heat : public density, virtual public increment
{

public:
    density_heat(lexer*,heat*&);
	virtual ~density_heat();

	virtual double roface(lexer*,fdm*,int,int,int);
	
	double H,roval,phival;
	int ii,jj,kk;
	const double epsi,eps;
    double psi;
    
    heat *pheat;
private:
    void material(lexer*);
    double material_ipol(double**,int,double);

    static int iocheck,iter;
    //--
    double visc_1,visc_2,ro_1,ro_2,alpha_air,alpha_water;
	double **water_density;
	double **water_viscosity;
	double **air_density;
	double **air_viscosity;

	int water_density_num;
	int water_viscosity_num;
	int air_density_num;
	int air_viscosity_num;

};

#endif

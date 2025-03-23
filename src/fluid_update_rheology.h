/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

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

#ifndef FLUID_UPDATE_RHEOLOGY_H_
#define FLUID_UPDATE_RHEOLOGY_H_

#include"fluid_update.h"
#include"increment.h"

class lexer;
class fdm;
class ghostcell;
class rheology;

class fluid_update_rheology : public fluid_update, increment
{
public:
<<<<<<< HEAD
    fluid_update_rheology(lexer*);
    virtual ~fluid_update_rheology();

    void start(lexer*, fdm*, ghostcell*) override;

private:
    rheology *prheo;
    int iter;
    int n;
    const double ro1,ro2;
    const double visc2;
    double visc1;
=======
    fluid_update_rheology(lexer*, fdm*);
	virtual ~fluid_update_rheology();

	virtual void start(lexer*, fdm*, ghostcell*);

private:
	rheology *prheo;
	static int iocheck,iter;
    int gcval_ro,gcval_visc;
	int n;
	const double dx,ro1,visc2,ro2;
	double visc1;
>>>>>>> parent of 516fad2a7 (Replaced \t with 4 spaces)
    double epsi;

    bool iocheck;

};

#endif

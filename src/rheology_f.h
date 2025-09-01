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
Author: Hans Bihs, Alexander Hanke
--------------------------------------------------------------------*/

#ifndef RHEOLOGY_F_H_
#define RHEOLOGY_F_H_

#include"rheology.h"
#include"strain.h"
#include"field4.h"

class rheology_f : public rheology, public strain
{

public:

    rheology_f(lexer*);
    virtual ~rheology_f()=default;

    double viscosity(lexer*,fdm*,ghostcell*) override;
    
    void u_source(lexer*,fdm*) override;
    void v_source(lexer*,fdm*) override;
    void w_source(lexer*,fdm*) override;
    
    void filltau(lexer*,fdm*,ghostcell*) override;

private:
    double Herschel_Bulkley(lexer*,fdm*,ghostcell*);
    double Mohr_Coulomb_and_Herschel_Bulkley(lexer*,fdm*,ghostcell*);
    double heaviside(int);
    double yield_stress(lexer*,fdm*);
    void yieldStressGradient(lexer*,fdm*,int,int,int);
    void pressurePhi(lexer*,fdm*,int,int,int,bool=false);
    void pressurePhiGradient(lexer*,fdm*,int,int,int);

    field4 tau_x,tau_y,tau_z;
    
    double gamma;
    double val,f,H,phival,pressureval;
    
    double tau0;
    double tau0_p,tau0_m;
    double tanphi;

    double pressureval1,pressureval2;
    double tau01,tau02;
    
    const double epsi;
    const double gravity;
    const double density_interstitial_fluid;
    
    int count;
};
#endif

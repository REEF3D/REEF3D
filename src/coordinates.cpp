/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2026 Hans Bihs

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

#include"coordinates.h"
#include"lexer.h"

coordinates::coordinates(lexer *pp) 
{	
    p=pp;
}

coordinates::~coordinates()
{
}

double coordinates::Xin(double xworld, double yworld)
{
    double xmodel=xworld;
    
    xmodel = (xworld-p->global_orig_x)*cos(-p->alpha_grid) - (yworld-p->global_orig_y)*sin(-p->alpha_grid);
    
    return xmodel;
}

double coordinates::Yin(double xworld, double yworld)
{
    double ymodel=yworld;
    
    ymodel = (xworld-p->global_orig_x)*sin(-p->alpha_grid) + (yworld-p->global_orig_y)*cos(-p->alpha_grid);
    
    return ymodel;
}

double coordinates::Xout(double xmodel, double ymodel)
{
    double xworld=xmodel;
    
    xworld = (xmodel)*cos(p->alpha_grid) - (ymodel)*sin(p->alpha_grid) + p->global_orig_x;
    
    return xworld;
}

double coordinates::Yout(double xmodel, double ymodel)
{
    double yworld=ymodel;
    
    yworld = (xmodel)*sin(p->alpha_grid) + (ymodel)*cos(p->alpha_grid) + p->global_orig_y;
    
    return yworld;
}

// 

double coordinates::Alpha_rad_in(double alpha_world)
{
    alpha_world -= p->alpha_grid;
    
    return alpha_world;
}

double coordinates::Alpha_rad_out(double alpha_model)
{
    alpha_model += p->alpha_grid;
    
    return alpha_model;
}

double coordinates::Alpha_deg_in(double alpha_world)
{
    alpha_world -= p->alpha_grid*180.0/PI;
    
    return alpha_world;
}

double coordinates::Alpha_deg_out(double alpha_model)
{
    alpha_model += p->alpha_grid*180.0/PI;
    
    return alpha_model;
}
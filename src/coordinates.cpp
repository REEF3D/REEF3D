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

void coordinates::XYin(double &Xcoor, double &Ycoor)
{
    if(p->cms_flag==1)
    {
    Xtemp = (Xcoor-p->global_orig_x)*cos(-p->alpha_grid) - (Ycoor-p->global_orig_y)*sin(-p->alpha_grid);
    Ytemp = (Xcoor-p->global_orig_x)*cos(-p->alpha_grid) - (Ycoor-p->global_orig_y)*sin(-p->alpha_grid);
    
    Xcoor = Xtemp;
    Ycoor = Ytemp;
    }
}

void coordinates::XYout(double &Xcoor, double &Ycoor)
{
    if(p->cms_flag==1)
    {
    Xtemp = (Xcoor)*cos(p->alpha_grid) - (Ycoor)*sin(p->alpha_grid) + p->global_orig_x;
    Ytemp = (Xcoor)*sin(p->alpha_grid) + (Ycoor)*cos(p->alpha_grid) + p->global_orig_y;
    
    Xcoor = Xtemp;
    Ycoor = Ytemp;
    }
}

// ----------------

double coordinates::Xin(double xworld, double yworld)
{
    double xmodel=xworld;

    if(p->cms_flag==1)
    xmodel = (xworld-p->global_orig_x)*cos(-p->alpha_grid) - (yworld-p->global_orig_y)*sin(-p->alpha_grid);
    
    return xmodel;
}

double coordinates::Yin(double xworld, double yworld)
{
    double ymodel=yworld;
    
    if(p->cms_flag==1)
    ymodel = (xworld-p->global_orig_x)*sin(-p->alpha_grid) + (yworld-p->global_orig_y)*cos(-p->alpha_grid);
    
    return ymodel;
}

double coordinates::Xout(double xmodel, double ymodel)
{
    double xworld=xmodel;
    
    if(p->cms_flag==1)
    xworld = (xmodel)*cos(p->alpha_grid) - (ymodel)*sin(p->alpha_grid) + p->global_orig_x;
    
    return xworld;
}

double coordinates::Yout(double xmodel, double ymodel)
{
    double yworld=ymodel;
    
    if(p->cms_flag==1)
    yworld = (xmodel)*sin(p->alpha_grid) + (ymodel)*cos(p->alpha_grid) + p->global_orig_y;
    
    return yworld;
}

// 

double coordinates::Alpha_rad_in(double alpha_world)
{
    if(p->cms_flag==1)
    alpha_world -= p->alpha_grid;
    
    return alpha_world;
}

double coordinates::Alpha_rad_out(double alpha_model)
{
    if(p->cms_flag==1)
    alpha_model += p->alpha_grid;
    
    return alpha_model;
}

double coordinates::Alpha_deg_in(double alpha_world)
{
    if(p->cms_flag==1)
    alpha_world -= p->alpha_grid*180.0/PI;
    
    return alpha_world;
}

double coordinates::Alpha_deg_out(double alpha_model)
{
    if(p->cms_flag==1)
    alpha_model += p->alpha_grid*180.0/PI;
    
    return alpha_model;
}
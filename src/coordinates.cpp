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
    
    
    
    return xmodel;
}

double coordinates::Yin(double xworld, double yworld)
{
    double ymodel=yworld;
    
    
    
    return ymodel;
}

double coordinates::Xout(double xmodel, double ymodel)
{
    double xworld=xmodel;
    
    
    
    return xworld;
}

double coordinates::Yout(double xmodel, double ymodel)
{
    double yworld=ymodel;
    
    
    
    return yworld;
}
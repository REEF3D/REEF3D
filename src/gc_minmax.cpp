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

#include"ghostcell.h"
#include"lexer.h"
#include"fdm.h"
#include"field.h"

double ghostcell::mini1(fdm* a,lexer* p, field& b)
{
double mini;

mini=1.0e18;

     LOOP
     {
      if(0.5*(a->u(i,j,k)+a->u(i-1,j,k))<mini)
      mini=0.5*(a->u(i,j,k)+a->u(i-1,j,k));
     }
    return mini;
}

double ghostcell::maxi1(fdm* a,lexer* p, field& b)
{
double maxi;

maxi=-1.0e18;

     LOOP
     {
      if(0.5*(a->u(i,j,k)+a->u(i-1,j,k))>maxi)
      maxi=0.5*(a->u(i,j,k)+a->u(i-1,j,k));
     }

   return maxi;
}

double ghostcell::mini2(fdm* a,lexer* p, field& b)
{
double mini;

mini=1.0e18;

     LOOP
     {
      if(0.5*(a->v(i,j,k)+a->v(i,j-1,k))<mini)
      mini=0.5*(a->v(i,j,k)+a->v(i,j-1,k));
     }
    return mini;
}

double ghostcell::maxi2(fdm* a,lexer* p, field& b)
{
double maxi;

maxi=-1.0e18;

     LOOP
     {
      if(0.5*(a->v(i,j,k)+a->v(i,j-1,k))>maxi)
      maxi=0.5*(a->v(i,j,k)+a->v(i,j-1,k));
     }

   return maxi;
}

double ghostcell::mini3(fdm* a,lexer* p, field& b)
{
double mini;

mini=1.0e18;

     LOOP
     {
      if(0.5*(a->w(i,j,k)+a->w(i,j,k-1))<mini)
      mini=0.5*(a->w(i,j,k)+a->w(i,j,k-1));
     }
    return mini;
}

double ghostcell::maxi3(fdm* a,lexer* p, field& b)
{
double maxi;

maxi=-1.0e18;

     LOOP
     {
      if(0.5*(a->w(i,j,k)+a->w(i,j,k-1))>maxi)
      maxi=0.5*(a->w(i,j,k)+a->w(i,j,k-1));
     }

   return maxi;
}


double ghostcell::mini4(fdm* a,lexer* p, field& b)
{
double mini;

mini=1.0e18;

     LOOP

     if(b(i,j,k)<mini)
     mini=b(i,j,k);

    return mini;
}

double ghostcell::maxi4(fdm* a,lexer* p, field& b)
{
double maxi;

maxi=-1.0e18;

    LOOP
    if(b(i,j,k)>maxi)
    maxi=b(i,j,k);


   return maxi;
}



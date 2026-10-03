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

#include"6DOF_obj.h"
#include"lexer.h"
#include"ghostcell.h"

void sixdof_obj::objects_create(lexer *p, ghostcell *pgc)
{
    // surface triangles of all input objects (geom), in the order box, cylinders, wedges,
    // hexahedra, wavemakers, STL
    geom.allocate(p);
    
    geom.primitives(p,pgc);
    
    if(p->X170==1)
    {
        tstart[entity_count]=tricount;
        piston(p,pgc,0);
        tricount+=12;
        tend[entity_count]=tricount;
        ++entity_count;
    }
    
    if(p->X172==1)
    {
        tstart[entity_count]=tricount;
        flap_double(p,pgc,0);
        tricount+=28;
        tend[entity_count]=tricount;
        ++entity_count;
    }
    
    if(p->X180==1)
    {
        geom.read_stl(p,pgc);
		++entity_count;
    }
	
    if(p->mpirank==0)
	cout<<"Surface triangles: "<<tricount<<endl;

    // Initialise STL geometric parameters
	geometry_stl(p,pgc);
	
    // Order Triangles for correct inside/outside orientation
    geom.orient(p,pgc);

    // Refine triangles: X 185 1-3 split, X 185 4 remesh
    if(p->X185>0 && p->X60!=2 && entity_count>0 && p->X170==0 && p->X171==0 && p->X172==0)
    {
        geom.refine(p,pgc);
        
        // holes closed by the remesher: volume and mass (X 21) or density (X 22) of the closed surface
        if(geom.holes_closed)
        {
            geometry_stl(p,pgc);
            
            if(p->mpirank==0)
            cout<<"  volume of the closed surface: "<<Vfb<<", mass: "<<Mass_fb<<", density: "<<Rfb<<endl<<endl;
        }
    }

    if(p->mpirank==0)
	cout<<"Refined surface triangles: "<<tricount<<endl;
}

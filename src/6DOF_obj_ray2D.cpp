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

#include"6DOF_obj_2D.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"fieldint.h"

void sixdof_obj_2D::ray_cast_2D(lexer *p, ghostcell *pgc)
{
    // curvilinear grid (G 50): no rays along the grid lines (6DOF_obj_ray2D_curv.cpp)
    if(p->geo_curv==1)
    {
    ray_cast_2D_curv(p,pgc);
    return;
    }

	SLICELOOP4
	{
        fsio(i,j)=1;
	    fs(i,j)=1.0e9;
        Rxmin(i,j)=1.0e9;
        Rxmax(i,j)=-1.0e9;
        Rymin(i,j)=1.0e9;
        Rymax(i,j)=-1.0e9;
        Ls(i,j) = 1.0e-6;
        Bs(i,j) = 1.0e-6;
        draft(i,j) = 0.0;
	}
	
    for(int rayiter=0; rayiter<2; ++rayiter)
    {
        if(rayiter==0)
        {
            ray_cast_2D_io_x(p,pgc,0,tricount);
            ray_cast_2D_io_ycorr(p,pgc,0,tricount);
        }
    
        if(rayiter==1)
        {
            pgc->gcslparax_int(p,fsio,1);
            
            ray_cast_2D_x(p,pgc,0,tricount);
            ray_cast_2D_y(p,pgc,0,tricount);
            ray_cast_2D_z(p,pgc,0,tricount);
        }
    }
    
	SLICELOOP4
    {
        if(fsio(i,j)==-1)
        fs(i,j) = -fabs(fs(i,j));
        
        if(fsio(i,j)==1)
        fs(i,j) = fabs(fs(i,j));
    }
    
	SLICELOOP4
	{
		if(fs(i,j) > 10.0*p->DXM)
		fs(i,j) = 10.0*p->DXM;
		
		if(fs(i,j) < -10.0*p->DXM)
		fs(i,j) = -10.0*p->DXM;
	}
    
    SLICELOOP4
	{
    if(Rxmax(i,j)>-1.0e9 && Rxmin(i,j)<1.0e9)
    Ls(i,j) = Rxmax(i,j)-Rxmin(i,j);
    
    if(Rymax(i,j)>-1.0e8 && Rymin(i,j)<1.0e8)
    Bs(i,j) = Rymax(i,j)-Rymin(i,j);
    }
    
	pgc->gcsl_start4(p,fs,50);
    pgc->gcsl_start4(p,draft,50);

    // X 411 1: the draft continued outside the hull outline over the smoothing width of the patch
    if(p->X411==1)
    draft_extension(p,pgc);
}

void sixdof_obj_2D::draft_extension(lexer *p, ghostcell *pgc)
{
    // X 411 1 (STL draft, X 400 10): the cells outside the hull outline (fs > 0) within the smoothing width of the
    // pressure patch (fs < psi, psi = X41 dx as in Hsolidface_2D) take the draft of their neighbour closest to the
    // hull. The patch H(-fs,psi)*draft then goes to zero without the jump of about draft/2 at the hull side, and the
    // smoothing adds outside what it removes inside, so the hull section stays the same on grids of different
    // resolution (with draft = 0 outside, the section a grid sees changes by O(dx) with the position of the cell
    // centres relative to the hull side). Cartesian grids without mesh refinement (not in sflow_amr_ship).
    const int sweeps = int(p->X41) + 2;
    double psi,fmin,val;
    int changed;

    for(int q=0; q<sweeps; ++q)
    {
        changed = 0;

        SLICELOOP4
        {
        psi = p->X41*(1.0/2.0)*(p->DXN[IP] + p->DYN[JP]);

            if(fs(i,j)>0.0 && fs(i,j)<psi && draft(i,j)<=0.0)
            {
            fmin = fs(i,j);
            val = 0.0;

            // the neighbour closest to the hull that has a draft
            if(draft(i-1,j)>0.0 && fs(i-1,j)<fmin)
            {
            fmin = fs(i-1,j);
            val = draft(i-1,j);
            }

            if(draft(i+1,j)>0.0 && fs(i+1,j)<fmin)
            {
            fmin = fs(i+1,j);
            val = draft(i+1,j);
            }

                if(p->j_dir==1)
                {
                if(draft(i,j-1)>0.0 && fs(i,j-1)<fmin)
                {
                fmin = fs(i,j-1);
                val = draft(i,j-1);
                }

                if(draft(i,j+1)>0.0 && fs(i,j+1)<fmin)
                {
                fmin = fs(i,j+1);
                val = draft(i,j+1);
                }
                }

                if(val>0.0)
                {
                draft(i,j) = val;
                ++changed;
                }
            }
        }

    pgc->gcsl_start4(p,draft,50);

    if(pgc->globalisum(changed)==0)
    break;
    }
}






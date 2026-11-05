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

#include"sediment_exner.h"
#include"sediment_curv.h"
#include"lexer.h"
#include"ghostcell.h"
#include"topo_relax.h"
#include"sediment_fdm.h"

void sediment_exner::topovel_curv(lexer* p, ghostcell *pgc, sediment_fdm *s, int form)
{
    // Exner on a curvilinear grid: (sum of the face fluxes times the face lengths)/area
    //   S 31 1: face flux qb_face * sg, sg = u.n/|u| the transport direction of the face
    //   S 31 2: face flux (qbx,qby)_face . n of the cell flux vector qb*u/|u|, closed faces none
    //   S 32 1: the upwind cell of the face (by sg), S 32 2: the mean of the two cells
    // On Cartesian nodes this is topovel1/topovel2 with sediment_fou/sediment_cds.
    double u,v,umag;

    bend_prep(p,pgc,s);

    if(form==2)
    {
        SLICEBASELOOP
        {
        qbx(i,j) = 0.0;
        qby(i,j) = 0.0;
        }
        
        SEDSLICELOOP
        {
        sediment_curv::cellvel(p,s,i,j,u,v);
        if(btd(i,j)!=0.0)
        {
        const double t = btd(i,j), c = 1.0/sqrt(1.0 + t*t), uu = u;
        u = (uu - t*v)*c;
        v = (v + t*uu)*c;
        }
        umag = sqrt(u*u + v*v);
        qbx(i,j) = (fabs(umag)>1.0e-10?u/fabs(umag):0.0)*s->qb(i,j);
        qby(i,j) = (fabs(umag)>1.0e-10?v/fabs(umag):0.0)*s->qb(i,j);
        }
        
        pgc->gcsl_start4(p,qbx,1);
        pgc->gcsl_start4(p,qby,1);
    }
    
    // flux through face (dir, ii,jj) between cell (ii,jj) and its +dir neighbour (ia,ja), sg its direction
    auto face_flux = [&](int dir, int ii, int jj, double sg) -> double
    {
        double nx,ny,len;
        sediment_curv::normal(p,dir,ii,jj,nx,ny,len);
        
        const int ia = ii + (dir==1?1:0);
        const int ja = jj + (dir==2?1:0);
        
        if(form==1)
        {
            double f;
            
            if(p->S32==2)
            f = 0.5*(s->qb(ii,jj) + s->qb(ia,ja));
            else
            f = (sg>=0.0) ? s->qb(ii,jj) : s->qb(ia,ja);
            
            return f*sg*len;
        }
        
        double fx,fy;
        
        if(p->S32==2)
        {
        fx = 0.5*(qbx(ii,jj) + qbx(ia,ja));
        fy = 0.5*(qby(ii,jj) + qby(ia,ja));
        }
        
        else
        {
        fx = (sg>=0.0) ? qbx(ii,jj) : qbx(ia,ja);
        fy = (sg>=0.0) ? qby(ii,jj) : qby(ia,ja);
        }
        
        return (fx*nx + fy*ny)*(sg!=0.0?1.0:0.0)*len;
    };
    
    double sgx1,sgx2,sgy1,sgy2;
    
    SEDSLICELOOP
	if(p->pos_x()>=p->S71 && p->pos_x()<=p->S72)
	{
        sgx1 = sediment_curv::face_sign(p,s,1,i-1,j,bend_td(p,1,i-1,j));
        sgx2 = sediment_curv::face_sign(p,s,1,i,j,bend_td(p,1,i,j));
        sgy1 = sediment_curv::face_sign(p,s,2,i,j-1,bend_td(p,2,i,j-1));
        sgy2 = sediment_curv::face_sign(p,s,2,i,j,bend_td(p,2,i,j));
        
        face_closure(p,sgx1,sgx2,sgy1,sgy2);
        
        const int ci=i, cj=j;
        const double div = (face_flux(1,ci,cj,sgx2) - face_flux(1,ci-1,cj,sgx1)
                          + face_flux(2,ci,cj,sgy2) - face_flux(2,ci,cj-1,sgy1))/p->AREA2D[IJ] + bend_div(p,s,ci,cj);
        
        s->vz(i,j) =  -p->S35*s->guard(i,j)*prelax->rf(p,pgc)*(1.0/(1.0-p->S24))*(div + susp_ED(p,pgc,s));

        if(p->DFBED[IJ]<0)
        s->vz(i,j) = 0.0;
	}
    
    pgc->gcsl_start4(p,s->vz,1);
}

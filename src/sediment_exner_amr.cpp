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

// The Exner step of start() in parts for the mesh refinement and the multi-block grids (sediment_exner.h).
// The face fluxes are those of topovel1/topovel2 with sediment_fou (S 32 1) and sediment_cds (S 32 2),
// one value per face (the two cells of a face see the same flux), on curvilinear grids those of
// topovel_curv.

void sediment_exner::amr_qb(lexer *p, ghostcell *pgc, sediment_fdm *s)
{
    if(p->S33==0)
    SEDSLICELOOP
    s->qb(i,j)=s->qbe(i,j);

    if(p->S33>0)
    non_equillibrium_solve(p,pgc,s);

    qb_clear(p,s);
    pgc->gcsl_start4(p,s->qb,1);

    bend_prep(p,pgc,s);

    if(p->S62==2)
    susp_qs(p,pgc,s);
}

double& sediment_exner::amr_field(sediment_fdm *s, int n, int ii, int jj)
{
    if(n==1)
    return bK(ii,jj);
    if(n==2)
    return btd(ii,jj);
    if(n==3)
    return bz(ii,jj);
    if(n==4)
    return bg(ii,jj);
    return s->qb(ii,jj);
}

void sediment_exner::amr_flux(lexer *p, ghostcell *pgc, sediment_fdm *s, slice &Fx, slice &Fy)
{
    const int curv = p->geo_curv;
    const int form = (p->S31==1) ? 1 : 2;
    double u,v,umag;

    // S 31 2: the cell flux vector qb u/|u|
    if(form==2)
    {
        SLICEBASELOOP
        {
        qbx(i,j) = 0.0;
        qby(i,j) = 0.0;
        }

        SEDSLICELOOP
        {
            if(curv==1)
            sediment_curv::cellvel(p,s,i,j,u,v);
            else
            {
            u = 0.5*(s->P(i,j)+s->P(i-1,j));
            v = 0.5*(s->Q(i,j)+s->Q(i,j-1));
            }

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

    // transport direction of a face (as topovel1/2 and sediment_curv::face_sign)
    auto sign_i = [&](int ii, int jj) -> double
    {
        if(curv==1)
        return sediment_curv::face_sign(p,s,1,ii,jj,bend_td(p,s,1,ii,jj));

        const double ux = s->P(ii,jj);
        const double vx = 0.25*(s->Q(ii+1,jj)+s->Q(ii,jj)+s->Q(ii+1,jj-1)+s->Q(ii,jj-1));
        const double ua = sqrt(ux*ux + vx*vx);
        return bend_turn(ux,vx,ua,bend_td(p,s,1,ii,jj));
    };

    auto sign_j = [&](int ii, int jj) -> double
    {
        if(curv==1)
        return sediment_curv::face_sign(p,s,2,ii,jj,bend_td(p,s,2,ii,jj));

        const double uy = 0.25*(s->P(ii,jj)+s->P(ii,jj+1)+s->P(ii-1,jj)+s->P(ii-1,jj+1));
        const double vy = s->Q(ii,jj);
        const double ua = sqrt(uy*uy + vy*vy);
        return bend_turn(vy,-uy,ua,bend_td(p,s,2,ii,jj));
    };

    auto dfbed = [&](int ii, int jj) { return s->DFBED[(ii-p->imin)*p->jmax + (jj-p->jmin)]; };

    // flux per unit length through the face between (ii,jj) and (ia,ja), direction sg
    auto flux = [&](int dir, int ii, int jj, int ia, int ja, double sg) -> double
    {
        if(dfbed(ii,jj)<0 || dfbed(ia,ja)<0)
        return 0.0;

        if(form==1)
        {
            double f;
            if(p->S32==2)
            f = 0.5*(s->qb(ii,jj) + s->qb(ia,ja));
            else
            f = (sg>=0.0) ? s->qb(ii,jj) : s->qb(ia,ja);

            return f*sg;
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

        const double open = (sg!=0.0) ? 1.0 : 0.0;

        if(curv==1)
        {
        double nx,ny,len;
        sediment_curv::normal(p,dir,ii,jj,nx,ny,len);
        return (fx*nx + fy*ny)*open;
        }

        return ((dir==1) ? fx : fy)*open;
    };

    SLICEBASELOOP
    {
    Fx(i,j) = 0.0;
    Fy(i,j) = 0.0;
    }

    // the faces of the cells of the grid: i-faces i = -1 .. knox-1, j-faces j = -1 .. knoy-1
    for(int ii=-1; ii<p->knox; ++ii)
    for(int jj=0; jj<p->knoy; ++jj)
    Fx(ii,jj) = flux(1,ii,jj,ii+1,jj,sign_i(ii,jj)) + bend_flux(p,s,1,ii,jj);

    if(p->j_dir==1)
    for(int ii=0; ii<p->knox; ++ii)
    for(int jj=-1; jj<p->knoy; ++jj)
    Fy(ii,jj) = flux(2,ii,jj,ii,jj+1,sign_j(ii,jj)) + bend_flux(p,s,2,ii,jj);
}

double sediment_exner::amr_vz(lexer *p, ghostcell *pgc, sediment_fdm *s, slice &Fx, slice &Fy, slice *mask)
{
    double div;

    SLICEBASELOOP
    s->vz(i,j) = 0.0;

    SEDSLICELOOP
	if(p->pos_x()>=p->S71 && p->pos_x()<=p->S72)
	{
        if(p->geo_curv==1)
        div = (Fx(i,j)*p->LEN1[IJ] - Fx(i-1,j)*p->LEN1[Im1J] + Fy(i,j)*p->LEN2[IJ] - Fy(i,j-1)*p->LEN2[IJm1])/p->AREA2D[IJ];
        else
        div = (Fx(i,j)-Fx(i-1,j))/p->DXN[IP] + (Fy(i,j)-Fy(i,j-1))/p->DYN[JP];

        s->vz(i,j) =  -p->S35*s->guard(i,j)*prelax->rf(p,pgc)*(1.0/(1.0-p->S24))*(div + susp_ED(p,pgc,s));

        if(s->DFBED[IJ]<0)
        s->vz(i,j) = 0.0;
	}

    pgc->gcsl_start4(p,s->vz,1);

    // the filter (S 100) runs over all grids in steps (sediment_xfilter, the driver), then amr_vzmax
    return 0.0;
}

double sediment_exner::amr_vzmax(lexer *p, sediment_fdm *s, slice *mask)
{
    double mx = 0.0;
    SEDSLICELOOP
    if(mask==nullptr || (*mask)(i,j)>0.5)
    mx = MAX(mx,fabs(s->vz(i,j)));

    return mx;
}

void sediment_exner::amr_bed(lexer *p, ghostcell *pgc, sediment_fdm *s)
{
    SEDSLICELOOP
    s->dh(i,j) = p->dtsed*s->vz(i,j);

    dry_deposit(p,s);

	SEDSLICELOOP
    s->bedzh(i,j) += s->dh(i,j);

	pgc->gcsl_start4(p,s->bedzh,1);
}

double sediment_exner::amr_dx(lexer *p, sediment_fdm *s)
{
    double dx = 1.0e8;

    SEDSLICELOOP
    {
        if(p->geo_curv==1)
        {
        dx = MIN(dx,sediment_curv::width(p,1,i,j));
        dx = MIN(dx,sediment_curv::width(p,2,i,j));
        }
        else
        {
        dx = MIN(dx,p->DXN[IP]);
        if(p->j_dir==1 && p->knoy>1)
        dx = MIN(dx,p->DYN[JP]);
        }
    }

    return dx;
}

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

#include"interface_width.h"
#include"cell_index.h"
#include"initialize.h"
#include"fdm.h"
#include"lexer.h"
#include"ghostcell.h"
#include<algorithm>
#include<cmath>
#include<vector>
#include"heaviside.h"

//  The level set of the start as a signed distance (positive in the water).
//
//  Every region of the water (F 50 - 59, F 70, F 72, F 112 - 115) is put in with its exact signed
//  distance, and the regions are combined as sets: a region of water by the union (max), a region of
//  air (F 71) by the subtraction (min with the negative distance), F 60 sets the whole field, F 72
//  replaces the column.  The zero level set is then exactly the boundary of the regions (not the cell
//  faces next to it) and the values next to it are the distance, which the fast sweeping of the start
//  (reini_fsm, from reini_RK3 at count 0) keeps and continues into the domain.  The max/min of
//  distances is a distance next to the faces of the regions; at their edges and corners it is lower
//  than the distance, the fast sweeping takes care of that.
//
//  In 2D (j_dir 0) a region is the cut of its 3D shape with the plane of the cells (the cells of the
//  region the same as before), the distance is taken in that plane (x, z).
namespace
{
const double no_water = -1.0e20;

// the cells whose density and viscosity of the start are smoothed over the interface (the level set of
// before was a distance there: F 60, F 62/F 63, F 72), sharp elsewhere (+1/-1 before); (i,j,k) of the
// partition and its ghost cells
std::vector<char> smooth;

// signed distance to the interval [x0,x1], positive outside (an empty interval x1 < x0: outside everywhere)
inline double outside1(double x, double x0, double x1)
{
    return std::max(x0-x, x-x1);
}

// signed distance of the box, positive inside; usey=false: without y (2D)
double sd_box(double x, double y, double z, double x0, double x1, double y0, double y1, double z0, double z1, bool usey)
{
    const double qx = std::max(x0-x, x-x1);
    const double qy = usey ? std::max(y0-y, y-y1) : -1.0e30;
    const double qz = std::max(z0-z, z-z1);

    const double ox = std::max(qx,0.0), oy = std::max(qy,0.0), oz = std::max(qz,0.0);

    return -(std::sqrt(ox*ox + oy*oy + oz*oz) + std::min(std::max(qx,std::max(qy,qz)),0.0));
}

// signed distance in the plane of two outside measures a, b (> 0 outside), positive inside: the
// intersection of two convex sets whose boundaries meet at a right angle (capped cylinder)
inline double sd_two(double a, double b)
{
    return -(std::min(std::max(a,b),0.0) + std::sqrt(std::max(a,0.0)*std::max(a,0.0) + std::max(b,0.0)*std::max(b,0.0)));
}

// union of the water with a region of the signed distance d
inline void add(double &phi, double d)
{
    phi = std::max(phi,d);
}

// subtraction of a region of the signed distance d (air)
inline void cut(double &phi, double d)
{
    phi = std::min(phi,-d);
}
}

void initialize::iniphi(lexer* p, fdm* a, ghostcell* pgc)
{
    double r;
    double phidiff, xdiff;
    p->phimean=p->F56;
    const bool jd = p->j_dir==1;

    // interface width used for the initial density and viscosity below; inipsi() refines it later
    // (p->psi was used here before it was set)
    if(p->j_dir==0)
    p->psi = p->F45*(1.0/2.0)*(p->DXM+p->DZM);

    if(p->j_dir==1)
    p->psi = p->F45*(1.0/3.0)*(p->DXM+p->DYM+p->DZM);


    smooth.assign(p->imax*p->jmax*p->kmax,0);

    // the ghost cells as before (from -1, for the density of the start below)
    LOOP
	a->phi(i,j,k)=-1.0;

	pgc->start4(p,a->phi,50);

    LOOP
	a->phi(i,j,k)=no_water;

    // box F 51 - 56; in 2D the cells of the box as before (the nodes in y)
    if(p->F50_flag==1)
	LOOP
    if(jd || (p->YN[JP]>=p->F52 && p->YN[JP]<p->F55))
    add(a->phi(i,j,k), sd_box(p->XP[IP],p->YP[JP],p->ZP[KP],p->F51,p->F54,p->F52,p->F55,p->F53,p->F56,jd));


    // half space a x + b y + c z < d
    if(p->F57_1>0||p->F57_2>0||p->F57_3>0||p->F57_4>0)
    {
        const double norm = std::sqrt(p->F57_1*p->F57_1 + (jd?p->F57_2*p->F57_2:0.0) + p->F57_3*p->F57_3);

        LOOP
        {
        const double s = p->F57_4 - (p->F57_1*p->XP[IP]+ p->F57_2*p->YP[JP]+ p->F57_3*p->ZP[KP]);

        if(norm>1.0e-20)
        add(a->phi(i,j,k), s/norm);

        // 2D with the plane parallel to the plane of the cells: all water or none
        else if(s>0.0)
        add(a->phi(i,j,k), 1.0);
        }
    }

    // sphere
    if(p->F58_4>0.0)
    {
        LOOP
        {
        if(jd)
        {
        r = sqrt( pow(p->XP[IP]-p->F58_1,2.0)+pow(p->YP[JP]-p->F58_2,2.0)+pow(p->ZP[KP]-p->F58_3,2.0));
        add(a->phi(i,j,k), p->F58_4 - r);
        }

        // 2D: the circle of the cut
        if(!jd)
        {
        const double rc2 = p->F58_4*p->F58_4 - pow(p->YP[JP]-p->F58_2,2.0);

            if(rc2>0.0)
            {
            r = sqrt( pow(p->XP[IP]-p->F58_1,2.0)+pow(p->ZP[KP]-p->F58_3,2.0));
            add(a->phi(i,j,k), sqrt(rc2) - r);
            }
        }
        }
    }

    // vertical cylinder from F59_zs to F59_ze
    if(p->F59_r>0.0)
    {
        LOOP
        {
        double rc = p->F59_r;
        r = sqrt( pow(p->XP[IP]-p->F59_xm,2.0)+pow(p->YP[JP]-p->F59_ym,2.0));

        // 2D: the strip of the cut
        if(!jd)
        {
        const double rc2 = p->F59_r*p->F59_r - pow(p->YP[JP]-p->F59_ym,2.0);
        rc = rc2>0.0 ? sqrt(rc2) : -1.0;
        r = fabs(p->XP[IP]-p->F59_xm);
        }

        if(rc>0.0)
        add(a->phi(i,j,k), sd_two(r-rc, outside1(p->ZP[KP],p->F59_zs,p->F59_ze)));
        }
    }

    if(p->F60>-1.0e20)
    {
        LOOP
        {
        a->phi(i,j,k)=p->F60-p->pos_z();
        smooth[cell_index(p,i,j,k)]=1;
        }

        p->phimean=p->F60;

    }

    // free surface sloping linearly from F63 to the outflow level F62, as a distance
    if((p->F60>-1.0e20 || p->F56>-1.0e20) && p->F62>-1.0e-20&& p->F63>-1.0e-20 )
    {
        phidiff=p->F62-p->phimean;
        xdiff=p->xcoormax-p->F63;
        const double slope = phidiff/xdiff;

        LOOP
        if(p->pos_x() > p->F63)
        {
        a->phi(i,j,k)=(slope*(p->pos_x()-p->F63) + p->phimean - p->pos_z())/sqrt(1.0+slope*slope);
        smooth[cell_index(p,i,j,k)]=1;
        }
    }

	iniphi_box(p,a,pgc);

    iniphi_wedge(p,a,pgc);

    // no water at all: the level set of the air as before
    LOOP
    if(a->phi(i,j,k)<=0.5*no_water)
    a->phi(i,j,k)=-1.0;

	double H=0.0;

    // the density and viscosity of the start as before: smoothed where the level set of before was a
    // distance (F 60, F 62/F 63, F 72), sharp from the sign elsewhere (+1/-1 before).  The hydrostatic
    // pressure of the start (I 12 1, (h - z) ro g) is hydrostatic only with the density it had: with the
    // density of the distance smoothed over the box of F 56 the still water of
    // cfd_3d_still_amr_geodat_mpi2 took 7e-5 m/s instead of 7e-7 m/s
	BASELOOP
	{
        const double f = a->phi(i,j,k);
        const double s = smooth[cell_index(p,i,j,k)] ? f : (f>0.0 ? 1.0 : (f<0.0 ? -1.0 : 0.0));
		H = heaviside(s,interface_width(p,a->phi,i,j,k));

		a->ro(i,j,k)=p->W1*H + p->W3*(1.0-H);
		a->visc(i,j,k)= p->W2*H + p->W4*(1.0-H);
	}

	pgc->start4(p,a->phi,50);
	pgc->start4(p,a->ro,1);
	pgc->start4(p,a->visc,1);


	p->phimean=p->F56;

    if(p->F60>-1.0e20)
    {
    p->phimean=p->F60;
    p->phiout=p->F60;
    p->fsfin=p->F60;
    p->fsfout=p->F60;
    }

    if(p->F61>-1.0e20)
    {
    p->phiin=p->F61;
    p->fsfin=p->F61;
    }

    if(p->F62>-1.0e20)
    {
    p->phiout=p->F62;
    p->fsfout=p->F62;
    }


	pgc->start4(p,a->phi,50);
}

void initialize::iniphi_io(fdm*a, lexer* p, ghostcell* pgc)
{
    if(p->F61>-1.0e20)
    GC4LOOP
    {
        if(p->gcb4[n][4]==1)
        {
        i=p->gcb4[n][0];
        j=p->gcb4[n][1];
        k=p->gcb4[n][2];

        a->phi(i-1,j,k)=p->F61-p->pos_z();
        a->phi(i-2,j,k)=p->F61-p->pos_z();
        a->phi(i-3,j,k)=p->F61-p->pos_z();
        }
    p->phiin=p->F61;
    }

}

void initialize::iniphi_box(lexer* p, fdm *a, ghostcell* pgc)
{
    int jstart, jend;
    int qn;
    const bool jd = p->j_dir==1;

    // boxes of water; in 2D the cells of the box as before (the index range in y)
    if(p->F70>0)
    for(qn=0;qn<p->F70;++qn)
    {
        jstart = p->posc_j(p->F70_ys[qn]);
        jend = p->posc_j(p->F70_ye[qn]);

        LOOP
        if(jd || (j>=jstart && j<jend))
        {
        const double d = sd_box(p->XP[IP],p->YP[JP],p->ZP[KP],p->F70_xs[qn],p->F70_xe[qn],
                                p->F70_ys[qn],p->F70_ye[qn],p->F70_zs[qn],p->F70_ze[qn],jd);
        add(a->phi(i,j,k), d);

        if(d>0.0)
        smooth[cell_index(p,i,j,k)]=0;
        }
    }
	
    // boxes of air
    if(p->F71>0)
	for(qn=0;qn<p->F71;++qn)
    {
        jstart = p->posc_j(p->F71_ys[qn]);
        jend = p->posc_j(p->F71_ye[qn]);

        LOOP
        if(jd || (j>=jstart && j<jend))
        {
        const double d = sd_box(p->XP[IP],p->YP[JP],p->ZP[KP],p->F71_xs[qn],p->F71_xe[qn],
                                p->F71_ys[qn],p->F71_ye[qn],p->F71_zs[qn],p->F71_ze[qn],jd);
        cut(a->phi(i,j,k), d);

        if(d>0.0)
        smooth[cell_index(p,i,j,k)]=0;
        }
    }
	
    // columns with the water level F72_h: inside the column the free surface h - z, outside the
    // level set as it is
    if(p->F72>0)
	for(qn=0;qn<p->F72;++qn)
    {
        jstart = p->posc_j(p->F72_ys[qn]);
        jend = p->posc_j(p->F72_ye[qn]);

        LOOP
        if(jd || (j>=jstart && j<jend))
        {
        const double dc = sd_box(p->XP[IP],p->YP[JP],0.0,p->F72_xs[qn],p->F72_xe[qn],
                                 p->F72_ys[qn],p->F72_ye[qn],-1.0e30,1.0e30,jd);

        a->phi(i,j,k) = std::max(std::min(a->phi(i,j,k),-dc), std::min(dc,p->F72_h[qn]-p->pos_z()));

        if(dc>0.0)
        smooth[cell_index(p,i,j,k)]=1;
        }
	}
}

// wedge of water: the box xs-xe, ys-ye, from the lower to the upper of zs, ze, below (F 112, 113) or
// above (F 114, 115) the plane through zs at xs (ys) and ze at xe (ye), sloping along x (dir 0) or y (dir 1)
static void iniphi_wedge_add(lexer *p, fdm *a, double xs, double xe, double ys, double ye, double zs, double ze, int dir, bool below)
{
    int i,j,k;
    const int marge = increment::marge;
    const bool jd = p->j_dir==1;

    const double slope = (ze-zs)/(dir==0 ? xe-xs : ye-ys);
    const double zlo = std::min(zs,ze);
    const double zhi = std::max(zs,ze);
    const double norm = std::sqrt(1.0 + ((dir==0 || jd) ? slope*slope : 0.0));

    LOOP
    {
    const double x = p->XP[IP];
    const double y = p->YP[JP];
    const double z = p->ZP[KP];

    // 2D: the cells of the wedge as before (the cell centre in y)
    if(!jd && !(y>=ys && y<ye))
    continue;

    const double zp = zs + slope*(dir==0 ? x-xs : y-ys);
    const double dh = (below ? zp-z : z-zp)/norm;

    const double d = std::min(sd_box(x,y,z,xs,xe,ys,ye,zlo,zhi,jd), dh);
    add(a->phi(i,j,k), d);

    if(d>0.0)
    smooth[cell_index(p,i,j,k)]=0;
    }
}

void initialize::iniphi_wedge(lexer* p, fdm *a, ghostcell* pgc)
{
    for(int qn=0; qn<p->F112; ++qn)
    iniphi_wedge_add(p,a,p->F112_xs[qn],p->F112_xe[qn],p->F112_ys[qn],p->F112_ye[qn],p->F112_zs[qn],p->F112_ze[qn],0,true);

    for(int qn=0; qn<p->F113; ++qn)
    iniphi_wedge_add(p,a,p->F113_xs[qn],p->F113_xe[qn],p->F113_ys[qn],p->F113_ye[qn],p->F113_zs[qn],p->F113_ze[qn],1,true);

    for(int qn=0; qn<p->F114; ++qn)
    iniphi_wedge_add(p,a,p->F114_xs[qn],p->F114_xe[qn],p->F114_ys[qn],p->F114_ye[qn],p->F114_zs[qn],p->F114_ze[qn],0,false);

    for(int qn=0; qn<p->F115; ++qn)
    iniphi_wedge_add(p,a,p->F115_xs[qn],p->F115_xe[qn],p->F115_ys[qn],p->F115_ye[qn],p->F115_zs[qn],p->F115_ze[qn],1,false);
}

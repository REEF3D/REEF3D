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

#include"iowave.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

double iowave::xgen_calc(lexer *p)
{
	double x1,y1;
	double x0,y0;
	double dist=0.0;
	
	x0=p->B105_2;
	y0=p->B105_3;
	
	x1 = p->pos_x();
	y1 = p->pos_y();
	
	dist = fabs(y1 - tan_alpha*x1 + tan_alpha*x0 - y0)/sqrt(pow(tan_alpha,2.0)+1.0);
	
	return dist;
}

double iowave::xgen1(lexer *p)
{
	double x1,y1;
	double x0,y0;
	double dist=0.0;
	
	x0=p->B105_2;
	y0=p->B105_3;
	
	x1 = p->pos1_x();
	y1 = p->pos1_y();
	
	dist = fabs(y1 - tan_alpha*x1 + tan_alpha*x0 - y0)/sqrt(pow(tan_alpha,2.0)+1.0);
	
	return dist;
}

double iowave::xgen2(lexer *p)
{
	double x1,y1;
	double x0,y0;
	double dist=0.0;
	
	x0=p->B105_2;
	y0=p->B105_3;
	
	x1 = p->pos2_x();
	y1 = p->pos2_y();
	
	dist = fabs(y1 - tan_alpha*x1 + tan_alpha*x0 - y0)/sqrt(pow(tan_alpha,2.0)+1.0);
	
	return dist;
}

double iowave::ygen_calc(lexer *p)
{
	double x1,y1;
	double x0,y0;
	double dist=0.0;
	
	x0=p->B105_2;
	y0=p->B105_3;
	
	x1 = p->pos_x();
	y1 = p->pos_y();
	
	dist = fabs(x1 - tan_alpha*y1 + tan_alpha*y0 - x0)/sqrt(pow(tan_alpha,2.0)+1.0);
    
	return dist;
}

double iowave::ygen1(lexer *p)
{
	double x1,y1;
	double x0,y0;
	double dist=0.0;
	
	x0=p->B105_2;
	y0=p->B105_3;
	
	x1 = p->pos1_x();
	y1 = p->pos1_y();
	
	dist = fabs(x1 - tan_alpha*y1 + tan_alpha*y0 - x0)/sqrt(pow(tan_alpha,2.0)+1.0);
	
	return dist;
}

double iowave::ygen2(lexer *p)
{
	double x1,y1;
	double x0,y0;
	double dist=0.0;
	
	x0=p->B105_2;
	y0=p->B105_3;
	
	x1 = p->pos2_x();
	y1 = p->pos2_y();
	
	dist = fabs(x1 - tan_alpha*y1 + tan_alpha*y0 - x0)/sqrt(pow(tan_alpha,2.0)+1.0);
	
	return dist;
}

double iowave::distgen_calc(lexer *p)
{
    double x0,y0,denom;
	double dist=1.0e20;
    int test1,test2;    
    
    x0 = p->pos_x();
    y0 = p->pos_y();
    
    for(int qn=0;qn<p->B108;++qn)
    {
    test1=0;
    test2=0;
    
    test1=intriangle(p,G1[qn][0],G1[qn][1],G3[qn][0],G3[qn][1],G2[qn][0],G2[qn][1],x0,y0);
    test2=intriangle(p,G3[qn][0],G3[qn][1],G4[qn][0],G4[qn][1],G2[qn][0],G2[qn][1],x0,y0);

        if(test1==1||test2==1)
        {
        denom = sqrt(pow(Ge[qn][1]-Gs[qn][1],2.0) + pow(Ge[qn][0]-Gs[qn][0],2.0));
        denom = denom>1.0e-20?denom:1.0e20;
        
        dist = MIN(fabs((Ge[qn][1]-Gs[qn][1])*x0 - (Ge[qn][0]-Gs[qn][0])*y0 
                  + Ge[qn][0]*Gs[qn][1] - Ge[qn][1]*Gs[qn][0])/denom,dist);
        
        }
    }
    
	return dist;
}

double iowave::distbeach_calc(lexer *p)
{
    double x0,y0,denom;
	double dist=1.0e20;
    int test1,test2;    
    
    x0 = p->pos_x();
    y0 = p->pos_y();
    
    for(int qn=0;qn<p->B107;++qn)
    {
    test1=0;
    test2=0;
    
    test1=intriangle(p,B1[qn][0],B1[qn][1],B3[qn][0],B3[qn][1],B2[qn][0],B2[qn][1],x0,y0);
    test2=intriangle(p,B3[qn][0],B3[qn][1],B4[qn][0],B4[qn][1],B2[qn][0],B2[qn][1],x0,y0);

        if(test1==1||test2==1)
        {
        denom = sqrt(pow(Be[qn][1]-Bs[qn][1],2.0) + pow(Be[qn][0]-Bs[qn][0],2.0));
        denom = denom>1.0e-20?denom:1.0e20;
        
        dist = MIN(fabs((Be[qn][1]-Bs[qn][1])*x0 - (Be[qn][0]-Bs[qn][0])*y0 
                  + Be[qn][0]*Bs[qn][1] - Be[qn][1]*Bs[qn][0])/denom,dist);
        
        }
    }
    
	return dist;
}

// distgen()/distbeach() test every slice cell against all B107/B108 zone
// polygons. The zones and the grid are fixed, so this is done once here and
// the per-step relaxation loops visit only the cells that lie in a zone.
void iowave::relaxzone4_build(lexer *p)
{
    rz4_i.clear();
    rz4_j.clear();
    rz4_xg.clear();
    rz4_yg.clear();
    rz4_dg.clear();
    rz4_db.clear();
    
    SLICELOOP4
    {
        const double dgv = distgen(p);
        const double dbv = distbeach(p);
        
        if(dgv<1.0e20 || dbv<1.0e20)
        {
        rz4_i.push_back(i);
        rz4_j.push_back(j);
        rz4_xg.push_back(xgen(p));
        rz4_yg.push_back(ygen(p));
        rz4_dg.push_back(dgv);
        rz4_db.push_back(dbv);
        }
    }
    
    rz4_built=true;
}

// distgen/distbeach depend only on (i,j) through XP[IP], YP[JP] and on the
// relaxation-zone polygons, which are fixed after the constructor.
// They are called for every 3D cell in every relax function and RK stage,
// so they are tabulated once per (i,j) (including ghost columns).
void iowave::dist_cache_build(lexer *p)
{
    const int is=i, js=j;

    if(dgcache==nullptr)
    {
    p->Darray(dgcache,p->imax*p->jmax);
    p->Darray(dbcache,p->imax*p->jmax);

    for(i=p->imin; i<p->imin+p->imax; ++i)
    for(j=p->jmin; j<p->jmin+p->jmax; ++j)
    {
    dgcache[IJ] = distgen_calc(p);
    dbcache[IJ] = distbeach_calc(p);
    }
    }

    i=is;
    j=js;
}

// xgen/ygen additionally depend on tan_alpha (set at the end of the constructor)
void iowave::xy_cache_build(lexer *p)
{
    const int is=i, js=j;

    p->Darray(xgcache,p->imax*p->jmax);
    p->Darray(ygcache,p->imax*p->jmax);

    for(i=p->imin; i<p->imin+p->imax; ++i)
    for(j=p->jmin; j<p->jmin+p->jmax; ++j)
    {
    xgcache[IJ] = xgen_calc(p);
    ygcache[IJ] = ygen_calc(p);
    }

    i=is;
    j=js;
}

double iowave::distgen(lexer *p)
{
    if(dgcache==nullptr)
    dist_cache_build(p);

    return dgcache[IJ];
}

double iowave::distbeach(lexer *p)
{
    if(dbcache==nullptr)
    dist_cache_build(p);

    return dbcache[IJ];
}

double iowave::xgen(lexer *p)
{
    if(xgcache==nullptr)
    xy_cache_build(p);

    return xgcache[IJ];
}

double iowave::ygen(lexer *p)
{
    if(ygcache==nullptr)
    xy_cache_build(p);

    return ygcache[IJ];
}
// Generation-zone columns (distgen<1e20) over all (i,j) in ILOOP/JLOOP order,
// registered with the wave library for cached-point evaluation. FNPF loops
// skip flagslice4<=0 columns (SLICELOOP4 order), NHFLOW loops add KLOOP with
// PCHECK (LOOP order), so both reproduce the count sequence of the relaxation
// functions that consume the value arrays.
void iowave::genzone4_build(lexer *p, ghostcell *pgc)
{
    gen_i.clear(); gen_j.clear();
    gen_idx.assign(size_t(p->imax)*size_t(p->jmax),-1);
    
    std::vector<double> xg_, yg_;
    
    ILOOP
    JLOOP
    {
        dg = distgen(p);
        
        if(dg<1.0e20)
        {
        gen_idx[IJ] = int(gen_i.size());
        gen_i.push_back(i);
        gen_j.push_back(j);
        xg_.push_back(xgen(p));
        yg_.push_back(ygen(p));
        }
    }
    
    wave_cache_points(p,pgc,xg_,yg_);
    
    gen_built=true;
}

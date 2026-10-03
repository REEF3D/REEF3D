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
#include<string>

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
    return zones.relax_dist(p->pos_x(),p->pos_y());
}

double iowave::distbeach_calc(lexer *p)
{
    return zones.beach_dist(p->pos_x(),p->pos_y());
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
    gen_i.clear(); gen_j.clear(); gen_src.clear();
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
        
        const bc_zone *z = zones.relax_zone_at(p->pos_x(),p->pos_y());
        gen_src.push_back((z!=nullptr && !z->sources.empty()) ? &z->sources : nullptr);
        xg_.push_back(xgen(p));
        yg_.push_back(ygen(p));
        }
    }
    
    wave_cache_points(p,pgc,xg_,yg_);
    
    gen_built=true;
}

// zone input (B 520-524) that iowave cannot honour yet
void iowave::zones_check(lexer *p)
{
    std::string err;
    
    if(zones.user_relax() && p->B98!=2)
    err = "relaxation zones (B 520 method 1) need relaxation wave generation (B 98 2)";
    
    if(zones.has_sources() && p->A10!=3 && p->A10!=5)
    err = "zone sources (B 524) are available for FNPF and NHFLOW only, so far";
    
    if(zones.has_sources() && p->B89==1)
    err = "zone sources (B 524) do not work with decomposed precalc (B 89 1) yet";
    
    for(const bc_zone &z : zones.relax)
    for(int s : z.sources)
    if(!source_exists(s))
    err = "zone "+std::to_string(z.id)+" uses source "+std::to_string(s)+", which is not defined (B 92 is 1, B 500 the others)";
    
    for(const bc_zone &z : zones.beach)
    if(!z.sources.empty())
    err = "beach zone "+std::to_string(z.id)+" cannot have sources (B 524)";
    
    if(!err.empty())
    {
        if(p->mpirank==0)
        cout<<endl<<"!!! iowave: "<<err<<" !!!"<<endl<<endl;
        exit(1);
    }
}

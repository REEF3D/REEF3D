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
#include"wave_lib.h"
#include<string>

// Wave-frame coordinates of a point: xgen along the wave direction B 105_1,
// ygen along the crest (90 deg counter-clockwise), both measured from the
// origin (B 105_2, B 105_3) and signed. Unsigned distances mirror the wave
// field at the B 105 line, so an origin inside the domain reversed the phase
// of everything upstream of it (generation zone included).
static inline double wave_xframe(double x1, double y1, double x0, double y0, double g)
{
    return (x1-x0)*cos(g) + (y1-y0)*sin(g);
}

static inline double wave_yframe(double x1, double y1, double x0, double y0, double g)
{
    return -(x1-x0)*sin(g) + (y1-y0)*cos(g);
}

double iowave::xgen_calc(lexer *p)
{
	return wave_xframe(p->pos_x(),p->pos_y(),p->B105_2,p->B105_3,gamma);
}

double iowave::xgen1(lexer *p)
{
	return wave_xframe(p->pos1_x(),p->pos1_y(),p->B105_2,p->B105_3,gamma);
}

double iowave::xgen2(lexer *p)
{
	return wave_xframe(p->pos2_x(),p->pos2_y(),p->B105_2,p->B105_3,gamma);
}

double iowave::ygen_calc(lexer *p)
{
	return wave_yframe(p->pos_x(),p->pos_y(),p->B105_2,p->B105_3,gamma);
}

double iowave::ygen1(lexer *p)
{
	return wave_yframe(p->pos1_x(),p->pos1_y(),p->B105_2,p->B105_3,gamma);
}

double iowave::ygen2(lexer *p)
{
	return wave_yframe(p->pos2_x(),p->pos2_y(),p->B105_2,p->B105_3,gamma);
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

// xgen/ygen additionally depend on gamma (set at the end of the constructor)
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
    
    for(const bc_zone &z : zones.edges)
    if(!z.sources.empty() && z.method!=bc_method::riemann)
    err = "Flather / clamped edge "+std::to_string(z.id)+" carries the background only (no B 524); waves come in through a Riemann edge or a relaxation zone";
    
    for(const bc_zone &z : zones.edges)
    for(int s : z.sources)
    if(!source_exists(s))
    err = "zone "+std::to_string(z.id)+" uses source "+std::to_string(s)+", which is not defined (B 92 is 1, B 500 the others)";
    
    if(zones.has_background() && p->A10!=5)
    err = "backgrounds (B 523) and Riemann / Flather edges are available for NHFLOW only, so far";
    
    // decomposed precalc (B 89 1) keeps the spatial parts of the waves from the start, so k
    // cannot follow the background
    if(p->B530>0 && p->B89==1)
    err = "waves on the background (B 530) do not work with decomposed precalc (B 89 1); B 530 0 or B 89 0";
    
    for(const std::vector<bc_zone> *v : {&zones.relax, &zones.beach, &zones.edges})
    for(const bc_zone &z : *v)
    if(z.bg>0 && bgs.index(z.bg)<0)
    err = "zone "+std::to_string(z.id)+" uses background "+std::to_string(z.bg)+", which has no B 510";
    
    // waves on the background (B 530)
    if(p->B530!=0)
    {
        if(p->B530<0 || p->B530>2)
        err = "B 530 mode is 1 (waves on h_eff) or 2 (h_eff + Doppler)";
        
        if(!zones.has_background() || p->A10!=5)
        err = "waves on the background (B 530) need a background (B 510, B 523) in NHFLOW";
        
        for(int n=0; n<wave_nsources(); ++n)
        {
            int id, type;
            double rot;
            const wave_lib *lib = wave_source_lib(n,id,type,rot);
            
            if((lib==nullptr || lib->wave_ncomp()==0) && !(n==0 && type==0))
            err = "waves on the background (B 530) work with linear waves (type 2) and spectral irregular waves (type 31) only, so far; source "+std::to_string(id)+" is type "+std::to_string(type);
        }
    }
    
    if(!err.empty())
    {
        if(p->mpirank==0)
        cout<<endl<<"!!! iowave: "<<err<<" !!!"<<endl<<endl;
        exit(1);
    }
}

// B 530 not given (-1): waves on the background with Doppler whenever a zone's background
// carries a current, on h_eff when it only sets a level; off without a background, without
// waves, with decomposed precalc (B 89 1), or for wave types that do not support it
void iowave::b530_auto(lexer *p)
{
    if(p->B530>=0)
    return;
    
    int mode = 0;
    
    if(zones.has_background() && p->A10==5)
    for(const std::vector<bc_zone> *v : {&zones.relax, &zones.beach, &zones.edges})
    for(const bc_zone &z : *v)
    {
        const int b = z.bg>0 ? bgs.index(z.bg) : -1;
        
        if(bgs.carries_current(b))
        mode = 2;
        
        if(mode==0 && bgs.carries_level(b))
        mode = 1;
    }
    
    bool waves=false, ok=true;
    
    for(int n=0; n<wave_nsources(); ++n)
    {
        int id, type;
        double rot;
        const wave_lib *lib = wave_source_lib(n,id,type,rot);
        
        if(n==0 && type==0)
        continue;
        
        waves = true;
        
        if(lib==nullptr || lib->wave_ncomp()==0)
        ok = false;
    }
    
    if(mode>0 && waves && !ok && p->mpirank==0)
    cout<<"iowave: waves on the background (B 530) stay off: the wave type does not support it (linear 2 and irregular 31 do)"<<endl;
    
    if(mode>0 && waves && ok && p->B89==1 && p->mpirank==0)
    cout<<"iowave: waves on the background (B 530) stay off with decomposed precalc (B 89 1)"<<endl;
    
    p->B530 = (waves && ok && p->B89==0) ? mode : 0;
    
    if(p->B530>0 && p->mpirank==0)
    cout<<"iowave: waves on the background (B 530 "<<p->B530<<" "<<p->B530_N<<", default with a background "<<(p->B530==2?"current":"level")<<")"<<endl;
}

// NHFLOW active beach (B 99 3 / 4): the old ghost-cell velocity eta sqrt(g/h) (or the linear-theory
// profile) with a zero-gradient water level reflected 0.6-0.8 of the waves (validation 10). The
// x+ outflow now becomes an absorbing Riemann edge with still water outside: h_g and the depth
// mean u_g from the outgoing characteristic of the interior and the incoming one of still water.
// The old condition left the ghost water level at still water and imposed the outflow velocity on
// top, so the HLL flux at the boundary face counted the outgoing wave twice. B 99 4 is treated as
// B 99 3: a linear-theory profile or the celerity omega / k in the characteristic gave more
// reflection for kh 2-3 (0.09-0.18 against 0.08-0.12).
// A Riemann / Flather / clamped edge of the user at x+ (B 520, B 521 edge 2) takes precedence.
void iowave::nhflow_active_beach_edge(lexer *p, ghostcell *pgc)
{
    nhf_active_edge = false;
    
    if(p->A10!=5 || (p->B99!=3 && p->B99!=4))
    return;
    
    // only with an outflow boundary at x+
    int nout = 0;
    for(int q=0; q<p->gcslout_count; ++q)
    if(p->gcslout[q][3]==4)
    ++nout;
    
    if(pgc->globalisum(nout)==0)
    return;
    
    if(zones.open_edge(2)!=nullptr)
    {
        if(p->mpirank==0)
        cout<<"iowave: B 99 "<<p->B99<<" ignored at x+: the edge of B 520 zone "<<zones.open_edge(2)->id<<" is used"<<endl;
        return;
    }
    
    bc_zone e(-99,bc_method::riemann,0.0,0.0,0.0,0.0,0.0,1.0);
    e.edge = 2;
    e.active = p->B99;
    zones.edges.push_back(e);
    nhf_active_edge = true;
    
    if(p->mpirank==0)
    cout<<"iowave: active beach B 99 "<<p->B99<<": absorbing Riemann edge at x+ with still water outside"<<endl;
}

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
Architect: Hans Bihs
--------------------------------------------------------------------*/

#include"pic_f.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"reini.h"
#include"fluid_update.h"
#include"levelset_solid_extend.h"
#include"ioflow.h"
#include<algorithm>
#include<cmath>
#include<iostream>

// The free surface from the particles (F 81 1, 2): at every cell centre X the particles within
// the kernel radius R = 2 h (h the cell size), kernel k = (1 - d^2/R^2)^3, give the water fraction
//     c(X) = sum_particles k / sum_fluid k
// where the sum below runs over the sub-cell positions of the seeding (F 83 per cell and
// direction) around X that are not in a solid: c = 1 in water seeded as at the start, also next
// to solids (they are not counted as air); at the walls of the domain the particles within 3 cells
// have mirror images beyond them, which make c there as away from walls. c = 0.5 at the
// surface of a half-space of particles (the edge of their last layer), whatever the position of
// the surface in the cell. Near the surface (0.01 < c < 0.99, within about 1.7 h of it)
//     phi = h d(c)
// with d(c) the signed distance at which a flat surface gives the water fraction c (a table,
// pic_f::fraction_table): a distance for a flat surface on the lattice of the particles.
// Deeper in the water (c >= 0.99) phi keeps its old value, at least the distance of the edge of
// the band, higher in the air (c <= 0.01) likewise, and cells without particles within R the
// smaller of the old value and -R. In the band the signed distance comes then directly from the zero
// of h d(c) (pic_f::redistance), the reinitialisation of the level set (F 40, F 44 iterations) continues
// it beyond the band, and the density and viscosity follow (the fluid update of the level set). The
// particles of the neighbours near the subdomain boundary come in as a halo (pic_f::halo), so phi
// does not depend on the decomposition.
// (The mean position of the particles, Zhu & Bridson 2005, Chen 2017 eq. 3.21, was tried first:
// it falls off nonlinearly across the surface, so its zero between the cell centres and after
// the reinitialisation moved by up to 0.2 h with the position of the surface in the cell; a linear
// map of c, (c - 0.5) h / slope, was off by up to 2 % of a distance next to the surface, which in
// still water started currents.)

namespace
{
const double KR = 2.0;      // kernel radius in cells
const double FULL = 0.99;   // water fraction from which a cell is in the water
const double EMPTY = 0.01;  // water fraction up to which a cell is in the air

double cellsize(lexer *p, int m, int i, int j, int k)
{
    if(p->j_dir==1)
    return std::cbrt(p->DXN[i+m]*p->DYN[j+m]*p->DZN[k+m]);

    return std::sqrt(p->DXN[i+m]*p->DZN[k+m]);
}
}

int pic_f::lid(lexer *p, int i, int j, int k) const
{
    return (i*p->knoy + j)*p->knoz + k;
}

// The water fraction across the surface of a half-space of particles (F 83 per cell and
// direction, below z = 0) at the height z0 (in cells) above the surface, the mean over the
// lateral position in the lattice of the particles: ctab_c falls from 1 to 0 as ctab_d = z0
// goes from -R to R
void pic_f::fraction_table(lexer *p)
{
    const double R2 = KR*KR;
    const double s = 1.0/double(npd);
    const int m = int(KR/s) + 3;
    const bool d3 = p->j_dir==1;
    const int nd = 161;

    ctab_d.assign(nd,0.0);
    ctab_c.assign(nd,0.0);

    for(int q=0; q<nd; ++q)
    {
        const double z0 = -KR + 2.0*KR*double(q)/double(nd-1);
        double sum=0.0;
        int count=0;

        for(int oa=0; oa<4; ++oa)
        for(int ob=0; ob<(d3?4:1); ++ob)
        {
            const double ox = 0.25*double(oa)*s;
            const double oy = 0.25*double(ob)*s;
            double w=0.0, wf=0.0;

            for(int ix=-m; ix<=m; ++ix)
            for(int iy=(d3?-m:0); iy<=(d3?m:0); ++iy)
            for(int iz=-2*m; iz<=2*m; ++iz)
            {
                const double x = (double(ix)+0.5)*s - ox;
                const double y = d3 ? (double(iy)+0.5)*s - oy : 0.0;
                const double z = (double(iz)+0.5)*s;
                const double kw = kernel(x*x+y*y+(z-z0)*(z-z0), R2);

                wf += kw;

                if(z<0.0)
                w += kw;
            }

            sum += w/wf;
            ++count;
        }

        ctab_d[q] = z0;
        ctab_c[q] = sum/double(count);
    }
}

// the signed distance (in cells, water > 0) at which a flat surface gives the water fraction c
double pic_f::fraction_distance(double c) const
{
    const int nd = int(ctab_c.size());

    if(c>=ctab_c[0])
    return -ctab_d[0];

    if(c<=ctab_c[nd-1])
    return -ctab_d[nd-1];

    int lo=0, hi=nd-1;
    while(hi-lo>1)
    {
        const int mid = (lo+hi)/2;
        if(ctab_c[mid]>c)
        lo = mid;
        else
        hi = mid;
    }

    const double t = (ctab_c[lo]-c)/(ctab_c[lo]-ctab_c[hi]);
    return -(ctab_d[lo] + t*(ctab_d[hi]-ctab_d[lo]));
}

// the table of the water fraction and the kernel sums
void pic_f::phi_ini(lexer *p, fdm *a, ghostcell *pgc)
{
    const double starttime = pgc->timer();

    fraction_table(p);
    dfull = fraction_distance(FULL);
    support(p,a);

    if(p->mpirank==0)
    cout<<"PIC free surface: water fraction of the particles, kernel radius "<<KR<<" cells ("
        <<pgc->timer()-starttime<<" s)"<<endl;
}

// per cell the kernel radius and the kernel sum over the fluid: the sub-cell positions of the
// seeding in the cells around that are not in a solid (pic_f::blocked; with floating bodies,
// X 10, again at every step); beyond the walls of the domain all count (the mirror images)
void pic_f::support(lexer *p, fdm *a)
{
    const int fm = fmarg;
    const int ny = p->j_dir==1 ? npd : 1;
    const int nx_g = p->knox + 2*fm;
    const int ny_g = p->j_dir==1 ? p->knoy + 2*fm : 1;
    const int nz_g = p->knoz + 2*fm;
    const int jlo = p->j_dir==1 ? -fm : 0;
    const int jhi = p->j_dir==1 ? p->knoy+fm-1 : 0;

    // the sub-cell positions in the fluid, cells with the margin
    auto gid = [&](int ci, int cj, int ck) { return ((ci+fm)*ny_g + (cj-jlo))*nz_g + (ck+fm); };
    std::vector<char> fluid(size_t(nx_g)*ny_g*nz_g*npc, 0);

    for(int ci=-fm; ci<p->knox+fm; ++ci)
    for(int cj=jlo; cj<=jhi; ++cj)
    for(int ck=-fm; ck<p->knoz+fm; ++ck)
    for(int qi=0; qi<npd; ++qi)
    for(int qj=0; qj<ny; ++qj)
    for(int qk=0; qk<npd; ++qk)
    {
        const double x = p->XN[ci+marge] + (double(qi)+0.5)/double(npd)*p->DXN[ci+marge];
        const double y = p->j_dir==1 ? p->YN[cj+marge] + (double(qj)+0.5)/double(npd)*p->DYN[cj+marge] : p->YP[cj+marge];
        const double z = p->ZN[ck+marge] + (double(qk)+0.5)/double(npd)*p->DZN[ck+marge];

        // outside the domain: the mirror images of the particles (pic_f::phi_build) fill it
        const bool outside = x<xglo_s || x>xglo_e || z<zglo_s || z>zglo_e || (p->j_dir==1 && (y<yglo_s || y>yglo_e));

        if(!outside && blocked(p,a,x,y,z))
        continue;

        fluid[size_t(gid(ci,cj,ck))*npc + (qi*ny + qj)*npd + qk] = 1;
    }

    const int n = p->knox*p->knoy*p->knoz;
    cnorm.assign(n,0.0);
    crad.assign(n,0.0);

    for(int ci=0; ci<p->knox; ++ci)
    for(int cj=0; cj<p->knoy; ++cj)
    for(int ck=0; ck<p->knoz; ++ck)
    {
        const double h = cellsize(p,marge,ci,cj,ck);
        const double R = KR*h;
        const double R2 = R*R;
        const double X = p->XP[ci+marge];
        const double Y = p->YP[cj+marge];
        const double Z = p->ZP[ck+marge];
        double sw=0.0;

        for(int ni=std::max(ci-3,-fm); ni<=std::min(ci+3,p->knox+fm-1); ++ni)
        for(int nj=(p->j_dir==1 ? std::max(cj-3,-fm) : 0); nj<=(p->j_dir==1 ? std::min(cj+3,p->knoy+fm-1) : 0); ++nj)
        for(int nk=std::max(ck-3,-fm); nk<=std::min(ck+3,p->knoz+fm-1); ++nk)
        {
            const size_t base = size_t(gid(ni,nj,nk))*npc;

            for(int qi=0; qi<npd; ++qi)
            for(int qj=0; qj<ny; ++qj)
            for(int qk=0; qk<npd; ++qk)
            {
                if(!fluid[base + (qi*ny + qj)*npd + qk])
                continue;

                const double x = p->XN[ni+marge] + (double(qi)+0.5)/double(npd)*p->DXN[ni+marge];
                const double y = p->j_dir==1 ? p->YN[nj+marge] + (double(qj)+0.5)/double(npd)*p->DYN[nj+marge] : Y;
                const double z = p->ZN[nk+marge] + (double(qk)+0.5)/double(npd)*p->DZN[nk+marge];

                sw += kernel((x-X)*(x-X) + (y-Y)*(y-Y) + (z-Z)*(z-Z), R2);
            }
        }

        const int l = lid(p,ci,cj,ck);
        crad[l] = R;
        cnorm[l] = sw;
    }
}

void pic_f::phi_build(lexer *p, fdm *a, ghostcell *pgc)
{
    const double starttime = pgc->timer();

    // floating bodies (X 10): the fluid around them changes
    if(p->X10>0 && p->count>0)
    support(p,a);

    const int n = p->knox*p->knoy*p->knoz;
    std::vector<double> sw(n,0.0);

    // the particles of this rank and the halo: their kernel weights at the cell centres around
    auto deposit = [&](double x, double y, double z)
    {
        const int ci = cell(p->XN,p->knox,x);
        const int cj = p->j_dir==1 ? cell(p->YN,p->knoy,y) : 0;
        const int ck = cell(p->ZN,p->knoz,z);

        for(int ti=std::max(ci-3,0); ti<=std::min(ci+3,p->knox-1); ++ti)
        for(int tj=(p->j_dir==1 ? std::max(cj-3,0) : 0); tj<=(p->j_dir==1 ? std::min(cj+3,p->knoy-1) : 0); ++tj)
        for(int tk=std::max(ck-3,0); tk<=std::min(ck+3,p->knoz-1); ++tk)
        {
            const int l = lid(p,ti,tj,tk);
            const double dx = x - p->XP[ti+marge];
            const double dy = p->j_dir==1 ? y - p->YP[tj+marge] : 0.0;
            const double dz = z - p->ZP[tk+marge];
            const double d2 = dx*dx + dy*dy + dz*dz;
            const double R2 = crad[l]*crad[l];

            if(d2<R2)
            sw[l] += kernel(d2,R2);
        }
    };

    // with the mirror images across the walls of the domain within 3 cells of them: the water
    // fraction of a flat surface is then the same next to a wall as away from it
    double band[3] = {0.0,0.0,0.0};
    for(int ci=0; ci<p->knox; ++ci)
    band[0] = std::max(band[0], 3.0*p->DXN[ci+marge]);
    for(int cj=0; cj<p->knoy; ++cj)
    band[1] = std::max(band[1], 3.0*p->DYN[cj+marge]);
    for(int ck=0; ck<p->knoz; ++ck)
    band[2] = std::max(band[2], 3.0*p->DZN[ck+marge]);

    auto images = [&](double x, double y, double z)
    {
        double xs[2] = {x,x}, ys[2] = {y,y}, zs[2] = {z,z};
        int nx=1, ny=1, nz=1;

        if(x-xglo_s<band[0])
        {
        // B 60: no images across the inflow, its particles are in the ghost cells (pic_f_inflow.cpp)
        if(!inflowsrc)
        xs[nx++] = 2.0*xglo_s - x;
        }
        else if(xglo_e-x<band[0])
        xs[nx++] = 2.0*xglo_e - x;

        if(p->j_dir==1)
        {
        if(y-yglo_s<band[1])
        ys[ny++] = 2.0*yglo_s - y;
        else if(yglo_e-y<band[1])
        ys[ny++] = 2.0*yglo_e - y;
        }

        if(z-zglo_s<band[2])
        zs[nz++] = 2.0*zglo_s - z;
        else if(zglo_e-z<band[2])
        zs[nz++] = 2.0*zglo_e - z;

        for(int a1=0; a1<nx; ++a1)
        for(int a2=0; a2<ny; ++a2)
        for(int a3=0; a3<nz; ++a3)
        deposit(xs[a1],ys[a2],zs[a3]);
    };

    for(const auto &q : P)
    images(q.x,q.y,q.z);

    std::vector<double> H;
    halo(p,pgc,H);

    for(size_t m=0; m+2<H.size(); m+=3)
    images(H[m],H[m+1],H[m+2]);

    // the new level set: near the surface from the water fraction, the old one extended elsewhere
    std::vector<char> inband(n,0);

    for(int ci=0; ci<p->knox; ++ci)
    for(int cj=0; cj<p->knoy; ++cj)
    for(int ck=0; ck<p->knoz; ++ck)
    {
        if(p->flag4[(ci-p->imin)*p->jmax*p->kmax + (cj-p->jmin)*p->kmax + ck-p->kmin]<0)
        continue;

        const int l = lid(p,ci,cj,ck);
        const double R = crad[l];
        const double h = R/KR;
        double &phi = a->phi(ci,cj,ck);

        if(sw[l]<=0.0 || cnorm[l]<=0.0)
        {
        phi = std::min(phi,-R);
        continue;
        }

        const double c = sw[l]/cnorm[l];
        const double val = h*fraction_distance(c);

        // beyond the band at least as deep in the water / high in the air as its edge (a fixed
        // value: the particles' value there would ratchet the old one up with their noise)
        if(c>=FULL)
        phi = std::max(phi,h*dfull);
        else if(c<=EMPTY)
        phi = std::min(phi,-h*dfull);
        else
        {
        phi = val;
        inband[l] = 1;
        }
    }

    pgc->start4(p,a->phi,gcval_phi);

    // the signed distance in the band from the zero of h d(c) (pic_f_redist.cpp)
    redistance(p,a,inband);
    pgc->start4(p,a->phi,gcval_phi);

    // the relaxation zones of the waves (B 98 2, B 99 1, 2): phi relaxed towards that of the waves,
    // as for the level set; the columns of the particles there follow it (pic_f::zone_stretch)
    if(relaxzones)
    {
    // the surface of the particles near the relaxed one of the step before, the relaxed one near it
    if(etar.empty())
    etar.assign(size_t(p->gknox)*p->gknoy, p->phimean);

    zone_surface(p,a,pgc,etap,etar);
    pflow->phi_relax(p,pgc,a->phi);
    pgc->start4(p,a->phi,gcval_phi);
    zone_surface(p,a,pgc,etar,etap);
    }

    // inside the solids of the direct forcing (and the 6DOF bodies with X 48 3) the level set is
    // continued from the fluid around, as for the level set (levelset_solid_extend.cpp)
    levelset_solid_extend(p,a,pgc,a->phi,gcval_phi);

    // the signed distance (at the start, count 0: the iterations of the level set's initialisation)
    p->reini_iter = p->F44;
    preini->start(a,p,a->phi,pgc,pflow);
    pgc->start4(p,a->phi,gcval_phi);

    // B 60: phi in the ghost cells of the inflow, as for the level set
    if(inflowsrc)
    pflow->fsfinflow(p,a,pgc);

    // density and viscosity
    pupdate->start(p,a,pgc,a->u,a->v,a->w);

    // F 81 1, single phase: those of the water everywhere (the air is not in the momentum; the
    // pressure is 0 at the surface, pic_gfm.cpp)
    if(single)
    {
        std::fill(a->ro.begin(),a->ro.end(), p->W1);
        std::fill(a->visc.begin(),a->visc.end(), p->W2);
    }

    phitime += pgc->timer() - starttime;
}

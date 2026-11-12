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
#include<algorithm>
#include<cmath>
#include<mpi.h>

// The particles in the relaxation zones of the waves (B 98 2, B 99 1, 2; F 81 1, 2). There phi
// of the particles is relaxed towards that of the waves (ioflow::phi_relax, pic_f::phi_build),
// and the particles follow: each column of cells in a zone is stretched or compressed between its
// bottom and the surface, from the surface of the particles (the lowest zero of phi above the water)
// to the relaxed one,
//     z -> zb + (z - zb) (eta_relaxed - zb) / (eta_particles - zb),
// so the surface of the particles is the relaxed one at the next step and the relaxation adds up
// over the steps as for the level set. No particles are added or removed for it (the density of
// the particles in a column changes a little; pic_f::reseed keeps the water filled): adding them
// in the gaps of the zones and removing the surplus (the first two versions) pumped water into the
// tank at the edge of the zone of the wave generation, where the flow leaves it every half period.
// The surface of a column is its lowest zero from the water to the air, and the stretch moves it by
// at most half a cell per step. The columns and their bottoms are those of the whole domain (one value
// per column, reduced over the ranks), so a split in z gives the same.

void pic_f::zone_surface(lexer *p, fdm *a, ghostcell *pgc, std::vector<double> &eta, const std::vector<double> &ref)
{
    const size_t nc = size_t(p->gknox)*p->gknoy;
    std::vector<double> dist(nc,1.0e30), loc(nc,1.0e30);

    // the zero of phi from the water to the air nearest to the reference height of the column (the
    // relaxed surface of the step before, at the start the still water level): the surface of the
    // water, not a drop above it or a gap in the particles below it. The highest zero (the first
    // version) took the spray of a steep wave for the surface and stretched the column up to it,
    // the lowest (the second) a gap in the particles and stretched it down.
    for(int i=0; i<p->knox; ++i)
    for(int j=0; j<p->knoy; ++j)
    {
        const size_t g = size_t(i+p->origin_i)*p->gknoy + (j+p->origin_j);

        for(int k=0; k<p->knoz; ++k)
        {
            const double f0 = a->phi(i,j,k);
            const double f1 = a->phi(i,j,k+1);

            if(f0>=0.0 && f1<0.0)
            {
            const double za = p->ZP[k+marge], zc = p->ZP[k+1+marge];
            const double z0 = za + f0/(f0-f1)*(zc-za);

                if(std::fabs(z0-ref[g])<dist[g])
                {
                dist[g] = std::fabs(z0-ref[g]);
                loc[g] = z0;
                }
            }
        }
    }

    // over the ranks: the crossing of the smallest distance
    std::vector<double> dmin(nc,1.0e30);
    MPI_Allreduce(dist.data(),dmin.data(),int(nc),MPI_DOUBLE,MPI_MIN,pgc->mpi_comm);

    for(size_t g=0; g<nc; ++g)
    if(dist[g]>dmin[g])
    loc[g] = 1.0e30;

    eta.assign(nc,1.0e30);
    MPI_Allreduce(loc.data(),eta.data(),int(nc),MPI_DOUBLE,MPI_MIN,pgc->mpi_comm);

    // the bottoms of the columns: the lowest face of a fluid cell above the topography
    if(zbot.empty())
    {
        std::vector<double> zl(nc,1.0e30);

        for(int i=0; i<p->knox; ++i)
        for(int j=0; j<p->knoy; ++j)
        for(int k=0; k<p->knoz; ++k)
        if(p->flag4[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin]>0 && a->topo(i,j,k)>=0.0)
        {
        const size_t g = size_t(i+p->origin_i)*p->gknoy + (j+p->origin_j);
        zl[g] = std::min(zl[g], p->ZN[k+marge]);
        break;
        }

        zbot.assign(nc,1.0e30);
        MPI_Allreduce(zl.data(),zbot.data(),int(nc),MPI_DOUBLE,MPI_MIN,pgc->mpi_comm);
    }
}

void pic_f::zone_stretch(lexer *p, fdm *a)
{
    if(etap.empty() || etar.empty())
    return;

    for(auto &q : P)
    {
        if(!inzone(p,q.x,q.y))
        continue;

        const int ci = std::min(std::max(cell(p->XN,p->knox,q.x),0),p->knox-1);
        const int cj = p->j_dir==1 ? std::min(std::max(cell(p->YN,p->knoy,q.y),0),p->knoy-1) : 0;
        const size_t g = size_t(ci+p->origin_i)*p->gknoy + (cj+p->origin_j);

        const double ep = etap[g], zb = zbot[g];

        if(ep>1.0e29 || etar[g]>1.0e29 || zb>1.0e29 || ep-zb<=0.0 || etar[g]-zb<=0.0)
        continue;

        // at most half a cell per step
        const double er = std::max(ep-0.5*p->DXM, std::min(ep+0.5*p->DXM, etar[g]));

        const double zn = zb + (q.z-zb)*(er-zb)/(ep-zb);

        if(!blocked(p,a,q.x,q.y,zn) && zn<zglo_e)
        q.z = zn;
    }
}

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

#include"levelset_solid_extend.h"
#include"solid_closed.h"
#include"solid_aperture.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"field.h"
#include"field4.h"
#include"increment.h"
#include<algorithm>
#include<cmath>
#include<vector>
#include<utility>

//  The level set inside the solids of the direct forcing.
//
//  With the direct forcing (G 505 0) the cells of a solid are fluid cells: their velocity is forced
//  to zero, and the level set there is not moved.  It keeps the free surface the start (F 60, or
//  the hot start I 30 from FNPF) put through the solid.  A pile through the free surface then holds
//  a water column of its own.  As soon as the free surface outside moves past it (waves, or the
//  surface of the hot start relaxing), the level set has a step at the wall.  The pressure of that
//  column against the water outside drives a jet along the wall through the forcing band, and the
//  air just above the free surface at the wall takes velocities of the order of 1 m/s (a pier with
//  a hot start from FNPF: 0.7 m/s down along the piles where the waves had 0.07 m/s, validation 55).
//
//  Here the level set of the cells inside a solid is extrapolated from the fluid, the way the free
//  surface at rest lies: constant across (x, y), linear along the gravity (z).  The cells of the
//  solid are put in layers by their distance in cells from the fluid (layer 1: a face next to a
//  fluid cell, layer 2: next to layer 1, ...); every cell takes the mean over the directions in
//  which its neighbour lies in a lower layer (towards the fluid): in x and y the value of that
//  neighbour, in z 2 f(+1) - f(+2), the change limited to the height of the cell (|grad f| = 1 of a
//  distance).  The layers come from the sign of the solid level set alone: inside some solids (S 61
//  a wedge) it is not a distance.
//  A free surface at rest (f = h - z) is continued exactly through every wall: at a pile the free
//  surface meets the wall at a right angle and the pile holds the level of the water around it, at a
//  gently sloping bed (a beach) the still water level stays a plane through the bed.  Constant along
//  the normal of the solid would bend the free surface along the shore line of a beach and drive a
//  flow there (0.75 m/s in still water, the case of cfd_2d_beach_still); linear across as well amplifies the slope of
//  the free surface at the wall of a pile a few cells wide into a column of its own (validation 55).
//  The layers are done in order, so that the neighbours towards the fluid are done first: one sweep,
//  a second one after the ghost cells for the solids across a partition edge.  The layers go through
//  the whole solid (up to 250 layers): its deep cells are not left to the reinitialisation.
//
//  The fluid below a solid does not take part: a cell next to the fluid only below it is not in
//  layer 1, and the fluid below a cell is not one of its directions.  The fluid under a solid roof (a
//  culvert in a lock wall) is enclosed water: its level set is the distance the reinitialisation
//  finds, to the zero crossing the extension put into the wall above it, not the free surface.
//  Extrapolated upwards into the wall it lowered that crossing, the reinitialisation then the level
//  set of the culvert, and so on: in a 3D culvert under a 0.36 m wall the level set at the culvert
//  roof fell from 0.08 to 0.0 in 4 s; in the full lock the wrong values reached the top of the walls
//  and the air above them (phi -0.04 instead of -0.28 at the wall top after 25 s), and from there
//  the air over the whole lock.  So the solid takes the level set of the free surface around it,
//  sideways and from above.  A solid that touches the fluid only from below (a roof over the whole
//  width at the top of the domain) takes it from below, in layers after the others.  (Until 2026-11
//  the fluid below took part and only 8 layers of level 0 were extended; X 48 4 was this as an
//  option, now the same as 3.)
//
//  Along the gravity the level set of a free surface at rest grows with the depth.  The step of the
//  extension along z is limited to that side (from 0 to the height of the cell downwards): the water
//  next to a body between two water levels (the foot of a lock gate on its sill) has the distance to
//  the air behind the body, which grows with z; continued linearly it fell with the depth and put a
//  zero crossing into the sill below the water.  The reinitialisation then took the water over the
//  sill for the distance to that crossing, the extension followed, and the crossing rose to the top
//  of the sill and spread along the bed (3 cm/s in the lock of validation 27): a layer of part-air
//  density under the water, which rose in plumes of 1 m/s.
//
//  X 48 2 (the default); X 48 0 keeps the level set inside the solids as it is, X 48 1 is the copy
//  into the ghost cells of the flagged solids (G 505 1, ghostcell::solid_forcing_lsm).
//
//  X 48 3: as 2, and inside the 6DOF bodies of the direct forcing as well (X 10, the level set a->fb
//  of all bodies at the pose of the stage).  A body between two water levels (a lock gate, its top
//  above both) holds the free surface of the start inside, a step at the face of the high side: the
//  same jet along the face, and the water level of the high side falls (the water goes into the air
//  above, the level set loses it), 0.12 m in 7 s at 0.32 m head (validation 26_cfd_ship_lock).  The
//  layers of the body follow the level set of their own side; with fewer than about four cells across
//  the body the two sides meet in the middle.
//
//  With the gravity along -z (W 22 < 0) the solids take the level of the nearest fluid that sees the
//  free surface (in levelset_solid_extend below, 2026-11): fluid that reaches the free surface only
//  around a solid (the air over a dry top, the water under a roof) no longer feeds the solid next to
//  it.  The layers above remain for the other directions of the gravity.

namespace
{
// the neighbour after the next one (l2) towards the fluid as well, from the next (l1)
inline bool towards(unsigned char l1, unsigned char l2)
{
    return l2<l1 || l2==0;
}

// linear from the two neighbours towards the fluid (f1 the next at the distance h1, f2 the one after,
// h2 from f1), the change limited to the distance h1 and to the side of the free surface at rest
// (lo, hi: 0 to h1 downwards, -h1 to 0 upwards): below water the level set of the solid does not
// fall with the depth
inline double extrapolate(double f1, double f2, double h1, double h2, double lo, double hi)
{
    return f1 + std::max(lo*h1, std::min(hi*h1, (f1-f2)*h1/h2));
}
}

void levelset_solid_extend(lexer *p, fdm *a, ghostcell *pgc, field &f, int gcval)
{
    // with a free surface (F 30): solids of the direct forcing (G 505 0; X 48 2, 3), the 6DOF bodies
    // (X 10; X 48 3) with the level set a->fb of all bodies
    const bool solids = (p->X48==2 || p->X48==3 || p->X48==4) && p->solidread>0 && p->G505==0;
    const bool bodies = (p->X48==3 || p->X48==4) && p->X10>0;
    
    if(p->F30==0 || (!solids && !bodies))
    return;

    int i,j,k;
    const int marge = increment::marge;
    const int jd = p->j_dir==1 ? 1 : 0;
    const unsigned char outside = 255;

    // D 1 2: the open fractions of the faces (before the first projection as well)
    solid_aperture_update(p,a,pgc);

    // the layers, on the cells of the partition and its ghost cells across a partition edge or a
    // periodic edge (the fluid there within the ghost cells); fluid 0.  The ghost cells beyond the
    // edges of the domain are left out (outside): the level set there is the boundary condition.
    bool xm = p->nb1>=0 || p->periodic1==1, xp = p->nb4>=0 || p->periodic1==1;
    bool ym = p->nb3>=0 || p->periodic2==1, yp = p->nb2>=0 || p->periodic2==1;
    bool zm = p->nb5>=0 || p->periodic3==1, zp = p->nb6>=0 || p->periodic3==1;

    // a patch of the mesh refinement (cfd_amr): the cells around it inside the domain hold the level
    // set of the patch next to it or of the coarser level (cfd_amr::fill_phi_out) and the level sets of
    // the solids and the bodies.  They take part as the ghost cells across a partition edge: the fluid
    // there starts the layers, the solid cells there pass their values on.  Without them a body cut
    // by the edge between two patches took the level set of the fluid of its own patch alone (a lock
    // gate: on the patch of the high side the water up to the edge, on the other the air), a step in
    // the level set at the edge, on the face of the gate, and a film of 1.3 m/s down that face
    if(p->amr_patch==1)
    {
        xm = p->origin_i>0;
        xp = p->origin_i+p->knox<p->gknox;
        ym = p->origin_j>0;
        yp = p->origin_j+p->knoy<p->gknoy;
        zm = p->origin_k>0;
        zp = p->origin_k+p->knoz<p->gknoz;
    }
    const int n = p->imax*p->jmax*p->kmax;
    std::vector<unsigned char> layer(n,outside);
    std::vector<int> solid;
    for(i=p->imin; i<p->imin+p->imax; ++i)
    for(j=p->jmin; j<p->jmin+p->jmax; ++j)
    for(k=p->kmin; k<p->kmin+p->kmax; ++k)
    {
        if((i<0 && !xm) || (i>=p->knox && !xp) || (k<0 && !zm) || (k>=p->knoz && !zp))
        continue;
        if(jd && ((j<0 && !ym) || (j>=p->knoy && !yp)))
        continue;

        // (D 1 2: the cells of a solid with an open face take part in the flow, their level set is
        // transported; the solid is its cells with all faces closed, solid_aperture.h)
        if(!((solids && a->solid(i,j,k)<0.0 && (p->D1!=2 || solid_cell_closed(p,a,i,j,k))) || (bodies && a->fb(i,j,k)<0.0)))
        layer[IJK]=0;
        else
        solid.push_back(IJK);
    }

    const int si = p->jmax*p->kmax, sj = p->kmax, sk = 1;

    // the side of the step along z (in units of the cell height): down along the gravity the level set
    // does not fall, up it does not rise; without gravity along z both
    const double dn_lo = p->W22<0.0 ? 0.0 : -1.0, dn_hi = p->W22>0.0 ? 0.0 : 1.0;
    const double up_lo = p->W22>0.0 ? 0.0 : -1.0, up_hi = p->W22<0.0 ? 0.0 : 1.0;

    // the depth of the layers: through the whole solid, as many layers as the marks hold (254, the
    // mark 255 is outside).  On a patch of the mesh refinement (cfd_amr) the solid cells at the edge of
    // the patch (filled from the coarser level) and those inside it then have the level set of the
    // same extension (8 cells of the patch against 8 of level 0 left a step in the density in a thick
    // wall at a patch edge, and the pressure correction blew the velocities up there)
    const int nband = 250;

    // layer 1: the solid cells next to the fluid; layer l: the solid cells next to layer l-1 (from
    // the cells of layer l-1, the same cells as a sweep over all solid cells per layer)
    std::vector<char> insolid(n,0);
    for(int q : solid)
    insolid[q]=1;

    // Gravity along -z: the level of the nearest fluid that sees the free surface.
    //
    // A fluid cell sees the free surface when its own column leads to it: water with air over it
    // before any solid, air with water under it before any solid.  The other fluid reaches the free
    // surface only around a solid: the water in a culvert and under a hull, the air over the crest
    // of a wall, over a deck, over a dry sill or bed.  Its level set is what the reinitialisation
    // makes of the level set of the solid next to it (the larger value for air, the smaller one for
    // water: the way to the free surface).  Taken back into the solid, as until 2026-11, the two
    // follow each other and nothing ties them to the free surface: over the walls of a lock (1.5
    // cells of air under the top of the domain) the level of the upper head crept along the wall
    // tops into the chamber part (0.18 m off after 4 s of still water), every transient added to
    // it, and during the filling the level set over the walls reached zero: water density in the
    // air there, 3.6 m/s (validation 49).  The water of a culvert and the wall over it did the same.
    //
    // So the solid cells and the fluid cells that do not see the free surface take the level
    // eta = f + z of the nearest fluid that does, in layers by their distance in cells from it: the
    // mean over the neighbours of a lower layer, f the same across and f + dz one cell down (the
    // plane of the free surface at rest; |grad f| = 1).  Only the solid cells are written, the
    // fluid keeps its level set and carries the level: the crest of a wall takes the level of the
    // chamber through the air over it, a dry sill in front of a gate the level in front of it and
    // not the one behind the gate, the wall of a culvert the level of the chamber through the ports.
    // What no such fluid reaches (a pocket closed all around) takes the fluid next to it.
    if(p->W22<0.0)
    {
        const int si = p->jmax*p->kmax, sj = p->kmax, sk = 1;
        const int unreached = 1<<30;
        field4 g(p);

        // sees 1: the fluid cells that see the free surface, air 1: the fluid cells with f < 0
        // (the end of the partition or of a patch in a column counts as open)
        std::vector<char> sees(n,0), air(n,0);
        for(i=p->imin; i<p->imin+p->imax; ++i)
        for(j=p->jmin; j<p->jmin+p->jmax; ++j)
        {
            char open = 1;
            for(k=p->kmin+p->kmax-1; k>=p->kmin; --k)
            {
                const int q = IJK;
                if(layer[q]!=0)
                {
                open = 0;
                continue;
                }

                if(f(i,j,k)<0.0)
                {
                air[q] = 1;
                open = 1;
                }
                else
                sees[q] = open;
            }

            open = 1;
            for(k=p->kmin; k<p->kmin+p->kmax; ++k)
            {
                const int q = IJK;
                if(layer[q]!=0)
                {
                open = 0;
                continue;
                }

                if(air[q])
                sees[q] = open;
                else
                open = 1;
            }
        }

        // lay 0: the fluid that sees the free surface; the cells to do (med): the other fluid, then
        // the solid.  Around a patch of the mesh refinement the fluid that does not see the free
        // surface is a source with its own level set (the level it carries on the coarser level is
        // not known on the patch; its level set follows the solids there, which are tied to the free
        // surface).  Across a partition edge the level comes with the exchange.
        std::vector<int> lay(n,unreached);
        std::vector<char> med(n,0);
        for(i=p->imin; i<p->imin+p->imax; ++i)
        for(j=p->jmin; j<p->jmin+p->jmax; ++j)
        for(k=p->kmin; k<p->kmin+p->kmax; ++k)
        {
            const int q = IJK;
            g(i,j,k) = f(i,j,k);

            if(layer[q]==0)
            {
                if(sees[q])
                lay[q] = 0;
                else if(p->amr_patch==0 || (i>=0 && i<p->knox && j>=0 && j<p->knoy && k>=0 && k<p->knoz))
                med[q] = 1;
                else
                lay[q] = 0;
            }
        }

        // the layers of the cells med from the cells src next to them
        std::vector<std::vector<int> > band(1);
        std::vector<char> src(n,0);
        auto layers = [&](const bool phase)
        {
            const int l0 = int(band.size());
            band.push_back(std::vector<int>());
            for(int q=0; q<n; ++q)
            if(med[q] && lay[q]==unreached)
            {
                const int ii = q/si, jj = (q/sj)%p->jmax, kk = q%p->kmax;
                auto from = [&](int r) { return src[r] && (!phase || air[r]==air[q]); };
                if((ii>0 && from(q-si)) || (ii<p->imax-1 && from(q+si))
                || (jd && jj>0 && from(q-sj)) || (jd && jj<p->jmax-1 && from(q+sj))
                || (kk>0 && from(q-sk)) || (kk<p->kmax-1 && from(q+sk)))
                band[l0].push_back(q);
            }
            for(int q : band[l0])
            lay[q] = l0;

            for(int l=l0+1; !band[l-1].empty(); ++l)
            {
                band.push_back(std::vector<int>());
                for(int r : band[l-1])
                {
                    const int ii = r/si, jj = (r/sj)%p->jmax, kk = r%p->kmax;
                    auto add = [&](int q)
                    {
                        if(med[q] && lay[q]==unreached && (!phase || air[q]==air[r]))
                        {
                        lay[q]=l;
                        band[l].push_back(q);
                        }
                    };
                    if(ii<p->imax-1) add(r+si);
                    if(ii>0) add(r-si);
                    if(jd && jj<p->jmax-1) add(r+sj);
                    if(jd && jj>0) add(r-sj);
                    if(kk<p->kmax-1) add(r+sk);
                    if(kk>0) add(r-sk);
                }
            }
            band.pop_back();
        };

        // 1. the fluid that does not see the free surface, through the fluid of its own phase from
        // the one that does: the air over a dry sill in front of a gate takes the level in front of
        // the gate through the air, however near the water behind the gate is (through the leaf, or
        // in the gap at its heel)
        for(int q=0; q<n; ++q)
        src[q] = lay[q]==0;
        layers(true);

        // 2. the solid, from the fluid next to it: the level the fluid carries (its own level set in
        // a pocket no fluid that sees the free surface reaches)
        for(int q=0; q<n; ++q)
        {
            if(layer[q]==0 && lay[q]==unreached)
            lay[q] = 0;

            src[q] = layer[q]==0;
            med[q] = insolid[q];
        }
        layers(false);

        // (every rank takes part in the ghost cell exchange)
        for(int pass=0; pass<2; ++pass)
        {
            for(size_t l=1; l<band.size(); ++l)
            for(int q : band[l])
            {
                i = q/si + p->imin;
                j = (q/sj)%p->jmax + p->jmin;
                k = q%p->kmax + p->kmin;

                if(i<0 || i>=p->knox || j<0 || j>=p->knoy || k<0 || k>=p->knoz)
                continue;

                const int ll = int(l);
                double sum = 0.0;
                int count = 0;

                const bool fl = layer[q]==0;
                auto from = [&](int r) { return lay[r]<ll && (!fl || air[r]==air[q]); };

                if(from(q+si)) { sum += g(i+1,j,k); ++count; }
                if(from(q-si)) { sum += g(i-1,j,k); ++count; }
                if(jd && from(q+sj)) { sum += g(i,j+1,k); ++count; }
                if(jd && from(q-sj)) { sum += g(i,j-1,k); ++count; }
                if(from(q+sk)) { sum += g(i,j,k+1) + p->DZP[KP]; ++count; }
                if(from(q-sk)) { sum += g(i,j,k-1) - p->DZP[KM1]; ++count; }

                if(count>0)
                g(i,j,k) = sum/double(count);
            }

            pgc->start4(p,g,gcval);
        }

        for(int q : solid)
        if(lay[q]!=unreached)
        {
            i = q/si + p->imin;
            j = (q/sj)%p->jmax + p->jmin;
            k = q%p->kmax + p->kmin;

            if(i<0 || i>=p->knox || j<0 || j>=p->knoy || k<0 || k>=p->knoz)
            continue;

            f(i,j,k) = g(i,j,k);
        }

        // The solid next to the fluid.  The reinitialisation takes a fluid cell from the solid cell
        // next to it where that is nearer to the free surface.  Where the solid holds another level
        // than the fluid at its face, it puts a free surface there, and the fluid takes the distance
        // to the face: a film of the other density on it (a dry bed ahead of a dam break, where the
        // bed holds the level of the front; a pipe under pressure above the level outside).  So
        // within three cells of the fluid along the gravity and two across, the solid does not lie
        // nearer to the free surface than the fluid continued into it (across: its value, along the
        // gravity: the plane through it).  With the solid on that plane the fluid follows the level
        // of the solid, which is tied to the free surface: the air over a dry sill comes back to the
        // level in front of it after a disturbance.  Where the solid is more than one cell beyond the
        // plane, the solid cell takes the value of the fluid cell (blended over a second cell): the
        // fluid then has no way to the free surface through the solid, and the air over a dry bed
        // ahead of a front keeps its distance to the front.  The cells with fluid as near on the other
        // side (a thin wall between two levels) keep their value.
        {
            std::vector<std::pair<int,double> > lim;
            const int st[6] = {sk,-sk,si,-si,sj,-sj};

            for(int q : solid)
            {
                i = q/si + p->imin;
                j = (q/sj)%p->jmax + p->jmin;
                k = q%p->kmax + p->kmin;

                if(i<0 || i>=p->knox || j<0 || j>=p->knoy || k<0 || k>=p->knoz)
                continue;

                const double s0 = f(i,j,k);
                double best = 0.0, val = s0;

                for(int d=0; d<6; ++d)
                {
                    if(!jd && d>=4)
                    continue;

                    const int mmax = d<2 ? 3 : 2;

                    // the fluid cell m cells away, solid in between
                    int m = 0;
                    for(int r=1; r<=mmax; ++r)
                    {
                        if(layer[q+r*st[d]]==0)
                        {
                        m = r;
                        break;
                        }
                        if(!insolid[q+r*st[d]])
                        break;
                    }
                    if(m==0)
                    continue;

                    const int qf = q+m*st[d];
                    const int fi = qf/si + p->imin, fj = (qf/sj)%p->jmax + p->jmin, fk = qf%p->kmax + p->kmin;

                    // fluid as near on the other side (a thin wall): the cell belongs to that side
                    bool other = false;
                    for(int r=1; r<=m; ++r)
                    {
                        const int qo = q-r*st[d];
                        if(layer[qo]==0)
                        {
                        other = true;
                        break;
                        }
                        if(!insolid[qo])
                        break;
                    }
                    if(other)
                    continue;

                    const double ff = f(fi,fj,fk);
                    const double c = ff + (d<2 ? p->ZP[fk+marge]-p->ZP[KP] : 0.0);
                    const double h = d<2 ? p->DZN[KP] : (d<4 ? p->DXN[IP] : p->DYN[JP]);
                    const double ex = (air[qf] ? s0-c : c-s0)/h;

                    if(ex>best)
                    {
                    best = ex;
                    val = c + std::max(0.0,std::min(1.0,ex-1.0))*(ff-c);
                    }
                }

                if(best>0.0)
                lim.push_back(std::make_pair(q,val));
            }

            for(const auto &c : lim)
            {
                i = c.first/si + p->imin;
                j = (c.first/sj)%p->jmax + p->jmin;
                k = c.first%p->kmax + p->kmin;
                f(i,j,k) = c.second;
            }
        }

        pgc->start4(p,f,gcval);
        return;
    }

    std::vector<std::vector<int> > band(nband+1);
    for(int q : solid)
    {
        const int ii = q/si, jj = (q/sj)%p->jmax, kk = q%p->kmax;
        bool next = false;
        if(ii>0 && layer[q-si]==0) next = true;
        if(ii<p->imax-1 && layer[q+si]==0) next = true;
        if(jd && jj>0 && layer[q-sj]==0) next = true;
        if(jd && jj<p->jmax-1 && layer[q+sj]==0) next = true;
        if(kk<p->kmax-1 && layer[q+sk]==0) next = true;
        if(next)
        band[1].push_back(q);
    }
    for(int q : band[1])
    layer[q]=1;

    for(int l=2; l<=nband; ++l)
    {
        for(int r : band[l-1])
        {
            const int ii = r/si, jj = (r/sj)%p->jmax, kk = r%p->kmax;
            auto add = [&](int q)
            {
                if(insolid[q] && layer[q]==outside)
                {
                layer[q]=l;
                band[l].push_back(q);
                }
            };
            if(ii<p->imax-1) add(r+si);
            if(ii>0) add(r-si);
            if(jd && jj<p->jmax-1) add(r+sj);
            if(jd && jj>0) add(r-sj);
            if(kk<p->kmax-1) add(r+sk);
            if(kk>0) add(r-sk);
        }
        if(band[l].empty())
        break;
    }

    // a solid that touches the fluid only from below (a roof over the whole width at the top of the
    // domain): its cells next to that fluid start layers after the others, and they take the fluid
    // below them (belowok), so that no solid cell is left to the reinitialisation
    std::vector<char> belowok;
    {
        int lend = 1;
        while(lend<nband && !band[lend+1].empty())
        ++lend;

        std::vector<int> seed;
        for(int q : solid)
        if(layer[q]==outside)
        {
            const int kk = q%p->kmax;
            if(kk>0 && layer[q-sk]==0)
            seed.push_back(q);
        }

        if(!seed.empty() && lend+1<=nband)
        {
            belowok.assign(n,0);
            band[lend+1] = seed;
            for(int q : seed)
            {
            layer[q] = lend+1;
            belowok[q] = 1;
            }

            for(int l=lend+2; l<=nband; ++l)
            {
                for(int r : band[l-1])
                {
                    const int ii = r/si, jj = (r/sj)%p->jmax, kk = r%p->kmax;
                    auto add = [&](int q)
                    {
                        if(insolid[q] && layer[q]==outside)
                        {
                        layer[q]=l;
                        band[l].push_back(q);
                        }
                    };
                    if(ii<p->imax-1) add(r+si);
                    if(ii>0) add(r-si);
                    if(jd && jj<p->jmax-1) add(r+sj);
                    if(jd && jj>0) add(r-sj);
                    if(kk<p->kmax-1) add(r+sk);
                    if(kk>0) add(r-sk);
                }
                if(band[l].empty())
                break;
            }
        }
    }

    // (every rank takes part in the ghost cell exchange)
    for(int pass=0; pass<2; ++pass)
    {
        for(int l=1; l<=nband; ++l)
        for(int q : band[l])
        {
            i = q/si + p->imin;
            j = (q/sj)%p->jmax + p->jmin;
            k = q%p->kmax + p->kmin;

            // the cells of the partition (the ghost cells come from the exchange)
            if(i<0 || i>=p->knox || j<0 || j>=p->knoy || k<0 || k>=p->knoz)
            continue;

            double sum = 0.0;
            int count = 0;

            if(layer[q+si]<l)
            {
                sum += f(i+1,j,k);
                ++count;
            }
            if(layer[q-si]<l)
            {
                sum += f(i-1,j,k);
                ++count;
            }
            if(jd && layer[q+sj]<l)
            {
                sum += f(i,j+1,k);
                ++count;
            }
            if(jd && layer[q-sj]<l)
            {
                sum += f(i,j-1,k);
                ++count;
            }
            if(layer[q+sk]<l)
            {
                sum += towards(layer[q+sk],layer[q+2*sk]) ? extrapolate(f(i,j,k+1),f(i,j,k+2),p->DZP[KP],p->DZP[KP1],dn_lo,dn_hi) : f(i,j,k+1);
                ++count;
            }
            if(layer[q-sk]<l && !(layer[q-sk]==0 && (belowok.empty() || !belowok[q])))
            {
                sum += towards(layer[q-sk],layer[q-2*sk]) ? extrapolate(f(i,j,k-1),f(i,j,k-2),p->DZP[KM1],p->DZP[KM2],up_lo,up_hi) : f(i,j,k-1);
                ++count;
            }

            if(count>0)
            f(i,j,k) = sum/double(count);
        }

        pgc->start4(p,f,gcval);
    }
}

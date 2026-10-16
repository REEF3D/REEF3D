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

#include"coral/core/coral_model.h"

namespace coral
{

int model::add_body(const rigid_body &b)
{
    bodies.push_back(b);
    bodies.back().x0 = b.x;
    bodies.back().R0 = b.R;
    return nbody()-1;
}

int model::find_body(const std::string &name) const
{
    for(int i=0; i<nbody(); ++i)
        if(bodies[i].name==name)
            return i;
    return -1;
}

constraint &model::add_constraint(std::unique_ptr<constraint> c)
{
    constraint_list.push_back(std::move(c));
    return *constraint_list.back();
}

load &model::add_load(std::unique_ptr<load> l)
{
    load_list.push_back(std::move(l));
    return *load_list.back();
}

int model::ncon() const
{
    int n = 0;
    for(const auto &c : constraint_list)
        n += c->size();
    return n;
}

void model::forces(double t, VecX &f) const
{
    f.setZero(ndof());
    for(const auto &l : load_list)
        l->eval(*this, t, f);
}

void model::inertia_forces(VecX &r) const
{
    r.setZero(ndof());
    for(int i=0; i<nbody(); ++i)
    {
        const rigid_body &b = bodies[i];
        r.segment<3>(6*i) = b.mass*b.vd;
        r.segment<3>(6*i+3) = b.inertia*b.wd + b.w.cross(b.inertia*b.w);
    }
}

void model::constraints(VecX &g, MatX &G) const
{
    const int m = ncon();
    g.setZero(m);
    G.setZero(m, ndof());
    int row = 0;
    for(const auto &c : constraint_list)
    {
        const int k = c->size();
        c->eval(*this, g.segment(row,k), G.middleRows(row,k));
        row += k;
    }
}

void model::initialize(double t0)
{
    for(auto &l : load_list)
        l->initialize(*this, t0);
}

void model::commit(double t, double dt)
{
    for(auto &l : load_list)
        l->commit(*this, t, dt);
}

void model::add_force(VecX &f, int i, const Vec3 &p, const Vec3 &F) const
{
    if(i<0)
        return;
    const rigid_body &b = bodies[i];
    f.segment<3>(6*i) += F;
    f.segment<3>(6*i+3) += b.R.transpose()*((p - b.x).cross(F));
}

void model::add_moment(VecX &f, int i, const Vec3 &M) const
{
    if(i<0)
        return;
    f.segment<3>(6*i+3) += bodies[i].R.transpose()*M;
}

int model::coloring(std::vector<int> &color, std::vector<std::vector<int>> &nbr) const
{
    const int nb = nbody();
    std::vector<std::vector<char>> adj(nb, std::vector<char>(nb, 0));
    bool dense = false;
    auto clique = [&](const std::vector<int> &v)
    {
        for(int a : v)
            for(int b : v)
                if(a >= 0 && b >= 0)
                    adj[a][b] = 1;
    };
    std::vector<int> list;
    for(const auto &l : load_list)
    {
        const coupling_kind k = l->coupling();
        if(k==coupling_kind::global)
            dense = true;
        else if(k==coupling_kind::local)
        {
            list.clear();
            l->coupled_bodies(list);
            clique(list);
        }
    }
    for(const auto &c : constraint_list)
        clique({c->body_a, c->body_b});
    nbr.assign(nb, std::vector<int>());
    color.assign(nb, 0);
    for(int i=0; i<nb; ++i)
    {
        adj[i][i] = 1;
        for(int j=0; j<nb; ++j)
            if(dense || adj[i][j])
                nbr[i].push_back(j);
    }
    if(dense)
    {
        for(int i=0; i<nb; ++i) color[i] = i;
        return nb;
    }
    // greedy distance-2 colouring
    int nc = 0;
    std::vector<int> used;
    for(int i=0; i<nb; ++i)
    {
        used.assign(nc + 1, 0);
        for(int j : nbr[i])
        {
            if(j < i) used[color[j]] = 1;
            for(int k : nbr[j])
                if(k < i) used[color[k]] = 1;
        }
        int c = 0;
        while(c < nc && used[c]) ++c;
        color[i] = c;
        if(c==nc) ++nc;
    }
    return nc;
}

double model::kinetic_energy() const
{
    double e = 0.0;
    for(const auto &b : bodies)
        e += b.kinetic_energy();
    return e;
}

double model::potential_energy() const
{
    double e = 0.0;
    for(const auto &l : load_list)
        e += l->potential(*this);
    return e;
}

}

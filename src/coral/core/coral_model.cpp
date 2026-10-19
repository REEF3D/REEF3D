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
#include"coral/core/coral_state.h"
#include"coral/core/coral_rotation.h"
#include<algorithm>

namespace coral
{

int model::add_body(const rigid_body &b)
{
    bodies.push_back(b);
    bodies.back().x0 = b.x;
    bodies.back().R0 = b.R;
    return nbody()-1;
}

void model::add_point_mass(int i, const Vec3 &p, double mass, const Mat3 &J)
{
    if(i < 0)
        return;
    const rigid_body &b = bodies[i];
    point_mass q;
    q.body = i;
    q.s = b.to_body(p);
    q.mass = mass;
    q.J = b.R.transpose()*J*b.R;
    point_masses.push_back(q);
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

int model::add_coordinates(int n, const std::string &name)
{
    const int k = ncoord();
    cq.conservativeResize(k + n);
    cv.conservativeResize(k + n);
    ca.conservativeResize(k + n);
    cscale.conservativeResize(k + n);
    cq.tail(n).setZero();
    cv.tail(n).setZero();
    ca.tail(n).setZero();
    cscale.tail(n).setOnes();
    coord_blocks.emplace_back(name, k);
    return k;
}

void model::forces(double t, VecX &f, bool fd_only) const
{
    f.setZero(ndof());
    for(const auto &l : load_list)
        if(!fd_only || !l->has_tangent())
            l->eval(*this, t, f);
}

bool model::has_tangents() const
{
    for(const auto &l : load_list)
        if(l->has_tangent())
            return true;
    return false;
}

void model::tangents(double t, double cK, double cC, double cM, std::vector<triplet> &T) const
{
    for(const auto &l : load_list)
        if(l->has_tangent())
            l->tangent(*this, t, cK, cC, cM, T);
}

model::config model::save_config() const
{
    config c;
    c.x.reserve(bodies.size());
    c.R.reserve(bodies.size());
    for(const auto &b : bodies) {c.x.push_back(b.x); c.R.push_back(b.R);}
    c.cq = cq;
    return c;
}

void model::set_config(const config &c, const VecX &dq, double s)
{
    for(int i=0; i<nbody(); ++i)
    {
        rigid_body &b = bodies[i];
        b.x = c.x[i] + s*dq.segment<3>(6*i);
        b.R = c.R[i]*exp_so3(s*dq.segment<3>(6*i+3));
    }
    if(ncoord())
        cq = c.cq + s*dq.tail(ncoord());
}

double model::increment_norm(const VecX &d) const
{
    const int nb = 6*nbody();
    double r = nb ? d.head(nb).cwiseAbs().maxCoeff() : 0.0;
    for(int k=0; k<ncoord(); ++k)
        r = std::max(r, std::fabs(d(nb + k))*cscale(k));
    return r;
}

void model::zero_rates()
{
    for(auto &b : bodies)
    {
        b.v.setZero(); b.w.setZero(); b.vd.setZero(); b.wd.setZero();
    }
    cv.setZero(ncoord());
    ca.setZero(ncoord());
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
    // point masses at offsets: m a_p, a_p = vd + R (wd x s + w x (w x s)); moment s x (R^T m a_p)
    for(const auto &q : point_masses)
    {
        const rigid_body &b = bodies[q.body];
        const Vec3 ab = b.wd.cross(q.s) + b.w.cross(b.w.cross(q.s));
        const Vec3 Fb = q.mass*(b.R.transpose()*b.vd + ab);
        r.segment<3>(6*q.body) += b.R*Fb;
        r.segment<3>(6*q.body+3) += q.s.cross(Fb) + q.J*b.wd + b.w.cross(q.J*b.w);
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

int model::coloring(std::vector<int> &color, std::vector<std::vector<int>> &nbr, bool fd_only) const
{
    const int nb = nbody(), nn = nb + ncoord();
    std::vector<std::vector<int>> adj(nn);
    bool dense = false;
    std::vector<int> list;
    auto clique = [&]()
    {
        for(int a : list)
            for(int b : list)
                if(a >= 0 && b >= 0)
                    adj[a].push_back(b);
    };
    std::vector<int> cl;
    for(const auto &l : load_list)
    {
        if(fd_only && l->has_tangent())
            continue;
        const coupling_kind k = l->coupling();
        if(k==coupling_kind::global)
            dense = true;
        else if(k==coupling_kind::local)
        {
            list.clear();
            l->coupled_bodies(list);
            cl.clear();
            l->coupled_coordinates(cl);
            for(int c : cl) list.push_back(nb + c);
            clique();
        }
    }
    for(const auto &c : constraint_list)
    {
        list = {c->body_a, c->body_b};
        clique();
    }
    nbr.assign(nn, std::vector<int>());
    color.assign(nn, 0);
    for(int i=0; i<nn; ++i)
    {
        if(dense)
        {
            nbr[i].resize(nn);
            for(int j=0; j<nn; ++j) nbr[i][j] = j;
            color[i] = i;
            continue;
        }
        std::vector<int> &v = adj[i];
        v.push_back(i);
        std::sort(v.begin(), v.end());
        v.erase(std::unique(v.begin(), v.end()), v.end());
        nbr[i] = v;
    }
    if(dense)
        return nn;
    // greedy distance-2 colouring
    int nc = 0;
    std::vector<int> used;
    for(int i=0; i<nn; ++i)
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
    for(const auto &q : point_masses)
    {
        const rigid_body &b = bodies[q.body];
        e += 0.5*q.mass*(b.v + b.R*b.w.cross(q.s)).squaredNorm() + 0.5*b.w.dot(q.J*b.w);
    }
    return e;
}

double model::potential_energy() const
{
    double e = 0.0;
    for(const auto &l : load_list)
        e += l->potential(*this);
    return e;
}

void model::save_state(std::vector<double> &v) const
{
    state_out o{v};
    o.list(bodies, [&](const rigid_body &b){o.d(b.mass); o.m(b.inertia); o.m(b.x); o.m(b.R); o.m(b.v); o.m(b.w); o.m(b.vd); o.m(b.wd);});
    o.m(cq);
    o.m(cv);
    o.m(ca);
    o.i((long)load_list.size());
    for(const auto &l : load_list)
    {
        std::vector<double> s;
        l->save_state(s);
        o.vd(s);
    }
}

bool model::load_state(const std::vector<double> &v, size_t &k)
{
    state_in in{v, k};
    if(!in.has(1) || size_t(v[in.k]) != bodies.size())
        return false;
    in.list(bodies, [&](rigid_body &b){in.d(b.mass); in.m(b.inertia); in.m(b.x); in.m(b.R); in.m(b.v); in.m(b.w); in.m(b.vd); in.m(b.wd);});
    in.m(cq);
    in.m(cv);
    in.m(ca);
    int nl = 0;
    in.i(nl);
    if(!in.ok || nl != (int)load_list.size())
        return false;
    for(auto &l : load_list)
    {
        std::vector<double> s;
        in.vd(s);
        size_t j = 0;
        if(!in.ok || !l->load_state(s, j))
            return false;
    }
    k = in.k;
    return in.ok;
}

}

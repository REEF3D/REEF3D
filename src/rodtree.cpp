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

#include"rodtree.h"
#include<cmath>
#include<map>
#include<sstream>
#include<fstream>
#include<stdexcept>
#include<algorithm>

namespace
{
    const double rt_pi = 3.14159265358979323846;

    std::runtime_error input_error(int line, const std::string& msg)
    {
        std::ostringstream s;
        s<<"rodtree.dat line "<<line<<": "<<msg;
        return std::runtime_error(s.str());
    }
}

rodtree::rodtree() : g(0.0,0.0,-9.81), rhof(1000.0), t(0.0), integrator(0), substeps(1), nsub_last(1), reaction_mode(0), hsub_prev(0.0), pattern_ready(false)
{
}

// ---------------------------------------------------------------------------
// input
// ---------------------------------------------------------------------------

void rodtree::read(std::istream& in)
{
    std::string line;
    int lineno = 0;
    pcolony *cur = nullptr;

    while(std::getline(in,line))
    {
        ++lineno;
        size_t hash = line.find('#');
        if(hash!=std::string::npos)
        line = line.substr(0,hash);

        std::istringstream ls(line);
        std::string key;
        if(!(ls>>key))
        continue;

        // global keywords (outside colony blocks)
        if(!cur && key=="integrator")
        {
            std::string m;
            ls>>m;
            if(m=="implicit") integrator = 0;
            else if(m=="explicit") integrator = 1;
            else throw input_error(lineno,"integrator must be 'implicit' or 'explicit'");
            continue;
        }
        if(!cur && key=="substeps")
        {
            if(!(ls>>substeps) || substeps<1)
            throw input_error(lineno,"substeps needs an integer >= 1");
            continue;
        }
        if(!cur && key=="reaction")
        {
            std::string m;
            ls>>m;
            if(m=="full") reaction_mode = 0;
            else if(m=="drag") reaction_mode = 1;
            else if(m=="none") reaction_mode = 2;
            else throw input_error(lineno,"reaction must be 'full', 'drag' or 'none'");
            continue;
        }

        if(key=="colony")
        {
            if(cur)
            throw input_error(lineno,"'colony' inside an open colony (missing 'end')");
            input.push_back(pcolony());
            cur = &input.back();
            if(!(ls>>cur->name))
            cur->name = "colony" + std::to_string(input.size());
            continue;
        }

        if(!cur)
        throw input_error(lineno,"'" + key + "' outside a colony block");

        if(key=="end")
        {
            cur = nullptr;
        }
        else if(key=="material")
        {
            if(!(ls>>cur->mat.E>>cur->mat.nu>>cur->mat.rho))
            throw input_error(lineno,"material needs: E nu rho");
        }
        else if(key=="hydro")
        {
            if(!(ls>>cur->mat.Cdn>>cur->mat.Cdt>>cur->mat.Ca))
            throw input_error(lineno,"hydro needs: Cdn Cdt Ca [Cm]");
            if(!(ls>>cur->mat.Cm))
            cur->mat.Cm = 1.0 + cur->mat.Ca;
        }
        else if(key=="damping")
        {
            if(!(ls>>cur->mat.beta))
            throw input_error(lineno,"damping needs: beta [s]");
        }
        else if(key=="polyps")
        {
            if(!(ls>>cur->mat.hp>>cur->mat.phip>>cur->mat.Cdp) || cur->mat.hp<0.0 || cur->mat.phip<0.0 || cur->mat.phip>1.0)
            throw input_error(lineno,"polyps needs: h_p phi_p Cd_p [k_t]  (h_p>=0, 0<=phi_p<=1)");
            if(!(ls>>cur->mat.ktp))
            cur->mat.ktp = 1.0;
        }
        else if(key=="polyp_response")
        {
            if(!(ls>>cur->mat.Ur>>cur->mat.dUr>>cur->mat.taup) || cur->mat.Ur<0.0 || cur->mat.dUr<0.0 || cur->mat.taup<0.0)
            throw input_error(lineno,"polyp_response needs: U_retract dU tau  (all >= 0)");
        }
        else if(key=="refine")
        {
            if(!(ls>>cur->refine) || cur->refine<1)
            throw input_error(lineno,"refine needs an integer >= 1");
        }
        else if(key=="node")
        {
            pnode n;
            if(!(ls>>n.id>>n.x(0)>>n.x(1)>>n.x(2)>>n.r) || n.r<=0.0)
            throw input_error(lineno,"node needs: id x y z r (r>0)");
            cur->nodes.push_back(n);
        }
        else if(key=="edge")
        {
            pedge e;
            if(!(ls>>e.p>>e.c))
            throw input_error(lineno,"edge needs: parent_node child_node");
            cur->edges.push_back(e);
        }
        else if(key=="clamp")
        {
            int id;
            if(!(ls>>id))
            throw input_error(lineno,"clamp needs: node id");
            cur->clamps.push_back(id);
        }
        else if(key=="instance")
        {
            Vec3 d;
            if(!(ls>>d(0)>>d(1)>>d(2)))
            throw input_error(lineno,"instance needs: dx dy dz");
            cur->instances.push_back(d);
        }
        else
        throw input_error(lineno,"unknown keyword '" + key + "'");
    }

    if(cur)
    throw std::runtime_error("rodtree.dat: last colony block not closed with 'end'");
}

void rodtree::finalize_setup()
{
    el.clear();
    jt.clear();
    col.clear();

    for(const pcolony& pc : input)
    {
        std::map<int,int> idx;
        for(size_t i=0; i<pc.nodes.size(); ++i)
        {
            if(idx.count(pc.nodes[i].id))
            throw std::runtime_error("rodtree: duplicate node id " + std::to_string(pc.nodes[i].id) + " in colony " + pc.name);
            idx[pc.nodes[i].id] = (int)i;
        }

        const int nn = (int)pc.nodes.size();
        std::vector<int> parent_edge(nn,-1);
        for(size_t k=0; k<pc.edges.size(); ++k)
        {
            if(!idx.count(pc.edges[k].p) || !idx.count(pc.edges[k].c))
            throw std::runtime_error("rodtree: edge references unknown node in colony " + pc.name);
            int c = idx[pc.edges[k].c];
            if(parent_edge[c]>=0)
            throw std::runtime_error("rodtree: node " + std::to_string(pc.edges[k].c) + " has two parent edges (colony " + pc.name + " must be a tree)");
            parent_edge[c] = (int)k;
        }

        std::vector<char> clamped(nn,0);
        for(int id : pc.clamps)
        {
            if(!idx.count(id))
            throw std::runtime_error("rodtree: clamp references unknown node in colony " + pc.name);
            if(parent_edge[idx[id]]>=0)
            throw std::runtime_error("rodtree: clamped node " + std::to_string(id) + " has a parent edge (colony " + pc.name + ")");
            clamped[idx[id]] = 1;
        }

        // every root of an edge must be clamped; check for cycles by walking up
        for(int i=0; i<nn; ++i)
        {
            int n = i, steps = 0;
            while(parent_edge[n]>=0)
            {
                n = idx[pc.edges[parent_edge[n]].p];
                if(++steps>nn)
                throw std::runtime_error("rodtree: cycle in colony " + pc.name);
            }
            bool has_edges = parent_edge[i]>=0;
            for(const pedge& e : pc.edges)
            if(idx[e.p]==i) has_edges = true;
            if(has_edges && !clamped[n])
            throw std::runtime_error("rodtree: root node " + std::to_string(pc.nodes[n].id) + " is not clamped (colony " + pc.name + ")");
        }

        // process edges parents-first
        std::vector<int> order, depth(pc.edges.size(),0);
        for(size_t k=0; k<pc.edges.size(); ++k)
        {
            int n = idx[pc.edges[k].p];
            while(parent_edge[n]>=0) {++depth[k]; n = idx[pc.edges[parent_edge[n]].p];}
            order.push_back((int)k);
        }
        std::stable_sort(order.begin(),order.end(),[&](int a, int b){return depth[a]<depth[b];});

        std::vector<Vec3> offsets = pc.instances;
        if(offsets.empty())
        offsets.push_back(Vec3::Zero());

        for(size_t inst=0; inst<offsets.size(); ++inst)
        {
            colony_info ci;
            ci.name = pc.name + (offsets.size()>1 ? "_" + std::to_string(inst) : "");
            ci.mat = pc.mat;
            const int cid = (int)col.size();

            std::vector<int> last_elem_of_edge(pc.edges.size(),-1);
            const material& mat = pc.mat;
            const double G = mat.E/(2.0*(1.0+mat.nu));

            for(int k : order)
            {
                const pnode& np = pc.nodes[idx[pc.edges[k].p]];
                const pnode& nc = pc.nodes[idx[pc.edges[k].c]];
                Vec3 x0 = np.x + offsets[inst];
                Vec3 x1 = nc.x + offsets[inst];
                Vec3 d = x1 - x0;
                double L = d.norm();
                if(L<=0.0)
                throw std::runtime_error("rodtree: zero-length edge in colony " + pc.name);
                Vec3 dir = d/L;
                Quat q = Quat::FromTwoVectors(Vec3(0,0,1),dir);
                q.normalize();

                int parent_elem = -1;
                int pe = parent_edge[idx[pc.edges[k].p]];
                if(pe>=0)
                parent_elem = last_elem_of_edge[pe];

                const int ns = pc.refine;
                for(int s=0; s<ns; ++s)
                {
                    double s0 = double(s)/ns, s1 = double(s+1)/ns;
                    double r0 = np.r + (nc.r-np.r)*s0;
                    double r1 = np.r + (nc.r-np.r)*s1;
                    element e;
                    e.colony = cid;
                    e.n0 = pc.edges[k].p; e.n1 = pc.edges[k].c;
                    e.l = L/ns;
                    e.r = 0.5*(r0+r1);
                    double V = rt_pi*e.l*(r0*r0 + r0*r1 + r1*r1)/3.0;
                    e.m = mat.rho*V;
                    e.c = x0 + (0.5*(s0+s1)*L)*dir;
                    e.q = q;
                    e.c0 = e.c; e.q0 = e.q;
                    e.v.setZero(); e.w.setZero(); e.a.setZero();
                    e.uf.setZero(); e.af.setZero(); e.chi = 0.0;
                    e.Fh.setZero(); e.Fdrag.setZero(); e.Finert.setZero();
                    e.ext = (mat.hp>0.0 && mat.phip>0.0) ? 1.0 : 0.0;
                    const int eid = (int)el.size();
                    el.push_back(e);
                    ci.elements.push_back(eid);

                    joint j;
                    j.b = eid;
                    j.beta = mat.beta;
                    j.F.setZero(); j.T.setZero();
                    double rj = r0, lj;
                    j.pb = Vec3(0,0,-0.5*e.l);
                    if(parent_elem>=0)
                    {
                        const element& ea = el[parent_elem];
                        j.a = parent_elem;
                        j.pa = Vec3(0,0,0.5*ea.l);
                        j.qrel0 = ea.q.conjugate()*e.q;
                        lj = 0.5*(ea.l + e.l);
                    }
                    else
                    {
                        j.a = -1;
                        j.pa = x0;                   // world anchor point
                        j.qrel0 = e.q;
                        lj = 0.5*e.l;
                        ci.roots.push_back((int)jt.size());
                    }
                    j.qrel0.normalize();
                    const double A = rt_pi*rj*rj, I = 0.25*rt_pi*rj*rj*rj*rj;
                    j.k  = mat.E*A/lj;
                    j.kb = mat.E*I/lj;
                    j.kt = G*2.0*I/lj;
                    jt.push_back(j);

                    parent_elem = eid;
                }
                last_elem_of_edge[k] = parent_elem;
            }

            // tip: element end farthest from the first root anchor
            if(!ci.roots.empty())
            {
                Vec3 base = jt[ci.roots[0]].pa;
                double dmax = -1.0;
                for(int eid : ci.elements)
                {
                    double dd = (end1(eid)-base).norm();
                    if(dd>dmax) {dmax = dd; ci.tip_node = eid;}
                }
            }
            col.push_back(ci);
        }
    }

    pattern_ready = false;
}

// ---------------------------------------------------------------------------
// helpers
// ---------------------------------------------------------------------------

rodtree::Vec3 rodtree::logq(const Quat& qin)
{
    Quat q = qin;
    if(q.w()<0.0) q.coeffs() *= -1.0;
    Vec3 v = q.vec();
    double s = v.norm();
    if(s<1.0e-12)
    return 2.0*v/std::max(q.w(),1.0e-300);
    double th = 2.0*std::atan2(s,q.w());
    return (th/s)*v;
}

rodtree::Quat rodtree::expq(const Vec3& th)
{
    double a = th.norm();
    if(a<1.0e-12)
    {
        Quat q(1.0,0.5*th(0),0.5*th(1),0.5*th(2));
        q.normalize();
        return q;
    }
    Vec3 ax = th/a;
    double s = std::sin(0.5*a);
    return Quat(std::cos(0.5*a),s*ax(0),s*ax(1),s*ax(2));
}

double rodtree::volume(int e) const
{
    return el[e].m/col[el[e].colony].mat.rho;
}

rodtree::Mat3 rodtree::mass_trans(int e) const
{
    const element& E = el[e];
    Vec3 tt = axis(e);
    Mat3 P = Mat3::Identity() - tt*tt.transpose();
    double ma = E.chi*col[E.colony].mat.Ca*rhof*volume(e);
    return E.m*Mat3::Identity() + ma*P;
}

rodtree::Mat3 rodtree::inertia_world(int e, bool with_fluid) const
{
    const element& E = el[e];
    double It = E.m*(0.25*E.r*E.r + E.l*E.l/12.0);
    double Ia = 0.5*E.m*E.r*E.r;
    if(with_fluid)
    It += E.chi*col[E.colony].mat.Ca*rhof*volume(e)*E.l*E.l/12.0;
    Mat3 R = E.q.toRotationMatrix();
    return R*Vec3(It,It,Ia).asDiagonal()*R.transpose();
}

void rodtree::joint_force(const joint& j, const Vec3& ca, const Quat& qa, const Vec3& va, const Vec3& wa,
                          const Vec3& cb, const Quat& qb, const Vec3& vb, const Vec3& wb,
                          Vec3& Fa, Vec3& Ta, Vec3& Fb, Vec3& Tb) const
{
    // attachment points (world) and their velocities
    Vec3 ra = (j.a>=0) ? Vec3(qa*j.pa) : Vec3::Zero();
    Vec3 xa = (j.a>=0) ? Vec3(ca + ra) : j.pa;
    Vec3 xadot = (j.a>=0) ? Vec3(va + wa.cross(ra)) : Vec3::Zero();
    Vec3 rb = qb*j.pb;
    Vec3 xb = cb + rb;
    Vec3 xbdot = vb + wb.cross(rb);

    Vec3 gap = xb - xa;
    Vec3 gdot = xbdot - xadot;
    Fb = -j.k*(gap + j.beta*gdot);
    Fa = -Fb;

    // relative rotation, error vector in parent frame
    Quat qA = (j.a>=0) ? qa : Quat::Identity();
    Quat qrel = qA.conjugate()*qb;
    Vec3 th = logq(qrel*j.qrel0.conjugate());
    Vec3 tp = j.qrel0*Vec3(0,0,1);
    Mat3 K = j.kb*(Mat3::Identity() - tp*tp.transpose()) + j.kt*tp*tp.transpose();
    Vec3 wrel = qA.conjugate()*(wb - ((j.a>=0) ? wa : Vec3::Zero()));
    Vec3 Tp = -K*(th + j.beta*wrel);
    Vec3 Tw = qA*Tp;

    Tb = Tw + rb.cross(Fb);
    Ta = -Tw + ra.cross(Fa);
}

void rodtree::drag_coeffs(int e, double& cn, double& ct) const
{
    // drag per |u|u, normal and tangential, incl. the polyp layer:
    //   normal:     0.5 rho l chi (Cdn D + Cdp * 2 ext h_p phi_p)
    //   tangential: 0.5 rho l chi (Cdt pi D + k_t Cdp * 2 ext h_p phi_p)
    const element& E = el[e];
    const material& mat = col[E.colony].mat;
    const double D = 2.0*E.r;
    const double bp = 2.0*E.ext*mat.hp*mat.phip;          // extra frontal width of the polyps
    cn = 0.5*rhof*E.l*E.chi*(mat.Cdn*D + mat.Cdp*bp);
    ct = 0.5*rhof*E.l*E.chi*(mat.Cdt*rt_pi*D + mat.ktp*mat.Cdp*bp);
}

double rodtree::polyp_target(int e) const
{
    const element& E = el[e];
    const material& mat = col[E.colony].mat;
    if(mat.hp<=0.0 || mat.phip<=0.0) return 0.0;
    if(mat.Ur<0.0) return 1.0;                             // no response: always extended
    Vec3 tt = axis(e);
    Vec3 ur = E.uf - E.v;
    double un = (ur - tt*tt.dot(ur)).norm();
    if(un<=mat.Ur) return 1.0;
    if(mat.dUr<=0.0 || un>=mat.Ur+mat.dUr) return 0.0;
    return 1.0 - (un-mat.Ur)/mat.dUr;
}

void rodtree::update_polyps(double dt)
{
    // first-order response, implicit in time: d ext/dt = (target - ext)/tau
    for(int e=0; e<nelem(); ++e)
    {
        const material& mat = col[el[e].colony].mat;
        double tg = polyp_target(e);
        if(el[e].chi<=0.0 && mat.Ur>=0.0) tg = 0.0;           // emerged polyps retract
        if(mat.taup<=0.0) el[e].ext = tg;
        else el[e].ext = (el[e].ext + dt/mat.taup*tg)/(1.0 + dt/mat.taup);
    }
}

void rodtree::external_loads(int e, Vec3& F, Vec3& T, Mat3& Cdrag) const
{
    const element& E = el[e];
    const material& mat = col[E.colony].mat;
    const double V = volume(e);

    F = E.m*g - E.chi*rhof*V*g;
    T.setZero();
    Cdrag.setZero();

    if(E.chi<=0.0)
    return;

    Vec3 tt = axis(e);
    Mat3 P = Mat3::Identity() - tt*tt.transpose();
    Vec3 ur = E.uf - E.v;
    Vec3 un = P*ur;
    Vec3 ut = tt*tt.dot(ur);
    double cn, ct;
    drag_coeffs(e,cn,ct);
    double unm = un.norm(), utm = ut.norm();

    F += cn*unm*un + ct*utm*ut;
    F += E.chi*rhof*V*mat.Cm*(P*E.af);

    Cdrag = cn*unm*P + 2.0*ct*utm*tt*tt.transpose();
    if(unm>1.0e-12)
    Cdrag += cn*un*un.transpose()/unm;
}

void rodtree::compute_hydro()
{
    for(int e=0; e<nelem(); ++e)
    {
        element& E = el[e];
        E.Fdrag.setZero(); E.Finert.setZero(); E.Fh.setZero();
        if(E.chi<=0.0) continue;
        const material& mat = col[E.colony].mat;
        Vec3 tt = axis(e);
        Mat3 P = Mat3::Identity() - tt*tt.transpose();
        Vec3 ur = E.uf - E.v;
        Vec3 un = P*ur;
        Vec3 ut = tt*tt.dot(ur);
        double cn, ct;
        drag_coeffs(e,cn,ct);
        E.Fdrag = cn*un.norm()*un + ct*ut.norm()*ut;
        E.Finert = E.chi*rhof*volume(e)*(P*(mat.Cm*E.af - mat.Ca*E.a));
        E.Fh = E.Fdrag + E.Finert;
    }
}

rodtree::Vec3 rodtree::fluid_reaction(int e) const
{
    if(reaction_mode==1) return -el[e].Fdrag;
    if(reaction_mode==2) return Vec3::Zero();
    return -el[e].Fh;
}

double rodtree::drag_slope(int e) const
{
    const element& E = el[e];
    if(E.chi<=0.0) return 0.0;
    Vec3 tt = axis(e);
    Vec3 ur = E.uf - E.v;
    Vec3 un = ur - tt*tt.dot(ur);
    double cn, ct;
    drag_coeffs(e,cn,ct);
    return 2.0*cn*un.norm();
}

void rodtree::assemble_forces(std::vector<Vec3>& F, std::vector<Vec3>& T, bool store)
{
    const int ne = nelem();
    F.assign(ne,Vec3::Zero());
    T.assign(ne,Vec3::Zero());

    Mat3 Cd;
    for(int e=0; e<ne; ++e)
    {
        Vec3 f, tq;
        external_loads(e,f,tq,Cd);
        F[e] += f;
        T[e] += tq;
        Mat3 Iw = inertia_world(e,true);
        T[e] -= el[e].w.cross(Iw*el[e].w);
    }

    const Vec3 z = Vec3::Zero();
    const Quat qi = Quat::Identity();
    for(joint& j : jt)
    {
        Vec3 Fa,Ta,Fb,Tb;
        const element& B = el[j.b];
        if(j.a>=0)
        {
            const element& A = el[j.a];
            joint_force(j,A.c,A.q,A.v,A.w,B.c,B.q,B.v,B.w,Fa,Ta,Fb,Tb);
            F[j.a] += Fa; T[j.a] += Ta;
        }
        else
        joint_force(j,z,qi,z,z,B.c,B.q,B.v,B.w,Fa,Ta,Fb,Tb);

        F[j.b] += Fb; T[j.b] += Tb;
        if(store) {j.F = Fb; j.T = Tb;}
    }
}

// ---------------------------------------------------------------------------
// time integration
// ---------------------------------------------------------------------------

void rodtree::advance(double dt)
{
    if(nelem()==0 || dt<=0.0)
    return;

    update_polyps(dt);

    if(integrator==1)
    {
        double h = explicit_dt();
        double nd = std::ceil(dt/h);
        // Kelvin-Voigt damping is stiffness proportional: its explicit limit
        // h < 2/(beta*omega^2) collapses for stiff joints -> use the implicit
        // integrator whenever beta > 0 is needed
        if(nd>2.0e4)
        throw std::runtime_error("rodtree: explicit integrator would need more than 20000 sub-steps per time step "
                                 "(damping beta too large or joints too stiff); use the implicit integrator");
        int n = std::max(1,(int)nd);
        nsub_last = n;
        for(int s=0; s<n; ++s)
        step_explicit(dt/n);
    }
    else
    {
        // Linearly implicit Euler does one Newton step per sub-step; it is
        // unconditionally stable for the linearised system but diverges when
        // an element rotates/translates too much within one sub-step.
        // Adaptive sub-stepping with roll-back: reject a sub-step if any
        // element rotates by more than dtheta_max or moves by more than
        // dx_max*l, halve h and retry; grow h again after success.
        const double dtheta_max = 0.25, dx_max = 0.25;
        const double hnom = dt/substeps, hmin = dt/65536.0;
        double h = std::min(hnom,hsub_prev>0.0 ? 2.0*hsub_prev : hnom);
        double rem = dt;
        int nacc = 0;
        std::vector<element> save;

        while(rem>1.0e-12*dt)
        {
            h = std::min(h,rem);
            if(rem-h<0.01*h) h = rem;           // avoid a sliver step
            save = el;
            step_implicit(h);

            bool ok = true;
            for(const element& E : el)
            {
                double dth = E.w.norm()*h, dxl = E.v.norm()*h/E.l;
                if(!std::isfinite(dth) || !std::isfinite(dxl) || dth>dtheta_max || dxl>dx_max)
                {ok = false; break;}
            }

            if(!ok && h>hmin)
            {
                el = save;
                h *= 0.5;
                continue;
            }
            if(!ok)
            throw std::runtime_error("rodtree: implicit sub-step below dt/65536 still rejected (structure unstable)");

            rem -= h;
            ++nacc;
            hsub_prev = h;
            h = std::min(2.0*h,hnom);
        }
        nsub_last = nacc;
    }

    // store reaction forces of the new state
    std::vector<Vec3> F,T;
    assemble_forces(F,T,true);
    t += dt;
}

double rodtree::explicit_dt() const
{
    double hmin = 1.0e30;
    for(const joint& j : jt)
    {
        double m = el[j.b].m, It = el[j.b].m*(0.25*el[j.b].r*el[j.b].r + el[j.b].l*el[j.b].l/12.0), Ia = 0.5*el[j.b].m*el[j.b].r*el[j.b].r;
        if(j.a>=0)
        {
            const element& A = el[j.a];
            m = std::min(m,A.m);
            It = std::min(It,A.m*(0.25*A.r*A.r + A.l*A.l/12.0));
            Ia = std::min(Ia,0.5*A.m*A.r*A.r);
        }
        double w2 = std::max(std::max(3.0*j.k/m, j.kb/It), j.kt/Ia);
        double w = 2.0*std::sqrt(w2);
        double zeta = 0.5*j.beta*w;
        double h = 0.5*(2.0/w)*(std::sqrt(1.0+zeta*zeta) - zeta);
        hmin = std::min(hmin,h);
    }
    return hmin;
}

void rodtree::step_explicit(double h)
{
    std::vector<Vec3> F,T;
    assemble_forces(F,T,false);

    Mat3 Cd;
    for(int e=0; e<nelem(); ++e)
    {
        element& E = el[e];
        Vec3 f, tq;
        external_loads(e,f,tq,Cd);                 // for the drag Jacobian only
        Mat3 M = mass_trans(e) + h*Cd;
        Vec3 dv = M.ldlt().solve(h*F[e]);
        E.v += dv;
        E.a = dv/h;
        Mat3 Iw = inertia_world(e,true);
        E.w += Iw.ldlt().solve(h*T[e]);
        E.c += h*E.v;
        E.q = expq(h*E.w)*E.q;
        E.q.normalize();
    }
}

void rodtree::step_implicit(double h)
{
    const int ne = nelem();
    const int N = 6*ne;

    std::vector<Vec3> F,T;
    assemble_forces(F,T,false);

    Eigen::VectorXd rhs(N), u(N), Ku(N);
    for(int e=0; e<ne; ++e)
    {
        u.segment<3>(6*e) = el[e].v;
        u.segment<3>(6*e+3) = el[e].w;
        rhs.segment<3>(6*e) = F[e];
        rhs.segment<3>(6*e+3) = T[e];
    }
    Ku.setZero();

    trip.clear();
    trip.reserve(ne*18 + jt.size()*288);

    // mass and drag blocks
    Mat3 Cd;
    for(int e=0; e<ne; ++e)
    {
        Vec3 f, tq;
        external_loads(e,f,tq,Cd);
        Mat3 Mt = mass_trans(e) + h*Cd;
        Mat3 Iw = inertia_world(e,true);
        // gyroscopic term -w x (I w) linearised implicitly: explicit treatment is
        // unstable for slender elements (I_axial << I_transverse)
        {
            const Vec3& w = el[e].w;
            Vec3 Lw = Iw*w;
            Mat3 Sw, SL;
            Sw << 0.0,-w(2),w(1), w(2),0.0,-w(0), -w(1),w(0),0.0;
            SL << 0.0,-Lw(2),Lw(1), Lw(2),0.0,-Lw(0), -Lw(1),Lw(0),0.0;
            Iw += h*(Sw*Iw - SL);                // M - h d(-w x Iw)/dw
        }
        for(int r=0; r<3; ++r)
        for(int c=0; c<3; ++c)
        {
            trip.emplace_back(6*e+r,6*e+c,Mt(r,c));
            trip.emplace_back(6*e+3+r,6*e+3+c,Iw(r,c));
        }
    }

    // joint Jacobians by central differences
    const Vec3 z = Vec3::Zero();
    const Quat qi = Quat::Identity();
    for(const joint& j : jt)
    {
        const int nb = (j.a>=0) ? 2 : 1;
        int body[2] = {j.b, j.a};

        for(int sb=0; sb<nb; ++sb)
        {
            const int s = body[sb];
            const double epsx = 1.0e-6*std::max(el[s].l,1.0e-6);
            const double epsr = 1.0e-6;
            const double epsv = 1.0e-6;

            for(int kind=0; kind<2; ++kind)      // 0: position, 1: velocity
            for(int k=0; k<6; ++k)
            {
                Eigen::Matrix<double,12,1> col_d;
                double eps = kind==0 ? (k<3 ? epsx : epsr) : epsv;
                Eigen::Matrix<double,12,1> fp, fm;

                for(int sign=-1; sign<=1; sign+=2)
                {
                    Vec3 c[2], v[2], w[2]; Quat q[2];
                    for(int bb=0; bb<nb; ++bb)
                    {
                        const element& E = el[body[bb]];
                        c[bb]=E.c; q[bb]=E.q; v[bb]=E.v; w[bb]=E.w;
                    }
                    Vec3 d = Vec3::Zero(); d(k%3) = sign*eps;
                    if(kind==0) {if(k<3) c[sb] += d; else q[sb] = expq(d)*q[sb];}
                    else        {if(k<3) v[sb] += d; else w[sb] += d;}

                    Vec3 Fa,Ta,Fb,Tb;
                    if(nb==2)
                    joint_force(j,c[1],q[1],v[1],w[1],c[0],q[0],v[0],w[0],Fa,Ta,Fb,Tb);
                    else
                    joint_force(j,z,qi,z,z,c[0],q[0],v[0],w[0],Fa,Ta,Fb,Tb);

                    Eigen::Matrix<double,12,1>& f = sign<0 ? fm : fp;
                    f.segment<3>(0)=Fb; f.segment<3>(3)=Tb; f.segment<3>(6)=Fa; f.segment<3>(9)=Ta;
                }
                col_d = (fp-fm)/(2.0*eps);

                const int colidx = 6*s + k;
                for(int rb=0; rb<nb; ++rb)
                {
                    const int r = body[rb];
                    for(int i=0; i<6; ++i)
                    {
                        double dfd = col_d(6*rb+i);
                        if(dfd==0.0) continue;
                        if(kind==0)
                        {
                            trip.emplace_back(6*r+i,colidx,-h*h*dfd);
                            Ku(6*r+i) += dfd*u(colidx);
                        }
                        else
                        trip.emplace_back(6*r+i,colidx,-h*dfd);
                    }
                }
            }
        }
    }

    Eigen::SparseMatrix<double> A(N,N);
    A.setFromTriplets(trip.begin(),trip.end());
    A.makeCompressed();

    if(!pattern_ready)
    {
        lu.analyzePattern(A);
        pattern_ready = true;
    }
    lu.factorize(A);
    if(lu.info()!=Eigen::Success)
    {
        // pattern may have changed (exact zeros); redo the analysis once
        lu.analyzePattern(A);
        lu.factorize(A);
        if(lu.info()!=Eigen::Success)
        throw std::runtime_error("rodtree: sparse LU factorisation failed");
    }

    Eigen::VectorXd du = lu.solve(h*(rhs + h*Ku));

    for(int e=0; e<ne; ++e)
    {
        element& E = el[e];
        Vec3 dv = du.segment<3>(6*e);
        E.v += dv;
        E.a = dv/h;
        E.w += du.segment<3>(6*e+3);
        E.c += h*E.v;
        E.q = expq(h*E.w)*E.q;
        E.q.normalize();
    }
}

// ---------------------------------------------------------------------------
// diagnostics and output
// ---------------------------------------------------------------------------

double rodtree::kinetic_energy() const
{
    double ek = 0.0;
    for(int e=0; e<nelem(); ++e)
    {
        const element& E = el[e];
        ek += 0.5*E.m*E.v.squaredNorm() + 0.5*E.w.dot(inertia_world(e,false)*E.w);
    }
    return ek;
}

double rodtree::elastic_energy() const
{
    double ep = 0.0;
    for(const joint& j : jt)
    {
        const element& B = el[j.b];
        Vec3 xa = (j.a>=0) ? Vec3(el[j.a].c + el[j.a].q*j.pa) : j.pa;
        Vec3 xb = B.c + B.q*j.pb;
        ep += 0.5*j.k*(xb-xa).squaredNorm();
        Quat qA = (j.a>=0) ? el[j.a].q : Quat::Identity();
        Vec3 th = logq((qA.conjugate()*B.q)*j.qrel0.conjugate());
        Vec3 tp = j.qrel0*Vec3(0,0,1);
        double tz = th.dot(tp);
        ep += 0.5*j.kb*(th.squaredNorm()-tz*tz) + 0.5*j.kt*tz*tz;
    }
    return ep;
}

rodtree::Vec3 rodtree::tip_displacement(int c) const
{
    int e = col[c].tip_node;
    if(e<0) return Vec3::Zero();
    Vec3 x0 = el[e].c0 + 0.5*el[e].l*(el[e].q0*Vec3(0,0,1));
    return end1(e) - x0;
}

rodtree::Vec3 rodtree::base_force(int c) const
{
    Vec3 f = Vec3::Zero();
    for(int j : col[c].roots)
    f -= jt[j].F;                                   // force of the colony on the bed
    return f;
}

rodtree::Vec3 rodtree::hydro_force(int c) const
{
    Vec3 f = Vec3::Zero();
    for(int e : col[c].elements)
    f += el[e].Fh;
    return f;
}

double rodtree::polyp_extension(int c) const
{
    double s = 0.0;
    for(int e : col[c].elements) s += el[e].ext;
    return col[c].elements.empty() ? 0.0 : s/col[c].elements.size();
}

void rodtree::write_vtp(const std::string& filename) const
{
    std::ofstream out(filename.c_str());
    const int ne = nelem();
    out<<"<?xml version=\"1.0\"?>\n<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\">\n<PolyData>\n";
    out<<"<Piece NumberOfPoints=\""<<2*ne<<"\" NumberOfLines=\""<<ne<<"\">\n";
    out<<"<Points>\n<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for(int e=0; e<ne; ++e)
    {
        Vec3 a = end0(e), b = end1(e);
        out<<a(0)<<" "<<a(1)<<" "<<a(2)<<" "<<b(0)<<" "<<b(1)<<" "<<b(2)<<"\n";
    }
    out<<"</DataArray>\n</Points>\n<PointData Scalars=\"radius\">\n<DataArray type=\"Float64\" Name=\"radius\" format=\"ascii\">\n";
    for(int e=0; e<ne; ++e)
    out<<el[e].r<<" "<<el[e].r<<"\n";
    out<<"</DataArray>\n<DataArray type=\"Float64\" Name=\"velocity\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for(int e=0; e<ne; ++e)
    {
        Vec3 v0 = el[e].v + el[e].w.cross(end0(e)-el[e].c);
        Vec3 v1 = el[e].v + el[e].w.cross(end1(e)-el[e].c);
        out<<v0(0)<<" "<<v0(1)<<" "<<v0(2)<<" "<<v1(0)<<" "<<v1(1)<<" "<<v1(2)<<"\n";
    }
    out<<"</DataArray>\n</PointData>\n<CellData>\n<DataArray type=\"Int32\" Name=\"colony\" format=\"ascii\">\n";
    for(int e=0; e<ne; ++e)
    out<<el[e].colony<<"\n";
    out<<"</DataArray>\n<DataArray type=\"Float64\" Name=\"hydro_force\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for(int e=0; e<ne; ++e)
    out<<el[e].Fh(0)<<" "<<el[e].Fh(1)<<" "<<el[e].Fh(2)<<"\n";
    out<<"</DataArray>\n<DataArray type=\"Float64\" Name=\"submerged\" format=\"ascii\">\n";
    for(int e=0; e<ne; ++e)
    out<<el[e].chi<<"\n";
    out<<"</DataArray>\n<DataArray type=\"Float64\" Name=\"polyp_extension\" format=\"ascii\">\n";
    for(int e=0; e<ne; ++e)
    out<<el[e].ext<<"\n";
    out<<"</DataArray>\n</CellData>\n<Lines>\n<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    for(int e=0; e<ne; ++e)
    out<<2*e<<" "<<2*e+1<<"\n";
    out<<"</DataArray>\n<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    for(int e=0; e<ne; ++e)
    out<<2*e+2<<"\n";
    out<<"</DataArray>\n</Lines>\n</Piece>\n</PolyData>\n</VTKFile>\n";
}

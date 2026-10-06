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

#include"fem_solid.h"
#include<cmath>
#include<stdexcept>
#include<algorithm>
#include<sstream>

fem_solid::fem_solid()
{
}

fem_solid::~fem_solid()
{
}

void fem_solid::set_lattice(double ox_,double oy_,double oz_,double hx_,double hy_,double hz_)
{
    ox=ox_; oy=oy_; oz=oz_;
    hx=hx_; hy=hy_; hz=hz_;
}

int fem_solid::add_material(const material& mt)
{
    material q = mt;
    q.lambda = q.E*q.nu/((1.0+q.nu)*(1.0-2.0*q.nu));
    q.mu = q.E/(2.0*(1.0+q.nu));
    q.cp = std::sqrt((q.lambda+2.0*q.mu)/q.rho);

    for(size_t n=0; n<mats.size(); ++n)
    if(mats[n].id==q.id)
    {
        mats[n] = q;
        return (int)n;
    }

    mats.push_back(q);
    return (int)mats.size()-1;
}

void fem_solid::add_box(double x0,double x1,double y0,double y1,double z0,double z1,int matid)
{
    shape_box b = {x0,x1,y0,y1,z0,z1,matid,false};
    boxes.push_back(b);
    shapes.push_back({0,(int)boxes.size()-1});
}

void fem_solid::remove_box(double x0,double x1,double y0,double y1,double z0,double z1)
{
    shape_box b = {x0,x1,y0,y1,z0,z1,-1,true};
    boxes.push_back(b);
    shapes.push_back({0,(int)boxes.size()-1});
}

void fem_solid::add_fix(double x0,double x1,double y0,double y1,double z0,double z1,bool fx,bool fy,bool fz)
{
    fix_box f = {x0,x1,y0,y1,z0,z1,{fx,fy,fz}};
    fixes.push_back(f);
}

// ----------------------------------------------------------------------
// setup
// ----------------------------------------------------------------------

void fem_solid::build()
{
    if(mats.empty())
    throw std::runtime_error("FEM: no material defined");

    if(hx<=0.0 || hy<=0.0 || hz<=0.0)
    throw std::runtime_error("FEM: lattice spacing missing or not positive (keyword 'lattice')");

    voxelise();
    make_nodes_elements();

    if(elems.empty())
    throw std::runtime_error("FEM: the geometry does not contain any voxel");

    if(snap_on)
    snap_surface();

    shape_derivatives();
    compute_mass();

    // supports
    fixed.assign(nnode(),0);
    make_fixes();

    // monitors: nearest node, or the top of the structure
    mons.clear();
    for(auto& mp : mon_pts)
    {
        int best = 0;
        double dbest = 1.0e300;
        for(int i=0; i<nnode(); ++i)
        {
            double d = (X[i]-mp.second).squaredNorm();
            if(d<dbest) {dbest = d; best = i;}
        }
        mons.push_back({mp.first,best});
    }
    if(!mon_autos.empty())
    {
        Vec3 c = Vec3::Zero();
        double ms = 0.0;
        for(int i=0; i<nnode(); ++i) {c += m[i]*X[i]; ms += m[i];}
        c /= ms;
        double zmax = -1.0e300;
        for(int i=0; i<nnode(); ++i) zmax = std::max(zmax,X[i](2));
        for(const mon_auto& ma : mon_autos)
        {
            int best = 0;
            double dbest = 1.0e300;
            for(int i=0; i<nnode(); ++i)
            if(X[i](2)>zmax-0.01*hz)
            {
                const double d = std::pow(X[i](0)-c(0),2) + std::pow(X[i](1)-c(1),2);
                if(d<dbest) {dbest = d; best = i;}
            }
            mons.push_back({ma.name,best});
        }
    }

    // rigid bodies (free bodies of rigid materials)
    build_surface();
    count_bodies();
    setup_rigid();

    // critical time step of the explicit scheme (hex8, lumped mass); rigid
    // bodies have no internal forces and only need the contact time step
    dtcrit = 1.0e30;
    cp_contact = 0.0;
    double hmax = 0.0;
    for(const element& e : elems)
    {
        const egeom& G = geom(e);
        hmax = std::max(hmax,G.h);
        if(e.rigid)
        continue;
        const double c = mats[e.mat].cp;
        dtcrit = std::min(dtcrit, G.L/c);
        cp_contact = std::max(cp_contact,c);
    }

    // contact wave speed per node: the stiffest material at the node (deformable parts)
    cnode.assign(nnode(),0.0);
    for(const element& e : elems)
    if(!e.rigid)
    for(int a=0; a<8; ++a)
    cnode[e.n[a]] = std::max(cnode[e.n[a]],mats[e.mat].cp);

    // rigid bodies: the impact (duration pi sqrt(M/k)) is resolved with at least
    // 20 time steps when a body can touch a wall, the ground, the bed or another
    // rigid body within the step (against deformable parts the contact is limited
    // to the penalty of their surface nodes, stable with the element time step)
    dtcrit_rig = 1.0e30;
    if(!rbs.empty())
    {
        double Mmin = 1.0e300, kmax = 0.0;
        for(const rigid_body& rb : rbs)
        {
            Mmin = std::min(Mmin,rb.M);
            kmax = std::max(kmax,rb.k);
        }
        if(kmax>0.0)
        dtcrit_rig = 3.14159265358979*std::sqrt(Mmin/kmax)/20.0;
    }
    dtcrit_el = dtcrit;
    dtcrit = std::min(dtcrit_el,dtcrit_rig);

    // crack band check: the softening branch must not snap back
    for(const material& mt : mats)
    if(mt.type==MAT_CONCRETE)
    {
        double e0 = mt.ft/mt.E;
        if(mt.ft>0.0 && mt.Gf/(mt.ft*hmax) <= 0.5*e0)
        throw std::runtime_error("FEM: concrete material "+std::to_string(mt.id)+": elements too large for the crack band (h > 2 E Gf / ft^2), use a finer resolution or increase Gf");
        double e0c = mt.fc/mt.E;
        if(mt.fc>0.0 && mt.Gc/(mt.fc*hmax) <= 0.5*e0c)
        throw std::runtime_error("FEM: concrete material "+std::to_string(mt.id)+": elements too large for the compressive crack band (h > 2 E Gc / fc^2), use a finer resolution");
    }

    x = X;
    v.assign(nnode(),Vec3::Zero());
    fint.assign(nnode(),Vec3::Zero());
    fext.assign(nnode(),Vec3::Zero());
    fcon.assign(nnode(),Vec3::Zero());
    Mcpl.assign(nnode(),Vec3::Zero());
    mfl.assign(nnode(),0.0);
    Mu_cpl.assign(nnode(),Vec3::Zero());
    fcpl.assign(nnode(),Vec3::Zero());
    vbar.assign(nnode(),Vec3::Zero());

    build_surface();
    count_bodies();

    t = 0.0;
    built = true;
}

void fem_solid::compute_mass()
{
    m.assign(nnode(),0.0);
    vnode.assign(nnode(),0.0);
    Vel = hx*hy*hz;
    helem = std::cbrt(Vel);

    // lumped mass: integral of the shape functions (row sum)
    for(const element& e : elems)
    {
        const egeom& G = geom(e);
        for(int a=0; a<8; ++a)
        {
            m[e.n[a]] += mats[e.mat].rho*G.mw[a];
            vnode[e.n[a]] += G.mw[a];
        }
    }
}

// ----------------------------------------------------------------------
// time stepping
// ----------------------------------------------------------------------

void fem_solid::clear_loads()
{
    std::fill(fext.begin(),fext.end(),Vec3::Zero());
}

void fem_solid::clear_coupling()
{
    std::fill(Mcpl.begin(),Mcpl.end(),Vec3::Zero());
    std::fill(mfl.begin(),mfl.end(),0.0);
    std::fill(Mu_cpl.begin(),Mu_cpl.end(),Vec3::Zero());
}

double fem_solid::extent() const
{
    Vec3 lo = Vec3::Constant(1.0e300), hi = Vec3::Constant(-1.0e300);
    for(const Vec3& p : X)
    {
        lo = lo.cwiseMin(p);
        hi = hi.cwiseMax(p);
    }
    return X.empty() ? 0.0 : (hi-lo).maxCoeff();
}

void fem_solid::advance(double dt)
{
    if(!built)
    throw std::runtime_error("FEM: advance() before build()");

    if(dt<=0.0)
    return;
    dt_last = dt;

    // effective masses: the enclosed fluid reduces the mass that is accelerated
    // by the structural forces (at most to 10 %), the time step follows
    std::vector<double> mr(nnode());
    double rmin = 1.0;
    for(int i=0; i<nnode(); ++i)
    {
        mr[i] = std::max(m[i]-mfl[i], 0.1*m[i]);
        if(m[i]>0.0)
        rmin = std::min(rmin, mr[i]/m[i]);
    }

    double dts_max = cfl*dtcrit_el*std::sqrt(rmin);
    if(dtcrit_rig<1.0e30 && rigid_contact_near(dt))
    dts_max = std::min(dts_max,cfl*dtcrit_rig);
    if(dts_max>=1.0e29)
    dts_max = dt;           // rigid bodies only, free motion
    const int nsub = std::max(1,(int)std::ceil(dt/dts_max - 1.0e-9));
    const double dts = dt/double(nsub);
    dts_cur = dts;
    nsub_last = nsub;

    // attached fluid parcels: their momentum exchange with the nodes over the
    // step, M (u_f - v_0), is applied as a constant force over the substeps,
    // and the parcels move with the nodes (mass m - m_f + M). Over one fluid
    // step this equals an inelastic merge, without kicking the structure.
    std::vector<double> mt(nnode());
    std::vector<Vec3> fpar(nnode());
    for(int i=0; i<nnode(); ++i)
    {
        mt[i] = mr[i] + Mcpl[i].maxCoeff();
        // relative to the mean velocity of the last step: vibrations faster
        // than the fluid step are not resolved by the flow and must not be
        // driven by the coupling
        fpar[i] = (Mu_cpl[i] - Mcpl[i].cwiseProduct(vbar[i]))/(alpha_cpl*dt);
        if(plane_strain)
        fpar[i](1) = 0.0;
    }

    for(rigid_body& rb : rbs)
    {
        rb.V0 = rb.V;
        rb.w0 = rb.w;
        rb.fcstep = 0.0;
        rb.Jc.setZero();
        rb.Hc.setZero();
    }

    std::vector<Vec3> vsum(nnode(),Vec3::Zero());
    std::vector<Vec3> vrig;
    std::vector<Vec3> frig(rbs.empty() ? 0 : nnode(),Vec3::Zero());

    for(int s=0; s<nsub; ++s)
    {
        std::fill(fint.begin(),fint.end(),Vec3::Zero());
        std::fill(fcon.begin(),fcon.end(),Vec3::Zero());

        internal_forces(dts);

        rpairs.clear();
        if(ground_on || !planes.empty() || bed_on)
        contact_ground();

        if(contact_on && (nbodies>1 || surf_dirty || n_eroded()>0) && bodies_near())
        contact_nodes();

        if(!rbs.empty())
        rigid_contact_forces();

        const double alpha = alpha_damp + (t<relax_time ? relax_alpha : 0.0);

        // structural damping of the deformation velocity
        if(alpha_struct>0.0)
        rigid_velocity(vrig);

        for(int i=0; i<nnode(); ++i)
        {
            if(m[i]<=0.0)
            continue;

            const Vec3 F = fext[i] + fcon[i] + fpar[i] - fint[i] + mr[i]*grav;

            if(!rnode.empty() && rnode[i]>=0)
            {
                frig[i] = F;
                continue;
            }

            v[i] += dts*(F/mt[i] - alpha*v[i]);
            if(alpha_struct>0.0 && body[i]>=0)
            v[i] -= (dts*alpha_struct)*(v[i]-vrig[i]);

            if(fixed[i])
            for(int d=0; d<3; ++d)
            if(fixed[i] & (1<<d)) v[i](d) = 0.0;

            if(plane_strain)
            v[i](1) = 0.0;

            if(!v[i].allFinite())
            {
                std::ostringstream os;
                os<<"FEM: non-finite velocity at node "<<i<<" x "<<x[i].transpose()<<" m "<<m[i]<<" m_f "<<mfl[i]
                  <<" M_cpl "<<Mcpl[i].transpose()<<" M u_f "<<Mu_cpl[i].transpose()<<" f_ext "<<fext[i].transpose()<<" f_int "<<fint[i].transpose()
                  <<" f_con "<<fcon[i].transpose()<<" intact elements "<<nalive[i];
                throw std::runtime_error(os.str());
            }

            x[i] += dts*v[i];
            vsum[i] += v[i];
        }

        if(!rbs.empty())
        {
            rigid_step(dts,frig);
            for(const rigid_body& rb : rbs)
            for(int i : rb.nodes)
            vsum[i] += v[i];
        }

        t += dts;
    }

    for(int i=0; i<nnode(); ++i)
    vbar[i] = vsum[i]/double(nsub);

    // accelerations of the rigid bodies over the fluid step without the contact
    // (added-mass stabilisation of the fluid loads)
    for(rigid_body& rb : rbs)
    {
        rb.a_prev = (rb.V-rb.V0-rb.Jc/rb.M)/dt;
        const Eigen::Matrix3d I = rb.R*rb.I0*rb.R.transpose();
        Vec3 dwc = Vec3::Zero();
        if(plane_strain)
        {
            if(I(1,1)>0.0) dwc(1) = rb.Hc(1)/I(1,1);
        }
        else
        dwc = I.ldlt().solve(rb.Hc);
        rb.al_prev = (rb.w-rb.w0-dwc)/dt;
    }

    // force of the attached fluid on the nodes over the step
    fcpl = fpar;

    if(surf_dirty)
    {
        // parts that broke off a supported structure in this step become rigid bodies
        std::vector<unsigned char> sup(nnode(),0);
        for(int i=0; i<nnode(); ++i)
        sup[i] = body[i]>=0 && body_fixed[body[i]];
        build_surface();
        count_bodies();
        if(fragments_rigid)
        make_rigid_fragments(sup);
    }
}

// ----------------------------------------------------------------------
// diagnostics
// ----------------------------------------------------------------------

int fem_solid::n_alive() const
{
    int n = 0;
    for(const element& e : elems)
    if(e.alive) ++n;
    return n;
}

double fem_solid::kinetic_energy() const
{
    double ek = 0.0;
    for(int i=0; i<nnode(); ++i)
    ek += 0.5*m[i]*v[i].squaredNorm();
    return ek;
}

double fem_solid::strain_energy() const
{
    double es = 0.0;

    for(int e=0; e<nelem(); ++e)
    {
        const element& el = elems[e];
        if(!el.alive)
        continue;

        const material& mt = mats[el.mat];

        for(int g=0; g<ngp; ++g)
        {
            Mat3 F = Mat3::Zero();
            for(int a=0; a<8; ++a)
            for(int i=0; i<3; ++i)
            for(int J=0; J<3; ++J)
            F(i,J) += x[el.n[a]](i)*(ngp==1 ? geom(el).dN0[a][J] : geom(el).dNg[g][a][J]);

            Mat3 E = 0.5*(F.transpose()*F - Mat3::Identity());
            const gpstate& st = gps[e*ngp+g];

            if(mt.type==MAT_J2)
            {
                E(0,0)-=st.Ep[0]; E(1,1)-=st.Ep[1]; E(2,2)-=st.Ep[2];
                E(0,1)-=st.Ep[3]; E(1,0)-=st.Ep[3];
                E(1,2)-=st.Ep[4]; E(2,1)-=st.Ep[4];
                E(2,0)-=st.Ep[5]; E(0,2)-=st.Ep[5];
            }

            const double tr = E.trace();
            Mat3 S = mt.lambda*tr*Mat3::Identity() + 2.0*mt.mu*E;
            es += (1.0-st.d)*0.5*(S.array()*E.array()).sum()*(ngp==1 ? geom(el).V : geom(el).wg[g]);
        }
    }
    return es;
}

fem_solid::Vec3 fem_solid::support_force() const
{
    return support_force(-1.0e300,1.0e300,-1.0e300,1.0e300,-1.0e300,1.0e300);
}

fem_solid::Vec3 fem_solid::support_force(double x0,double x1,double y0,double y1,double z0,double z1) const
{
    // force of the structure on its supports = load arriving at the fixed dofs
    Vec3 F = Vec3::Zero();
    for(int i=0; i<nnode(); ++i)
    if(fixed[i])
    {
        const Vec3& p = X[i];
        if(p(0)<x0 || p(0)>x1 || p(1)<y0 || p(1)>y1 || p(2)<z0 || p(2)>z1)
        continue;

        Vec3 r = fext[i] + fcon[i] + fcpl[i] + (m[i]-mfl[i])*grav - fint[i];
        for(int d=0; d<3; ++d)
        if(fixed[i] & (1<<d)) F(d) += r(d);
    }
    return F;
}

fem_solid::Vec3 fem_solid::total_load() const
{
    Vec3 F = Vec3::Zero();
    // fluid load: pressure/drag loads, direct-forcing reaction and buoyancy
    // (deformable structure and debris; rigid bodies are reported on their own)
    for(int i=0; i<nnode(); ++i)
    if(!is_rigid_node(i))
    F += fext[i] + fcpl[i] - mfl[i]*grav;
    return F;
}

double fem_solid::max_vonmises() const
{
    double s = 0.0;
    for(const element& e : elems)
    if(e.alive) s = std::max(s,e.svm);
    return s;
}

double fem_solid::max_displacement() const
{
    // deformable structure and debris; rigid bodies are reported on their own
    double d = 0.0;
    for(int i=0; i<nnode(); ++i)
    if(!is_rigid_node(i))
    d = std::max(d,(x[i]-X[i]).norm());
    return d;
}

// ----------------------------------------------------------------------
// rigid bodies
// ----------------------------------------------------------------------

void fem_solid::setup_rigid()
{
    rbs.clear();
    rnode.assign(nnode(),-1);
    for(element& e : elems) e.rigid = false;

    bool any = false;
    for(const material& mt : mats) if(mt.rigid) any = true;
    if(!any)
    return;

    // bodies whose intact elements are all rigid and that have no supports
    const int nb = (int)body_fixed.size();
    std::vector<int> allrigid(nb,1);
    for(const element& e : elems)
    if(e.alive)
    {
        const int b = body[e.n[0]];
        if(b>=0 && !mats[e.mat].rigid) allrigid[b] = 0;
    }
    std::vector<int> rid(nb,-1);
    for(int b=0; b<nb; ++b)
    if(allrigid[b] && !body_fixed[b])
    {
        rid[b] = (int)rbs.size();
        rbs.push_back(rigid_body());
    }
    for(element& e : elems)
    if(e.alive && body[e.n[0]]>=0 && rid[body[e.n[0]]]>=0)
    {
        e.rigid = true;
        rbs[rid[body[e.n[0]]]].Vol += geom(e).V;
    }

    for(int i=0; i<nnode(); ++i)
    if(body[i]>=0 && rid[body[i]]>=0)
    {
        rnode[i] = rid[body[i]];
        rbs[rnode[i]].nodes.push_back(i);
    }

    for(rigid_body& rb : rbs)
    {
        rb.M = 0.0;
        rb.c0.setZero();
        for(int i : rb.nodes) {rb.M += m[i]; rb.c0 += m[i]*X[i];}
        rb.c0 /= rb.M;
        rb.c = rb.c0;
        rb.I0.setZero();
        rb.r0.clear();
        for(int i : rb.nodes)
        {
            const Vec3 r = X[i]-rb.c0;
            rb.r0.push_back(r);
            rb.I0 += m[i]*(r.squaredNorm()*Eigen::Matrix3d::Identity() - r*r.transpose());
        }
        rb.R.setIdentity();
        rb.V.setZero(); rb.L.setZero(); rb.w.setZero();

        rigid_props(rb,(int)(&rb-&rbs[0]));
    }
}

void fem_solid::rigid_props(rigid_body& rb,int k)
{
    // added mass per unit fluid density, upper estimate: plate of the two
    // largest dimensions moving normal to itself, rho pi/4 L1 L2^2
    // (2D: rho pi/4 L1^2 per slice width)
    Vec3 lo = Vec3::Constant(1.0e300), hi = Vec3::Constant(-1.0e300);
    for(const Vec3& r : rb.r0) {lo = lo.cwiseMin(r); hi = hi.cwiseMax(r);}
    const Vec3 L = hi-lo;
    if(plane_strain)
    {
        const double L1 = std::max(L(0),L(2));
        rb.Aunit = 0.25*3.14159265358979*L1*L1*L(1);
        rb.Aface = L1*L(1);
    }
    else
    {
        double d[3] = {L(0),L(1),L(2)};
        std::sort(d,d+3);
        rb.Aunit = 0.25*3.14159265358979*d[2]*d[1]*d[1];
        rb.Aface = d[2]*d[1];
    }

    // impact stiffness of the debris: given, or the axial stiffness of a
    // bar along the longest dimension, k = E A / L with A = V / L (the
    // contact-stiffness approach of ASCE 7 for logs and poles); materials
    // without E ('material rigid <rho>') use E = rho c^2 (rigid_contact_speed)
    rb.Lbar = plane_strain ? std::max(L(0),L(2)) : L.maxCoeff();
    const double wslice = plane_strain ? std::max(L(1),1.0e-30) : 1.0;
    std::map<int,double> vm;
    for(const element& e : elems)
    if(e.alive && e.rigid && rnode[e.n[0]]==k)
    vm[e.mat] += geom(e).V;
    int md = -1;
    double vmax = -1.0;
    for(const auto& q : vm)
    if(q.second>vmax) {vmax = q.second; md = q.first;}
    const material& mt = mats[md];
    if(mt.kdebris>0.0)
    {
        rb.k = mt.kdebris*wslice;
        rb.ksrc = 0;
    }
    else
    {
        const double Eb = mt.Egiven ? mt.E : mt.rho*c_rigid*c_rigid;
        rb.k = Eb*(rb.Vol/rb.Lbar)/rb.Lbar;
        rb.ksrc = mt.Egiven ? 1 : 2;
    }
    rb.Fcap = mt.fcrush*wslice;
    if(rb.fragment)
    {
        // a fragment of the structure: its stiffness is the bar of its own
        // material, and it crushes at the strength of its cross-section
        // (concrete f_c A, steel f_y A) unless a crushing force is given
        rb.ksrc = 3;
        if(mt.fcrush<=0.0)
        {
            const double A = rb.Vol/rb.Lbar;
            if(mt.type==MAT_CONCRETE) rb.Fcap = mt.fc*A;
            else if(mt.type==MAT_J2) rb.Fcap = mt.sigy*A;
        }
    }
}

void fem_solid::rigid_step(double dts,const std::vector<Vec3>& F)
{
    for(rigid_body& rb : rbs)
    {
        Vec3 Ft = Vec3::Zero(), T = Vec3::Zero(), Fc = Vec3::Zero(), Tc = Vec3::Zero();
        for(int i : rb.nodes)
        {
            Ft += F[i];
            T += (x[i]-rb.c).cross(F[i]);
            Fc += fcon[i];
            Tc += (x[i]-rb.c).cross(fcon[i]);
        }
        if(plane_strain)
        {
            Ft(1) = 0.0; Fc(1) = 0.0;
            T(0) = T(2) = 0.0; Tc(0) = Tc(2) = 0.0;
        }

        // added-mass stabilisation of the explicit pressure loads: the body
        // carries the added mass A, and A times its acceleration of the last
        // fluid step is added back (exact for steady accelerations, stable
        // for light bodies as long as A is not far below the true added mass).
        // The contact (impact, bed, walls) acts on the mass of the body alone.
        const double ra = rb.A/rb.M;
        rb.V += (dts/(rb.M+rb.A))*(Ft - Fc + rb.A*rb.a_prev) + (dts/rb.M)*Fc;
        rb.Jc += dts*Fc;
        rb.Hc += dts*Tc;

        const Eigen::Matrix3d I = rb.R*rb.I0*rb.R.transpose();
        if(plane_strain)
        {
            const double Iy = I(1,1);
            if(Iy>0.0)
            rb.w(1) += dts*(T(1) - Tc(1) + ra*Iy*rb.al_prev(1))/((1.0+ra)*Iy) + dts*Tc(1)/Iy;
            rb.w(0) = rb.w(2) = 0.0;
        }
        else
        {
            const Vec3 rhs = T - Tc + ra*(I*rb.al_prev) - rb.w.cross(I*rb.w);
            const Eigen::LDLT<Eigen::Matrix3d> Id = I.ldlt();
            rb.w += dts*Id.solve(rhs)/(1.0+ra) + dts*Id.solve(Tc);
        }
        rb.L = I*rb.w;

        // rotation over the substep (Rodrigues), then re-orthonormalised
        const double ang = rb.w.norm()*dts;
        if(ang>0.0)
        {
            const Eigen::Matrix3d Q = Eigen::AngleAxisd(ang,rb.w/rb.w.norm()).toRotationMatrix();
            rb.R = Q*rb.R;
            Eigen::HouseholderQR<Eigen::Matrix3d> qr(rb.R);
            Eigen::Matrix3d Qo = qr.householderQ();
            // keep the orientation of the columns
            for(int k=0; k<3; ++k) if(Qo.col(k).dot(rb.R.col(k))<0.0) Qo.col(k) *= -1.0;
            rb.R = Qo;
        }
        rb.c += dts*rb.V;

        for(size_t q=0; q<rb.nodes.size(); ++q)
        {
            const int i = rb.nodes[q];
            const Vec3 r = rb.R*rb.r0[q];
            x[i] = rb.c + r;
            v[i] = rb.V + rb.w.cross(r);
            if(plane_strain)
            {
                v[i](1) = 0.0;
                x[i](1) = X[i](1);
            }
        }

        if(!rb.c.allFinite() || !rb.V.allFinite())
        throw std::runtime_error("FEM: non-finite rigid body motion");

        rb.vmax = std::max(rb.vmax,rb.V.norm());
        rb.dmax = std::max(rb.dmax,(rb.c-rb.c0).norm());
        rb.fcmax = std::max(rb.fcmax,Fc.norm());
        rb.fcstep = std::max(rb.fcstep,Fc.norm());
    }
}

bool fem_solid::bodies_near() const
{
    // part contact is only needed when two bodies come closer than the
    // contact distance; with debris particles or eroded elements always
    if(!orphan.empty() || n_eroded()>0 || surf_dirty)
    return true;
    const int nb = (int)body_fixed.size();
    if(nb<2)
    return false;
    std::vector<Vec3> lo(nb,Vec3::Constant(1.0e300)), hi(nb,Vec3::Constant(-1.0e300));
    for(int i=0; i<nnode(); ++i)
    if(body[i]>=0)
    {
        lo[body[i]] = lo[body[i]].cwiseMin(x[i]);
        hi[body[i]] = hi[body[i]].cwiseMax(x[i]);
    }
    const double d0 = contact_dist*hmin();
    for(int a=0; a<nb; ++a)
    for(int b=a+1; b<nb; ++b)
    if((lo[a].array()-d0 <= hi[b].array()).all() && (lo[b].array()-d0 <= hi[a].array()).all())
    return true;
    return false;
}

double fem_solid::rigid_terminal_speed(int k,double rho_f,double g) const
{
    // speed at which the drag on the largest face (cd 1) balances buoyancy minus weight
    const rigid_body& rb = rbs[k];
    if(rb.Aface<=0.0 || rho_f<=0.0)
    return 1.0;
    return std::max(1.0, std::sqrt(2.0*g*std::fabs(rho_f*rb.Vol-rb.M)/(rho_f*rb.Aface)));
}

void fem_solid::limit_rigid_speed(int k,const Vec3& uf,double vmax)
{
    // speed of rigid body k relative to the water uf after the fluid step: the
    // velocity change without the contact impulse may not raise it above
    // max(vmax, speed at the start of the step), and above vmax if it reverses
    rigid_body& rb = rbs[k];
    const Vec3 vc = rb.Jc/rb.M;
    const Vec3 vf = rb.V - vc;
    const Vec3 w0 = rb.V0 - uf, w1 = vf - uf;
    const double n1 = w1.norm();
    const double lim = w1.dot(w0)<0.0 ? vmax : std::max(vmax,w0.norm());
    if(n1<=lim)
    return;
    Vec3 dV = uf + (lim/n1)*w1 - vf;
    if(plane_strain)
    dV(1) = 0.0;
    rb.V += dV;
    rb.a_prev += dV/std::max(dt_last,1.0e-30);
    for(int i : rb.nodes)
    v[i] += dV;
    ++rb.nlimit;
}

bool fem_solid::rigid_contact_near(double dt) const
{
    // can a rigid body reach a wall, the ground, the bed or another rigid body within dt?
    const double d0 = contact_dist*hmin();
    const int nr = (int)rbs.size();
    std::vector<double> marg(nr);
    std::vector<Vec3> lo(nr,Vec3::Constant(1.0e300)), hi(nr,Vec3::Constant(-1.0e300));
    for(int k=0; k<nr; ++k)
    {
        const rigid_body& rb = rbs[k];
        double rmax = 0.0;
        for(const Vec3& r : rb.r0) rmax = std::max(rmax,r.norm());
        marg[k] = 2.0*(rb.V.norm() + rb.w.norm()*rmax)*dt + grav.norm()*dt*dt + d0 + 0.5*hmin();
        for(int i : rb.nodes)
        {
            lo[k] = lo[k].cwiseMin(x[i]);
            hi[k] = hi[k].cwiseMax(x[i]);
            if(ground_on && x[i](2)-zground<marg[k])
            return true;
            for(const cplane& pl : planes)
            if(pl.n.dot(x[i])-pl.d<marg[k])
            return true;
            if(bed_on && !bed_ok.empty() && bed_ok[i] && bed_phi[i] + bed_n[i].dot(x[i]-bed_x[i])<marg[k])
            return true;
        }
    }
    for(int a=0; a<nr; ++a)
    for(int b=a+1; b<nr; ++b)
    {
        const double g = marg[a]+marg[b];
        if((lo[a].array()-g <= hi[b].array()).all() && (lo[b].array()-g <= hi[a].array()).all())
        return true;
    }
    return false;
}

void fem_solid::make_rigid_fragments(const std::vector<unsigned char>& was_supported)
{
    // a free body (no supports) with nodes that belonged to a supported body
    // before the last failure has broken off the structure: it becomes a rigid
    // body with the momentum and angular momentum of its nodes, its current
    // (cracked, deformed) shape as reference and the bar stiffness of its
    // material for impacts. It does not crack further; its elements no longer
    // limit the time step, and its loads are the probed pressure with the
    // added-mass stabilisation like other floating debris.
    const int nb = (int)body_fixed.size();
    std::vector<int> frag(nb,0);
    for(int i=0; i<nnode(); ++i)
    if(body[i]>=0 && !body_fixed[body[i]] && was_supported[i] && (rnode.empty() || rnode[i]<0))
    frag[body[i]] = 1;

    bool any = false;
    for(int b=0; b<nb; ++b) any = any || frag[b];
    if(!any)
    return;

    if(rnode.empty())
    rnode.assign(nnode(),-1);

    std::vector<int> rid(nb,-1);
    for(int b=0; b<nb; ++b)
    if(frag[b])
    {
        rid[b] = (int)rbs.size();
        rbs.push_back(rigid_body());
        rbs.back().fragment = true;
    }
    for(element& e : elems)
    if(e.alive && body[e.n[0]]>=0 && rid[body[e.n[0]]]>=0)
    {
        e.rigid = true;
        rbs[rid[body[e.n[0]]]].Vol += geom(e).V;
    }
    for(int i=0; i<nnode(); ++i)
    if(body[i]>=0 && rid[body[i]]>=0)
    {
        rnode[i] = rid[body[i]];
        rbs[rnode[i]].nodes.push_back(i);
    }

    for(int b=0; b<nb; ++b)
    if(rid[b]>=0)
    {
        rigid_body& rb = rbs[rid[b]];
        rb.M = 0.0;
        Vec3 c = Vec3::Zero(), P = Vec3::Zero();
        for(int i : rb.nodes) {rb.M += m[i]; c += m[i]*x[i]; P += m[i]*v[i];}
        c /= rb.M;
        rb.c0 = rb.c = c;
        rb.V = P/rb.M;
        rb.I0.setZero();
        rb.r0.clear();
        Vec3 Lm = Vec3::Zero();
        for(int i : rb.nodes)
        {
            const Vec3 r = x[i]-c;
            rb.r0.push_back(r);
            rb.I0 += m[i]*(r.squaredNorm()*Eigen::Matrix3d::Identity() - r*r.transpose());
            Lm += m[i]*r.cross(v[i]-rb.V);
        }
        rb.R.setIdentity();
        rb.w.setZero();
        if(plane_strain)
        {
            rb.V(1) = 0.0;
            if(rb.I0(1,1)>0.0) rb.w(1) = Lm(1)/rb.I0(1,1);
        }
        else
        rb.w = rb.I0.ldlt().solve(Lm);
        rb.L = rb.I0*rb.w;
        rb.V0 = rb.V; rb.w0 = rb.w;
        for(size_t q=0; q<rb.nodes.size(); ++q)
        v[rb.nodes[q]] = rb.V + rb.w.cross(rb.r0[q]);
        rigid_props(rb,rid[b]);
    }
    update_time_steps();
}

void fem_solid::update_time_steps()
{
    // deformable elements (failed ones too: their nodes are debris particles whose
    // contact penalty follows the element frequency) and the impacts of the rigid bodies
    dtcrit_el = 1.0e30;
    for(const element& e : elems)
    if(!e.rigid)
    dtcrit_el = std::min(dtcrit_el, geom(e).L/mats[e.mat].cp);
    dtcrit_rig = 1.0e30;
    double Mmin = 1.0e300, kmax = 0.0;
    for(const rigid_body& rb : rbs)
    {
        Mmin = std::min(Mmin,rb.M);
        kmax = std::max(kmax,rb.k);
    }
    if(kmax>0.0)
    dtcrit_rig = 3.14159265358979*std::sqrt(Mmin/kmax)/20.0;
    dtcrit = std::min(dtcrit_el,dtcrit_rig);
}

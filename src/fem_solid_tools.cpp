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

#include"fem_solid.h"
#include<cmath>
#include<sstream>
#include<ostream>
#include<iomanip>
#include<stdexcept>
#include<algorithm>

// ----------------------------------------------------------------------
// material presets (mean values, Eurocode-type)
// ----------------------------------------------------------------------

namespace
{
    std::string upper(std::string s)
    {
        for(char& c : s) c = (char)std::toupper((unsigned char)c);
        return s;
    }
}

bool fem_solid::preset(const std::string& type,const std::string& name_in,material& mt)
{
    const std::string name = upper(name_in);
    mt.name = type+" "+name_in;

    if(type=="concrete")
    {
        // EN 1992-1-1 Table 3.1: fck, fctm, Ecm; fcm = fck + 8 MPa
        // Gf = 73 fcm^0.18 N/m (fib Model Code 2010), Gc = 8.8 sqrt(fcm) N/mm (Nakamura & Higai 2001)
        struct row {const char* n; double fck, fctm, Ecm;};
        static const row tab[] = {{"C20",20,2.2,30},{"C25",25,2.6,31},{"C30",30,2.9,33},{"C35",35,3.2,34},
                                  {"C40",40,3.5,35},{"C45",45,3.8,36},{"C50",50,4.1,37}};
        const std::string key = name.substr(0,name.find('/'));
        for(const row& r : tab)
        if(key==r.n)
        {
            const double fcm = r.fck + 8.0;
            mt.type = MAT_CONCRETE;
            mt.rho = 2400.0; mt.nu = 0.2; mt.E = r.Ecm*1.0e9;
            mt.ft = r.fctm*1.0e6; mt.fc = fcm*1.0e6;
            mt.Gf = 73.0*std::pow(fcm,0.18);
            mt.Gc = 8.8*std::sqrt(fcm)*1000.0;
            return true;
        }
        return false;
    }

    if(type=="steel")
    {
        // EN 10025, bilinear: hardening modulus 1 % of E, erosion at 15 % plastic strain
        struct row {const char* n; double fy, ef;};
        static const row tab[] = {{"S235",235,0.15},{"S275",275,0.15},{"S355",355,0.15},{"S420",420,0.12},{"S460",460,0.12}};
        for(const row& r : tab)
        if(name==r.n)
        {
            mt.type = MAT_J2;
            mt.rho = 7850.0; mt.E = 2.1e11; mt.nu = 0.3;
            mt.sigy = r.fy*1.0e6; mt.H = 0.01*mt.E; mt.epsfail = r.ef;
            return true;
        }
        return false;
    }

    if(type=="aluminium" || type=="aluminum")
    {
        if(name=="6061" || name=="6061-T6")
        {
            mt.type = MAT_J2;
            mt.rho = 2700.0; mt.E = 6.9e10; mt.nu = 0.33;
            mt.sigy = 2.76e8; mt.H = 0.01*mt.E; mt.epsfail = 0.10;
            return true;
        }
        return false;
    }

    if(type=="timber" || type=="wood")
    {
        // EN 338 / EN 14080 mean values along the grain, isotropic elastic simplification
        struct row {const char* n; double E, rho;};
        static const row tab[] = {{"C16",8.0,370},{"C24",11.0,420},{"C30",12.0,460},{"GL24H",11.5,420},{"GL28H",12.6,460}};
        for(const row& r : tab)
        if(name==r.n)
        {
            mt.type = MAT_ELASTIC;
            mt.rho = r.rho; mt.E = r.E*1.0e9; mt.nu = 0.3;
            return true;
        }
        return false;
    }

    if(type=="rock" || type=="stone")
    {
        struct row {const char* n; double rho, E;};
        static const row tab[] = {{"GRANITE",2650,50},{"BASALT",2900,70},{"LIMESTONE",2600,40},{"SANDSTONE",2300,15}};
        for(const row& r : tab)
        if(name==r.n)
        {
            mt.type = MAT_ELASTIC;
            mt.rho = r.rho; mt.E = r.E*1.0e9; mt.nu = 0.25;
            return true;
        }
        return false;
    }

    if(type=="masonry")
    {
        // rough values, calibrate for the actual masonry
        if(name=="BRICK")
        {
            mt.type = MAT_CONCRETE;
            mt.rho = 1800.0; mt.E = 5.0e9; mt.nu = 0.2;
            mt.ft = 0.3e6; mt.Gf = 20.0; mt.fc = 8.0e6; mt.Gc = 15000.0;
            return true;
        }
        return false;
    }

    if(type=="rubber")
    {
        if(name=="DEFAULT" || name=="RUBBER" || name.empty())
        {
            mt.type = MAT_ELASTIC;
            mt.rho = 1100.0; mt.E = 1.2e7; mt.nu = 0.45;
            return true;
        }
        return false;
    }

    return false;
}

std::string fem_solid::preset_list()
{
    return "concrete C20 C25 C30 C35 C40 C45 C50 | steel S235 S275 S355 S420 S460 | aluminium 6061 | "
           "timber C16 C24 C30 GL24h GL28h | rock granite basalt limestone sandstone | masonry brick | rubber";
}

// ----------------------------------------------------------------------
// lattice spacing, supports
// ----------------------------------------------------------------------

void fem_solid::set_default_spacing(double h,double hy_)
{
    if(lattice_given())
    return;
    hx = hy = hz = h;
    if(hy_>0.0)
    hy = hy_;
}

void fem_solid::make_fixes()
{
    // boxes
    for(const fix_box& f : fixes)
    if(f.mode==0)
    for(int i=0; i<nnode(); ++i)
    {
        const Vec3& p = X[i];
        const double eps = 1.0e-9*hmin();
        if(p(0)>=f.x0-eps && p(0)<=f.x1+eps && p(1)>=f.y0-eps && p(1)<=f.y1+eps && p(2)>=f.z0-eps && p(2)<=f.z1+eps)
        for(int d=0; d<3; ++d)
        if(f.f[d]) fixed[i] |= (unsigned char)(1<<d);
    }

    // base / top: lowest / highest node level of every connected body
    bool any = false;
    for(const fix_box& f : fixes) if(f.mode>0) any = true;
    if(any)
    {
        std::vector<int> comp(nnode(),-1);
        std::vector<int> stack;
        int nc = 0;
        for(int s=0; s<nnode(); ++s)
        {
            if(comp[s]>=0 || node_elem_start[s]==node_elem_start[s+1])
            continue;
            comp[s] = nc;
            stack.push_back(s);
            while(!stack.empty())
            {
                const int i = stack.back();
                stack.pop_back();
                for(int q=node_elem_start[i]; q<node_elem_start[i+1]; ++q)
                for(int a=0; a<8; ++a)
                {
                    const int j = elems[node_elem[q]].n[a];
                    if(comp[j]<0) {comp[j] = nc; stack.push_back(j);}
                }
            }
            ++nc;
        }

        std::vector<double> zlo(nc,1.0e300), zhi(nc,-1.0e300);
        for(int i=0; i<nnode(); ++i)
        if(comp[i]>=0)
        {
            zlo[comp[i]] = std::min(zlo[comp[i]],X[i](2));
            zhi[comp[i]] = std::max(zhi[comp[i]],X[i](2));
        }

        for(const fix_box& f : fixes)
        if(f.mode>0)
        for(int i=0; i<nnode(); ++i)
        if(comp[i]>=0)
        {
            const bool hit = (f.mode==1) ? X[i](2)<=zlo[comp[i]]+0.01*hz : X[i](2)>=zhi[comp[i]]-0.01*hz;
            if(hit)
            for(int d=0; d<3; ++d)
            if(f.f[d]) fixed[i] |= (unsigned char)(1<<d);
        }
    }

    // base centre: centroid of the supports at their lowest level
    base_c.setZero();
    int n = 0;
    double zmin = 1.0e300;
    for(int i=0; i<nnode(); ++i)
    if(fixed[i])
    {
        base_c += X[i];
        zmin = std::min(zmin,X[i](2));
        ++n;
    }
    if(n>0)
    {
        base_c /= double(n);
        base_c(2) = zmin;
    }
}

void fem_solid::fix_nodes(const std::vector<int>& nodes)
{
    for(int i : nodes)
    if(i>=0 && i<nnode())
    fixed[i] = 7;

    const std::vector<fix_box> keep = fixes;
    fixes.clear();
    make_fixes();       // only recomputes the base centre (no boxes)
    fixes = keep;
}

bool fem_solid::has_supports() const
{
    for(unsigned char f : fixed)
    if(f) return true;
    return false;
}

double fem_solid::min_density() const
{
    double r = 1.0e300;
    for(const element& e : elems)
    r = std::min(r,mats[e.mat].rho);
    return r;
}

// ----------------------------------------------------------------------
// gravity settling: kinetic damping (velocities reset at every peak of
// the kinetic energy) converges to the static equilibrium
// ----------------------------------------------------------------------

bool fem_solid::settle(int maxsteps,double tol,double* residual,bool supported_only)
{
    if(!built)
    throw std::runtime_error("FEM: settle() before build()");

    const double dts = cfl*dtcrit;
    double ke_prev = 0.0;
    double wref = 0.0;
    for(int i=0; i<nnode(); ++i)
    wref += m[i]*grav.norm() + fext[i].norm();
    if(wref<=0.0)
    wref = 1.0;

    std::fill(v.begin(),v.end(),Vec3::Zero());
    double res = 1.0e30;
    bool conv = false;

    for(int n=0; n<maxsteps; ++n)
    {
        std::fill(fint.begin(),fint.end(),Vec3::Zero());
        std::fill(fcon.begin(),fcon.end(),Vec3::Zero());
        internal_forces(dts);
        if(ground_on)
        contact_ground();
        if(contact_on && (nbodies>1 || surf_dirty || n_eroded()>0))
        contact_nodes();

        double ke = 0.0, r2 = 0.0;
        for(int i=0; i<nnode(); ++i)
        {
            if(m[i]<=0.0)
            continue;
            // free parts without ground would fall for ever: kept in place
            if(supported_only && (body[i]<0 || !body_fixed[body[i]]))
            continue;
            Vec3 F = fext[i] + fcon[i] - fint[i] + m[i]*grav;
            for(int d=0; d<3; ++d)
            if(fixed[i] & (1<<d)) F(d) = 0.0;
            if(plane_strain) F(1) = 0.0;
            r2 += F.squaredNorm();

            v[i] += dts*F/m[i];
            for(int d=0; d<3; ++d)
            if(fixed[i] & (1<<d)) v[i](d) = 0.0;
            if(plane_strain) v[i](1) = 0.0;
            x[i] += dts*v[i];
            ke += 0.5*m[i]*v[i].squaredNorm();
        }

        if(ke<ke_prev)
        {
            std::fill(v.begin(),v.end(),Vec3::Zero());
            ke_prev = 0.0;
        }
        else
        ke_prev = ke;

        res = std::sqrt(r2)/wref;
        if(n>10 && res<tol)
        {
            conv = true;
            break;
        }
    }

    std::fill(v.begin(),v.end(),Vec3::Zero());
    std::fill(vbar.begin(),vbar.end(),Vec3::Zero());
    if(surf_dirty)
    {
        build_surface();
        count_bodies();
    }
    if(residual)
    *residual = res;
    return conv;
}

// ----------------------------------------------------------------------
// utilisation, diagnostics
// ----------------------------------------------------------------------

void fem_solid::update_utilisation()
{
    const Mat3 I = Mat3::Identity();

    for(int e=0; e<nelem(); ++e)
    {
        element& el = elems[e];
        el.util = -1.0;
        if(!el.alive)
        continue;

        const material& mt = mats[el.mat];
        if(mt.type==MAT_ELASTIC)
        continue;

        const egeom& G = geom(el);
        double u = 0.0;

        for(int g=0; g<ngp; ++g)
        {
            const double (*dN)[3] = (ngp==1) ? G.dN0 : G.dNg[g];
            Mat3 F = Mat3::Zero();
            for(int a=0; a<8; ++a)
            for(int i=0; i<3; ++i)
            for(int J=0; J<3; ++J)
            F(i,J) += x[el.n[a]](i)*dN[a][J];

            const double Jd = F.determinant();
            if(Jd<=0.0)
            continue;

            Mat3 E = 0.5*(F.transpose()*F - I);
            const gpstate& st = gps[e*ngp+g];
            Mat3 S;

            if(mt.type==MAT_J2)
            {
                Mat3 Ep;
                Ep << st.Ep[0], st.Ep[3], st.Ep[5],
                      st.Ep[3], st.Ep[1], st.Ep[4],
                      st.Ep[5], st.Ep[4], st.Ep[2];
                const Mat3 Ee = E - Ep;
                S = mt.lambda*Ee.trace()*I + 2.0*mt.mu*Ee;
            }
            else
            S = (1.0-st.d)*(mt.lambda*E.trace()*I + 2.0*mt.mu*E);

            const Mat3 sig = F*S*F.transpose()/Jd;

            if(mt.type==MAT_J2)
            {
                const Mat3 dev = sig - sig.trace()/3.0*I;
                const double svm = std::sqrt(1.5*(dev.array()*dev.array()).sum());
                if(mt.sigy>0.0)
                u = std::max(u, svm/mt.sigy);
            }
            else
            {
                Eigen::SelfAdjointEigenSolver<Mat3> es;
                es.computeDirect(sig,Eigen::EigenvaluesOnly);
                const Eigen::Vector3d sp = es.eigenvalues();
                if(mt.ft>0.0) u = std::max(u, sp(2)/mt.ft);
                if(mt.fc>0.0) u = std::max(u, -sp(0)/mt.fc);
            }
        }
        el.util = u;
    }
}

double fem_solid::max_utilisation(int* elem) const
{
    double u = -1.0;
    int best = -1;
    for(int e=0; e<nelem(); ++e)
    if(elems[e].alive && elems[e].util>u)
    {
        u = elems[e].util;
        best = e;
    }
    if(elem) *elem = best;
    return u;
}

double fem_solid::max_damage(int* elem) const
{
    double d = 0.0;
    int best = -1;
    for(int e=0; e<nelem(); ++e)
    for(int g=0; g<ngp; ++g)
    if(gps[e*ngp+g].d>d)
    {
        d = gps[e*ngp+g].d;
        best = e;
    }
    if(elem) *elem = best;
    return d;
}

double fem_solid::max_plastic_strain(int* elem) const
{
    double p = 0.0;
    int best = -1;
    for(int e=0; e<nelem(); ++e)
    for(int g=0; g<ngp; ++g)
    if(gps[e*ngp+g].ep>p)
    {
        p = gps[e*ngp+g].ep;
        best = e;
    }
    if(elem) *elem = best;
    return p;
}

fem_solid::Vec3 fem_solid::elem_centre(int e) const
{
    Vec3 c = Vec3::Zero();
    for(int a=0; a<8; ++a)
    c += x[elems[e].n[a]];
    return c/8.0;
}

fem_solid::Vec3 fem_solid::support_moment() const
{
    Vec3 M = Vec3::Zero();
    for(int i=0; i<nnode(); ++i)
    if(fixed[i])
    {
        Vec3 r = fext[i] + fcon[i] + fcpl[i] + (m[i]-mfl[i])*grav - fint[i];
        for(int d=0; d<3; ++d)
        if(!(fixed[i] & (1<<d))) r(d) = 0.0;
        M += (X[i]-base_c).cross(r);
    }
    return M;
}

double fem_solid::eroded_mass_fraction() const
{
    double me = 0.0, mt = 0.0;
    for(const element& e : elems)
    {
        const double mm = mats[e.mat].rho*geom(e).V;
        mt += mm;
        if(!e.alive) me += mm;
    }
    return mt>0.0 ? me/mt : 0.0;
}

// ----------------------------------------------------------------------
// check mode
// ----------------------------------------------------------------------

fem_solid::check_info fem_solid::check()
{
    check_info ci;

    ci.lo = Vec3::Constant(1.0e300);
    ci.hi = Vec3::Constant(-1.0e300);
    for(int i=0; i<nnode(); ++i)
    {
        ci.mass += m[i];
        ci.cog += m[i]*X[i];
        ci.lo = ci.lo.cwiseMin(X[i]);
        ci.hi = ci.hi.cwiseMax(X[i]);
        if(fixed[i]) ++ci.nfixed;
    }
    ci.cog /= ci.mass;
    for(const element& e : elems)
    ci.volume += geom(e).V;

    for(const face& f : faces)
    ci.surface_area += 0.5*((X[f.n[2]]-X[f.n[0]]).cross(X[f.n[3]]-X[f.n[1]])).norm();

    // one element thick parts: exposed faces on opposite sides
    {
        static const int off[6][3] = {{-1,0,0},{1,0,0},{0,-1,0},{0,1,0},{0,0,-1},{0,0,1}};
        int nsurf = 0, nthin = 0;
        for(const element& el : elems)
        {
            bool ex[6];
            bool any = false;
            for(int f=0; f<6; ++f)
            {
                const int nb = voxel(el.ix+off[f][0],el.iy+off[f][1],el.iz+off[f][2]);
                ex[f] = (nb<0 || vox_elem[nb]<0);
                any = any || ex[f];
            }
            if(!any) continue;
            ++nsurf;
            const bool thinx = ex[0] && ex[1], thiny = !plane_strain && ex[2] && ex[3], thinz = ex[4] && ex[5];
            if(thinx || thiny || thinz) ++nthin;
        }
        ci.thin_fraction = nsurf>0 ? double(nthin)/double(nsurf) : 0.0;
    }

    // warnings that need no simulation
    if(ci.thin_fraction>0.1)
    {
        std::ostringstream os;
        os<<"about "<<int(100.0*ci.thin_fraction+0.5)<<" % of the surface elements belong to parts that are only one element thick: "
            "bending stiffness and stresses there are inaccurate, use 'resolution fine' or a smaller lattice spacing";
        ci.warnings.push_back(os.str());
    }
    if(ci.nfixed==0 && !ground_on)
    ci.warnings.push_back("the structure has no supports ('fix base', 'fix bed' or a fix box) and no ground: it will fall or drift");
    for(const material& mt : mats)
    if(mt.type==MAT_CONCRETE && mt.ft>0.0)
    {
        double hmax = 0.0;
        for(const element& e : elems) if(e.mat>=0 && &mats[e.mat]==&mt) hmax = std::max(hmax,geom(e).h);
        const double hlim = 2.0*mt.E*mt.Gf/(mt.ft*mt.ft);
        if(hmax>0.5*hlim)
        {
            std::ostringstream os;
            os<<"material "<<mt.id<<" ("<<mt.name<<"): elements of "<<hmax<<" m are more than half of the crack-band limit "<<hlim<<" m, cracking will be brittle";
            ci.warnings.push_back(os.str());
        }
    }
    if(nelem()>200000)
    ci.warnings.push_back("more than 200000 elements: the solid is computed on every rank, expect it to dominate the run time");

    // keep the initial state
    const std::vector<Vec3> x0 = x;
    const std::vector<gpstate> g0 = gps;
    const std::vector<element> e0 = elems;
    const Vec3 grav0 = grav;
    const double wdiss0 = wdiss;

    // self weight
    if(ci.nfixed>0 || ground_on)
    {
        double res;
        ci.settled = settle(200000,1.0e-4,&res,!ground_on);
        std::fill(fint.begin(),fint.end(),Vec3::Zero());
        std::fill(fcon.begin(),fcon.end(),Vec3::Zero());
        internal_forces(cfl*dtcrit);
        ci.sw_maxdisp = max_displacement();
        ci.sw_maxvm = max_vonmises();
        update_utilisation();
        ci.sw_maxutil = max_utilisation();
        ci.sw_support = support_force();
        if(!ci.settled)
        ci.warnings.push_back("the structure did not reach equilibrium under its own weight within 200000 steps (very soft or not stable)");
        if(n_eroded()>0)
        ci.warnings.push_back("elements fail under the self weight alone");

        x = x0; gps = g0; elems = e0; wdiss = wdiss0;
        build_surface();
        count_bodies();
    }

    // first natural frequencies (Rayleigh, static deflection under 1 g in each direction)
    if(ci.nfixed>0)
    rayleigh_frequencies(ci.freq,ci.freq_ok);
    grav = grav0;
    std::fill(v.begin(),v.end(),Vec3::Zero());

    return ci;
}

void fem_solid::write_check(std::ostream& os,const check_info& ci) const
{
    os<<std::setprecision(4);
    os<<"FEM check\n";
    os<<"  mesh:        "<<nelem()<<" hex8 elements ("<<(ngp==8 ? "full" : "reduced")<<" integration), "<<nnode()<<" nodes, "
      <<geos.size()<<" elements snapped to the surface, element size "<<hx<<" x "<<hy<<" x "<<hz<<" m\n";
    os<<"  extent:      x "<<ci.lo(0)<<" .. "<<ci.hi(0)<<"   y "<<ci.lo(1)<<" .. "<<ci.hi(1)<<"   z "<<ci.lo(2)<<" .. "<<ci.hi(2)<<" m\n";
    os<<"  volume:      "<<ci.volume<<" m3,  mass "<<ci.mass<<" kg,  weight "<<ci.mass*grav.norm()/1000.0<<" kN\n";
    os<<"  centre of gravity: "<<ci.cog(0)<<" "<<ci.cog(1)<<" "<<ci.cog(2)<<" m\n";
    os<<"  surface:     "<<ci.surface_area<<" m2\n";
    for(const material& mt : mats)
    {
        double vol = 0.0;
        for(const element& e : elems) if(&mats[e.mat]==&mt) vol += geom(e).V;
        os<<"  material "<<mt.id<<":  "<<(mt.name.empty() ? std::string("custom") : mt.name)<<", "<<vol<<" m3, rho "<<mt.rho<<" kg/m3, E "<<mt.E/1.0e9<<" GPa";
        if(mt.type==MAT_J2) os<<", yield "<<mt.sigy/1.0e6<<" MPa";
        if(mt.type==MAT_CONCRETE) os<<", tension "<<mt.ft/1.0e6<<" MPa, compression "<<mt.fc/1.0e6<<" MPa";
        os<<"\n";
    }
    os<<"  supports:    "<<ci.nfixed<<" nodes";
    if(ci.nfixed>0) os<<", base centre "<<base_c(0)<<" "<<base_c(1)<<" "<<base_c(2)<<" m";
    os<<(ground_on ? ", ground contact on" : "")<<"\n";
    os<<"  time step:   "<<dtcrit*cfl<<" s (solid)\n";
    if(ci.nfixed>0 || ground_on)
    {
        os<<"  self weight: "<<(ci.settled ? "in equilibrium" : "NOT in equilibrium")<<", max displacement "<<ci.sw_maxdisp*1000.0<<" mm, max von Mises "<<ci.sw_maxvm/1.0e6<<" MPa";
        if(ci.sw_maxutil>=0.0) os<<", max utilisation "<<ci.sw_maxutil;
        os<<"\n               support force "<<ci.sw_support.transpose()/1000.0<<" kN (weight "<<ci.mass*grav.norm()/1000.0<<" kN)\n";
    }
    const char* dn[3] = {"x","y","z"};
    for(int d=0; d<3; ++d)
    if(ci.freq_ok[d])
    os<<"  first natural frequency in "<<dn[d]<<" (Rayleigh, dry): "<<ci.freq[d]<<" Hz,  period "<<1.0/ci.freq[d]<<" s\n";
    {
        double f1 = 0.0;
        for(int d=0; d<3; ++d) if(ci.freq_ok[d] && (f1==0.0 || ci.freq[d]<f1)) f1 = ci.freq[d];
        if(zeta>0.0 && f1>0.0)
        os<<"  structural damping: "<<100.0*zeta<<" % of critical at "<<f1<<" Hz (less for higher modes, rigid motion undamped)\n";
        else if(zeta>0.0)
        os<<"  structural damping: "<<100.0*zeta<<" % requested, but no natural frequency (no supports): not applied\n";
        else
        os<<"  structural damping: off\n";
    }
    if(ci.warnings.empty())
    os<<"  no warnings\n";
    for(const std::string& w : ci.warnings)
    os<<"  WARNING: "<<w<<"\n";
}

// ----------------------------------------------------------------------
// first natural frequencies: Rayleigh quotient of the static deflection
// under 1 g in each direction (supports needed), state restored
// ----------------------------------------------------------------------

void fem_solid::rayleigh_frequencies(double f[3],bool ok[3])
{
    const std::vector<Vec3> x0 = x, v0 = v, fe0 = fext;
    const std::vector<gpstate> g0 = gps;
    const std::vector<element> e0 = elems;
    const Vec3 grav0 = grav;
    const double wdiss0 = wdiss;
    const int ne0 = n_eroded();
    const std::vector<material> m0 = mats;

    // linear-elastic stiffness: no cracking, yielding or erosion under the test load
    std::fill(fext.begin(),fext.end(),Vec3::Zero());
    for(material& mt : mats)
    mt.type = MAT_ELASTIC;

    for(int d=0; d<3; ++d)
    {
        f[d] = 0.0;
        ok[d] = false;
        if(plane_strain && d==1)
        continue;

        // from the undeformed state: deflection under 1 g only
        x = X;
        grav = Vec3::Zero();
        grav(d) = 9.81;
        const bool conv = settle(200000,1.0e-5,nullptr,true);
        double num = 0.0, den = 0.0;
        for(int i=0; i<nnode(); ++i)
        {
            const Vec3 u = x[i]-X[i];
            num += m[i]*9.81*u(d);
            den += m[i]*u.squaredNorm();
        }
        if(conv && den>0.0 && num>0.0 && n_eroded()==ne0)
        {
            f[d] = std::sqrt(num/den)/(2.0*3.14159265358979);
            ok[d] = true;
        }
        x = x0; gps = g0; elems = e0; wdiss = wdiss0;
        build_surface();
        count_bodies();
    }
    grav = grav0;
    v = v0;
    fext = fe0;
    mats = m0;
}

void fem_solid::prepare_damping()
{
    alpha_struct = 0.0;
    f_damp = 0.0;
    if(zeta<=0.0 || !has_supports())
    return;

    double f[3];
    bool ok[3];
    rayleigh_frequencies(f,ok);
    for(int d=0; d<3; ++d)
    if(ok[d] && (f_damp==0.0 || f[d]<f_damp))
    f_damp = f[d];

    // mass-proportional on the deformation: ratio zeta at the first frequency,
    // less for the higher modes (bulk viscosity takes care of those)
    if(f_damp>0.0)
    alpha_struct = 2.0*zeta*2.0*3.14159265358979*f_damp;
}

void fem_solid::rigid_velocity(std::vector<Vec3>& vr) const
{
    // supported bodies: rest, free bodies: translation + rotation fitted to
    // their momentum and angular momentum
    vr.assign(nnode(),Vec3::Zero());
    const int nb = (int)body_fixed.size();
    std::vector<double> M(nb,0.0);
    std::vector<Vec3> P(nb,Vec3::Zero()), C(nb,Vec3::Zero());
    for(int i=0; i<nnode(); ++i)
    {
        const int b = body[i];
        if(b<0 || body_fixed[b])
        continue;
        M[b] += m[i];
        P[b] += m[i]*v[i];
        C[b] += m[i]*x[i];
    }
    std::vector<Vec3> L(nb,Vec3::Zero());
    std::vector<Eigen::Matrix3d> J(nb,Eigen::Matrix3d::Zero());
    for(int b=0; b<nb; ++b)
    if(M[b]>0.0)
    C[b] /= M[b];
    for(int i=0; i<nnode(); ++i)
    {
        const int b = body[i];
        if(b<0 || body_fixed[b] || M[b]<=0.0)
        continue;
        const Vec3 r = x[i]-C[b];
        L[b] += m[i]*r.cross(v[i]);
        J[b] += m[i]*(r.squaredNorm()*Eigen::Matrix3d::Identity() - r*r.transpose());
    }
    std::vector<Vec3> W(nb,Vec3::Zero());
    for(int b=0; b<nb; ++b)
    if(M[b]>0.0 && !body_fixed[b])
    {
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es(J[b]);
        const Eigen::Vector3d lam = es.eigenvalues();
        const double lmax = lam.maxCoeff();
        Eigen::Vector3d inv;
        for(int k=0; k<3; ++k)
        inv(k) = lam(k)>1.0e-9*lmax ? 1.0/lam(k) : 0.0;
        W[b] = es.eigenvectors()*inv.asDiagonal()*es.eigenvectors().transpose()*L[b];
        if(plane_strain)
        {
            W[b](0) = 0.0;
            W[b](2) = 0.0;
        }
    }
    for(int i=0; i<nnode(); ++i)
    {
        const int b = body[i];
        if(b<0 || body_fixed[b] || M[b]<=0.0)
        continue;
        vr[i] = P[b]/M[b] + W[b].cross(x[i]-C[b]);
    }
}

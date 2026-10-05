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
#include<sstream>
#include<istream>
#include<stdexcept>
#include<cstdlib>
#include<vector>
#include<algorithm>

// fem.dat: one keyword per line, '#' starts a comment, SI units.
// See FEM.md (FEM documentation) for the full list.

void fem_solid::read(std::istream& is)
{
    std::string line;
    int lineno = 0;

    auto fail = [&](const std::string& msg)
    {
        throw std::runtime_error("FEM: fem.dat line "+std::to_string(lineno)+": "+msg);
    };

    while(std::getline(is,line))
    {
        ++lineno;

        const size_t c = line.find('#');
        if(c!=std::string::npos)
        line.erase(c);

        std::istringstream ls(line);
        std::string kw;
        if(!(ls>>kw))
        continue;

        auto need = [&](bool ok)
        {
            if(!ok) fail("missing or invalid values for '"+kw+"'");
        };

        if(kw=="lattice")
        {
            need(static_cast<bool>(ls>>hx>>hy>>hz));
            double a,b,cc;
            if(ls>>a>>b>>cc) {ox=a; oy=b; oz=cc;}
        }
        else if(kw=="material")
        {
            // material [id] <type> <preset> [key value ...]
            // material [id] <type> <numbers ...>          (elastic rho E nu | plastic rho E nu sigy H epsfail |
            //                                                concrete rho E nu ft Gf fc Gc [derode])
            std::vector<std::string> tok;
            std::string w;
            while(ls>>w) tok.push_back(w);
            if(tok.empty()) fail("material: type missing");

            material mt;
            size_t q = 0;
            bool hasid = false;
            {
                char* end = nullptr;
                const long id = std::strtol(tok[0].c_str(),&end,10);
                if(end && *end=='\0') {mt.id = (int)id; hasid = true; q = 1;}
            }
            if(!hasid)
            {
                int idmax = 0;
                for(const material& m0 : mats) idmax = std::max(idmax,m0.id);
                mt.id = idmax+1;
            }
            if(q>=tok.size()) fail("material: type missing");
            std::string type = tok[q++];
            if(type=="plastic" || type=="j2") type = "steel_custom";

            auto isnum = [](const std::string& t){char* e=nullptr; std::strtod(t.c_str(),&e); return e && *e=='\0';};

            if(q<tok.size() && !isnum(tok[q]))
            {
                // preset
                if(!preset(type,tok[q],mt))
                fail("unknown material preset '"+type+" "+tok[q]+"', available: "+preset_list());
                ++q;
            }
            else if(q==tok.size() && type=="rubber")
            preset("rubber","default",mt);
            else if(type=="rigid")
            {
                // material rigid <rho>: rigid body of the given density (debris, floating objects)
                if(q>=tok.size() || !isnum(tok[q])) fail("material rigid: density missing (material rigid 500)");
                mt.type = MAT_ELASTIC; mt.rho = std::strtod(tok[q++].c_str(),nullptr); mt.E = 1.0e9; mt.nu = 0.25;
                mt.rigid = true; mt.name = "rigid"; mt.Egiven = false;
            }
            else
            {
                std::vector<double> num;
                while(q<tok.size() && isnum(tok[q])) num.push_back(std::strtod(tok[q++].c_str(),nullptr));
                mt.name = type;
                if(type=="elastic" && num.size()>=3)
                {mt.type = MAT_ELASTIC; mt.rho=num[0]; mt.E=num[1]; mt.nu=num[2];}
                else if((type=="steel_custom" || type=="steel") && num.size()>=6)
                {mt.type = MAT_J2; mt.rho=num[0]; mt.E=num[1]; mt.nu=num[2]; mt.sigy=num[3]; mt.H=num[4]; mt.epsfail=num[5]; mt.name="plastic";}
                else if(type=="concrete" && num.size()>=7)
                {mt.type = MAT_CONCRETE; mt.rho=num[0]; mt.E=num[1]; mt.nu=num[2]; mt.ft=num[3]; mt.Gf=num[4]; mt.fc=num[5]; mt.Gc=num[6]; if(num.size()>=8) mt.derode=num[7];}
                else
                fail("material: give a preset (e.g. 'material concrete C30') or the numbers of the type, presets: "+preset_list());
            }

            // key value overrides
            while(q<tok.size())
            {
                const std::string k = tok[q];
                if(k=="rigid") {mt.rigid = true; ++q; continue;}
                if(q+1>=tok.size()) break;
                if(!isnum(tok[q+1])) fail("material: value missing for '"+k+"'");
                const double val = std::strtod(tok[q+1].c_str(),nullptr);
                if(k=="rho") mt.rho = val;
                else if(k=="E") {mt.E = val; mt.Egiven = true;}
                else if(k=="nu") mt.nu = val;
                else if(k=="fy" || k=="sigma_y") mt.sigy = val;
                else if(k=="H") mt.H = val;
                else if(k=="eps_fail") mt.epsfail = val;
                else if(k=="ft") mt.ft = val;
                else if(k=="Gf") mt.Gf = val;
                else if(k=="fc") mt.fc = val;
                else if(k=="Gc") mt.Gc = val;
                else if(k=="derode") mt.derode = val;
                else if(k=="stiffness") {if(val<=0.0) fail("material: stiffness must be positive [N/m]"); mt.kdebris = val;}
                else if(k=="crush") {if(val<=0.0) fail("material: crush must be positive [N]"); mt.fcrush = val;}
                else fail("material: unknown parameter '"+k+"' (rho E nu fy H eps_fail ft Gf fc Gc derode stiffness crush, or 'rigid')");
                q += 2;
            }
            if(q<tok.size()) fail("material: cannot read '"+tok[q]+"'");

            if(mt.rho<=0.0 || mt.E<=0.0 || mt.nu<=-1.0 || mt.nu>=0.5)
            fail("invalid material constants");

            add_material(mt);
            cur_mat = mt.id;
        }
        else if(kw=="box")
        {
            double a[6];
            need(static_cast<bool>(ls>>a[0]>>a[1]>>a[2]>>a[3]>>a[4]>>a[5]));
            int id;
            if(!(ls>>id))
            {
                if(cur_mat<0) fail("box: no material defined before");
                id = cur_mat;
            }
            add_box(a[0],a[1],a[2],a[3],a[4],a[5],id);
        }
        else if(kw=="remove")
        {
            double a[6];
            need(static_cast<bool>(ls>>a[0]>>a[1]>>a[2]>>a[3]>>a[4]>>a[5]));
            remove_box(a[0],a[1],a[2],a[3],a[4],a[5]);
        }
        else if(kw=="stl" || kw=="remove_stl")
        {
            // stl file [id] [move dx dy dz] [rotate deg] [scale s]
            shape_stl st;
            st.remove = (kw=="remove_stl");
            st.mat = -1;
            need(static_cast<bool>(ls>>st.file));
            std::string w;
            while(ls>>w)
            {
                if(w=="move") need(static_cast<bool>(ls>>st.move(0)>>st.move(1)>>st.move(2)));
                else if(w=="rotate") need(static_cast<bool>(ls>>st.rot));
                else if(w=="scale") need(static_cast<bool>(ls>>st.scale));
                else
                {
                    char* e = nullptr;
                    const long id = std::strtol(w.c_str(),&e,10);
                    if(!(e && *e=='\0') || st.remove) fail("stl: cannot read '"+w+"'");
                    st.mat = (int)id;
                }
            }
            if(!st.remove && st.mat<0)
            {
                if(cur_mat<0) fail("stl: no material defined before");
                st.mat = cur_mat;
            }
            stls.push_back(st);
            shapes.push_back({1,(int)stls.size()-1});
        }
        else if(kw=="fix")
        {
            // fix base [xyz] | fix top [xyz] | fix bed | fix x0 x1 y0 y1 z0 z1 [xyz]
            std::string first;
            need(static_cast<bool>(ls>>first));
            if(first=="bed")
            copt.fix_bed = 1;
            else if(first=="base" || first=="top")
            {
                std::string dofs = "xyz";
                ls>>dofs;
                fix_box f = {0,0,0,0,0,0,{dofs.find('x')!=std::string::npos,dofs.find('y')!=std::string::npos,dofs.find('z')!=std::string::npos}};
                f.mode = (first=="base") ? 1 : 2;
                fixes.push_back(f);
            }
            else
            {
                double a[6];
                a[0] = std::strtod(first.c_str(),nullptr);
                need(static_cast<bool>(ls>>a[1]>>a[2]>>a[3]>>a[4]>>a[5]));
                std::string dofs = "xyz";
                ls>>dofs;
                add_fix(a[0],a[1],a[2],a[3],a[4],a[5],
                        dofs.find('x')!=std::string::npos,dofs.find('y')!=std::string::npos,dofs.find('z')!=std::string::npos);
            }
        }
        else if(kw=="element")
        {
            std::string type;
            need(static_cast<bool>(ls>>type));
            if(type=="full") full_int = 1;
            else if(type=="reduced") full_int = 0;
            else fail("element type must be 'full' or 'reduced'");
            double hg;
            if(ls>>hg) hg_coef = hg;
        }
        else if(kw=="cfl")           need(static_cast<bool>(ls>>cfl) && cfl>0.0 && cfl<=1.0);
        else if(kw=="damping")
        {
            // damping 2%  |  damping 0.02  |  damping off  |  damping mass alpha[1/s]
            std::string t;
            need(static_cast<bool>(ls>>t));
            if(t=="off" || t=="none")
            zeta = 0.0;
            else if(t=="mass")
            {
                need(static_cast<bool>(ls>>alpha_damp) && alpha_damp>=0.0);
                zeta = 0.0;
            }
            else
            {
                bool pct = false;
                if(!t.empty() && t.back()=='%')
                {
                    pct = true;
                    t.pop_back();
                }
                std::string u;
                if(ls>>u && u=="%")
                pct = true;
                double z = 0.0;
                try {z = std::stod(t);} catch(...) {fail("damping: give a ratio (0.02), a percentage (2%), 'off' or 'mass alpha'");}
                if(pct) z *= 0.01;
                if(z<0.0 || z>=1.0)
                fail("damping: the ratio of critical damping must be between 0 and 1 (2% = 0.02)");
                zeta = z;
            }
        }
        else if(kw=="relax")         need(static_cast<bool>(ls>>relax_time>>relax_alpha));
        else if(kw=="bulk_viscosity")need(static_cast<bool>(ls>>bulkq1>>bulkq2));
        else if(kw=="erode_J")       need(static_cast<bool>(ls>>erode_J));
        else if(kw=="ground")
        {
            // ground z [kfac mu]  |  ground bed [kfac mu]: bed and solids of the fluid grid
            std::string w;
            need(static_cast<bool>(ls>>w));
            if(w=="bed")
            bed_on = true;
            else
            {
                char* e = nullptr;
                zground = std::strtod(w.c_str(),&e);
                if(!(e && *e=='\0')) fail("ground: give a height or 'bed'");
                ground_on = true;
            }
            double k,mu;
            if(ls>>k) kground = k;
            if(ls>>mu) mu_ground = mu;
        }
        else if(kw=="contact")
        {
            std::string onoff;
            need(static_cast<bool>(ls>>onoff));
            contact_on = (onoff=="on" || onoff=="1");
            double k,mu,d;
            if(ls>>k) kcontact = k;
            if(ls>>mu) mu_contact = mu;
            if(ls>>d) contact_dist = d;
        }
        else if(kw=="plane_strain")
        {
            int b;
            need(static_cast<bool>(ls>>b));
            plane_strain = (b!=0);
        }
        else if(kw=="gravity")
        {
            double g[3];
            need(static_cast<bool>(ls>>g[0]>>g[1]>>g[2]));
            grav = Vec3(g[0],g[1],g[2]);
        }
        else if(kw=="monitor")
        {
            // monitor name x y z | monitor name auto|top
            std::string name, w;
            need(static_cast<bool>(ls>>name>>w));
            if(w=="auto" || w=="top")
            mon_autos.push_back({name,0});
            else
            {
                double p[3];
                p[0] = std::strtod(w.c_str(),nullptr);
                need(static_cast<bool>(ls>>p[1]>>p[2]));
                mon_pts.push_back(std::make_pair(name,Vec3(p[0],p[1],p[2])));
            }
        }
        else if(kw=="resolution")
        {
            std::string w;
            need(static_cast<bool>(ls>>w));
            if(w=="coarse") copt.resolution = 2.0;
            else if(w=="normal") copt.resolution = 1.0;
            else if(w=="fine") copt.resolution = 0.5;
            else
            {
                char* e = nullptr;
                copt.resolution = std::strtod(w.c_str(),&e);
                if(!(e && *e=='\0') || copt.resolution<=0.0) fail("resolution: coarse, normal, fine or a factor of the fluid cell size");
            }
        }
        else if(kw=="walls")
        {
            std::string w;
            need(static_cast<bool>(ls>>w) && (w=="on" || w=="off"));
            copt.walls = (w=="on") ? 1 : 0;
        }
        else if(kw=="shear" || kw=="air_forcing")
        {
            std::string w;
            need(static_cast<bool>(ls>>w) && (w=="on" || w=="off"));
            if(kw=="shear") copt.shear = (w=="on") ? 1 : 0;
            else copt.air_forcing = (w=="on") ? 1 : 0;
        }
        else if(kw=="settle" || kw=="snap")
        {
            std::string w;
            need(static_cast<bool>(ls>>w));
            const bool on = (w=="on" || w=="1");
            if(kw=="settle") copt.settle = on ? 1 : 0;
            else snap_on = on;
        }
        else if(kw=="fragments")
        {
            // fragments rigid|deformable: parts that break off a supported structure
            std::string w;
            need(static_cast<bool>(ls>>w) && (w=="rigid" || w=="deformable"));
            fragments_rigid = (w=="rigid");
        }
        else if(kw=="check")
        copt.check = 1;
        // coupling options
        else if(kw=="pressure_offset")  need(static_cast<bool>(ls>>copt.pressure_offset));
        else if(kw=="points_per_face")  need(static_cast<bool>(ls>>copt.points_per_face));
        else if(kw=="debris_cd")        need(static_cast<bool>(ls>>copt.debris_cd));
        else if(kw=="debris_reaction")  need(static_cast<bool>(ls>>copt.debris_reaction));
        else if(kw=="forcing")          need(static_cast<bool>(ls>>copt.forcing));
        else if(kw=="print")            need(static_cast<bool>(ls>>copt.print_dt));
        else if(kw=="hybrid_tau")       need(static_cast<bool>(ls>>copt.hybrid_tau));
        else if(kw=="added_mass")       need(static_cast<bool>(ls>>copt.added_mass) && copt.added_mass>=0.0);
        else if(kw=="rigid_contact_speed") need(static_cast<bool>(ls>>c_rigid) && c_rigid>0.0);
        else if(kw=="debris_damping")
        {
            // debris_damping 5%  |  debris_damping 0.05
            std::string t;
            need(static_cast<bool>(ls>>t));
            bool pct = false;
            if(!t.empty() && t.back()=='%') {pct = true; t.pop_back();}
            double z = -1.0;
            try {z = std::stod(t);} catch(...) {fail("debris_damping: give a ratio (0.05) or a percentage (5%)");}
            if(pct) z *= 0.01;
            if(z<0.0 || z>=1.0) fail("debris_damping: the ratio of critical damping must be between 0 and 1");
            debris_zeta = z;
        }
        else if(kw=="loads")
        {
            std::string t;
            need(static_cast<bool>(ls>>t));
            if(t=="hybrid") copt.loads = 0;
            else if(t=="pressure") copt.loads = 1;
            else if(t=="reaction") copt.loads = 2;
            else fail("loads must be 'hybrid', 'pressure' or 'reaction'");
        }
        else
        fail("unknown keyword '"+kw+"'");
    }
}

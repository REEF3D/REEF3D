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
#include<sstream>
#include<istream>
#include<stdexcept>

// fem.dat: one keyword per line, '#' starts a comment, SI units.
// See docs/FEM.md for the full list.

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
            material mt;
            std::string type;
            need(static_cast<bool>(ls>>mt.id>>type));

            if(type=="elastic")
            {
                mt.type = MAT_ELASTIC;
                need(static_cast<bool>(ls>>mt.rho>>mt.E>>mt.nu));
            }
            else if(type=="plastic" || type=="steel" || type=="j2")
            {
                mt.type = MAT_J2;
                need(static_cast<bool>(ls>>mt.rho>>mt.E>>mt.nu>>mt.sigy>>mt.H>>mt.epsfail));
            }
            else if(type=="concrete")
            {
                mt.type = MAT_CONCRETE;
                need(static_cast<bool>(ls>>mt.rho>>mt.E>>mt.nu>>mt.ft>>mt.Gf>>mt.fc>>mt.Gc));
                double de;
                if(ls>>de) mt.derode = de;
            }
            else
            fail("unknown material type '"+type+"' (elastic, plastic, concrete)");

            if(mt.rho<=0.0 || mt.E<=0.0 || mt.nu<=-1.0 || mt.nu>=0.5)
            fail("invalid material constants");

            add_material(mt);
        }
        else if(kw=="box")
        {
            double a[6]; int id;
            need(static_cast<bool>(ls>>a[0]>>a[1]>>a[2]>>a[3]>>a[4]>>a[5]>>id));
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
            shape_stl s;
            s.remove = (kw=="remove_stl");
            s.mat = -1;
            need(static_cast<bool>(ls>>s.file));
            if(!s.remove)
            need(static_cast<bool>(ls>>s.mat));
            stls.push_back(s);
            shapes.push_back({1,(int)stls.size()-1});
        }
        else if(kw=="fix")
        {
            double a[6];
            need(static_cast<bool>(ls>>a[0]>>a[1]>>a[2]>>a[3]>>a[4]>>a[5]));
            std::string dofs = "xyz";
            ls>>dofs;
            add_fix(a[0],a[1],a[2],a[3],a[4],a[5],
                    dofs.find('x')!=std::string::npos,dofs.find('y')!=std::string::npos,dofs.find('z')!=std::string::npos);
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
        else if(kw=="damping")       need(static_cast<bool>(ls>>alpha_damp));
        else if(kw=="relax")         need(static_cast<bool>(ls>>relax_time>>relax_alpha));
        else if(kw=="bulk_viscosity")need(static_cast<bool>(ls>>bulkq1>>bulkq2));
        else if(kw=="erode_J")       need(static_cast<bool>(ls>>erode_J));
        else if(kw=="ground")
        {
            need(static_cast<bool>(ls>>zground));
            ground_on = true;
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
            std::string name; double p[3];
            need(static_cast<bool>(ls>>name>>p[0]>>p[1]>>p[2]));
            mon_pts.push_back(std::make_pair(name,Vec3(p[0],p[1],p[2])));
        }
        // coupling options
        else if(kw=="pressure_offset")  need(static_cast<bool>(ls>>copt.pressure_offset));
        else if(kw=="points_per_face")  need(static_cast<bool>(ls>>copt.points_per_face));
        else if(kw=="debris_cd")        need(static_cast<bool>(ls>>copt.debris_cd));
        else if(kw=="debris_reaction")  need(static_cast<bool>(ls>>copt.debris_reaction));
        else if(kw=="forcing")          need(static_cast<bool>(ls>>copt.forcing));
        else if(kw=="print")            need(static_cast<bool>(ls>>copt.print_dt));
        else if(kw=="loads")
        {
            std::string t;
            need(static_cast<bool>(ls>>t));
            if(t=="reaction") copt.loads = 0;
            else if(t=="pressure") copt.loads = 1;
            else fail("loads must be 'reaction' or 'pressure'");
        }
        else
        fail("unknown keyword '"+kw+"'");
    }
}

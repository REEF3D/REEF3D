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

//  Resolved bodies of REEF3D::FNPF in checkpoints (P 48, I 45).  The blob: the bodies (the kernel
//  state of sixdof_obj and the FNPF part of sixdof_obj_fnpf) and the footprint state of level 0.
//  The arrays of the body on a grid (unit-mode and psi_0 solutions, the start values of the next psi
//  solves; footprint, its last extension of eta and Fifsf) with those of the grid: level 0 in the
//  fields of FNPF, a patch with the arrays of the patch (fnpf_amr_ckpt.cpp).  The body nodes and the
//  Neumann data of the body faces (FBF, FBu, FBv, FBw) are fields of fdm_fnpf (fnpf_ckpt_fields.h).
//  On a restart the body grids of the patches are made for the patches of the checkpoint before
//  their arrays are read; a patch without arrays (made by the regrid at the end of the step that
//  wrote the checkpoint) is set up after the bodies, as the next stage would have done.

#include"fnpf_6DOF.h"
#include"6DOF_obj_fnpf.h"
#include"lexer.h"
#include"fdm_fnpf.h"
#include"ghostcell.h"
#include"fnpf_amr.h"
#include"reefamr.h"
#include<cmath>
#include<limits>
#include<string>

void fnpf_6DOF::grid_fields(fnpf_6DOF_grid *G, const ckpt_fd &fd, const ckpt_fi &fi)
{
    static const char *pn[6] = {"fb_psi_0","fb_psi_1","fb_psi_2","fb_psi_3","fb_psi_4","fb_psi_5"};
    
    fd("fb_psi0",'F',G!=nullptr ? G->psi0 : nullptr);
    for(int m=0; m<6; ++m)
    fd(pn[m],'F',G!=nullptr ? G->psi[m] : nullptr);
    // several bodies: the unit modes of the bodies 1, 2, ... (fb_psi_<body>_<mode>)
    for(int nb=1; nb<nbody; ++nb)
    for(int m=0; m<6; ++m)
    {
        const size_t e = 6*(nb-1)+m;
        const std::string name = "fb_psi_"+std::to_string(nb)+"_"+std::to_string(m);
        fd(name.c_str(),'F',(G!=nullptr && G->psib.size()>e) ? G->psib[e] : nullptr);
    }
    fi("fb_mark",'F',G!=nullptr ? G->mark : nullptr);
    fd("fb_foot",'S',G!=nullptr ? G->foot->V : nullptr);
    fd("fb_psiD",'S',G!=nullptr ? G->psiD->V : nullptr);
    fd("fb_eta_ext",'S',G!=nullptr ? G->eta_ext->V : nullptr);
    fd("fb_fi_ext",'S',G!=nullptr ? G->fi_ext->V : nullptr);
}

void fnpf_6DOF::ckpt_fields(const ckpt_fd &fd, const ckpt_fi &fi)
{
    grid_fields(&g0,fd,fi);
}

void fnpf_6DOF::amr_ckpt_fields(lexer *pp, const ckpt_fd &fd, const ckpt_fi &fi)
{
    // the grid of the patch by its serial number (the grids may still be those of the layout before
    // the last regrid); on a restart the grids made for the patches of the checkpoint
    fnpf_6DOF_grid *G = nullptr;
    
    if(amr!=nullptr)
    for(int n=0; n<amr->patches(); ++n)
    if(amr->patch_lexer(n)==pp)
    {
        const int serial = amr->patch_serial(n);
        for(auto &g : gp)
        if(g.serial==serial && (ckpt_grids || !g.fresh))
        G = &g;
    }
    
    grid_fields(G,fd,fi);
}

void fnpf_6DOF::amr_ckpt_grids(lexer *p, ghostcell *pgc)
{
    // the grids of the patches of the checkpoint, without set-up; a NaN in the extension of eta
    // marks the grids whose arrays do not come
    ckpt_grids = true;
    amr_layout = -1;
    amr_grids(p,pgc);

    // X 36: the unit modes of the bodies 1, 2, ... live on the patches of their windows only (allocated
    // at their first solve); their arrays for the checkpoint's values on the patches that reach into
    // the zone of the body grown by its margin once more (a superset of the windows of the run; the
    // others stay unallocated)
    if(p->X36==1 && p->G12>0 && p->G7!=1 && nbody>1 && amr!=nullptr)
    for(auto &G : gp)
    {
        double x0,x1,y0,y1;
        amr->patch_box(G.id,x0,x1,y0,y1);
        for(int nb=1; nb<nbody; ++nb)
        {
            double a0,a1,b0,b1;
            if(!amr->zone_box(nb,a0,a1,b0,b1))
            continue;
            const double g = p->G12_r;
            if(x0<a1+g && x1>a0-g && (p->j_dir==0 || (y0<b1+g && y1>b0-g)))
            for(int m=0; m<6; ++m)
            if(psi[m]!=nullptr && !fb_obj[nb]->dof_fixed_fnpf(p,m))
            mode_array(G,nb,m);
        }
    }
    
    for(auto &G : gp)
    {
        lexer *pp = G.p;
        const double nan = std::numeric_limits<double>::quiet_NaN();
        const int n = pp->imax*pp->jmax;
        for(int q=0; q<n; ++q)
        G.eta_ext->V[q] = nan;
    }
}

void fnpf_6DOF::ckpt_state(lexer *p, std::vector<double> &v)
{
    v.push_back(2.0);   // version (2: with the added mass of X 28)
    v.push_back(double(nbody));
    
    for(int nb=0; nb<nbody; ++nb)
    {
        std::vector<double> s;
        fb_obj[nb]->ckpt_state(p,s);
        v.push_back(double(s.size()));
        v.insert(v.end(),s.begin(),s.end());
        
        s.clear();
        fb_obj[nb]->ckpt_fnpf_state(s);
        v.push_back(double(s.size()));
        v.insert(v.end(),s.begin(),s.end());
    }
    
    v.push_back(double(g0.footcount));
    v.push_back(g0.ext_ini ? 1.0 : 0.0);

    // version 2, X 28: the step and the poses of the last added mass, the added mass of every body and
    // the cross terms (a restart continues with them until the next one is due, as the run would)
    const bool am = (am_count>=0 && (int)am_c.size()==nbody);
    v.push_back(am ? double(am_count) : -1.0);
    if(am)
    {
        for(int nb=0; nb<nbody; ++nb)
        {
            for(int d=0; d<3; ++d)
            v.push_back(am_c[nb](d));
            for(int d=0; d<3; ++d)
            v.push_back(am_ang[nb](d));
            const Eigen::Matrix<double,6,6> &A = fb_obj[nb]->external_added_mass();
            v.push_back(fb_obj[nb]->external_added_mass_on() ? 1.0 : 0.0);
            for(int a=0; a<6; ++a)
            for(int b=0; b<6; ++b)
            v.push_back(A(a,b));
        }
        for(const auto &A : Across)
        for(int a=0; a<6; ++a)
        for(int b=0; b<6; ++b)
        v.push_back(A(a,b));
    }
}

bool fnpf_6DOF::ckpt_restore(lexer *p, ghostcell *pgc, const std::vector<double> &v)
{
    size_t k=0;
    auto next = [&](){ return k<v.size() ? v[k++] : -1.0; };
    
    const int version = int(next());
    if((version!=1 && version!=2) || int(next())!=nbody)
    return false;
    
    for(int nb=0; nb<nbody; ++nb)
    {
        size_t n = size_t(next());
        if(k+n>v.size() || !fb_obj[nb]->ckpt_restore(p,pgc,std::vector<double>(v.begin()+k,v.begin()+k+n)))
        return false;
        k += n;
        
        n = size_t(next());
        if(k+n>v.size() || !fb_obj[nb]->ckpt_fnpf_restore(std::vector<double>(v.begin()+k,v.begin()+k+n)))
        return false;
        k += n;
    }
    
    if(k+2>v.size())
    return false;
    g0.footcount = int(next());
    g0.ext_ini = next()>0.5;

    // version 2: the added mass of X 28 (version 1: computed at the first step)
    am_count = -1;
    am_c.clear();
    am_ang.clear();
    if(version>=2)
    {
        if(k+1>v.size())
        return false;
        const int cnt = int(next());
        if(cnt>=0)
        {
            const size_t need = size_t(nbody)*(6+1+36) + Across.size()*36;
            if(k+need>v.size())
            return false;
            am_c.assign(nbody,Eigen::Vector3d::Zero());
            am_ang.assign(nbody,Eigen::Vector3d::Zero());
            for(int nb=0; nb<nbody; ++nb)
            {
                for(int d=0; d<3; ++d)
                am_c[nb](d) = next();
                for(int d=0; d<3; ++d)
                am_ang[nb](d) = next();
                const bool on = next()>0.5;
                Eigen::Matrix<double,6,6> A;
                for(int a=0; a<6; ++a)
                for(int b=0; b<6; ++b)
                A(a,b) = next();
                if(on)
                fb_obj[nb]->fnpf_added_mass_restore(A);
            }
            for(auto &A : Across)
            for(int a=0; a<6; ++a)
            for(int b=0; b<6; ++b)
            A(a,b) = next();
            am_count = cnt;
        }
    }

    initialized = true;
    
    // the body grids of the patches: those with the arrays of the checkpoint as they were, the
    // others set up at the restored pose
    if(ckpt_grids)
    {
        for(auto &G : gp)
        {
            lexer *p = G.p;
            slice4 &foot = *G.foot;
            slice4 &eta_ext = *G.eta_ext;
            bool got = false;
            int count = 0;
            
            SLICELOOP4
            {
            if(!std::isnan(eta_ext(i,j)))
            got = true;
            
            if(foot(i,j)>0.5)
            ++count;
            }
            
            if(got)
            {
                G.fresh = false;
                G.ext_ini = true;
                G.footcount = count;
            }
            
            const int n = p->imax*p->jmax;
            for(int q=0; q<n; ++q)
            if(std::isnan(eta_ext.V[q]))
            eta_ext.V[q] = 0.0;
        }
        
        ckpt_grids = false;
        amr_fresh(pgc);
    }
    
    return true;
}

// the state of the rank: per body (X 18 running means of chi)
void fnpf_6DOF::ckpt_rank_state(std::vector<double> &v)
{
    v.push_back(double(nbody));
    for(int nb=0; nb<nbody; ++nb)
    {
        std::vector<double> s;
        fb_obj[nb]->ckpt_fnpf_rank_state(s);
        v.push_back(double(s.size()));
        v.insert(v.end(),s.begin(),s.end());
    }
}

void fnpf_6DOF::ckpt_rank_restore(const std::vector<std::vector<double>> &all)
{
    for(int nb=0; nb<nbody; ++nb)
    {
        std::vector<std::vector<double>> parts;
        for(const std::vector<double> &v : all)
        {
            size_t k=0;
            if(v.empty() || int(v[k++])!=nbody)
            continue;
            for(int q=0; q<=nb && k<v.size(); ++q)
            {
                const size_t n = size_t(v[k++]);
                if(k+n>v.size())
                break;
                if(q==nb)
                parts.push_back(std::vector<double>(v.begin()+k,v.begin()+k+n));
                k += n;
            }
        }
        fb_obj[nb]->ckpt_fnpf_rank_restore(parts);
    }
}

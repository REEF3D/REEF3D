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
Authors: Hans Bihs, Tobias Martin
--------------------------------------------------------------------*/

#include"6DOF_sflow.h"
#include"lexer.h"
#include"fdm2D.h"
#include"ghostcell.h"
#include<cstdio>
#include<string>
   
sixdof_sflow::sixdof_sflow(lexer *p, ghostcell *pgc) : press(p)
{
    if(p->mpirank==0)
    cout<<"6DOF startup ..."<<endl;
    
    number6DOF = 1;
    
    for (int nb = 0; nb < number6DOF; nb++)
    fb_obj.push_back(new sixdof_obj_2D(p,pgc,nb));
}

sixdof_sflow::~sixdof_sflow()
{
}

void sixdof_sflow::start_sflow(lexer *p, fdm2D *b, ghostcell *pgc, int iter, slice &fsglobal, slice &P, slice &Q, slice &w, slice &fx, slice &fy, slice &eta, bool finalize)
{
    starttime = pgc->timer();
    
    if(p->X10==2)
    start_oneway(p,pgc,iter,fsglobal,P,Q,w,fx,fy,eta,finalize);
    
    if(p->X10==3)
    start_shipwave(p,pgc,iter,fsglobal,P,Q,fx,fy,eta,finalize);
    
    p->fbtime+=pgc->timer()-starttime;
}

void sixdof_sflow::start_oneway(lexer *p, ghostcell *pgc, int iter, slice &fsglobal, slice &P, slice &Q, slice &w, slice &fx, slice &fy, slice &eta, bool finalize)
{
    for (int nb=0; nb<number6DOF;++nb)
    {
        // Advance body in time
        fb_obj[nb]->solve_eqmotion_oneway_sflow(p,pgc,iter,finalize);
        
        // Update transformation matrices
        fb_obj[nb]->quat_matrices(p);
        
        // Update position and trimesh
        fb_obj[nb]->update_position_2D(p,pgc,fsglobal);  
        
        // Save
        fb_obj[nb]->update_fbvel(p,pgc);
        
        // Update forcing terms
        fb_obj[nb]->update_forcing_sflow(p,pgc,P,Q,w,fx,fy,eta,iter);
        
            // Print
        if(finalize==true)
        {
            fb_obj[nb]->saveTimeStep(p,iter);
            
            if(p->X50==1)
            fb_obj[nb]->print_vtp(p,pgc);
            
            if(p->X50==2)
            fb_obj[nb]->print_stl(p,pgc);
            
            fb_obj[nb]->print_parameter(p,pgc);
        }
    }
}

void sixdof_sflow::start_shipwave(lexer *p, ghostcell *pgc, int iter, slice &fsglobal, slice &P, slice&Q, slice &fx, slice &fy, slice &eta, bool finalize)
{
    if(finalize==1)
    for (int nb=0; nb<number6DOF;++nb)
    {
        // Advance body in time
        fb_obj[nb]->solve_eqmotion_oneway_onestep(p,pgc,finalize);
        
        // Update transformation matrices
        fb_obj[nb]->quat_matrices(p);
        
        // Update position and trimesh
        fb_obj[nb]->update_position_2D(p,pgc,fsglobal);  
        
        // Save
        fb_obj[nb]->update_fbvel(p,pgc);
        
        // Update forcing terms
        if(p->X400==2)
        fb_obj[nb]->updateForcing_box(p,pgc,press);
        
        if(p->X400==3)
        fb_obj[nb]->updateForcing_oned(p,pgc,press);
        
        if(p->X400==10)
        fb_obj[nb]->updateForcing_stl(p,pgc,press,eta);
        
            // Print
            if(p->X50==1)
            fb_obj[nb]->print_vtp(p,pgc);
            
            if(p->X50==2)
            fb_obj[nb]->print_stl(p,pgc);
            
            fb_obj[nb]->print_parameter(p,pgc);
    }
}

// checkpoints (P 48, I 45)
void sixdof_sflow::ckpt_state(lexer *p, std::vector<double> &v)
{
    v.push_back(double(number6DOF));
    for(int nb=0; nb<number6DOF; ++nb)
    {
        std::vector<double> s;
        fb_obj[nb]->ckpt_state(p,s);
        v.push_back(double(s.size()));
        v.insert(v.end(),s.begin(),s.end());
    }
}

bool sixdof_sflow::ckpt_restore(lexer *p, ghostcell *pgc, const std::vector<double> &v)
{
    size_t k=0;
    if(v.empty() || int(v[k++])!=number6DOF)
    return false;

    for(int nb=0; nb<number6DOF; ++nb)
    {
        const size_t n = (k<v.size()) ? size_t(v[k++]) : 0;
        if(n==0 || k+n>v.size())
        return false;
        if(!fb_obj[nb]->ckpt_restore(p,pgc,std::vector<double>(v.begin()+k,v.begin()+k+n)))
        return false;
        k += n;
    }
    return true;
}

void sixdof_sflow::ckpt_fields(const std::function<void(const char*,char,double*)> &fd, const std::function<void(const char*,char,int*)> &fi)
{
    for(int nb=0; nb<number6DOF; ++nb)
    {
        char pre[32];
        snprintf(pre,sizeof(pre),"sixdof%d_",nb);
        fb_obj[nb]->ckpt_2D_fields(pre,fd,fi);
    }
}

void sixdof_obj_2D::ckpt_2D_fields(const char *pre, const std::function<void(const char*,char,double*)> &fd, const std::function<void(const char*,char,int*)> &fi)
{
    struct { const char *n; slice4 *s; } sl[] = {{"press",&press},{"lrk1",&lrk1},{"lrk2",&lrk2},{"K",&K},{"dts",&dts},{"fs",&fs},
        {"Ls",&Ls},{"Bs",&Bs},{"Rxmin",&Rxmin},{"Rxmax",&Rxmax},{"Rymin",&Rymin},{"Rymax",&Rymax},{"draft",&draft}};
    for(auto &e : sl)
    {
        const std::string n = std::string(pre) + e.n;
        fd(n.c_str(),'S',e.s->V);
    }
    struct { const char *n; sliceint5 *s; } si[] = {{"cl",&cl},{"cr",&cr},{"fsio",&fsio}};
    for(auto &e : si)
    {
        const std::string n = std::string(pre) + e.n;
        fi(n.c_str(),'S',e.s->V);
    }
}

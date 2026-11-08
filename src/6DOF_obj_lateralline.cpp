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

#include"6DOF_obj_cfd.h"
#include"6DOF_deform.h"
#include"6DOF_output_dir.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"reefckpt.h"
#include"6DOF_brain.h"
#include<algorithm>

//  Lateral line of a 6DOF body (X 194 id n d): the senses of a fish along its flanks. Probes sit a
//  distance d (0: two cells) outside the hull on both sides of n stations along the body axis
//  (s/L = (k + 1/2)/n from the head), in the plane of the undulation through the widest point of
//  the section. They are material points of the body frame: with a deforming body (X 195) they
//  follow the flank. At each probe the fluid pressure and the fluid velocity relative to the hull
//  (rigid-body motion and deformation), in the body frame (axial, lateral, third axis), each the
//  mean of 3 x 3 points (one cell apart along the axis, half a cell normal to the flank): the
//  pressure next to the forced hull varies from cell to cell. Sampled each step when the body has
//  a brain (X 193), else at the output steps (X 19).

void sixdof_obj::lateral_line_setup(lexer *p)
{
    ll_pts_.clear();
    ll_s_.clear();
    ll_side_.clear();
    
    int n = 0;
    double d = 0.0;
    
    for(int q=0; q<p->X194; ++q)
    if(p->X194_id[q]==n6DOF)
    {
    n = p->X194_n[q];
    d = p->X194_d[q];
    }
    
    if(n<=0 || tricount<=0)
    return;
    
    if(d<=0.0)
    d = 2.0*p->DXM;
    
    ll_dx_ = p->DXM;
    
    // axes of the body frame: x along the body, l the lateral axis (2D: z), v the third axis
    const int la = (p->j_dir==0) ? 2 : 1, va = 3 - la;
    double **c[3] = {tri_x0, tri_y0, tri_z0};
    
    double xh = -1.0e20, xt = 1.0e20;
    
    for(int t=0; t<tricount; ++t)
    for(int q=0; q<3; ++q)
    {
    xh = std::max(xh, tri_x0[t][q]);
    xt = std::min(xt, tri_x0[t][q]);
    }
    
    const double L = xh - xt;
    
    for(int k=0; k<n; ++k)
    {
        const double s = (k + 0.5)/double(n);
        const double xs = xh - s*L;
        
        // the section x = xs: points of the triangle edges that cross it, the extremes of l
        double lmax = -1.0e20, lmin = 1.0e20, vmax = 0.0, vmin = 0.0;
        
        for(int t=0; t<tricount; ++t)
        for(int q=0; q<3; ++q)
        {
            const int r = (q+1)%3;
            const double xa = tri_x0[t][q], xb = tri_x0[t][r];
            
            if((xa-xs)*(xb-xs)>0.0 || xa==xb)
            continue;
            
            const double w = (xs - xa)/(xb - xa);
            const double l = c[la][t][q] + w*(c[la][t][r] - c[la][t][q]);
            const double v = c[va][t][q] + w*(c[va][t][r] - c[va][t][q]);
            
            if(l>lmax) {lmax = l; vmax = v;}
            if(l<lmin) {lmin = l; vmin = v;}
        }
        
        if(lmax<lmin)
        continue;
        
        // in 2D the third axis is the width of the slab: the probe at its centre
        if(p->j_dir==0)
        vmax = vmin = 0.0;
        
        ll_pts_.push_back(Eigen::Vector3d(xs, lmax + d, vmax));
        ll_s_.push_back(s);
        ll_side_.push_back(1);
        
        ll_pts_.push_back(Eigen::Vector3d(xs, lmin - d, vmin));
        ll_s_.push_back(s);
        ll_side_.push_back(-1);
    }
    
    if(p->mpirank==0 && !ll_pts_.empty())
    {
        cout<<"6DOF body "<<n6DOF<<": lateral line (X 194), "<<ll_pts_.size()<<" probes, "<<d<<" m outside the hull"<<endl;
        
        char str[1000];
        snprintf(str,sizeof(str),"%s/REEF3D_6DOF_lateralline_%i.dat",sixdof_output_dir(p),n6DOF);
        
        if(reefckpt_open(printll,str,p->I45))
        {
            printll<<"# lateral line of body "<<n6DOF<<": per probe p [Pa] and the fluid velocity relative to the hull in the body frame"
                   <<" (axial, lateral, third axis) [m/s]; side +1: towards +"<<(la==1 ? "y" : "z")<<"_B"<<endl;
            printll<<"# probe \t s/L \t side"<<endl;
            
            for(size_t q=0; q<ll_pts_.size(); ++q)
            printll<<"# "<<q<<" \t "<<ll_s_[q]<<" \t "<<ll_side_[q]<<endl;
            
            printll<<"time";
            
            for(size_t q=0; q<ll_pts_.size(); ++q)
            printll<<" \t p_"<<q<<" \t ua_"<<q<<" \t ul_"<<q<<" \t uv_"<<q;
            
            printll<<endl;
        }
    }
}

Eigen::Vector3d sixdof_obj::lateral_line_point(int q, double da, double dn) const
{
    // the material point (x, l, v) of the body frame, moved by da along the axis and dn away from the
    // flank, moved by the deformation, to the inertial frame
    const Eigen::Vector3d m(ll_pts_[q](0) + da, ll_pts_[q](1) + ll_side_[q]*dn, ll_pts_[q](2));
    double dx = 0.0, dl = 0.0;
    
    if(fish_!=nullptr)
    fish_->displacement(m(0),m(1),dx,dl);
    
    const int la = (fish_!=nullptr) ? fish_->lateral() : (rb.twoD ? 2 : 1), va = 3 - la;
    
    Eigen::Vector3d xb = Eigen::Vector3d::Zero();
    xb(0)  = m(0) + dx;
    xb(la) = m(1) + dl;
    xb(va) = m(2);
    
    return c_ + R_*xb;
}

void sixdof_obj_cfd::lateral_line(lexer *p, fdm *a, ghostcell *pgc)
{
    const bool out = (p->count%p->X19==0);
    
    // vision of the brain at its decision steps
    std::vector<double> vis;
    if(brain_due(p) && brain_->vision_range()>0.0)
    vis = vision(p,a,pgc);
    const std::vector<double> *pv = vis.empty() ? nullptr : &vis;
    
    if(ll_pts_.empty() || (!out && brain_==nullptr))
    {
        brain_update(p,nullptr,pv);
        return;
    }
    
    const int n = int(ll_pts_.size());
    std::vector<double> val(4*n, 0.0);
    std::vector<Eigen::Vector3d> X(n);
    const Eigen::Vector3d U(p->ufbi, p->vfbi, p->wfbi), W(p->pfbi, p->qfbi, p->rfbi);
    
    for(int q=0; q<n; ++q)
    {
        X[q] = lateral_line_point(q,0.0,0.0);
        
        // the mean over 3 x 3 points; the rank that owns a point interpolates, the others contribute 0
        for(int ia=-1; ia<=1; ++ia)
        for(int in=-1; in<=1; ++in)
        {
            const Eigen::Vector3d Y = (ia==0 && in==0) ? X[q] : lateral_line_point(q, ia*ll_dx_, 0.5*in*ll_dx_);
            const double x = Y(0), y = Y(1), z = Y(2);
            
            const bool own = (x>=p->originx && x<p->endx) && (p->j_dir==0 || (y>=p->originy && y<p->endy))
                          && (z>=p->originz && z<p->endz);
            
            if(own)
            {
                // the fluid velocity relative to the hull at the point: rigid-body motion and deformation
                Eigen::Vector3d ub = U + W.cross(Y - c_);
                
                if(fish_!=nullptr)
                ub += deform_vel(x,y,z);
                
                val[4*q]   += p->ccipol4(a->press, x, y, z)/9.0;
                val[4*q+1] += (p->ccipol1(a->u, x, y, z) - ub(0))/9.0;
                val[4*q+2] += (p->ccipol2(a->v, x, y, z) - ub(1))/9.0;
                val[4*q+3] += (p->ccipol3(a->w, x, y, z) - ub(2))/9.0;
            }
        }
    }
    
    for(int q=0; q<4*n; ++q)
    val[q] = pgc->globalsum(val[q]);
    
    // the brain: the pressures of the probes
    std::vector<double> pr(n);
    for(int q=0; q<n; ++q)
    pr[q] = val[4*q];
    
    brain_update(p,&pr,pv);
    
    if(!out || p->mpirank!=0 || !printll.is_open())
    return;
    
    const int la = (fish_!=nullptr) ? fish_->lateral() : (rb.twoD ? 2 : 1), va = 3 - la;
    
    printll<<p->simtime;
    
    for(int q=0; q<n; ++q)
    {
        const Eigen::Vector3d ur = R_.transpose()*Eigen::Vector3d(val[4*q+1], val[4*q+2], val[4*q+3]);
        printll<<" \t "<<val[4*q]<<" \t "<<ur(0)<<" \t "<<ur(la)<<" \t "<<ur(va);
    }
    
    printll<<endl;
}

std::vector<double> sixdof_obj_cfd::vision(lexer *p, fdm *a, ghostcell *pgc)
{
    // rays from the head in the plane of the undulation: the distance to the walls of the domain
    // (exact) and to the solids (their level set, sampled every half cell), up to the range D
    const std::vector<double> &ang = sixdof_brain::rays();
    const double D = brain_->vision_range();
    const int la = fish_->lateral(), nrm = (la==1) ? 2 : 1;
    
    Eigen::Vector3d ex = R_.col(0), el = R_.col(la);
    ex(nrm) = 0.0;
    el(nrm) = 0.0;
    ex.normalize();
    el = (el - el.dot(ex)*ex).normalized();
    const Eigen::Vector3d O = c_ + fish_->xhead()*R_.col(0);
    
    const double lo[3] = {p->global_xmin, p->global_ymin, p->global_zmin};
    const double hi[3] = {p->global_xmax, p->global_ymax, p->global_zmax};
    
    std::vector<double> dist(ang.size(), D);
    const double ds = 0.5*p->DXM;
    const int ns = int(D/ds);
    
    for(size_t r=0; r<ang.size(); ++r)
    {
        const Eigen::Vector3d e = std::cos(ang[r])*ex + std::sin(ang[r])*el;
        
        // the walls of the domain
        for(int c=0; c<3; ++c)
        if(c!=nrm && std::fabs(e(c))>1.0e-12)
        {
            const double s = ((e(c)>0.0 ? hi[c] : lo[c]) - O(c))/e(c);
            if(s>=0.0)
            dist[r] = std::min(dist[r], s);
        }
        
        // the solids: the first sample inside (level set < 0), the rank that owns it
        if(p->solidread>0)
        {
            double hit = 1.0e20;
            
            for(int k=1; k<=ns; ++k)
            {
                const Eigen::Vector3d X = O + (k*ds)*e;
                const bool own = (X(0)>=p->originx && X(0)<p->endx) && (p->j_dir==0 || (X(1)>=p->originy && X(1)<p->endy))
                              && (X(2)>=p->originz && X(2)<p->endz);
                
                if(own && p->ccipol4a(a->solid, X(0), X(1), X(2))<0.0)
                {
                    hit = k*ds;
                    break;
                }
            }
            
            hit = pgc->globalmin(hit);
            dist[r] = std::min(dist[r], hit);
        }
    }
    
    return dist;
}

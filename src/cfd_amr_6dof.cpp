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

#include"cfd_amr_6dof.h"
#include"6DOF_cfd.h"
#include"6DOF_obj_cfd.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"reinidisc_f.h"
#include"field1.h"
#include"field2.h"
#include"field3.h"
#include<algorithm>

cfd_amr_6dof::cfd_amr_6dof(lexer *p, sixdof_cfd *bb) : b(bb), cutl(p), cutr(p), fbio(p), frk1(p), L(p), dt(p)
{
    prdisc = new reinidisc_f(p,true,p->X14==2 ? 1.0e-6 : 0.0);
    
    if(b->objects()>1)
    {
    fball = new field4a(p);
    fbc.resize(b->objects());
    lastb.resize(b->objects());
    }
}

cfd_amr_6dof::~cfd_amr_6dof()
{
    delete prdisc;
    delete fball;
    
    for(int q=0; q<3; ++q)
    delete sf[q];
}

// body nb on the grid: ray cast, reinitialisation
void cfd_amr_6dof::cast(lexer *p, fdm *a, ghostcell *pgc, int nb)
{
    sixdof_obj_cfd *o = b->object(nb);
    o->ray_cast_amr(p,a,pgc,fbio,cutl,cutr,box);
    o->reini_amr(p,a,pgc,a->fb,prdisc,frk1,L,dt,box);
}

// the solid Heaviside functions of the stage (sixdof_cfd::setup)
void cfd_amr_6dof::reset_fbh(lexer *p, fdm *a, ghostcell *pgc)
{
    // the weights of the solids and the topography (ghostcell::solid_forcing of this stage); the
    // fields at the first stage (the grid of a patch gets its solids after this object)
    if(sf[0]==nullptr && (p->solidread>0 || p->topoforcing>0))
    {
    sf[0] = new field1(p);
    sf[1] = new field2(p);
    sf[2] = new field3(p);
    }
    
    if(sf[0]!=nullptr)
    {
    std::copy(a->fbh1.begin(), a->fbh1.end(), sf[0]->data());
    std::copy(a->fbh2.begin(), a->fbh2.end(), sf[1]->data());
    std::copy(a->fbh3.begin(), a->fbh3.end(), sf[2]->data());
    }
    
    ULOOP
    a->fbh1(i,j,k) = 0.0;

    VLOOP
    a->fbh2(i,j,k) = 0.0;

    WLOOP
    a->fbh3(i,j,k) = 0.0;

    LOOP
    a->fbh4(i,j,k) = 0.0;

    pgc->start1(p,a->fbh1,10);
    pgc->start2(p,a->fbh2,11);
    pgc->start3(p,a->fbh3,12);
    pgc->start4(p,a->fbh4,40);
}

// true if no vertex coordinate of the bodies is further than tol from the triangles of the last cast
// on this grid (tol 0: the same triangles); otherwise the current triangles are kept for the next
// check (tol < 0: always)
bool cfd_amr_6dof::unchanged(double tol)
{
    const int nob = b->objects();
    bool same = ((int)last.size()==nob);

    for(int nb=0; nb<nob && same; ++nb)
    {
        sixdof_obj_cfd *o = b->object(nb);
        const int n = o->amr_tricount();
        const vector<double> &v = last[nb];
        if((int)v.size()!=9*n)
        {
            same = false;
            break;
        }
        size_t c = 0;
        for(int d=0; d<3 && same; ++d)
        {
            double **t = o->amr_tri(d);
            for(int q=0; q<n && same; ++q)
            for(int e=0; e<3; ++e)
            {
                if(!(fabs(v[c]-t[q][e])<=tol))
                same = false;
                ++c;
            }
        }
    }

    if(!same)
    {
        last.assign(nob,vector<double>());
        for(int nb=0; nb<nob; ++nb)
        {
            sixdof_obj_cfd *o = b->object(nb);
            const int n = o->amr_tricount();
            vector<double> &v = last[nb];
            v.resize(9*n);
            size_t c = 0;
            for(int d=0; d<3; ++d)
            {
                double **t = o->amr_tri(d);
                for(int q=0; q<n; ++q)
                for(int e=0; e<3; ++e)
                v[c++] = t[q][e];
            }
        }
    }
    return same;
}

void cfd_amr_6dof::body(lexer *p, fdm *a, ghostcell *pgc)
{
    // several bodies: fb of all bodies
    if(fball!=nullptr)
    {
        for(int nb=0; nb<b->objects(); ++nb)
        {
        moved(nb,-1.0);
        cast(p,a,pgc,nb);
        store(p,a,nb);
        combine(p,a,nb);
        }
        
        finish(p,a,pgc);
        pgc->solid_forcing_flag_update(p,a);
        pgc->gcdf_update(p,a);
        return;
    }
    
    unchanged(-1.0);

    for(int nb=0; nb<b->objects(); ++nb)
    cast(p,a,pgc,nb);

    pgc->solid_forcing_flag_update(p,a);
    pgc->gcdf_update(p,a);
}

// stage iter: the bodies are at the time of the stage (moved on level 0)
void cfd_amr_6dof::start_cfd(lexer *p, fdm *a, ghostcell *pgc, int iter, field &uvel, field &vvel, field &wvel,
                             field &fx, field &fy, field &fz, bool finalize)
{
    reset_fbh(p,a,pgc);
    
    // several bodies: each with its own level set (start_bodies)
    if(fball!=nullptr)
    {
    start_bodies(p,a,pgc,iter,uvel,vvel,wvel,fx,fy,fz,finalize);
    return;
    }

    // the bodies with the triangles of the last cast keep fb (a cast of one body overwrites fb, so
    // all are cast again when one of them has moved); G 42: up to that many cells off the last cast
    const bool keep = unchanged(MAX(p->G42,0.0)*p->DXM);

    for(int nb=0; nb<b->objects(); ++nb)
    {
        if(!keep)
        cast(p,a,pgc,nb);
        force(p,a,pgc,nb,iter,uvel,vvel,wvel,fx,fy,fz);
    }

    pgc->solid_forcing_flag_update(p,a);
    pgc->gcdf_update(p,a);
}

// several bodies: a body that has not moved more than G 42 cells since its last cast on this grid
// takes the level set of that cast (fbc), the others are cast; the forcing of a body with its own
// level set; fb of all bodies after the stage
void cfd_amr_6dof::start_bodies(lexer *p, fdm *a, ghostcell *pgc, int iter, field &uvel, field &vvel, field &wvel,
                                field &fx, field &fy, field &fz, bool finalize)
{
    const double tol = MAX(p->G42,0.0)*p->DXM;
    
    for(int nb=0; nb<b->objects(); ++nb)
    {
        if(!moved(nb,tol) && fbc[nb].valid)
        restore(p,a,pgc,nb);
        
        else
        {
        cast(p,a,pgc,nb);
        store(p,a,nb);
        }
        
        force(p,a,pgc,nb,iter,uvel,vvel,wvel,fx,fy,fz);
        
        combine(p,a,nb);
    }
    
    finish(p,a,pgc);
    
    pgc->solid_forcing_flag_update(p,a);
    pgc->gcdf_update(p,a);
}

// true if a vertex of body nb is further than tol from its triangles at the last cast (tol < 0:
// always); then the current triangles are kept for the next check
bool cfd_amr_6dof::moved(int nb, double tol)
{
    sixdof_obj_cfd *o = b->object(nb);
    const int n = o->amr_tricount();
    vector<double> &v = lastb[nb];
    
    bool same = (tol>=0.0 && (int)v.size()==9*n);
    
    size_t c = 0;
    for(int d=0; d<3 && same; ++d)
    {
        double **t = o->amr_tri(d);
        for(int q=0; q<n && same; ++q)
        for(int e=0; e<3; ++e)
        {
            if(!(fabs(v[c]-t[q][e])<=tol))
            same = false;
            ++c;
        }
    }
    
    if(!same)
    {
        v.resize(9*n);
        c = 0;
        for(int d=0; d<3; ++d)
        {
            double **t = o->amr_tri(d);
            for(int q=0; q<n; ++q)
            for(int e=0; e<3; ++e)
            v[c++] = t[q][e];
        }
    }
    
    return !same;
}

// the level set of body nb in its box after its cast
void cfd_amr_6dof::store(lexer *p, fdm *a, int nb)
{
    fbcache &f = fbc[nb];
    
    for(int q=0; q<6; ++q)
    f.box[q] = box[q];
    
    f.v.clear();
    for(i=box[0]; i<box[1]; ++i)
    for(j=box[2]; j<box[3]; ++j)
    for(k=box[4]; k<box[5]; ++k)
    f.v.push_back(a->fb(i,j,k));
    
    f.valid = true;
}

// fb of body nb from its last cast: 1.0e8 outside its box (as ray_cast_amr)
void cfd_amr_6dof::restore(lexer *p, fdm *a, ghostcell *pgc, int nb)
{
    const int m = p->margin;
    const int jm = (p->j_dir==1) ? m : 0;
    const fbcache &f = fbc[nb];
    
    for(i=-m; i<p->knox+m; ++i)
    for(j=-jm; j<p->knoy+jm; ++j)
    for(k=-m; k<p->knoz+m; ++k)
    a->fb(i,j,k)=1.0e8;
    
    size_t c = 0;
    for(i=f.box[0]; i<f.box[1]; ++i)
    for(j=f.box[2]; j<f.box[3]; ++j)
    for(k=f.box[4]; k<f.box[5]; ++k)
    a->fb(i,j,k) = f.v[c++];
    
    pgc->start4a(p,a->fb,50);
}

// fb of all bodies: the minimum of their level sets
void cfd_amr_6dof::combine(lexer *p, fdm *a, int nb)
{
    if(nb==0)
    std::copy(a->fb.begin(), a->fb.end(), fball->data());
    
    if(nb>0)
    for(int q=0; q<(p->imax*p->jmax*p->kmax); ++q)
    fball->data()[q] = MIN(fball->data()[q], a->fb.data()[q]);
}

void cfd_amr_6dof::finish(lexer *p, fdm *a, ghostcell *pgc)
{
    std::copy(fball->begin(), fball->end(), a->fb.data());
    
    pgc->start4a(p,a->fb,50);
}

// the forcing of body nb with the weights of the solids of this grid
void cfd_amr_6dof::force(lexer *p, fdm *a, ghostcell *pgc, int nb, int iter, field &uvel, field &vvel, field &wvel,
                         field &fx, field &fy, field &fz)
{
    sixdof_obj_cfd *o = b->object(nb);
    
    for(int q=0; q<3; ++q)
    o->sfh[q] = sf[q];
    
    o->update_forcing(p,a,pgc,uvel,vvel,wvel,fx,fy,fz,iter);
    
    for(int q=0; q<3; ++q)
    o->sfh[q] = nullptr;
}

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

#include "fnpf_fsfbc_wd.h"
#include "lexer.h"
#include "fdm_fnpf.h"
#include "ghostcell.h"
#include "fnpf_coastline.h"

void fnpf_fsfbc_wd::coast_cache(lexer *p, fdm_fnpf *c)
{
    cz3.clear(); cz4.clear(); cz5.clear(); czdry.clear();
    
    // fac = 1 here (only used when !(I30==1 && count==0)), so the bands and
    // factors are exactly those of the loops below
    SLICELOOP4
    {
        const double cl = c->coastline(i,j);
        
        if(cl>=0.0)
        {
            if(cl<dist3) cz3.push_back({i,j,rb3(p,cl)});
            if(cl<dist4) cz4.push_back({i,j,rb4(p,cl)});
            if(cl<dist5) cz5.push_back({i,j,rb5(p,cl)});
        }
        
        if(cl<0.0)
        czdry.push_back({i,j,0.0});
    }
    
    cz_built=true;
}

void fnpf_fsfbc_wd::coastline_eta(lexer *p, fdm_fnpf *c, ghostcell *pgc, slice &f)
{
    // A343 2 (full runup/rundown): no coastline damping. Dry cells are held
    // by wetdry_dynamic, dry columns get Fi=0 from the Laplace solve, Fz=0
    // from fsfwvel and U=V=W=0 from velcalc_sig
    if(p->A343==2)
    return;
    
    double fac=1.0;

    if((p->A347==1 || p->A347==2) && coast_cached(p))
    {
        if(!cz_built) coast_cache(p,c);
        
        for(const coast_cell &q : cz3)
        {
            i=q.i; j=q.j;
            
            // A343 3: temporarily dry cells are held by wetdry, not relaxed
            if(p->A343==3 && p->wet[IJ]==0)
            continue;
            
            f(i,j) = q.r*f(i,j);
        }
        
        return;
    }

    if(p->A347==1 || p->A347==2)
    SLICELOOP4
    {
        if(p->I30==1 && p->count==0)
            fac=p->A349;

        if(c->coastline(i,j)>=0.0)
        {
            db = c->coastline(i,j);

            if(db<fac*dist3)
            {
                f(i,j) = rb3(p,db)*f(i,j);
            }
        }
    }
}

void fnpf_fsfbc_wd::coastline_fi(lexer *p, fdm_fnpf *c, ghostcell *pgc, slice &f)
{
    // A343 2 (full runup/rundown): no coastline damping. Dry cells are held
    // by wetdry_dynamic, dry columns get Fi=0 from the Laplace solve, Fz=0
    // from fsfwvel and U=V=W=0 from velcalc_sig
    if(p->A343==2)
    return;
    
    double fac=1.0;

    if((p->A347==1 || p->A347==3) && coast_cached(p))
    {
        if(!cz_built) coast_cache(p,c);
        
        for(const coast_cell &q : cz4)
        {
            i=q.i; j=q.j;
            
            if(p->A343==3 && p->wet[IJ]==0)
            continue;
            
            f(i,j) = q.r*f(i,j);
        }
        
        if(p->A343>=1)
        for(const coast_cell &q : czdry)
        f(q.i,q.j) = 0.0;
        
        return;
    }

    if(p->A347==1 || p->A347==3 || (p->I30==1 && p->count==0))
    SLICELOOP4
    {
        if(p->I30==1 && p->count==0)
        fac=p->A349;

        if(c->coastline(i,j)>=0.0)
        {
            db = c->coastline(i,j);

            if(db<fac*dist4)
            {
                f(i,j) = rb4(p,db)*f(i,j);
            }
        }

        if(c->coastline(i,j)<0.0 && p->A343>=1)
            f(i,j)=0.0;
    }
}

void fnpf_fsfbc_wd::coastline_vel(lexer *p, fdm_fnpf *c, ghostcell *pgc, double *F)
{
    // A343 2 (full runup/rundown): no coastline damping. Dry cells are held
    // by wetdry_dynamic, dry columns get Fi=0 from the Laplace solve, Fz=0
    // from fsfwvel and U=V=W=0 from velcalc_sig
    if(p->A343==2)
    return;
    
    double fac=1.0;

    if(coast_cached(p))
    {
        if(!cz_built) coast_cache(p,c);
        
        for(const coast_cell &q : cz5)
        {
            i=q.i; j=q.j;
            FKLOOP
            F[FIJK] = q.r*F[FIJK];
        }
        
        if(p->A343>=1)
        for(const coast_cell &q : czdry)
        {
            i=q.i; j=q.j;
            FKLOOP
            F[FIJK]=0.0;
        }
        
        return;
    }

    SLICELOOP4
    {
        if(p->I30==1 && p->count==0)
        fac=p->A349;

        if(c->coastline(i,j)>=0.0)
        {
            db = c->coastline(i,j);

            if(db<fac*dist5)
            FKLOOP
            F[FIJK] = rb5(p,db)*F[FIJK];
        }

        if(c->coastline(i,j)<0.0 && p->A343>=1)
        FKLOOP
        F[FIJK]=0.0;
    }
}

void fnpf_fsfbc_wd::coastline_fi_ini(lexer *p, fdm_fnpf *c, ghostcell *pgc, slice &f)
{
    // A343 2 (full runup/rundown): no coastline damping. Dry cells are held
    // by wetdry_dynamic, dry columns get Fi=0 from the Laplace solve, Fz=0
    // from fsfwvel and U=V=W=0 from velcalc_sig
    if(p->A343==2)
    return;
    
    double fac=1.0;

    if(p->A347==1 || p->A347==3 || (p->I30==1 && p->count==0))
    FLOOP
    {
        if(p->I30==1 && p->count==0)
        fac=p->A349;

        if(c->coastline(i,j)>=0.0)
        {
            db = c->coastline(i,j);

            if(db<fac*dist4)
            {
                c->Fi[FIJK] = rb4(p,db)*c->Fi[FIJK];
            }
        }

        if(c->coastline(i,j)<0.0 && p->A343>=1)
        c->Fi[FIJK]=0.0;
    }
}

void fnpf_fsfbc_wd::coastline_Fz(lexer *p, fdm_fnpf *c, ghostcell *pgc, slice &f)
{
    // A343 2 (full runup/rundown): no coastline damping. Dry cells are held
    // by wetdry_dynamic, dry columns get Fi=0 from the Laplace solve, Fz=0
    // from fsfwvel and U=V=W=0 from velcalc_sig
    if(p->A343==2)
    return;
    
    double fac=1.0;

    if(coast_cached(p))
    {
        if(!cz_built) coast_cache(p,c);
        
        for(const coast_cell &q : cz5)
        f(q.i,q.j) = q.r*f(q.i,q.j);
        
        if(p->A343>=1)
        for(const coast_cell &q : czdry)
        f(q.i,q.j) = 0.0;
        
        return;
    }

    SLICELOOP4
    {
        if(p->I30==1 && p->count==0)
        fac=p->A349;

        if(c->coastline(i,j)>=0.0)
        {
            db = c->coastline(i,j);

            if(db<fac*dist5)
            {
                f(i,j) = rb5(p,db)*f(i,j);
            }
        }

        if(c->coastline(i,j)<0.0 && p->A343>=1)
        f(i,j)=0.0;
    }
}

double fnpf_fsfbc_wd::rb3(lexer *p, double x)
{
    double r=0.0;
    double fac=1.0;

    if(p->I30==1 && p->count==0)
    fac=p->A349;

    x=(fac*dist3-fabs(x))/(fac*dist3);
    x=MAX(x,0.0);

    r = 1.0 - (exp(pow(x,3.5))-1.0)/(EE-1.0);

    return r;
}

double fnpf_fsfbc_wd::rb4(lexer *p, double x)
{
    double r=0.0;
    double fac=1.0;

    if(p->I30==1 && p->count==0)
    fac=p->A349;

    x=(fac*dist4-fabs(x))/(fac*dist4);
    x=MAX(x,0.0);

    r = 1.0 - (exp(pow(x,3.5))-1.0)/(EE-1.0);

    return r;
}

double fnpf_fsfbc_wd::rb5(lexer *p, double x)
{
    double r=0.0;
    double fac=1.0;

    if(p->I30==1 && p->count==0)
    fac=p->A349;

    x=(fac*dist5-fabs(x))/(fac*dist5);
    x=MAX(x,0.0);

    r = 1.0 - (exp(pow(x,1.5))-1.0)/(EE-1.0);

    return r;
}

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

#include"patchBC.h"
#include"patchBC_codes.h"
#include"patch_obj.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include<vector>
#include<iomanip>

patchBC::patchBC(lexer *p, ghostcell *pgc) 
{
}

patchBC::~patchBC()
{
}

void patchBC::patchBC_ini(lexer *p, ghostcell *pgc)
{
    patch_setup(p,pgc);
    
    // faces: the wall faces (21, 22) of gcb4 inside a patch geometry, each face in one patch
    vector<vector<int>> list(obj_count);
    
    for(n=0;n<p->gcb4_count;++n)
    {
    i  = p->gcb4[n][0];
    j  = p->gcb4[n][1];
    k  = p->gcb4[n][2];
    int cs = p->gcb4[n][3];
    int bc = p->gcb4[n][4];
    
        if(cs<1 || cs>6 || (bc!=21 && bc!=22))
        continue;
    
    int qq = patch_find(p,i,j,k,cs,false);
    
        if(qq>=0)
        list[qq].insert(list[qq].end(),{i,j,k,cs,n});
    }
    
    for(int qq=0;qq<obj_count;++qq)
    patch_faces(p,qq,list[qq].data(),int(list[qq].size())/5);
    
    // boundary codes: gcb1-3 are copies of gcb4 with the same index (grid_helper::fillgcb1-3)
    int misaligned = (p->gcb1_count!=p->gcb4_count || p->gcb2_count!=p->gcb4_count || p->gcb3_count!=p->gcb4_count) ? 1 : 0;
    
    if(pgc->globalimax(misaligned)>0)
    {
    patch_error(p,pgc,"gcb1-3 are not aligned with gcb4",-1);
    pgc->final(true);
    }
    
    for(int qq=0;qq<obj_count;++qq)
    for(n=0;n<patch[qq]->gcb_count;++n)
    {
    int q = patch[qq]->gcb[n][4];
    int flag = patch[qq]->gcb_flag;
    
    p->gcb4[q][4] = flag;
    
        if(p->gcb1[q][4]==21 || p->gcb1[q][4]==22)
        p->gcb1[q][4] = flag;
        
        if(p->gcb2[q][4]==21 || p->gcb2[q][4]==22)
        p->gcb2[q][4] = flag;
        
        if(p->gcb3[q][4]==21 || p->gcb3[q][4]==22)
        p->gcb3[q][4] = flag;
    }
    
    patch_check(p,pgc);
} 

void patchBC::ghost(int cs, int q, int comp, int &ii, int &jj, int &kk)
{
    int axis = (cs==1||cs==4) ? 0 : ((cs==2||cs==3) ? 1 : 2);
    int sgn  = (cs==1||cs==3||cs==5) ? -1 : 1;
    int off  = (comp==axis && sgn>0) ? q-1 : q;
    
    ii=i;
    jj=j;
    kk=k;
    
    if(axis==0)
    ii += sgn*off;
    
    if(axis==1)
    jj += sgn*off;
    
    if(axis==2)
    kk += sgn*off;
}

void patchBC::patchBC_discharge(lexer *p, fdm *a, ghostcell *pgc)
{
    patch_hydrograph(p);
    
    for(int qq=0;qq<obj_count;++qq)
    {
    patch_obj *pt = patch[qq];
    double A=0.0, Qm=0.0, zsum=0.0;
    int hcount=0;
    
        for(n=0;n<pt->gcb_count;++n)
        {
        i  = pt->gcb[n][0];
        j  = pt->gcb[n][1];
        k  = pt->gcb[n][2];
        int cs = pt->gcb[n][3];
        
        double area=0.0, un=0.0;
        
            if(cs==1 || cs==4)
            area = p->DYN[JP]*p->DZN[KP];
            
            if(cs==2 || cs==3)
            area = p->DXN[IP]*p->DZN[KP];
            
            if(cs==5 || cs==6)
            area = p->DXN[IP]*p->DYN[JP];
            
            // velocity at the boundary face
            if(cs==1)
            un = a->u(i-1,j,k);
            
            if(cs==4)
            un = a->u(i,j,k);
            
            if(cs==3)
            un = a->v(i,j-1,k);
            
            if(cs==2)
            un = a->v(i,j,k);
            
            if(cs==5)
            un = a->w(i,j,k-1);
            
            if(cs==6)
            un = a->w(i,j,k);
            
        // wetted part of the face
        area *= wetfrac(a->phi(i,j,k),p->DZN[KP]);
        
        A  += area;
        Qm += area*un;
            
            // water level at the side faces
            if(cs<=4 && a->phi(i,j,k)>=0.0 && a->phi(i,j,k+1)<0.0)
            {
            zsum += p->ZP[KP] + a->phi(i,j,k)*p->DZP[KP]/(a->phi(i,j,k)-a->phi(i,j,k+1));
            ++hcount;
            }
        }
        
    A  = pgc->globalsum(A);
    Qm = pgc->globalsum(Qm);
    zsum = pgc->globalsum(zsum);
    hcount = pgc->globalisum(hcount);
    
    pt->A0 = A;
    pt->Q0 = Qm;
    pt->U0 = Qm/(A>1.0e-20?A:1.0e20);
    pt->h0 = hcount>0 ? zsum/double(hcount) : 0.0;
    
        // inlet: normal velocity from the discharge over the wetted area
        if(pt->Q_flag==1)
        pt->Uq = pt->Q/(A>1.0e-20?A:1.0e20);
    
        if(p->mpirank==0 && (p->count%p->P12==0))
        {
        cout<<"patchBC ID: "<<pt->ID<<(pt->kind==PATCH_INLET?" in ":" out")<<" | Q: "<<setprecision(5)<<Qm;
        
        if(pt->Q_flag==1)
        cout<<" Qset: "<<pt->Q<<" Uq: "<<pt->Uq;
        
        cout<<" U: "<<pt->U0<<" A: "<<A<<" h: "<<pt->h0<<endl;
        }
    }
}

void patchBC::patchBC_ioflow(lexer *p, fdm *a, ghostcell *pgc, field &u, field &v, field &w)
{
    int ii,jj,kk;
    double val[3];
    
    // inlets: velocity at the boundary face and in the ghost cells, the interior is left alone
    for(int qq=0;qq<obj_count;++qq)
    if(patch[qq]->kind==PATCH_INLET)
    {
    patch_obj *pt = patch[qq];
    
    double Un = (pt->Q_flag==1) ? pt->Uq : ((pt->Uio_flag==1) ? pt->Uio : 0.0);
    
        for(n=0;n<pt->gcb_count;++n)
        {
        i  = pt->gcb[n][0];
        j  = pt->gcb[n][1];
        k  = pt->gcb[n][2];
        int cs = pt->gcb[n][3];
        
        pt->velocity(cs,Un,val[0],val[1],val[2]);
        
            // discharge: the flow enters through the wetted part of the patch, closed for the gas
            if(pt->Q_flag==1 && wetfrac(a->phi(i,j,k),p->DZN[KP])<=0.0)
            val[0]=val[1]=val[2]=0.0;
            
            for(int q=1;q<=3;++q)
            {
            ghost(cs,q,0,ii,jj,kk);
            u(ii,jj,kk) = val[0];
            
            ghost(cs,q,1,ii,jj,kk);
            v(ii,jj,kk) = val[1];
            
            ghost(cs,q,2,ii,jj,kk);
            w(ii,jj,kk) = val[2];
            }
        }
    }
} 

void patchBC::patchBC_rkioflow(lexer *p, fdm *a, ghostcell *pgc, field &u, field &v, field &w)
{
    int ii,jj,kk;
    
    // inlets: the RK stage fields take the boundary values
    for(int qq=0;qq<obj_count;++qq)
    if(patch[qq]->kind==PATCH_INLET)
    for(n=0;n<patch[qq]->gcb_count;++n)
    {
    i  = patch[qq]->gcb[n][0];
    j  = patch[qq]->gcb[n][1];
    k  = patch[qq]->gcb[n][2];
    int cs = patch[qq]->gcb[n][3];
    
        for(int q=1;q<=3;++q)
        {
        ghost(cs,q,0,ii,jj,kk);
        u(ii,jj,kk) = a->u(ii,jj,kk);
        
        ghost(cs,q,1,ii,jj,kk);
        v(ii,jj,kk) = a->v(ii,jj,kk);
        
        ghost(cs,q,2,ii,jj,kk);
        w(ii,jj,kk) = a->w(ii,jj,kk);
        }
    }
}

void patchBC::patchBC_pressure(lexer *p, fdm *a, ghostcell *pgc, field &press)
{
    int ii,jj,kk;
    
    // hydrostatic part with a free surface method (level set or VOF), also when the free surface
    // lies above the domain; single phase without free surface: the patch pressure only
    if(fsf_domain<0)
    fsf_domain = (p->F30>0 || p->F80>0) ? 1 : 0;
    
    const double g = fabs(p->W22);
    
    // outlets: patch pressure (B 412, default 0) at the free surface plus the hydrostatic pressure
    // below the water level (B 413/B 422, otherwise the local level from the level set)
    for(int qq=0;qq<obj_count;++qq)
    if(patch[qq]->kind==PATCH_OUTLET)
    {
    patch_obj *pt = patch[qq];
    
        for(n=0;n<pt->gcb_count;++n)
        {
        i  = pt->gcb[n][0];
        j  = pt->gcb[n][1];
        k  = pt->gcb[n][2];
        int cs = pt->gcb[n][3];
        
        double eta = (pt->waterlevel_flag==1) ? pt->waterlevel : p->ZP[KP] + a->phi(i,j,k);
        
            for(int q=1;q<=3;++q)
            {
            ghost(cs,q,-1,ii,jj,kk);
            
            double pval = pt->pressure;
            
                if(fsf_domain==1)
                pval += a->ro(i,j,k)*g*(eta - p->ZP[kk+marge]);
            
            press(ii,jj,kk) = pval;
            }
        }
    }
} 

void patchBC::patchBC_waterlevel(lexer *p, fdm *a, ghostcell *pgc, field &phi)
{
    int ii,jj,kk;
    
    patch_hydrograph(p);
    
    // B 413/B 422: level set in the ghost cells from the water level
    for(int qq=0;qq<obj_count;++qq)
    if(patch[qq]->waterlevel_flag==1)
    for(n=0;n<patch[qq]->gcb_count;++n)
    {
    i  = patch[qq]->gcb[n][0];
    j  = patch[qq]->gcb[n][1];
    k  = patch[qq]->gcb[n][2];
    int cs = patch[qq]->gcb[n][3];
    
        for(int q=1;q<=3;++q)
        {
        ghost(cs,q,-1,ii,jj,kk);
        phi(ii,jj,kk) = patch[qq]->waterlevel - p->ZP[kk+marge];
        }
    }
} 

void patchBC::patchBC_ioflow2D(lexer*, ghostcell*, slice&, slice&, slice&, slice&)
{
}

void patchBC::patchBC_discharge2D(lexer*, fdm2D*, ghostcell*, slice&, slice&, slice&, slice&)
{
}

void patchBC::patchBC_waterlevel2D(lexer*, fdm2D*, ghostcell*, slice&)
{
}

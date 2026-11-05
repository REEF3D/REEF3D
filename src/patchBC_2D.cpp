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

#include"patchBC_2D.h"
#include"patchBC_codes.h"
#include"patch_obj.h"
#include"lexer.h"
#include"fdm2D.h"
#include"ghostcell.h"
#include"slice.h"
#include<vector>
#include<iomanip>
#include<cstdlib>

patchBC_2D::patchBC_2D(lexer *p, ghostcell *pgc) 
{
}

patchBC_2D::~patchBC_2D()
{
}

void patchBC_2D::patchBC_ini(lexer *p, ghostcell *pgc)
{
    patch_setup(p,pgc);
    
    // faces: the wall faces (21, 22) of gcbsl4 inside a patch geometry, each face in one patch
    vector<vector<int>> list(obj_count);
    
    for(n=0;n<p->gcbsl4_count;++n)
    {
    i  = p->gcbsl4[n][0];
    j  = p->gcbsl4[n][1];
    int cs = p->gcbsl4[n][3];
    int bc = p->gcbsl4[n][4];
    
        if(cs<1 || cs>4 || (bc!=21 && bc!=22))
        continue;
    
    int qq = patch_find(p,i,j,0,cs,true);
    
        if(qq>=0)
        list[qq].insert(list[qq].end(),{i,j,0,cs,n});
    }
    
    for(int qq=0;qq<obj_count;++qq)
    patch_faces(p,qq,list[qq].data(),int(list[qq].size())/5);
    
    // boundary codes
    for(int qq=0;qq<obj_count;++qq)
    for(n=0;n<patch[qq]->gcb_count;++n)
    p->gcbsl4[patch[qq]->gcb[n][4]][4] = patch[qq]->gcb_flag;
    
    // staggered slices, generated independently of gcbsl4
    for(n=0;n<p->gcbsl1_count;++n)
    if(p->gcbsl1[n][4]==21 || p->gcbsl1[n][4]==22)
    {
    int flag = staggered_flag(p->gcbsl1[n][0],p->gcbsl1[n][1],p->gcbsl1[n][3]);
    
        if(flag>0)
        p->gcbsl1[n][4] = flag;
    }
    
    for(n=0;n<p->gcbsl2_count;++n)
    if(p->gcbsl2[n][4]==21 || p->gcbsl2[n][4]==22)
    {
    int flag = staggered_flag(p->gcbsl2[n][0],p->gcbsl2[n][1],p->gcbsl2[n][3]);
    
        if(flag>0)
        p->gcbsl2[n][4] = flag;
    }
    
    patch_check(p,pgc);
} 

int patchBC_2D::staggered_flag(int ii, int jj, int cs)
{
    for(int qq=0;qq<obj_count;++qq)
    for(int m=0;m<patch[qq]->gcb_count;++m)
    if(patch[qq]->gcb[m][3]==cs)
    {
        if((cs==1 || cs==4) && patch[qq]->gcb[m][1]==jj && abs(patch[qq]->gcb[m][0]-ii)<=1)
        return patch[qq]->gcb_flag;
        
        if((cs==2 || cs==3) && patch[qq]->gcb[m][0]==ii && abs(patch[qq]->gcb[m][1]-jj)<=1)
        return patch[qq]->gcb_flag;
    }
    
    return 0;
}

void patchBC_2D::ghost(int cs, int q, int &ii, int &jj)
{
    ii=i;
    jj=j;
    
    if(cs==1)
    ii-=q;
    
    if(cs==4)
    ii+=q;
    
    if(cs==3)
    jj-=q;
    
    if(cs==2)
    jj+=q;
}

void patchBC_2D::patchBC_discharge2D(lexer *p, fdm2D *b, ghostcell *pgc, slice &P, slice &Q, slice &eta, slice &bed)
{
    int ii,jj;
    
    patch_hydrograph(p);
    
    for(int qq=0;qq<obj_count;++qq)
    {
    patch_obj *pt = patch[qq];
    double A=0.0, Qm=0.0, hsum=0.0;
    int hcount=0;
    
        for(n=0;n<pt->gcb_count;++n)
        {
        i  = pt->gcb[n][0];
        j  = pt->gcb[n][1];
        int cs = pt->gcb[n][3];
        
            if(p->wet[IJ]==0)
            continue;
        
        ghost(cs,1,ii,jj);
        
        // face width, water depth of the cell, normal velocity at the boundary face
        double dl = (cs==1 || cs==4) ? p->DYN[JP] : p->DXN[IP];
        double un = (cs==1 || cs==4) ? 0.5*(b->U(i,j)+b->U(ii,jj)) : 0.5*(b->V(i,j)+b->V(ii,jj));
        
        A  += dl*b->hp(i,j);
        Qm += dl*b->hp(i,j)*un;
        
        hsum += b->eta(i,j) + p->wd;
        ++hcount;
        }
        
    A  = pgc->globalsum(A);
    Qm = pgc->globalsum(Qm);
    hsum = pgc->globalsum(hsum);
    hcount = pgc->globalisum(hcount);
    
    pt->A0 = A;
    pt->Q0 = Qm;
    pt->U0 = Qm/(A>1.0e-20?A:1.0e20);
    pt->h0 = hcount>0 ? hsum/double(hcount) : 0.0;
    
        // inlet: normal velocity from the discharge over the wetted width x depth
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

void patchBC_2D::patchBC_ioflow2D(lexer *p, ghostcell *pgc, slice &U, slice &V, slice &bed, slice &eta)
{
    int ii,jj;
    double val[3];
    
    // inlets: cell centred velocities in the ghost cells, the interior is left alone;
    // outlets: zero gradient (sflow_momentum_func::vel_bc)
    for(int qq=0;qq<obj_count;++qq)
    if(patch[qq]->kind==PATCH_INLET)
    {
    patch_obj *pt = patch[qq];
    
    double Un = (pt->Q_flag==1) ? pt->Uq : ((pt->Uio_flag==1) ? pt->Uio : 0.0);
    
        for(n=0;n<pt->gcb_count;++n)
        {
        i  = pt->gcb[n][0];
        j  = pt->gcb[n][1];
        int cs = pt->gcb[n][3];
        
        pt->velocity(cs,Un,val[0],val[1],val[2]);
        
            if(pt->Q_flag==1 && p->wet[IJ]==0)
            val[0]=val[1]=0.0;
        
            for(int q=1;q<=3;++q)
            {
            ghost(cs,q,ii,jj);
            
            U(ii,jj) = val[0];
            V(ii,jj) = val[1];
            }
        }
    }
}

void patchBC_2D::patchBC_waterlevel2D(lexer *p, fdm2D *b, ghostcell *pgc, slice &eta)
{
    int ii,jj;
    
    patch_hydrograph(p);
    
    // B 413/B 422: water level in the ghost cells
    for(int qq=0;qq<obj_count;++qq)
    if(patch[qq]->waterlevel_flag==1)
    for(n=0;n<patch[qq]->gcb_count;++n)
    {
    i  = patch[qq]->gcb[n][0];
    j  = patch[qq]->gcb[n][1];
    int cs = patch[qq]->gcb[n][3];
    
        for(int q=1;q<=3;++q)
        {
        ghost(cs,q,ii,jj);
        
        eta(ii,jj)   = patch[qq]->waterlevel - p->wd;
        b->hp(ii,jj) = MAX(patch[qq]->waterlevel - b->bed(i,j), 0.0);
        }
    }
}

void patchBC_2D::patchBC_ioflow(lexer*, fdm*, ghostcell*, field&, field&, field&)
{
} 

void patchBC_2D::patchBC_rkioflow(lexer*, fdm*, ghostcell*, field&, field&, field&)
{
}

void patchBC_2D::patchBC_discharge(lexer*, fdm*, ghostcell*)
{
}

void patchBC_2D::patchBC_pressure(lexer*, fdm*, ghostcell*, field&)
{
} 

void patchBC_2D::patchBC_waterlevel(lexer*, fdm*, ghostcell*, field&)
{
}

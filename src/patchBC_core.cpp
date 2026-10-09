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

#include"patchBC_core.h"
#include"patchBC_codes.h"
#include"patch_obj.h"
#include"lexer.h"
#include"ghostcell.h"
#include<vector>
#include<algorithm>
#include<cmath>
#include<cstdio>
#include<fstream>
#include<iostream>

void patchBC_core::patch_setup(lexer *p, ghostcell *pgc)
{
    // patches: the IDs of the geometries B 440-442, in input order
    vector<int> ids;
    
    auto add = [&](int id)
    {
        if(std::find(ids.begin(),ids.end(),id)==ids.end())
        ids.push_back(id);
    };
    
    for(int n=0;n<p->B440;++n)
    add(p->B440_ID[n]);
    
    for(int n=0;n<p->B441;++n)
    add(p->B441_ID[n]);
    
    for(int n=0;n<p->B442;++n)
    add(p->B442_ID[n]);
    
    obj_count = int(ids.size());
    
    patch = new patch_obj*[obj_count];
    
    for(int qq=0; qq<obj_count;++qq)
    patch[qq] = new patch_obj(p,ids[qq]);
    
    
    // settings
    vector<int> nQ(obj_count,0), nhQ(obj_count,0), nUio(obj_count,0), nvel(obj_count,0);
    vector<int> nh(obj_count,0), nhFSF(obj_count,0);
    
    auto index = [&](int id, const char *key)
    {
        int qq = patch_index(id);
        
        if(qq<0 && p->mpirank==0)
        cout<<"patchBC warning: "<<key<<" for ID "<<id<<" without a patch geometry (B 440-442), ignored"<<endl;
        
        return qq;
    };
    
    // discharge
    for(int n=0;n<p->B411;++n)
    {
    int qq = index(p->B411_ID[n],"B 411");
    if(qq<0) continue;
    
    patch[qq]->Q_flag=1;
    patch[qq]->Q=p->B411_Q[n];
    ++nQ[qq];
    }
    
    // pressure
    for(int n=0;n<p->B412;++n)
    {
    int qq = index(p->B412_ID[n],"B 412");
    if(qq<0) continue;
    
    patch[qq]->pressure_flag=1;
    patch[qq]->pressure=p->B412_pressBC[n];
    }
    
    // water level
    for(int n=0;n<p->B413;++n)
    {
    int qq = index(p->B413_ID[n],"B 413");
    if(qq<0) continue;
    
    patch[qq]->waterlevel_flag=1;
    patch[qq]->waterlevel=p->B413_h[n];
    ++nh[qq];
    }
    
    // velocity normal to the face
    for(int n=0;n<p->B414;++n)
    {
    int qq = index(p->B414_ID[n],"B 414");
    if(qq<0) continue;
    
    patch[qq]->Uio_flag=1;
    patch[qq]->Uio=p->B414_Uio[n];
    ++nUio[qq];
    }
    
    // velocity components
    for(int n=0;n<p->B415;++n)
    {
    int qq = index(p->B415_ID[n],"B 415");
    if(qq<0) continue;
    
    patch[qq]->velcomp_flag=1;
    patch[qq]->U=p->B415_U[n];
    patch[qq]->V=p->B415_V[n];
    patch[qq]->W=p->B415_W[n];
    ++nvel[qq];
    }
    
    // horizontal inflow angle, counter-clockwise from the face normal
    for(int n=0;n<p->B416;++n)
    {
    int qq = index(p->B416_ID[n],"B 416");
    if(qq<0) continue;
    
    patch[qq]->angle_flag=1;
    patch[qq]->alpha=(PI/180.0)*p->B416_alpha[n];
    
        if(fabs(p->B416_alpha[n])>=84.0)
        patch_error(p,pgc,"B 416: the inflow angle is measured from the face normal (0 deg = perpendicular inflow), |alpha| must be < 84 deg",patch[qq]->ID);
    }
    
    // flow direction vector
    for(int n=0;n<p->B417;++n)
    {
    int qq = index(p->B417_ID[n],"B 417");
    if(qq<0) continue;
    
    double norm = sqrt(p->B417_Nx[n]*p->B417_Nx[n] + p->B417_Ny[n]*p->B417_Ny[n] + p->B417_Nz[n]*p->B417_Nz[n]);
    
        if(norm<1.0e-12)
        {
        patch_error(p,pgc,"B 417: zero direction vector",patch[qq]->ID);
        continue;
        }
    
    patch[qq]->dir_flag=1;
    patch[qq]->dirx=p->B417_Nx[n]/norm;
    patch[qq]->diry=p->B417_Ny[n]/norm;
    patch[qq]->dirz=p->B417_Nz[n]/norm;
    }
    
    // SFLOW free stream outflow
    for(int n=0;n<p->B418;++n)
    {
    int qq = index(p->B418_ID[n],"B 418");
    if(qq<0) continue;
    
    patch[qq]->pio_flag=1;
    }
    
    // discharge hydrograph
    for(int n=0;n<p->B421;++n)
    {
    int qq = index(p->B421_ID[n],"B 421");
    if(qq<0) continue;
    
    char name[100];
    snprintf(name,sizeof(name),"hydrograph_Q_%i.dat",patch[qq]->ID);
    
    patch[qq]->hydroQ_flag=1;
    patch[qq]->Q_flag=1;
    hydrograph_read(p,pgc,name,patch[qq]->ID,patch[qq]->hydroQ,patch[qq]->hydroQ_count);
    ++nhQ[qq];
    }
    
    // water level hydrograph
    for(int n=0;n<p->B422;++n)
    {
    int qq = index(p->B422_ID[n],"B 422");
    if(qq<0) continue;
    
    char name[100];
    snprintf(name,sizeof(name),"hydrograph_FSF_%i.dat",patch[qq]->ID);
    
    patch[qq]->hydroFSF_flag=1;
    patch[qq]->waterlevel_flag=1;
    hydrograph_read(p,pgc,name,patch[qq]->ID,patch[qq]->hydroFSF,patch[qq]->hydroFSF_count);
    ++nhFSF[qq];
    }
    
    
    // kind: inlet (velocity prescribed) or outlet (pressure / water level / free), never both
    for(int qq=0;qq<obj_count;++qq)
    {
    patch_obj *pt = patch[qq];
    
        if(nQ[qq]>0 && nhQ[qq]>0)
        patch_error(p,pgc,"B 411 and B 421 for the same patch",pt->ID);
        
        if(nh[qq]>0 && nhFSF[qq]>0)
        patch_error(p,pgc,"B 413 and B 422 for the same patch",pt->ID);
    
    int nin = (pt->Q_flag==1) + (pt->Uio_flag==1) + (pt->velcomp_flag==1);
    
        if(nin>1)
        patch_error(p,pgc,"more than one inflow velocity (B 411/421, B 414, B 415) for the same patch",pt->ID);
        
        if(nin>0 && (pt->pressure_flag==1 || pt->pio_flag==1))
        patch_error(p,pgc,"velocity (B 411/414/415/421) and pressure (B 412/418) for the same patch: over-specified",pt->ID);
        
        if(pt->angle_flag==1 && pt->dir_flag==1)
        patch_error(p,pgc,"B 416 and B 417 for the same patch",pt->ID);
        
        if((pt->angle_flag==1 || pt->dir_flag==1) && pt->Q_flag==0 && pt->Uio_flag==0)
        patch_error(p,pgc,"B 416/B 417 set the flow direction of B 411/421 or B 414, neither is given",pt->ID);
    
    pt->kind = (nin>0) ? PATCH_INLET : PATCH_OUTLET;
    
        if(pt->kind==PATCH_INLET)
        pt->gcb_flag = (pt->waterlevel_flag==1) ? PATCH_INLET_FSF : PATCH_INLET;
        
        if(pt->kind==PATCH_OUTLET)
        pt->gcb_flag = (pt->waterlevel_flag==1) ? PATCH_OUTLET_FSF : PATCH_OUTLET;
    }
    
    // collective: a missing hydrograph file may be seen by some ranks only
    if(pgc->globalimax(error_count)>0)
    pgc->final(true);
}

int patchBC_core::patch_find(lexer *p, int ii, int jj, int kk, int cs, bool twoD)
{
    double xf,yf,zf;
    const double tol = 1.0e-6*p->DXM;
    
    face_center(p,ii,jj,kk,cs,twoD,xf,yf,zf);
    
    auto in = [&](double val, double s, double e)
    {
        return val>=std::min(s,e)-tol && val<=std::max(s,e)+tol;
    };
    
    // line: x-y rectangle, all z
    for(int n=0;n<p->B440;++n)
    if(p->B440_face[n]==cs)
    if(in(xf,p->B440_xs[n],p->B440_xe[n]) && in(yf,p->B440_ys[n],p->B440_ye[n]))
    return patch_index(p->B440_ID[n]);
    
    // box
    for(int n=0;n<p->B441;++n)
    if(p->B441_face[n]==cs)
    if(in(xf,p->B441_xs[n],p->B441_xe[n]) && in(yf,p->B441_ys[n],p->B441_ye[n]) 
    && (twoD || in(zf,p->B441_zs[n],p->B441_ze[n])))
    return patch_index(p->B441_ID[n]);
    
    // circle: the face plane within one cell of the center, the face center within the radius
    for(int n=0;n<p->B442;++n)
    if(p->B442_face[n]==cs)
    {
    double dn,rt;
    const double zm = twoD ? 0.0 : p->B442_zm[n];
    
        if(cs==1 || cs==4)
        {
        dn = fabs(xf-p->B442_xm[n]) - p->DXN[ii+marge];
        rt = sqrt(pow(yf-p->B442_ym[n],2.0) + pow(zf-zm,2.0));
        }
        
        else
        if(cs==2 || cs==3)
        {
        dn = fabs(yf-p->B442_ym[n]) - p->DYN[jj+marge];
        rt = sqrt(pow(xf-p->B442_xm[n],2.0) + pow(zf-zm,2.0));
        }
        
        else
        {
        dn = fabs(zf-zm) - p->DZN[kk+marge];
        rt = sqrt(pow(xf-p->B442_xm[n],2.0) + pow(yf-p->B442_ym[n],2.0));
        }
        
        if(dn<=tol && rt<=p->B442_r[n]+tol)
        return patch_index(p->B442_ID[n]);
    }
    
    return -1;
}

void patchBC_core::patch_faces(lexer *p, int qq, const int *list, int count)
{
    patch[qq]->gcb_count = count;
    p->Iarray(patch[qq]->gcb,count>0?count:1,5);
    
    for(int n=0;n<count;++n)
    for(int q=0;q<5;++q)
    patch[qq]->gcb[n][q] = list[5*n+q];
}

void patchBC_core::patch_check(lexer *p, ghostcell *pgc)
{
    int err=0;
    
    for(int qq=0;qq<obj_count;++qq)
    for(int n=0;n<patch[qq]->gcb_count;++n)
    {
    int cs = patch[qq]->gcb[n][3];
    
        // horizontal angle on top / bottom faces
        if(patch[qq]->angle_flag==1 && (cs==5 || cs==6))
        err=1;
        
        // direction (nearly) parallel to the face
        if(patch[qq]->dir_flag==1)
        {
        double dn = (cs==1||cs==4) ? patch[qq]->dirx : ((cs==2||cs==3) ? patch[qq]->diry : patch[qq]->dirz);
        
            if(fabs(dn)<0.1)
            err=2;
        }
    }
    
    err = pgc->globalimax(err);
    
    if(err==1)
    patch_error(p,pgc,"B 416: horizontal inflow angle on a bottom or top face (5, 6), use B 417",-1);
    
    if(err==2)
    patch_error(p,pgc,"B 417: flow direction (nearly) parallel to a patch face",-1);
    
    if(pgc->globalimax(error_count)>0)
    pgc->final(true);
    
    for(int qq=0;qq<obj_count;++qq)
    {
    patch_obj *pt = patch[qq];
    int nf = pgc->globalisum(pt->gcb_count);
    
        if(p->mpirank==0)
        {
        cout<<"patchBC ID "<<pt->ID<<" | ";
        
            if(pt->kind==PATCH_INLET)
            {
            cout<<"inlet";
            
            if(pt->Q_flag==1 && pt->hydroQ_flag==0)
            cout<<", Q: "<<pt->Q<<" m3/s";
            
            if(pt->hydroQ_flag==1)
            cout<<", Q: hydrograph";
            
            if(pt->Uio_flag==1)
            cout<<", Un: "<<pt->Uio<<" m/s";
            
            if(pt->velcomp_flag==1)
            cout<<", U V W: "<<pt->U<<" "<<pt->V<<" "<<pt->W<<" m/s";
            
            if(pt->angle_flag==1)
            cout<<", angle: "<<pt->alpha*180.0/PI<<" deg";
            
            if(pt->dir_flag==1)
            cout<<", direction: "<<pt->dirx<<" "<<pt->diry<<" "<<pt->dirz;
            }
            
            if(pt->kind==PATCH_OUTLET)
            {
            cout<<"outlet";
            
            if(pt->pressure_flag==1)
            cout<<", p: "<<pt->pressure<<" Pa";
            }
            
            if(pt->waterlevel_flag==1 && pt->hydroFSF_flag==0)
            cout<<", water level: "<<pt->waterlevel<<" m";
            
            if(pt->hydroFSF_flag==1)
            cout<<", water level: hydrograph";
            
        cout<<" | faces: "<<nf<<endl;
        
            if(nf==0)
            cout<<"patchBC warning: ID "<<pt->ID<<" has no faces, check the B 440-442 geometry and face"<<endl;
        }
    }
}

void patchBC_core::patch_hydrograph(lexer *p)
{
    for(int qq=0;qq<obj_count;++qq)
    {
        if(patch[qq]->hydroQ_flag==1)
        patch[qq]->Q = hydrograph_ipol(p,patch[qq]->hydroQ,patch[qq]->hydroQ_count);
        
        if(patch[qq]->hydroFSF_flag==1)
        patch[qq]->waterlevel = hydrograph_ipol(p,patch[qq]->hydroFSF,patch[qq]->hydroFSF_count);
    }
}

double patchBC_core::wetfrac(double phival, double dz)
{
    return std::max(0.0, std::min(1.0, 0.5 + phival/dz));
}

void patchBC_core::patch_error(lexer *p, ghostcell *pgc, const char *msg, int ID)
{
    if(p->mpirank==0)
    {
        if(ID>=0)
        cout<<"patchBC error, ID "<<ID<<": "<<msg<<endl;
        
        if(ID<0)
        cout<<"patchBC error: "<<msg<<endl;
    }
    
    ++error_count;
}

void patchBC_core::hydrograph_read(lexer *p, ghostcell *pgc, const char *name, int ID, double **&hg, int &count)
{
    vector<double> t,val;
    double tval,vval;
    
    ifstream file(name, ios_base::in);
	
	if(!file)
	{
    patch_error(p,pgc,"hydrograph file not found",ID);
    
        if(p->mpirank==0)
        cout<<"   expected: "<<name<<endl;
    
    count=0;
    return;
	}
    
    while(file>>tval>>vval)
    {
    t.push_back(tval);
    val.push_back(vval);
    }
    
    file.close();
    
    count = int(t.size());
    
    if(count==0)
    {
    patch_error(p,pgc,"empty hydrograph file",ID);
    return;
    }
	
	p->Darray(hg,count,2);
	
	for(int n=0;n<count;++n)
	{
    hg[n][0]=t[n];
    hg[n][1]=val[n];
    }
}

double patchBC_core::hydrograph_ipol(lexer *p, double **hg, int count)
{
    if(count==0)
    return 0.0;
    
    if(p->simtime<=hg[0][0])
    return hg[0][1];
    
	if(p->simtime>=hg[count-1][0])
	return hg[count-1][1];
    
    for(int n=0;n<count-1;++n)
    if(p->simtime>=hg[n][0] && p->simtime<hg[n+1][0])
	return hg[n][1] + (hg[n+1][1]-hg[n][1])*(p->simtime-hg[n][0])/(hg[n+1][0]-hg[n][0]);
    
	return hg[count-1][1];
}

void patchBC_core::face_center(lexer *p, int ii, int jj, int kk, int cs, bool twoD, double &xf, double &yf, double &zf)
{
    xf = p->XP[ii+marge];
    yf = p->YP[jj+marge];
    zf = twoD ? 0.0 : p->ZP[kk+marge];
    
    if(cs==1)
    xf = p->XN[ii+marge];
    
    if(cs==4)
    xf = p->XN[ii+1+marge];
    
    if(cs==3)
    yf = p->YN[jj+marge];
    
    if(cs==2)
    yf = p->YN[jj+1+marge];
    
    if(cs==5 && !twoD)
    zf = p->ZN[kk+marge];
    
    if(cs==6 && !twoD)
    zf = p->ZN[kk+1+marge];
}

int patchBC_core::patch_index(int ID)
{
    for(int qq=0;qq<obj_count;++qq)
    if(patch[qq]->ID==ID)
    return qq;
    
    return -1;
}

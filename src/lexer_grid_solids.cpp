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

#include"lexer.h"
#include"ghostcell.h"
#include"geo_mesh.h"
#include"geo_raycast.h"
#include<vector>

// Solids (S) and topography (T) of the grid file, built by the geometry core for all
// hydrodynamic modules (grid format v2). Before v2 DIVEMesh wrote the signed distance
// fields solid_dist/topo_dist and the bed levels; now the triangles are ray cast here,
// on the Cartesian DIVEMesh grid of the subdomain including its ghost cells:
//
//   flag_solid, flag_topo    signed distance (negative inside), CFD and PTF
//   solidbed, topobed, bed   bed levels of the columns, all modules
//   *_gcb*_est               ghost cell estimates for the vector sizes (lexer::vecsize)
//
// The geodat bed level (G 10) is united with the entities of its role (G 9).
//
// NHFLOW with A 580 1 treats the solids as immersed solids (direct forcing); the bed
// level then only follows the topography.

namespace
{
    // signed distance field of one role
    bool build_field(lexer *p, geo_mesh &G, geo_raycast &ray, const geo_cart &gx, int role,
                     double *bedrole, vector<double> &phi)
    {
        int i,j,k,n;
        const int marge = increment::marge;
        
        const int nf = p->imax*p->jmax*p->kmax;
        const double dxm = (G.dxm>0.0) ? G.dxm : p->DXM;
        
        bool objects=false;
        
        for(const geo_object &ob : G.obj)
        if(ob.role==role && ob.te>ob.ts)
        objects=true;
        
        const bool geodat = (G.geodat==role);
        
        if(!objects && !geodat)
        return false;
        
        phi.assign(nf,0.0);
        
        if(objects)
        {
            vector<int> IO(nf,1), CL(nf,0), CR(nf,0);
            
            for(n=0; n<nf; ++n)
            phi[n] = 1.0e9;
            
            // solid: ghost cells outside the global domain count as solid for the crossing check,
            // entity faces on the domain boundary are no fluid-solid interface; topo: they count
            // as fluid, the faces on the domain boundary are interfaces (both as in DIVEMesh)
            vector<char> outside(nf,0);
            
            if(role==geo_mesh::role_solid)
            for(i=gx.is; i<gx.ie; ++i)
            for(j=gx.js; j<gx.je; ++j)
            for(k=gx.ks; k<gx.ke; ++k)
            if(i+p->origin_i<0 || i+p->origin_i>=p->gknox
            || j+p->origin_j<0 || j+p->origin_j>=p->gknoy
            || k+p->origin_k<0 || k+p->origin_k>=p->gknoz)
            {
            outside[IJK] = 1;
            IO[IJK] = -1;
            }
            
            geo_cart g = gx;
            
            // inside/outside, entity by entity
            for(const geo_object &ob : G.obj)
            if(ob.role==role && ob.te>ob.ts)
            {
                g.raymode = ob.raymode;
                
                for(int dir=0; dir<3; ++dir)
                ray.cart_io(p,dir,G.tri_x,G.tri_y,G.tri_z,ob.ts,ob.te,g,IO.data(),CL.data(),CR.data());
                
                if(ob.invert==1)
                for(n=0; n<nf; ++n)
                if(outside[n]==0)
                IO[n] = -IO[n];
            }
            
            // distance to the crossings along the grid lines, bed level of the columns
            for(const geo_object &ob : G.obj)
            if(ob.role==role && ob.te>ob.ts)
            for(int dir=0; dir<3; ++dir)
            ray.cart_dist(p,dir,G.tri_x,G.tri_y,G.tri_z,ob.ts,ob.te,g,IO.data(),phi.data(),(dir==2)?bedrole:nullptr);
            
            for(n=0; n<nf; ++n)
            {
                // outside the global domain: fluid for the ghost cell estimates (as in DIVEMesh)
                if(outside[n]==1)
                IO[n] = 1;
                
                if(IO[n]==-1)
                phi[n] = -fabs(phi[n]);
                
                else if(IO[n]==1)
                phi[n] = fabs(phi[n]);
                
                if(phi[n]>10.0*dxm)
                phi[n] = 10.0*dxm;
                
                else if(phi[n]<-10.0*dxm)
                phi[n] = -10.0*dxm;
            }
        }
        
        // geodat bed level: vertical distance, united with the entities
        if(geodat)
        for(i=gx.is; i<gx.ie; ++i)
        for(j=gx.js; j<gx.je; ++j)
        {
            const int ic = MAX(0,MIN(i,p->knox-1));
            const int jc = MAX(0,MIN(j,p->knoy-1));
            
            const double zb = p->geobed[(ic-p->imin)*p->jmax + (jc-p->jmin)];
            
            for(k=gx.ks; k<gx.ke; ++k)
            {
                const double dg = p->ZP[KP] - zb;
                
                phi[IJK] = objects ? MIN(phi[IJK],dg) : dg;
            }
        }
        
        return true;
    }
    
    // ghost cell estimates (formerly DIVEMesh solid/topo/surface::gcb_estimate)
    void estimate(lexer *p, const geo_cart &gx, const vector<double> &A, const vector<double> &B,
                  bool useB, int &gcb, int &gcbextra)
    {
        int i,j,k;
        
        const int nf = p->imax*p->jmax*p->kmax;
        
        auto fluid = [&](int q, bool strict)
        {
            if(useB)
            return A[q]>0.0 && B[q]>0.0;
            
            return strict ? A[q]>0.0 : A[q]>=0.0;
        };
        
        auto solidc = [&](int q)
        {
            if(useB)
            return A[q]<0.0 || B[q]<0.0;
            
            return A[q]<0.0;
        };
        
        // surface cells
        gcb=0;
        
        if(!useB)
        for(i=0; i<p->knox; ++i)
        for(j=0; j<p->knoy; ++j)
        for(k=0; k<p->knoz; ++k)
        if(fluid(IJK,false))
        if(solidc(Im1JK) || solidc(Ip1JK) || solidc(IJm1K) || solidc(IJp1K) || solidc(IJKm1) || solidc(IJKp1))
        ++gcb;
        
        // cells reached by up to three ghost cells
        vector<int> fgc(nf,0);
        
        for(i=gx.is; i<gx.ie; ++i)
        for(j=gx.js; j<gx.je; ++j)
        for(k=gx.ks; k<gx.ke; ++k)
        if(useB ? fluid(IJK,true) : fluid(IJK,false))
        ++fgc[IJK];
        
        auto add = [&](int ii, int jj, int kk)
        {
            if(ii>=gx.is && ii<gx.ie && jj>=gx.js && jj<gx.je && kk>=gx.ks && kk<gx.ke)
            ++fgc[(ii-p->imin)*p->jmax*p->kmax + (jj-p->jmin)*p->kmax + kk-p->kmin];
        };
        
        // sources: fluid cells of the global domain next to a solid cell; a source up to three
        // cells outside the subdomain reaches its cells
        for(i=gx.is; i<gx.ie; ++i)
        for(j=gx.js; j<gx.je; ++j)
        for(k=gx.ks; k<gx.ke; ++k)
        if(i+p->origin_i>=0 && i+p->origin_i<p->gknox
        && j+p->origin_j>=0 && j+p->origin_j<p->gknoy
        && k+p->origin_k>=0 && k+p->origin_k<p->gknoz)
        if(fluid(IJK,true))
        {
            if(i>gx.is   && solidc(Im1JK)) {add(i-1,j,k); add(i-2,j,k); add(i-3,j,k);}
            if(i<gx.ie-1 && solidc(Ip1JK)) {add(i+1,j,k); add(i+2,j,k); add(i+3,j,k);}
            if(j>gx.js   && solidc(IJm1K)) {add(i,j-1,k); add(i,j-2,k); add(i,j-3,k);}
            if(j<gx.je-1 && solidc(IJp1K)) {add(i,j+1,k); add(i,j+2,k); add(i,j+3,k);}
            if(k>gx.ks   && solidc(IJKm1)) {add(i,j,k-1); add(i,j,k-2); add(i,j,k-3);}
            if(k<gx.ke-1 && solidc(IJKp1)) {add(i,j,k+1); add(i,j,k+2); add(i,j,k+3);}
        }
        
        gcbextra=0;
        
        for(i=0; i<p->knox; ++i)
        for(j=0; j<p->knoy; ++j)
        for(k=0; k<p->knoz; ++k)
        if(fgc[IJK]>=2)
        ++gcbextra;
    }
}

void lexer::grid_solids(ghostcell *pgc)
{
    int i,j,n;
    lexer *p = this;
    
    if(gridgeo==nullptr)
    return;
    
    geo_mesh &G = *gridgeo;
    geo_raycast ray(this);
    
    const double dxm = (G.dxm>0.0) ? G.dxm : DXM;
    const geo_cart gx = geo_raycast::extended(this,dxm);
    
    // bed levels as initialised by DIVEMesh, geodat bed of its role
    for(n=0; n<imax*jmax; ++n)
    {
    solidbed[n] = -1.0e10;
    topobed[n] = global_zmin;
    }
    
    for(i=0; i<knox; ++i)
    for(j=0; j<knoy; ++j)
    {
        if(G.geodat==geo_mesh::role_solid)
        solidbed[IJ] = geobed[IJ];
        
        if(G.geodat==geo_mesh::role_topo)
        topobed[IJ] = geobed[IJ];
    }
    
    vector<double> phis, phit;
    
    const bool solid = build_field(this,G,ray,gx,geo_mesh::role_solid,solidbed,phis);
    const bool topo  = build_field(this,G,ray,gx,geo_mesh::role_topo,topobed,phit);
    
    // signed distance fields of the CFD solver
    if(A10==4 || A10==6)
    {
        for(i=0; i<knox; ++i)
        for(j=0; j<knoy; ++j)
        for(int k=0; k<knoz; ++k)
        {
            if(solid)
            flag_solid[IJK] = phis[IJK];
            
            if(topo)
            flag_topo[IJK] = phit[IJK];
        }
    }
    
    // the signed distance fields are only used by CFD and PTF
    if(A10!=4 && A10!=6)
    {
    del_Darray(flag_solid,imax*jmax*kmax);
    del_Darray(flag_topo,imax*jmax*kmax);
    flag_solid = flag_topo = nullptr;
    }
    
    // combined bed level
    // NHFLOW with immersed solids (A 580 1): the solids are not part of the bed
    for(i=0; i<knox; ++i)
    for(j=0; j<knoy; ++j)
    {
        if(A10==5 && A580==1)
        bed[IJ] = MAX(global_zmin,topobed[IJ]);
        
        else
        bed[IJ] = MAX(global_zmin,MAX(solidbed[IJ],topobed[IJ]));
    }
    
    // ghost cell estimates; absent fields are zero as in DIVEMesh
    const int nf = imax*jmax*kmax;
    
    if(!solid)
    phis.assign(nf,0.0);
    
    if(!topo)
    phit.assign(nf,0.0);
    
    int gcb_dummy=0;
    
    if(solid)
    estimate(this,gx,phis,phis,false,solid_gcb_est,solid_gcbextra_est);
    
    if(topo)
    {
    estimate(this,gx,phit,phit,false,topo_gcb_est,topo_gcbextra_est);
    topo_gcb_est *= 4;
    }
    
    estimate(this,gx,phis,phit,true,gcb_dummy,tot_gcbextra_est);
    
    // summary
    const int nsolid = pgc->globalimax(G.count(geo_mesh::role_solid));
    const int ntopo = pgc->globalimax(G.count(geo_mesh::role_topo));
    
    if(mpirank==0 && (nsolid>0 || ntopo>0 || G.geodat>0))
    {
    cout<<"grid geometry: "<<nsolid<<" solid and "<<ntopo<<" topo entities, "<<G.ntri<<" triangles";
    
    if(G.geodat>0)
    cout<<", geodat bed level ("<<(G.geodat==geo_mesh::role_topo?"topo":"solid")<<")";
    
    cout<<endl;
    }
}

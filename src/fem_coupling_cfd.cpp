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

#include"fem_coupling.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"field.h"
#include<mpi.h>
#include<cmath>
#include<algorithm>
#include<iostream>
#include<fstream>
#include<iomanip>

typedef fem_solid::Vec3 Vec3;

void fem_coupling::start_cfd(lexer *p, fdm *a, ghostcell *pgc, double alpha,
                             field &u, field &v, field &w, field &fx, field &fy, field &fz, bool finalize)
{
    starttime = pgc->timer();

    if(surf_version!=fs.surface_version())
    ini_points(p,pgc);

    // halo of the stage velocity must be current for the kernel interpolation
    pgc->start1(p,u,10);
    pgc->start2(p,v,11);
    pgc->start3(p,w,12);

    const int np = (int)pts.size();
    const std::vector<int>& deb = fs.debris();
    const int nd = (int)deb.size();
    const int nn = fs.nnode();
    const bool pressure_loads = (fs.coupling().loads==1);
    const double off = fs.coupling().pressure_offset;

    // buffer: 8 per point, 8 per debris particle, in the final stage 2 per node
    const size_t nbuf = 8*size_t(np+nd) + (finalize ? 2*size_t(nn) : 0);
    if(buf.size()<nbuf)
    buf.resize(nbuf);
    std::fill(buf.begin(),buf.begin()+nbuf,0.0);

    auto owns = [&](const Vec3& x)->bool
    {
        if(x(0)<p->originx || x(0)>=p->endx)
        return false;
        if(p->j_dir==1 && (x(1)<p->originy || x(1)>=p->endy))
        return false;
        if(x(2)<p->originz || x(2)>=p->endz)
        return false;
        return true;
    };

    // ------------------------------------------------------------------
    // 1. sample at the Lagrangian points (owner rank only):
    //    [0-3] fluid velocity, count  [4-6] pressure, count, water at probe
    //    [7] density times the cell size normal to the surface
    // ------------------------------------------------------------------
    Vec3 xp, vp, n;
    double A;

    for(int q=0; q<np; ++q)
    {
        point_state(q,xp,vp,n,A);
        double *b = &buf[8*q];

        if(owns(xp))
        {
            b[0] = interpolate_kernel(p,u,xp(0),xp(1),xp(2),1);
            b[1] = (p->j_dir==1) ? interpolate_kernel(p,v,xp(0),xp(1),xp(2),2) : 0.0;
            b[2] = interpolate_kernel(p,w,xp(0),xp(1),xp(2),3);
            b[3] = 1.0;

            if(finalize)
            {
                const int ii = p->posc_i(xp(0)), jj = (p->j_dir==1) ? p->posc_j(xp(1)) : 0, kk = p->posc_k(xp(2));
                const double hn = std::fabs(n(0))*p->DXN[ii+marge] + std::fabs(n(1))*(p->j_dir==1 ? p->DYN[jj+marge] : 0.0) + std::fabs(n(2))*p->DZN[kk+marge];
                b[7] = p->ccipol4(a->ro,xp(0),xp(1),xp(2))*hn;
            }
        }

        if(finalize && pressure_loads)
        {
            Vec3 pr = xp + off*dxmin*n;
            if(p->j_dir==0)
            pr(1) = p->YP[marge];

            if(owns(pr))
            {
                // probe distance in local cells
                const int ii = p->posc_i(pr(0)), jj = (p->j_dir==1) ? p->posc_j(pr(1)) : 0, kk = p->posc_k(pr(2));
                pr(0) = xp(0) + off*n(0)*p->DXN[std::max(0,std::min(ii,p->knox-1))+marge];
                if(p->j_dir==1)
                pr(1) = xp(1) + off*n(1)*p->DYN[std::max(0,std::min(jj,p->knoy-1))+marge];
                pr(2) = xp(2) + off*n(2)*p->DZN[std::max(0,std::min(kk,p->knoz-1))+marge];

                // probes inside the bed or a solid body carry no pressure
                if(p->ccipol4a(a->solid,pr(0),pr(1),pr(2))>=0.0 && p->ccipol4a(a->topo,pr(0),pr(1),pr(2))>=0.0)
                {
                    b[4] = p->ccipol4a(a->press,pr(0),pr(1),pr(2)) - p->pressgage;
                    b[5] = 1.0;
                    b[6] = p->ccipol4(a->phi,pr(0),pr(1),pr(2))>=0.0 ? 1.0 : 0.0;
                }
            }
        }
    }

    // debris particles: velocity, count, water, density
    for(int d=0; d<nd; ++d)
    {
        const Vec3& x = fs.pos(deb[d]);
        double *b = &buf[8*(np+d)];

        if(owns(x))
        {
            b[0] = p->ccipol1(u,x(0),x(1),x(2));
            b[1] = (p->j_dir==1) ? p->ccipol2(v,x(0),x(1),x(2)) : 0.0;
            b[2] = p->ccipol3(w,x(0),x(1),x(2));
            b[3] = 1.0;
            b[4] = p->ccipol4(a->phi,x(0),x(1),x(2))>=0.0 ? 1.0 : 0.0;
            b[5] = p->ccipol4(a->ro,x(0),x(1),x(2));
        }
    }

    // fluid density at the nodes (enclosed fluid of the immersed boundary)
    if(finalize && !pressure_loads)
    for(int i=0; i<nn; ++i)
    {
        Vec3 x = fs.pos(i);
        if(p->j_dir==0)
        x(1) = p->YP[marge];
        if(owns(x))
        {
            double *b = &buf[8*size_t(np+nd)+2*size_t(i)];
            b[0] = p->ccipol4(a->ro,x(0),x(1),x(2));
            b[1] = 1.0;
        }
    }

    if(nbuf>0)
    MPI_Allreduce(MPI_IN_PLACE,buf.data(),(int)nbuf,MPI_DOUBLE,MPI_SUM,pgc->mpi_comm);

    // ------------------------------------------------------------------
    // 2. direct forcing of the structural surface velocity
    // ------------------------------------------------------------------
    const double hy = fs.lattice_h(1);

    if(fs.coupling().forcing)
    for(int q=0; q<np; ++q)
    {
        const double *b = &buf[8*q];
        if(b[3]<0.5)
        continue;

        point_state(q,xp,vp,n,A);

        Vec3 uf(b[0]/b[3],b[1]/b[3],b[2]/b[3]);
        if(p->j_dir==0)
        uf(1) = vp(1);

        const Vec3 f = (vp - uf)/(alpha*p->dt);

        // in 2D the face area is taken per unit width of the fluid cell
        const double Aeff = (p->j_dir==0) ? A/hy : A;
        spread(p,fx,fy,fz,xp,f,Aeff,&n);
    }

    // ------------------------------------------------------------------
    // 3. debris drag (quadratic), reaction onto the fluid
    // ------------------------------------------------------------------
    const double cd = fs.coupling().debris_cd;

    for(int d=0; d<nd; ++d)
    {
        const double *b = &buf[8*(np+d)];
        fdeb[d].setZero();

        if(b[3]<0.5 || b[4]/b[3]<0.5)
        continue;

        const int i = deb[d];
        Vec3 ur = Vec3(b[0],b[1],b[2])/b[3] - fs.vel(i);
        if(p->j_dir==0)
        ur(1) = 0.0;

        const double rho = b[5]/b[3];
        const double Ad = std::pow(fs.node_volume(i),2.0/3.0);
        fdeb[d] = 0.5*rho*cd*Ad*ur.norm()*ur;

        if(fs.coupling().debris_reaction)
        {
            Vec3 fr = -fdeb[d]/rho;
            if(p->j_dir==0)
            fr /= hy;
            spread(p,fx,fy,fz,fs.pos(i),fr,1.0,nullptr);
        }
    }

    pgc->start1(p,fx,10);
    pgc->start2(p,fy,11);
    pgc->start3(p,fz,12);

    // ------------------------------------------------------------------
    // 4. loads and solid time step (final stage)
    // ------------------------------------------------------------------
    if(finalize)
    finish_step(p,pgc,alpha);
}

void fem_coupling::finish_step(lexer *p, ghostcell *pgc, double alpha)
{
    (void)alpha;
    const int np = (int)pts.size();
    const std::vector<int>& deb = fs.debris();
    const int nd = (int)deb.size();
    const int nn = fs.nnode();

    fs.clear_loads();
    fs.clear_coupling();

    Vec3 xp, vp, n;
    double A;

    if(fs.coupling().loads==0)
    {
        // reaction: the fluid parcel in the forcing volume of every surface
        // point (mass rho dV, velocity u_f) is attached to the face nodes for
        // the step; its momentum exchange is the load on the solid
        for(int q=0; q<np; ++q)
        {
            const double *b = &buf[8*q];
            if(b[3]<0.5)
            continue;

            point_state(q,xp,vp,n,A);

            Vec3 uf(b[0]/b[3],b[1]/b[3],b[2]/b[3]);
            if(p->j_dir==0)
            uf(1) = 0.0;

            // fluid mass in the forcing volume of the point
            const double c = (b[7]/b[3])*A;

            const lpoint& L = pts[q];
            const fem_solid::face& fc = fs.surface()[L.face];
            const double wgt[4] = {(1.0-L.s)*(1.0-L.t), L.s*(1.0-L.t), L.s*L.t, (1.0-L.s)*L.t};
            for(int k=0; k<4; ++k)
            fs.add_coupling(fc.n[k],wgt[k]*c,wgt[k]*c*uf);
        }

        // fluid enclosed by the immersed boundary: buoyancy and inertia
        const double *bn = &buf[8*size_t(np+nd)];
        for(int i=0; i<nn; ++i)
        if(!fs.is_debris(i) && bn[2*i+1]>0.5)
        fs.set_fluid_mass(i,(bn[2*i]/bn[2*i+1])*fs.node_volume(i));
    }
    else
    {
        // pressure on the surface faces (explicit)
        for(int q=0; q<np; ++q)
        {
            const double *b = &buf[8*q];
            if(b[5]<0.5)
            continue;

            point_state(q,xp,vp,n,A);

            const double pr = b[4]/b[5];
            const Vec3 F = -pr*A*n;

            const lpoint& L = pts[q];
            const fem_solid::face& fc = fs.surface()[L.face];
            const double wgt[4] = {(1.0-L.s)*(1.0-L.t), L.s*(1.0-L.t), L.s*L.t, (1.0-L.s)*L.t};
            for(int k=0; k<4; ++k)
            fs.add_load(fc.n[k],wgt[k]*F);
        }
    }

    // debris: drag and buoyancy
    const Vec3 g(p->W20,p->W21,p->W22);
    for(int d=0; d<nd; ++d)
    {
        const double *b = &buf[8*(np+d)];
        if(b[3]<0.5 || b[4]/b[3]<0.5)
        continue;

        const int i = deb[d];
        const double rho = b[5]/b[3];
        fs.add_load(i, fdeb[d] - rho*fs.node_volume(i)*g);
    }

    try
    {
        fs.advance(p->dt);
    }
    catch(std::exception& e)
    {
        if(p->mpirank==0)
        std::cout<<"\n!!! "<<e.what()<<" !!!\n"<<std::endl;
        MPI_Abort(pgc->mpi_comm,1);
    }

    // log
    if(p->mpirank==0)
    {
        const Vec3 Ff = fs.total_load();
        const Vec3 Fs = fs.support_force();

        std::ofstream lg((outdir+"/REEF3D_FEM_log.dat").c_str(),std::ios::app);
        lg<<std::setprecision(9)<<fs.time()<<" "<<fs.last_substeps()<<" "<<fs.n_alive()<<" "<<fs.n_eroded()<<" "<<fs.debris().size()<<" "
          <<Ff(0)<<" "<<Ff(1)<<" "<<Ff(2)<<" "<<Fs(0)<<" "<<Fs(1)<<" "<<Fs(2)<<" "
          <<fs.max_vonmises()<<" "<<fs.max_displacement()<<" "<<fs.kinetic_energy()<<" "<<fs.dissipated_energy()<<"\n";

        for(const fem_solid::monitor& mo : fs.monitors())
        {
            std::ofstream ts((outdir+"/REEF3D_FEM_monitor_"+mo.name+".dat").c_str(),std::ios::app);
            const Vec3 dx = fs.pos(mo.node)-fs.ref_pos(mo.node);
            const Vec3& vv = fs.vel(mo.node);
            ts<<std::setprecision(9)<<fs.time()<<" "<<dx(0)<<" "<<dx(1)<<" "<<dx(2)<<" "<<vv(0)<<" "<<vv(1)<<" "<<vv(2)<<"\n";
        }

        std::cout<<"FEM time: "<<pgc->timer()-starttime<<"  substeps: "<<fs.last_substeps()
                 <<"  eroded: "<<fs.n_eroded()<<"  debris: "<<fs.debris().size()<<std::endl;
    }

    print(p);
}

double fem_coupling::interpolate_kernel(lexer *p, field &f, double xs, double ys, double zs, int comp)
{
    // owner rank only: the point lies in this subdomain, the 5-point stencil
    // reaches 2 ghost cells
    const int ii = p->posc_i(xs);
    const int jj = (p->j_dir==1) ? p->posc_j(ys) : 0;
    const int kk = p->posc_k(zs);

    const double dx = p->DXN[ii+marge];
    const double dy = p->DYN[jj+marge];
    const double dz = p->DZN[kk+marge];

    double wx[5], wy[5], wz[5];
    for(int c=0; c<5; ++c)
    {
        const int i = ii+c-2, j = jj+c-2, k = kk+c-2;
        wx[c] = kernel(((comp==1 ? p->XN[i+1+marge] : p->XP[i+marge]) - xs)/dx);
        wy[c] = (p->j_dir==1) ? kernel(((comp==2 ? p->YN[j+1+marge] : p->YP[j+marge]) - ys)/dy) : (c==2 ? 1.0 : 0.0);
        wz[c] = kernel(((comp==3 ? p->ZN[k+1+marge] : p->ZP[k+marge]) - zs)/dz);
    }

    double s = 0.0, ws = 0.0;
    for(int ci=0; ci<5; ++ci)
    {
        if(wx[ci]==0.0) continue;
        for(int cj=0; cj<5; ++cj)
        {
            if(wy[cj]==0.0) continue;
            for(int ck=0; ck<5; ++ck)
            {
                const double D = wx[ci]*wy[cj]*wz[ck];
                if(D==0.0) continue;
                s += D*f(ii+ci-2,jj+cj-2,kk+ck-2);
                ws += D;
            }
        }
    }

    return ws>0.0 ? s/ws : 0.0;
}

void fem_coupling::spread(lexer *p, field &fx, field &fy, field &fz, const Vec3& xp, const Vec3& f, double A, const Vec3* n)
{
    // every rank spreads onto its own cells only (ghost cells are exchanged
    // afterwards). With a normal n, A is the point area and the volume is
    // A times the cell size along n; without, A is the volume itself.
    if(xp(0)<p->originx-2.0*p->DXN[marge] || xp(0)>p->endx+2.0*p->DXN[p->knox-1+marge])
    return;
    if(xp(2)<p->originz-2.0*p->DZN[marge] || xp(2)>p->endz+2.0*p->DZN[p->knoz-1+marge])
    return;
    if(p->j_dir==1 && (xp(1)<p->originy-2.0*p->DYN[marge] || xp(1)>p->endy+2.0*p->DYN[p->knoy-1+marge]))
    return;

    const int ii = p->posc_i(xp(0));
    const int jj = (p->j_dir==1) ? p->posc_j(xp(1)) : 0;
    const int kk = p->posc_k(xp(2));

    const int is = std::max(ii-2,0), ie = std::min(ii+2,p->knox-1);
    const int js = (p->j_dir==1) ? std::max(jj-2,0) : 0;
    const int je = (p->j_dir==1) ? std::min(jj+2,p->knoy-1) : 0;
    const int ks = std::max(kk-2,0), ke = std::min(kk+2,p->knoz-1);

    if(is>ie || js>je || ks>ke)
    return;

    for(int i=is; i<=ie; ++i)
    {
        const double dx = p->DXN[i+marge];
        const double wxC = kernel((p->XP[i+marge]-xp(0))/dx);
        const double wxF = kernel((p->XN[i+1+marge]-xp(0))/dx);
        if(wxC==0.0 && wxF==0.0) continue;

        for(int j=js; j<=je; ++j)
        {
            const double dy = p->DYN[j+marge];
            const double wyC = (p->j_dir==1) ? kernel((p->YP[j+marge]-xp(1))/dy) : 1.0;
            const double wyF = (p->j_dir==1) ? kernel((p->YN[j+1+marge]-xp(1))/dy) : 1.0;
            if(wyC==0.0 && wyF==0.0) continue;

            for(int k=ks; k<=ke; ++k)
            {
                const double dz = p->DZN[k+marge];
                const double wzC = kernel((p->ZP[k+marge]-xp(2))/dz);
                const double wzF = kernel((p->ZN[k+1+marge]-xp(2))/dz);
                if(wzC==0.0 && wzF==0.0) continue;

                double dV = A;
                if(n)
                dV = A*(std::fabs((*n)(0))*dx + std::fabs((*n)(1))*(p->j_dir==1 ? dy : 0.0) + std::fabs((*n)(2))*dz);

                // 2D: forces are per unit width, the fluid cell has width dy
                if(p->j_dir==0)
                dV *= dy;

                const double fac = dV/(dx*dy*dz);

                fx(i,j,k) += f(0)*wxF*wyC*wzC*fac;
                if(p->j_dir==1)
                fy(i,j,k) += f(1)*wxC*wyF*wzC*fac;
                fz(i,j,k) += f(2)*wxC*wyC*wzF*fac;
            }
        }
    }
}

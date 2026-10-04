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
static const size_t BP = 15;     // buffer entries per Lagrangian point / debris particle

void fem_coupling::start_cfd(lexer *p, fdm *a, ghostcell *pgc, double alpha,
                             field &u, field &v, field &w, field &fx, field &fy, field &fz, bool finalize)
{
    starttime = pgc->timer();

    if(!initialised)
    first_call(p,a,pgc);

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
    const bool parcels = (fs.coupling().loads!=1);
    const bool probes = (fs.coupling().loads!=2) || fs.n_rigid()>0;

    // buffer: 8 per point, 8 per debris particle, in the final stage 2 per node
    const size_t nbuf = BP*size_t(np+nd) + (finalize ? 2*size_t(nn) : 0);
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
    //    [7] cell size normal to the surface  [8-10] momentum rho u  [11-13] density  [14] level set
    //    (kernel weighted on the staggered faces: the momentum the forcing acts on)
    // ------------------------------------------------------------------
    Vec3 xp, vp, n;
    double A;

    for(int q=0; q<np; ++q)
    {
        point_state(q,xp,vp,n,A);
        double *b = &buf[BP*q];

        if(owns(xp))
        {
            b[0] = interpolate_kernel(p,u,xp(0),xp(1),xp(2),1);
            b[1] = (p->j_dir==1) ? interpolate_kernel(p,v,xp(0),xp(1),xp(2),2) : 0.0;
            b[2] = interpolate_kernel(p,w,xp(0),xp(1),xp(2),3);
            b[3] = 1.0;
            b[14] = p->ccipol4(a->phi,xp(0),(p->j_dir==1) ? xp(1) : p->YP[marge],xp(2));

            if(finalize && parcels)
            {
                const int ii = p->posc_i(xp(0)), jj = (p->j_dir==1) ? p->posc_j(xp(1)) : 0, kk = p->posc_k(xp(2));
                b[7] = std::fabs(n(0))*p->DXN[ii+marge] + std::fabs(n(1))*(p->j_dir==1 ? p->DYN[jj+marge] : 0.0) + std::fabs(n(2))*p->DZN[kk+marge];
                double r;
                b[8] = interpolate_kernel(p,u,xp(0),xp(1),xp(2),1,&a->ro,&r);  b[11] = r;
                if(p->j_dir==1)
                {
                    b[9] = interpolate_kernel(p,v,xp(0),xp(1),xp(2),2,&a->ro,&r);  b[12] = r;
                }
                b[10] = interpolate_kernel(p,w,xp(0),xp(1),xp(2),3,&a->ro,&r);  b[13] = r;
            }
        }

        if(finalize && probes)
        probe_pressure(p,a,xp,n,b);
    }

    // debris particles: velocity, count, water, density
    for(int d=0; d<nd; ++d)
    {
        const Vec3& x = fs.pos(deb[d]);
        double *b = &buf[BP*(np+d)];

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
    if(finalize && parcels)
    for(int i=0; i<nn; ++i)
    {
        Vec3 x = fs.pos(i);
        if(p->j_dir==0)
        x(1) = p->YP[marge];
        if(owns(x))
        {
            double *b = &buf[BP*size_t(np+nd)+2*size_t(i)];
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
        const double *b = &buf[BP*q];
        if(b[3]<0.5)
        continue;

        // no forcing of the air: points more than 1.5 cells away from the
        // water stay passive (forcing light air next to water around moving
        // fragments produced air velocities of 200 m/s in the seawall collapse)
        if(!fs.coupling().air_forcing && b[14]/b[3] < -1.5*dxmin)
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
        const double *b = &buf[BP*(np+d)];
        fdeb[d].setZero();

        if(b[3]<0.5 || b[4]/b[3]<0.5)
        continue;

        const int i = deb[d];
        Vec3 ur = Vec3(b[0],b[1],b[2])/b[3] - fs.vel(i);
        if(p->j_dir==0)
        ur(1) = 0.0;

        const double rho = b[5]/b[3];
        if(!ur.allFinite() || !std::isfinite(rho) || rho<=0.0)
        continue;
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
    finish_step(p,a,pgc,alpha);
}

void fem_coupling::finish_step(lexer *p, fdm *a, ghostcell *pgc, double alpha)
{
    fs.set_coupling_alpha(alpha);

    // contact with the bed / solids of the fluid grid: distance and normal at the nodes
    if(fs.bed_contact())
    sample_bed(p,a,pgc);
    const int np = (int)pts.size();
    const std::vector<int>& deb = fs.debris();
    const int nd = (int)deb.size();
    const int nn = fs.nnode();

    fs.clear_loads();
    fs.clear_coupling();

    Vec3 xp, vp, n;
    double A;

    const int mode = fs.coupling().loads;

    if(mode!=1)
    {
        // reaction: the fluid parcel in the forcing volume of every surface
        // point (mass rho dV, velocity u_f) is attached to the face nodes for
        // the step; its momentum exchange is the load on the solid
        for(int q=0; q<np; ++q)
        {
            const double *b = &buf[BP*q];
            if(b[3]<0.5)
            continue;

            point_state(q,xp,vp,n,A);

            // fluid mass and momentum in the forcing volume of the point
            const double dV = (b[7]/b[3])*A;
            Vec3 M(b[11],b[12],b[13]), Mu(b[8],b[9],b[10]);
            M *= dV/b[3];
            Mu *= dV/b[3];
            if(p->j_dir==0)
            {
                M(1) = 0.5*(M(0)+M(2));
                Mu(1) = 0.0;
            }

            // tangential momentum of the fluid is not handed to the solid
            // (default): at the usual grid sizes the tangential momentum the
            // direct forcing removes in the forcing layer is far above the real
            // skin friction (a run-up jet lifted wall pieces). The attached mass
            // stays complete, so a solid moving tangentially still loses the
            // momentum it gives to the fluid (no momentum is created).
            if(!fs.coupling().shear)
            {
                Mu = Mu.dot(n)*n;
                if(p->j_dir==0)
                Mu(1) = 0.0;
            }

            const lpoint& L = pts[q];
            const fem_solid::face& fc = fs.surface()[L.face];
            if(fs.is_rigid_node(fc.n[0]))
            continue;
            const double wgt[4] = {(1.0-L.s)*(1.0-L.t), L.s*(1.0-L.t), L.s*L.t, (1.0-L.s)*L.t};
            for(int k=0; k<4; ++k)
            fs.add_coupling(fc.n[k],wgt[k]*M,wgt[k]*Mu);
        }

        // fluid enclosed by the immersed boundary: buoyancy and inertia
        const double *bn = &buf[BP*size_t(np+nd)];
        for(int i=0; i<nn; ++i)
        if(!fs.is_debris(i) && !fs.is_rigid_node(i) && bn[2*i+1]>0.5)
        fs.set_fluid_mass(i,(bn[2*i]/bn[2*i+1])*fs.node_volume(i));
    }

    if(mode!=2 || fs.n_rigid()>0)
    {
        // pressure on the surface faces, probed outside the immersed boundary
        fprb.assign(nn,Vec3::Zero());
        for(int q=0; q<np; ++q)
        {
            const double *b = &buf[BP*q];
            if(b[5]<0.5)
            continue;

            point_state(q,xp,vp,n,A);

            const double pr = b[4]/b[5];
            const Vec3 F = -pr*A*n;

            const lpoint& L = pts[q];
            const fem_solid::face& fc = fs.surface()[L.face];
            const double wgt[4] = {(1.0-L.s)*(1.0-L.t), L.s*(1.0-L.t), L.s*L.t, (1.0-L.s)*L.t};
            for(int k=0; k<4; ++k)
            fprb[fc.n[k]] += wgt[k]*F;
        }

        if(mode==1)
        {
            // explicit pressure loads
            for(int i=0; i<nn; ++i)
            fs.add_load(i,fprb[i]);
        }
        else
        {
            // rigid bodies: explicit pressure loads on the whole body (full
            // mass, no parcels): the probes see the outer flow only, so
            // buoyancy and the restoring force of floating bodies are right
            for(int i=0; i<nn; ++i)
            if(fs.is_rigid_node(i))
            fs.add_load(i,fprb[i]);
        }

        // rigid bodies: added mass for the stabilisation, by the wetted fraction
        if(fs.n_rigid()>0)
        {
            std::vector<double> nw(fs.n_rigid(),0.0), nt(fs.n_rigid(),0.0);
            for(int q=0; q<np; ++q)
            {
                const int k = fs.rigid_of_node(fs.surface()[pts[q].face].n[0]);
                if(k<0)
                continue;
                const double *b = &buf[BP*q];
                nt[k] += 1.0;
                if(b[5]>0.5)
                nw[k] += b[6]/b[5];
            }
            for(int k=0; k<fs.n_rigid(); ++k)
            fs.set_rigid_added_mass(k, fs.coupling().added_mass*rho_w*fs.rigid(k).Aunit*(nt[k]>0.0 ? nw[k]/nt[k] : 0.0));
        }
        if(mode==0)
        {
            // hybrid: the parcels carry the fast (added mass, impact) part of
            // the loads and keep the coupling stable; a slow correction moves
            // their mean onto the probed pressure, which is free of the fluid
            // inside the immersed boundary (hydrostatics, one-sided walls)
            const double c = 1.0/std::max(1.0,fs.coupling().hybrid_tau);
            if(!hybrid_ini)
            fp_bar = fprb;
            else
            for(int i=0; i<nn; ++i)
            fp_bar[i] += c*(fprb[i]-fp_bar[i]);

            // only on the supported structure: free fragments and debris are
            // light and fast, the explicit pressure part would destabilise them
            if(hybrid_ini)
            for(int i=0; i<nn; ++i)
            if(fs.on_supported_body(i))
            fs.add_load(i,fp_bar[i]-fr_bar[i]);
        }
    }

    // debris: drag and buoyancy
    const Vec3 g(p->W20,p->W21,p->W22);
    for(int d=0; d<nd; ++d)
    {
        const double *b = &buf[BP*(np+d)];
        if(b[3]<0.5 || b[4]/b[3]<0.5)
        continue;

        const int i = deb[d];
        const double rho = b[5]/b[3];
        if(!std::isfinite(rho) || !fdeb[d].allFinite())
        continue;
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

    // hybrid: filtered parcel loads of the step
    if(fs.coupling().loads==0)
    {
        const double c = 1.0/std::max(1.0,fs.coupling().hybrid_tau);
        if(!hybrid_ini)
        {
            fr_bar.resize(nn);
            for(int i=0; i<nn; ++i)
            fr_bar[i] = fs.parcel_load(i);
            hybrid_ini = true;
        }
        else
        for(int i=0; i<nn; ++i)
        fr_bar[i] += c*(fs.parcel_load(i)-fr_bar[i]);
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

        // rigid bodies: trajectory and contact force
        for(int k=0; k<std::min(fs.n_rigid(),50); ++k)
        {
            const fem_solid::rigid_body& rb = fs.rigid(k);
            std::ofstream rf((outdir+"/REEF3D_FEM_rigid_"+std::to_string(k+1)+".dat").c_str(),nstep==0 ? std::ios::out : std::ios::app);
            if(nstep==0)
            rf<<"# rigid body "<<k+1<<": mass "<<rb.M<<" kg, volume "<<rb.Vol<<" m3"<<(p->j_dir==0 ? " (2D: of the FEM slice)" : "")
              <<"\n# time  centre x y z  velocity x y z  angular velocity x y z  max contact force in the step [N]\n";
            rf<<std::setprecision(9)<<fs.time()<<" "<<rb.c(0)<<" "<<rb.c(1)<<" "<<rb.c(2)<<" "<<rb.V(0)<<" "<<rb.V(1)<<" "<<rb.V(2)
              <<" "<<rb.w(0)<<" "<<rb.w(1)<<" "<<rb.w(2)<<" "<<rb.fcstep<<"\n";
        }

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

    ++nstep;
    update_summary(p,false);
    print(p);
}

double fem_coupling::interpolate_kernel(lexer *p, field &f, double xs, double ys, double zs, int comp, field *ro, double *rho)
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

    double s = 0.0, ws = 0.0, wr = 0.0;
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
                const int i = ii+ci-2, j = jj+cj-2, k = kk+ck-2;
                if(ro)
                {
                    // density on the staggered face of the component
                    const double rf = 0.5*((*ro)(i,j,k) + (comp==1 ? (*ro)(i+1,j,k) : comp==2 ? (*ro)(i,j+1,k) : (*ro)(i,j,k+1)));
                    s += D*rf*f(i,j,k);
                    wr += D*rf;
                }
                else
                s += D*f(i,j,k);
                ws += D;
            }
        }
    }

    if(rho)
    *rho = ws>0.0 ? wr/ws : 0.0;
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

void fem_coupling::probe_pressure(lexer *p, fdm *a, const Vec3& xp, const Vec3& n, double *b, bool hydrostatic)
{
    // pressure probe outside the surface point (owner rank of the probe):
    // b[4] pressure, b[5] count, b[6] water at the probe
    // the probe offset is shortened when the probe would leave the fluid domain
    // (a body close to the bottom or a side wall of the domain)
    double off = fs.coupling().pressure_offset;
    auto in_domain = [&](const Vec3& q)
    {
        const double eps = 1.0e-9;
        return q(0)>p->global_xmin+eps && q(0)<p->global_xmax-eps && q(2)>p->global_zmin+eps && q(2)<p->global_zmax-eps
            && (p->j_dir==0 || (q(1)>p->global_ymin+eps && q(1)<p->global_ymax-eps));
    };
    for(double f : {1.0, 2.0/3.0, 1.0/3.0})
    {
        off = f*fs.coupling().pressure_offset;
        if(in_domain(xp + 1.01*off*dxmin*n))
        break;
    }
    Vec3 pr = xp + off*dxmin*n;
    if(p->j_dir==0)
    pr(1) = p->YP[marge];

    if(pr(0)<p->originx || pr(0)>=p->endx || pr(2)<p->originz || pr(2)>=p->endz)
    return;
    if(p->j_dir==1 && (pr(1)<p->originy || pr(1)>=p->endy))
    return;

    // probe distance in local cells
    const int ii = p->posc_i(pr(0)), jj = (p->j_dir==1) ? p->posc_j(pr(1)) : 0, kk = p->posc_k(pr(2));
    pr(0) = xp(0) + off*n(0)*p->DXN[std::max(0,std::min(ii,p->knox-1))+marge];
    if(p->j_dir==1)
    pr(1) = xp(1) + off*n(1)*p->DYN[std::max(0,std::min(jj,p->knoy-1))+marge];
    pr(2) = xp(2) + off*n(2)*p->DZN[std::max(0,std::min(kk,p->knoz-1))+marge];

    // a probe inside the bed or a solid body is moved closer to the surface
    // (a body near the bed), or out horizontally (faces near the bed whose
    // normal points slightly downwards), else it carries no pressure
    auto inside = [&](const Vec3& q){return p->ccipol4a(a->solid,q(0),q(1),q(2))<0.0 || p->ccipol4a(a->topo,q(0),q(1),q(2))<0.0;};
    // shorter offsets first (a body close to the bed), then horizontal
    for(double f : {2.0/3.0, 1.0/3.0})
    if(inside(pr))
    {
        Vec3 q = xp;
        q(0) += f*off*n(0)*p->DXN[std::max(0,std::min(ii,p->knox-1))+marge];
        if(p->j_dir==1)
        q(1) += f*off*n(1)*p->DYN[std::max(0,std::min(jj,p->knoy-1))+marge];
        else
        q(1) = pr(1);
        q(2) += f*off*n(2)*p->DZN[std::max(0,std::min(kk,p->knoz-1))+marge];
        if(!inside(q))
        pr = q;
    }
    if(inside(pr))
    {
        Vec3 nh(n(0),n(1),0.0);
        if(nh.norm()>0.3)
        {
            nh.normalize();
            pr(0) = xp(0) + off*nh(0)*p->DXN[std::max(0,std::min(ii,p->knox-1))+marge];
            if(p->j_dir==1)
            pr(1) = xp(1) + off*nh(1)*p->DYN[std::max(0,std::min(jj,p->knoy-1))+marge];
            pr(2) = xp(2);
        }
        if(pr(0)<p->originx || pr(0)>=p->endx || pr(2)<p->originz || pr(2)>=p->endz)
        return;
        if(p->j_dir==1 && (pr(1)<p->originy || pr(1)>=p->endy))
        return;
    }

    if(!inside(pr))
    {
        const double phi = p->ccipol4(a->phi,pr(0),pr(1),pr(2));
        // the probe lies off*dx outside: its pressure is taken back to the
        // surface with the hydrostatic gradient of the fluid at the probe
        // (otherwise a floating body gets rho g A off*dx too much buoyancy)
        if(hydrostatic)
        b[4] = (phi>=0.0) ? p->W1*std::fabs(p->W22)*std::max(p->phimean-xp(2),0.0) : 0.0;
        else
        {
            const double rho = p->ccipol4(a->ro,pr(0),pr(1),pr(2));
            const double gdx = p->W20*(pr(0)-xp(0)) + (p->j_dir==1 ? p->W21*(pr(1)-xp(1)) : 0.0) + p->W22*(pr(2)-xp(2));
            b[4] = p->ccipol4a(a->press,pr(0),pr(1),pr(2)) - p->pressgage - rho*gdx;
        }
        b[5] = 1.0;
        b[6] = phi>=0.0 ? 1.0 : 0.0;
    }
}

void fem_coupling::pressure_loads(lexer *p, fdm *a, ghostcell *pgc, std::vector<Vec3>& F, bool hydrostatic)
{
    // nodal loads -p n dA of the current pressure field, or of the hydrostatic
    // pressure below the still water level (as I 12 1, all ranks)
    const int np = (int)pts.size();
    std::vector<double> b(BP*size_t(np),0.0);
    Vec3 xp, vp, n;
    double A;

    for(int q=0; q<np; ++q)
    {
        point_state(q,xp,vp,n,A);
        probe_pressure(p,a,xp,n,&b[BP*q],hydrostatic);
    }
    if(!b.empty())
    MPI_Allreduce(MPI_IN_PLACE,b.data(),(int)b.size(),MPI_DOUBLE,MPI_SUM,pgc->mpi_comm);

    F.assign(fs.nnode(),Vec3::Zero());
    for(int q=0; q<np; ++q)
    {
        const double *bq = &b[BP*q];
        if(bq[5]<0.5)
        continue;
        point_state(q,xp,vp,n,A);
        const Vec3 f = -(bq[4]/bq[5])*A*n;
        const lpoint& L = pts[q];
        const fem_solid::face& fc = fs.surface()[L.face];
        const double wgt[4] = {(1.0-L.s)*(1.0-L.t), L.s*(1.0-L.t), L.s*L.t, (1.0-L.s)*L.t};
        for(int k=0; k<4; ++k)
        F[fc.n[k]] += wgt[k]*f;
    }
}

void fem_coupling::sample_bed(lexer *p, fdm *a, ghostcell *pgc)
{
    // signed distance to the bed and to the solids of the fluid grid (the
    // smaller of topo and solid) and its gradient, at the nodes of free bodies
    // and debris particles near the bed; owner rank samples, one reduction
    const int nn = fs.nnode();
    std::vector<double> b(5*size_t(nn),0.0);
    const double hb = 2.0*std::max(fs.hmin(),dxmin);

    auto lsv = [&](double x,double y,double z)->double
    {
        return std::min(p->ccipol4a(a->topo,x,y,z),p->ccipol4a(a->solid,x,y,z));
    };

    for(int i=0; i<nn; ++i)
    {
        if(!fs.free_node(i))
        continue;
        Vec3 x = fs.pos(i);
        if(p->j_dir==0)
        x(1) = p->YP[marge];
        if(x(0)<p->originx || x(0)>=p->endx || x(2)<p->originz || x(2)>=p->endz)
        continue;
        if(p->j_dir==1 && (x(1)<p->originy || x(1)>=p->endy))
        continue;

        const double phi = lsv(x(0),x(1),x(2));
        if(phi>hb)
        continue;

        const double d = 0.5*dxmin;
        Vec3 g;
        g(0) = (lsv(x(0)+d,x(1),x(2)) - lsv(x(0)-d,x(1),x(2)))/(2.0*d);
        g(1) = (p->j_dir==1) ? (lsv(x(0),x(1)+d,x(2)) - lsv(x(0),x(1)-d,x(2)))/(2.0*d) : 0.0;
        g(2) = (lsv(x(0),x(1),x(2)+d) - lsv(x(0),x(1),x(2)-d))/(2.0*d);
        double* q = &b[5*size_t(i)];
        q[0] = phi; q[1] = g(0); q[2] = g(1); q[3] = g(2); q[4] = 1.0;
    }

    if(!b.empty())
    MPI_Allreduce(MPI_IN_PLACE,b.data(),(int)b.size(),MPI_DOUBLE,MPI_SUM,pgc->mpi_comm);

    fs.clear_bed_samples();
    for(int i=0; i<nn; ++i)
    {
        const double* q = &b[5*size_t(i)];
        if(q[4]<0.5)
        continue;
        Vec3 n(q[1]/q[4],q[2]/q[4],q[3]/q[4]);
        const double gn = n.norm();
        if(gn<1.0e-6)
        continue;
        fs.set_bed_sample(i,q[0]/q[4],n/gn);
    }
}

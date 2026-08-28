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
Author: Elyas Larkermani
--------------------------------------------------------------------*/

#include"6DOF_chrono.h"
#include"6DOF_chrono_backend.h"
#include"6DOF_obj.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

#include<algorithm>
#include<cmath>
#include<cstdlib>
#include<iostream>
#include<vector>

struct sixdof_chrono::impl
{
    void *backend;
    bool ready;
    impl() : backend(0), ready(false) {}
};

sixdof_chrono::sixdof_chrono()
{
    pimpl = new impl();
}

sixdof_chrono::~sixdof_chrono()
{
    if(pimpl->backend)
    reef3d_chrono_destroy(pimpl->backend);
    delete pimpl;
}

void sixdof_chrono::initialize(lexer *p, std::vector<sixdof_obj*> &fb_obj)
{
#ifndef REEF3D_CHRONO
    if(p->mpirank==0)
    {
        cout<<"ERROR: X 13 1 requested Chrono motion, but REEF3D was built without Chrono."<<endl;
        cout<<"Rebuild with Chrono enabled (CHRONO_DIR)."<<endl;
    }
    exit(1);
#else
    pimpl->backend = reef3d_chrono_create(p->X16==1 ? 1 : 0);
    reef3d_chrono_add_floor(pimpl->backend, p->originx, p->originy, p->originz,
                            p->endx, p->endy, p->endz);

    for(size_t nb=0; nb<fb_obj.size(); ++nb)
    {
        Eigen::Vector3d c, v, w;
        Eigen::Vector4d e;
        fb_obj[nb]->get_pose_vel(c, e, v, w);

        const Eigen::Matrix3d &I = fb_obj[nb]->inertia_body();
        double cc[3] = {c(0), c(1), c(2)};
        double ee[4] = {e(0), e(1), e(2), e(3)};
        double vv[3] = {v(0), v(1), v(2)};
        double ww[3] = {w(0), w(1), w(2)};

        double dx, dy, dz;
        fb_obj[nb]->collision_extents(dx, dy, dz);

        std::vector<double> xyz;
        const int ntri = fb_obj[nb]->copy_collision_mesh(xyz);
        int added = -1;
        if(ntri > 0)
        added = reef3d_chrono_add_mesh(pimpl->backend, ntri, xyz.data(), fb_obj[nb]->Mass_fb,
                                       I(0,0), I(1,1), I(2,2), I(0,1), I(0,2), I(1,2),
                                       cc, ee, vv, ww);

        if(added < 0)
        {
            added = reef3d_chrono_add_box(pimpl->backend, dx, dy, dz, fb_obj[nb]->Mass_fb,
                                          I(0,0), I(1,1), I(2,2), I(0,1), I(0,2), I(1,2),
                                          cc, ee, vv, ww);
            if(p->mpirank==0)
            cout<<"Chrono body "<<nb<<"  AABB contact  size "<<dx<<" "<<dy<<" "<<dz
                <<"  mass "<<fb_obj[nb]->Mass_fb
                <<"  CoG "<<c(0)<<" "<<c(1)<<" "<<c(2)
                <<"  U0 "<<v(0)<<" "<<v(1)<<" "<<v(2)<<endl;
        }
        else if(p->mpirank==0)
        cout<<"Chrono body "<<nb<<"  triangle mesh  ntri "<<ntri
            <<"  AABB "<<dx<<" "<<dy<<" "<<dz
            <<"  mass "<<fb_obj[nb]->Mass_fb
            <<"  CoG "<<c(0)<<" "<<c(1)<<" "<<c(2)
            <<"  U0 "<<v(0)<<" "<<v(1)<<" "<<v(2)<<endl;
    }

    reef3d_chrono_set_locks(pimpl->backend, p->X11_u, p->X11_v, p->X11_w,
                            p->X11_p, p->X11_q, p->X11_r);

    pimpl->ready = true;
    if(p->mpirank==0)
    {
        if(p->X16==1)
        cout<<"Chrono SMC soft contact (Hooke, kn=2e5 N/m)  bodies: "<<fb_obj.size()<<endl;
        else
        cout<<"Chrono NSC hard contact initialized  bodies: "<<fb_obj.size()<<endl;
    }
#endif
}

void sixdof_chrono::advance(lexer *p, fdm *a, ghostcell *pgc, std::vector<sixdof_obj*> &fb_obj, int iter, bool finalize)
{
    if(!pimpl->ready)
    return;

    for(size_t nb=0; nb<fb_obj.size(); ++nb)
    {
        fb_obj[nb]->prepare_external_loads(p, a, pgc, finalize);

        Eigen::Vector3d F, M;
        fb_obj[nb]->get_wrench(F, M);
        double FF[3] = {F(0), F(1), F(2)};
        double MM[3] = {M(0), M(1), M(2)};
        reef3d_chrono_set_wrench(pimpl->backend, int(nb), FF, MM);
    }

    if(finalize==true)
    {
        const double dt = p->dt;
        const double cdt = std::min(dt, 1.0e-4);
        int nsub = int(dt/cdt + 0.5);
        if(nsub<1)
        nsub=1;

        reef3d_chrono_step(pimpl->backend, dt, nsub);

        const int nc = reef3d_chrono_ncontacts(pimpl->backend);
        if(p->mpirank==0 && (p->count%20==0 || nc>0))
        cout<<"Chrono contacts: "<<nc<<"  substeps: "<<nsub<<endl;
    }

    for(size_t nb=0; nb<fb_obj.size(); ++nb)
    {
        double c[3], e[4], v[3], w[3];
        reef3d_chrono_get_state(pimpl->backend, int(nb), c, e, v, w);
        fb_obj[nb]->apply_kinematics(p,
            Eigen::Vector3d(c[0], c[1], c[2]),
            Eigen::Vector4d(e[0], e[1], e[2], e[3]),
            Eigen::Vector3d(v[0], v[1], v[2]),
            Eigen::Vector3d(w[0], w[1], w[2]));
    }
}

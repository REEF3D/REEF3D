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

/* Isolated Chrono backend: do not include REEF3D headers here. */

#ifdef REEF3D_CHRONO

#include "6DOF_chrono_backend.h"

#include "chrono/collision/ChCollisionModel.h"
#include "chrono/collision/ChCollisionShapeBox.h"
#include "chrono/collision/ChCollisionShapeTriangleMesh.h"
#include "chrono/core/ChFrame.h"
#include "chrono/geometry/ChTriangleMeshConnected.h"
#include "chrono/physics/ChBody.h"
#include "chrono/physics/ChBodyEasy.h"
#include "chrono/physics/ChContactMaterialNSC.h"
#include "chrono/physics/ChContactMaterialSMC.h"
#include "chrono/physics/ChSystem.h"
#include "chrono/physics/ChSystemNSC.h"
#include "chrono/physics/ChSystemSMC.h"

#include <algorithm>
#include <cmath>
#include <memory>
#include <string>
#include <vector>

namespace {

struct Backend
{
    std::unique_ptr<::chrono::ChSystem> sys;
    std::vector<std::shared_ptr<::chrono::ChBody> > bodies;
    std::vector<unsigned int> acc;
    std::vector<::chrono::ChVector3d> c0;
    std::shared_ptr<::chrono::ChContactMaterial> mat;
    int free_u, free_v, free_w, free_p, free_q, free_r;

    Backend() : free_u(1), free_v(1), free_w(1), free_p(1), free_q(1), free_r(1)
    {
    }

    void apply_locks()
    {
        for(size_t nb=0; nb<bodies.size(); ++nb)
        {
            auto& body = bodies[nb];
            ::chrono::ChVector3d p = body->GetPos();
            ::chrono::ChVector3d v = body->GetPosDt();
            if(!free_u) { p.x() = c0[nb].x(); v.x() = 0.0; }
            if(!free_v) { p.y() = c0[nb].y(); v.y() = 0.0; }
            if(!free_w) { p.z() = c0[nb].z(); v.z() = 0.0; }
            body->SetPos(p);
            body->SetPosDt(v);

            ::chrono::ChVector3d w = body->GetAngVelParent();
            if(!free_p) w.x() = 0.0;
            if(!free_q) w.y() = 0.0;
            if(!free_r) w.z() = 0.0;
            body->SetAngVelParent(w);

            if(!free_p && !free_q)
            {
                const ::chrono::ChQuaterniond q = body->GetRot();
                const double e0=q.e0(), e1=q.e1(), e2=q.e2(), e3=q.e3();
                const double psi = std::atan2(2.0*(e1*e2 + e3*e0),
                                             1.0 - 2.0*(e2*e2 + e3*e3));
                if(!free_r)
                    body->SetRot(::chrono::ChQuaterniond(1.0, 0.0, 0.0, 0.0));
                else
                    body->SetRot(::chrono::ChQuaterniond(std::cos(0.5*psi), 0.0, 0.0, std::sin(0.5*psi)));
            }
        }
    }
};

}  // namespace

extern "C" {

void* reef3d_chrono_create(int smc)
{
    Backend* b = new Backend();

    if(smc)
    {
        auto s = std::make_unique<::chrono::ChSystemSMC>("reef3d_chrono");
        s->UseMaterialProperties(false);
        s->SetContactForceModel(::chrono::ChSystemSMC::ContactForceModel::Hooke);
        auto m = chrono_types::make_shared<::chrono::ChContactMaterialSMC>();
        m->SetFriction(0.2f);
        m->SetRestitution(0.0f);
        m->SetKn(2.0e5f);
        m->SetKt(8.0e4f);
        m->SetGn(3.0e3f);
        m->SetGt(1.2e3f);
        b->mat = m;
        b->sys = std::move(s);
    }
    else
    {
        auto s = std::make_unique<::chrono::ChSystemNSC>("reef3d_chrono");
        auto m = chrono_types::make_shared<::chrono::ChContactMaterialNSC>();
        m->SetFriction(0.2f);
        m->SetRestitution(0.0f);
        b->mat = m;
        b->sys = std::move(s);
    }

    b->sys->SetGravitationalAcceleration(::chrono::ChVector3d(0.0, 0.0, 0.0));
    b->sys->SetCollisionSystemType(::chrono::ChCollisionSystem::Type::BULLET);
    ::chrono::ChCollisionModel::SetDefaultSuggestedEnvelope(0.008);
    ::chrono::ChCollisionModel::SetDefaultSuggestedMargin(0.008);
    return b;
}

void reef3d_chrono_destroy(void* ptr)
{
    delete static_cast<Backend*>(ptr);
}

void reef3d_chrono_add_floor(void* ptr, double ox, double oy, double oz,
                             double ex, double ey, double ez)
{
    Backend* b = static_cast<Backend*>(ptr);
    const double Lx = ex - ox;
    const double Ly = ey - oy;
    const double thick = 0.2;
    auto floor = chrono_types::make_shared<::chrono::ChBodyEasyBox>(
        Lx + 2.0, Ly + 2.0, thick, 1000.0, false, true, b->mat);
    floor->SetPos(::chrono::ChVector3d(0.5 * (ox + ex), 0.5 * (oy + ey), oz - 0.5 * thick));
    floor->SetFixed(true);
    floor->SetName("floor");
    b->sys->Add(floor);
}

int reef3d_chrono_add_box(void* ptr,
                          double dx, double dy, double dz, double mass,
                          double Ixx, double Iyy, double Izz,
                          double Ixy, double Ixz, double Iyz,
                          const double c[3], const double e[4],
                          const double v[3], const double w[3])
{
    Backend* b = static_cast<Backend*>(ptr);
    const double vol = dx * dy * dz;
    const double density = (vol > 1.0e-12) ? mass / vol : 500.0;

    auto body = chrono_types::make_shared<::chrono::ChBodyEasyBox>(
        dx, dy, dz, density, false, true, b->mat);
    body->SetName(std::string("body") + std::to_string(int(b->bodies.size())));
    body->SetSleepingAllowed(false);
    body->SetMass(mass);
    body->SetInertiaXX(::chrono::ChVector3d(Ixx, Iyy, Izz));
    body->SetInertiaXY(::chrono::ChVector3d(Ixy, Ixz, Iyz));
    body->SetPos(::chrono::ChVector3d(c[0], c[1], c[2]));
    body->SetRot(::chrono::ChQuaterniond(e[0], e[1], e[2], e[3]));
    body->SetPosDt(::chrono::ChVector3d(v[0], v[1], v[2]));
    body->SetAngVelParent(::chrono::ChVector3d(w[0], w[1], w[2]));

    unsigned int idx = body->AddAccumulator();
    b->sys->Add(body);
    b->bodies.push_back(body);
    b->acc.push_back(idx);
    b->c0.push_back(::chrono::ChVector3d(c[0], c[1], c[2]));
    return int(b->bodies.size() - 1);
}

int reef3d_chrono_add_mesh(void* ptr,
                           int ntri, const double* xyz,
                           double mass,
                           double Ixx, double Iyy, double Izz,
                           double Ixy, double Ixz, double Iyz,
                           const double c[3], const double e[4],
                           const double v[3], const double w[3])
{
    if(ntri < 1 || xyz == 0)
        return -1;

    Backend* b = static_cast<Backend*>(ptr);
    auto mesh = chrono_types::make_shared<::chrono::ChTriangleMeshConnected>();

    double xmin = 1.0e20, xmax = -1.0e20;
    double ymin = 1.0e20, ymax = -1.0e20;
    double zmin = 1.0e20, zmax = -1.0e20;

    for(int n=0; n<ntri; ++n)
    {
        const double* t = xyz + n*9;
        const ::chrono::ChVector3d v0(t[0], t[1], t[2]);
        const ::chrono::ChVector3d v1(t[3], t[4], t[5]);
        const ::chrono::ChVector3d v2(t[6], t[7], t[8]);
        mesh->AddTriangle(v0, v1, v2);

        xmin = std::min({xmin, v0.x(), v1.x(), v2.x()});
        xmax = std::max({xmax, v0.x(), v1.x(), v2.x()});
        ymin = std::min({ymin, v0.y(), v1.y(), v2.y()});
        ymax = std::max({ymax, v0.y(), v1.y(), v2.y()});
        zmin = std::min({zmin, v0.z(), v1.z(), v2.z()});
        zmax = std::max({zmax, v0.z(), v1.z(), v2.z()});
    }

    mesh->RepairDuplicateVertices(1.0e-9);

    const double charlen = std::min({std::max(xmax-xmin, 1.0e-4),
                                     std::max(ymax-ymin, 1.0e-4),
                                     std::max(zmax-zmin, 1.0e-4)});
    const double swept = std::min(8.0e-3, std::max(1.0e-5, 0.02*charlen));

    auto body = chrono_types::make_shared<::chrono::ChBody>();
    body->SetName(std::string("body") + std::to_string(int(b->bodies.size())));
    body->SetSleepingAllowed(false);
    body->SetMass(mass);
    body->SetInertiaXX(::chrono::ChVector3d(Ixx, Iyy, Izz));
    body->SetInertiaXY(::chrono::ChVector3d(Ixy, Ixz, Iyz));
    body->SetPos(::chrono::ChVector3d(c[0], c[1], c[2]));
    body->SetRot(::chrono::ChQuaterniond(e[0], e[1], e[2], e[3]));
    body->SetPosDt(::chrono::ChVector3d(v[0], v[1], v[2]));
    body->SetAngVelParent(::chrono::ChVector3d(w[0], w[1], w[2]));

    auto cshape = chrono_types::make_shared<::chrono::ChCollisionShapeTriangleMesh>(
        b->mat, mesh, false, false, swept);
    body->AddCollisionShape(cshape);
    body->EnableCollision(true);

    unsigned int idx = body->AddAccumulator();
    b->sys->Add(body);
    b->bodies.push_back(body);
    b->acc.push_back(idx);
    b->c0.push_back(::chrono::ChVector3d(c[0], c[1], c[2]));
    return int(b->bodies.size() - 1);
}

void reef3d_chrono_set_wrench(void* ptr, int nb, const double F[3], const double M[3])
{
    Backend* b = static_cast<Backend*>(ptr);
    auto& body = b->bodies[size_t(nb)];
    const unsigned int idx = b->acc[size_t(nb)];
    body->EmptyAccumulator(idx);
    body->AccumulateForce(idx, ::chrono::ChVector3d(F[0], F[1], F[2]), body->GetPos(), false);
    body->AccumulateTorque(idx, ::chrono::ChVector3d(M[0], M[1], M[2]), false);
}

void reef3d_chrono_set_locks(void* ptr, int free_u, int free_v, int free_w,
                            int free_p, int free_q, int free_r)
{
    Backend* b = static_cast<Backend*>(ptr);
    b->free_u = free_u;
    b->free_v = free_v;
    b->free_w = free_w;
    b->free_p = free_p;
    b->free_q = free_q;
    b->free_r = free_r;
}

void reef3d_chrono_step(void* ptr, double dt, int nsub)
{
    Backend* b = static_cast<Backend*>(ptr);
    if(nsub < 1)
        nsub = 1;
    const double h = dt / double(nsub);
    for(int s = 0; s < nsub; ++s)
    {
        b->sys->DoStepDynamics(h);
        b->apply_locks();
    }
}

void reef3d_chrono_get_state(void* ptr, int nb, double c[3], double e[4], double v[3], double w[3])
{
    Backend* b = static_cast<Backend*>(ptr);
    auto& body = b->bodies[size_t(nb)];
    const ::chrono::ChVector3d pos = body->GetPos();
    const ::chrono::ChQuaterniond rot = body->GetRot();
    const ::chrono::ChVector3d vel = body->GetPosDt();
    const ::chrono::ChVector3d omg = body->GetAngVelParent();
    c[0] = pos.x(); c[1] = pos.y(); c[2] = pos.z();
    e[0] = rot.e0(); e[1] = rot.e1(); e[2] = rot.e2(); e[3] = rot.e3();
    v[0] = vel.x(); v[1] = vel.y(); v[2] = vel.z();
    w[0] = omg.x(); w[1] = omg.y(); w[2] = omg.z();
}

int reef3d_chrono_ncontacts(void* ptr)
{
    Backend* b = static_cast<Backend*>(ptr);
    return int(b->sys->GetNumContacts());
}

}

#else

#include "6DOF_chrono_backend.h"

extern "C" {

void* reef3d_chrono_create(int) { return 0; }
void  reef3d_chrono_destroy(void*) {}
void  reef3d_chrono_add_floor(void*, double, double, double, double, double, double) {}
int   reef3d_chrono_add_box(void*, double, double, double, double,
                            double, double, double, double, double, double,
                            const double*, const double*, const double*, const double*) { return -1; }
int   reef3d_chrono_add_mesh(void*, int, const double*, double,
                             double, double, double, double, double, double,
                             const double*, const double*, const double*, const double*) { return -1; }
void  reef3d_chrono_set_wrench(void*, int, const double*, const double*) {}
void  reef3d_chrono_set_locks(void*, int, int, int, int, int, int) {}
void  reef3d_chrono_step(void*, double, int) {}
void  reef3d_chrono_get_state(void*, int, double*, double*, double*, double*) {}
int   reef3d_chrono_ncontacts(void*) { return 0; }

}

#endif

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

#ifndef SIXDOF_CHRONO_BACKEND_H_
#define SIXDOF_CHRONO_BACKEND_H_

#ifdef __cplusplus
extern "C" {
#endif

void* reef3d_chrono_create(int smc);
void  reef3d_chrono_destroy(void*);
void  reef3d_chrono_add_floor(void*, double ox, double oy, double oz,
                              double ex, double ey, double ez);
int   reef3d_chrono_add_box(void*,
                            double dx, double dy, double dz, double mass,
                            double Ixx, double Iyy, double Izz,
                            double Ixy, double Ixz, double Iyz,
                            const double c[3], const double e[4],
                            const double v[3], const double w[3]);
int   reef3d_chrono_add_mesh(void*,
                             int ntri, const double* xyz,
                             double mass,
                             double Ixx, double Iyy, double Izz,
                             double Ixy, double Ixz, double Iyz,
                             const double c[3], const double e[4],
                             const double v[3], const double w[3]);
void  reef3d_chrono_set_wrench(void*, int nb, const double F[3], const double M[3]);
void  reef3d_chrono_set_locks(void*, int free_u, int free_v, int free_w,
                             int free_p, int free_q, int free_r);
void  reef3d_chrono_step(void*, double dt, int nsub);
void  reef3d_chrono_get_state(void*, int nb, double c[3], double e[4], double v[3], double w[3]);
int   reef3d_chrono_ncontacts(void*);

#ifdef __cplusplus
}
#endif

#endif

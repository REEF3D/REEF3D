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
// REEF3D::CORAL: the seabed under mooring lines and points that is not flat: a plane (a slope) or a
// regular grid of heights (bilinear), with the unit normal pointing up. The resolved fluid's bed at the
// nodes (bathymetry fluid) comes from the sampled kinematics instead (wave_kinematics::bed_at).
//
// Grid file: first the header 'nx ny x0 y0 dx dy', then nx*ny heights z(x0 + i dx, y0 + j dy), i fastest
// ('#' comments); outside the grid the height of its edge.

#ifndef CORAL_SEABED_H_
#define CORAL_SEABED_H_

#include"coral/core/coral_types.h"
#include<string>
#include<vector>

namespace coral
{
    class seabed_surface
    {
    public:
        virtual ~seabed_surface() = default;
        virtual double z(double x, double y) const = 0;
        // unit normal (up) at (x, y)
        virtual Vec3 normal(double x, double y) const = 0;
    };

    // z = z0 + sx x + sy y
    class seabed_plane : public seabed_surface
    {
    public:
        seabed_plane(double z0, double sx, double sy) : z0_(z0), sx_(sx), sy_(sy) {}
        double z(double x, double y) const override {return z0_ + sx_*x + sy_*y;}
        Vec3 normal(double, double) const override {return Vec3(-sx_, -sy_, 1.0).normalized();}
    private:
        double z0_, sx_, sy_;
    };

    class seabed_grid : public seabed_surface
    {
    public:
        seabed_grid(int nx, int ny, double x0, double y0, double dx, double dy, std::vector<double> z);
        double z(double x, double y) const override;
        Vec3 normal(double x, double y) const override;
        // the grid file (see above); throws std::runtime_error
        static seabed_grid read(const std::string &file);
    private:
        // cell and local coordinates of (x, y), clamped to the grid
        void cell(double x, double y, int &i, int &j, double &fx, double &fy, bool &inx, bool &iny) const;
        int nx_, ny_;
        double x0_, y0_, dx_, dy_;
        std::vector<double> z_;
    };
}

#endif

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

#include"coral/mooring/coral_seabed.h"
#include<algorithm>
#include<cmath>
#include<fstream>
#include<sstream>
#include<stdexcept>

namespace coral
{

seabed_grid::seabed_grid(int nx, int ny, double x0, double y0, double dx, double dy, std::vector<double> z)
    : nx_(nx), ny_(ny), x0_(x0), y0_(y0), dx_(dx), dy_(dy), z_(std::move(z))
{
    if(nx_ < 1 || ny_ < 1 || (int)z_.size() != nx_*ny_ || (nx_ > 1 && !(dx_ > 0.0)) || (ny_ > 1 && !(dy_ > 0.0)))
        throw std::runtime_error("coral::seabed_grid: nx, ny >= 1, dx, dy > 0 and nx*ny heights");
}

void seabed_grid::cell(double x, double y, int &i, int &j, double &fx, double &fy, bool &inx, bool &iny) const
{
    auto one = [](double u, double u0, double du, int n, int &k, double &f, bool &in)
    {
        if(n < 2) {k = 0; f = 0.0; in = false; return;}
        const double s = (u - u0)/du;
        in = s >= 0.0 && s <= n - 1;
        const double c = std::clamp(s, 0.0, double(n - 1));
        k = std::min((int)std::floor(c), n - 2);
        f = c - k;
    };
    one(x, x0_, dx_, nx_, i, fx, inx);
    one(y, y0_, dy_, ny_, j, fy, iny);
}

double seabed_grid::z(double x, double y) const
{
    int i, j;
    double fx, fy;
    bool inx, iny;
    cell(x, y, i, j, fx, fy, inx, iny);
    auto Z = [&](int a, int b) {return z_[std::min(b, ny_-1)*nx_ + std::min(a, nx_-1)];};
    return (1.0 - fx)*(1.0 - fy)*Z(i, j) + fx*(1.0 - fy)*Z(i+1, j) + (1.0 - fx)*fy*Z(i, j+1) + fx*fy*Z(i+1, j+1);
}

Vec3 seabed_grid::normal(double x, double y) const
{
    int i, j;
    double fx, fy;
    bool inx, iny;
    cell(x, y, i, j, fx, fy, inx, iny);
    auto Z = [&](int a, int b) {return z_[std::min(b, ny_-1)*nx_ + std::min(a, nx_-1)];};
    // gradient of the bilinear height (flat outside the grid, beyond its edges)
    const double gx = (inx && nx_ > 1) ? ((1.0 - fy)*(Z(i+1, j) - Z(i, j)) + fy*(Z(i+1, j+1) - Z(i, j+1)))/dx_ : 0.0;
    const double gy = (iny && ny_ > 1) ? ((1.0 - fx)*(Z(i, j+1) - Z(i, j)) + fx*(Z(i+1, j+1) - Z(i+1, j)))/dy_ : 0.0;
    return Vec3(-gx, -gy, 1.0).normalized();
}

seabed_grid seabed_grid::read(const std::string &file)
{
    std::ifstream in(file);
    if(!in)
        throw std::runtime_error("coral::seabed_grid: cannot open " + file);
    std::stringstream ss;
    std::string line;
    while(std::getline(in, line))
    {
        const size_t c = line.find('#');
        ss<<(c==std::string::npos ? line : line.substr(0, c))<<"\n";
    }
    int nx = 0, ny = 0;
    double x0, y0, dx, dy;
    if(!(ss>>nx>>ny>>x0>>y0>>dx>>dy))
        throw std::runtime_error("coral::seabed_grid: " + file + ": header 'nx ny x0 y0 dx dy' missing");
    std::vector<double> z(std::max(nx, 0)*std::max(ny, 0));
    for(double &v : z)
        if(!(ss>>v))
            throw std::runtime_error("coral::seabed_grid: " + file + ": fewer than nx*ny heights");
    return seabed_grid(nx, ny, x0, y0, dx, dy, std::move(z));
}

}

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

#include"coral/wind/coral_wind_grid.h"
#include<unsupported/Eigen/FFT>
#include<algorithm>
#include<cmath>
#include<complex>
#include<cstdint>
#include<fstream>
#include<random>

namespace coral
{

Vec3 wind_grid::level(int it, double eta, double z) const
{
    const size_t base = (size_t)it*nz*ny*3;
    double fy = ny > 1 ? (eta - eta0())/dy : 0.0;
    fy = std::clamp(fy, 0.0, (double)(ny - 1));
    const int j0 = std::min((int)fy, std::max(ny - 2, 0)), j1 = std::min(j0 + 1, ny - 1);
    const double sy = fy - j0;
    auto at = [&](int i, int j, int c) {return (double)data[base + ((size_t)i*ny + j)*3 + c];};
    if(z < z_bottom && ntower > 0)
    {
        // centre line: the bottom row at the hub's eta (interpolated), then the tower points dz apart, then zero at z = 0
        const double fc = std::clamp(ny > 1 ? -eta0()/dy : 0.0, 0.0, (double)(ny - 1));
        const int c0 = std::min((int)fc, std::max(ny - 2, 0)), c1 = std::min(c0 + 1, ny - 1);
        const double sc = fc - c0;
        auto zt = [&](int k) {return z_bottom - k*dz;};      // k = 0: the grid's bottom row
        auto vt = [&](int k)
        {
            Vec3 v;
            for(int c=0; c<3; ++c)
                v(c) = k==0 ? (1.0 - sc)*at(0, c0, c) + sc*at(0, c1, c) : (double)tower[((size_t)it*ntower + (k - 1))*3 + c];
            return v;
        };
        const double fk = (z_bottom - z)/dz;
        const int k0 = (int)fk;
        if(k0 >= ntower)
        {
            const double zl = zt(ntower);
            return zl > 0.0 ? Vec3(vt(ntower)*std::max(z, 0.0)/zl) : vt(ntower);
        }
        const double s = (zt(k0) - z)/dz;
        return (1.0 - s)*vt(k0) + s*vt(k0 + 1);
    }
    double fz = nz > 1 ? (z - z_bottom)/dz : 0.0;
    fz = std::clamp(fz, 0.0, (double)(nz - 1));
    const int i0 = std::min((int)fz, std::max(nz - 2, 0)), i1 = std::min(i0 + 1, nz - 1);
    const double sz = fz - i0;
    Vec3 v;
    for(int c=0; c<3; ++c)
        v(c) = (1.0 - sz)*((1.0 - sy)*at(i0, j0, c) + sy*at(i0, j1, c)) + sz*((1.0 - sy)*at(i1, j0, c) + sy*at(i1, j1, c));
    return v;
}

Vec3 wind_grid::local(double eta, double z, double tau) const
{
    if(nt <= 0)
        return Vec3::Zero();
    double ft = tau/dt;
    if(periodic)
    {
        ft = std::fmod(ft, (double)nt);
        if(ft < 0.0) ft += nt;
        const int a = std::min((int)ft, nt - 1), b = (a + 1) % nt;
        const double s = ft - a;
        return (1.0 - s)*level(a, eta, z) + s*level(b, eta, z);
    }
    ft = std::clamp(ft, 0.0, (double)(nt - 1));
    const int a = std::min((int)ft, std::max(nt - 2, 0)), b = std::min(a + 1, nt - 1);
    const double s = ft - a;
    return (1.0 - s)*level(a, eta, z) + s*level(b, eta, z);
}

double wind_grid::tau(const Vec3 &p, double t) const
{
    const double xi = std::cos(direction)*(p(0) - origin(0)) + std::sin(direction)*(p(1) - origin(1));
    return t + (xi0 - xi)/std::max(u_hub, 1.0e-6);
}

Vec3 wind_grid::velocity(const Vec3 &p, double t) const
{
    const double c = std::cos(direction), s = std::sin(direction);
    const double eta = -s*(p(0) - origin(0)) + c*(p(1) - origin(1));
    const Vec3 f = local(eta, p(2), tau(p, t));
    return Vec3(f(0)*c - f(1)*s, f(0)*s + f(1)*c, f(2));
}

namespace
{
    template<class T> bool get(std::ifstream &in, T &x) {return (bool)in.read(reinterpret_cast<char*>(&x), sizeof(T));}
    template<class T> void put(std::ofstream &out, T x) {out.write(reinterpret_cast<const char*>(&x), sizeof(T));}
}

bool wind_grid::read_bts(const std::string &file, std::string &msg)
{
    std::ifstream in(file, std::ios::binary);
    if(!in)
    {
        msg = "cannot open '" + file + "'";
        return false;
    }
    int16_t id;
    int32_t n_z, n_y, n_tw, n_t;
    float f_dz, f_dy, f_dt, f_u, f_zh, f_zb, so[6];
    int32_t nc;
    bool ok = get(in, id) && get(in, n_z) && get(in, n_y) && get(in, n_tw) && get(in, n_t) && get(in, f_dz) && get(in, f_dy)
              && get(in, f_dt) && get(in, f_u) && get(in, f_zh) && get(in, f_zb);
    for(int i=0; ok && i<6; ++i) ok = get(in, so[i]);
    ok = ok && get(in, nc);
    if(!ok || (id != 7 && id != 8) || n_z < 1 || n_y < 1 || n_tw < 0 || n_t < 1 || nc < 0 || f_dt <= 0.0f)
    {
        msg = file + ": not a TurbSim full-field file (.bts)";
        return false;
    }
    in.seekg(nc, std::ios::cur);
    nz = n_z; ny = n_y; ntower = n_tw; nt = n_t;
    dz = f_dz; dy = f_dy; dt = f_dt; u_hub = f_u; z_hub = f_zh; z_bottom = f_zb;
    periodic = id==8;
    has_mean = true;
    data.assign((size_t)nt*nz*ny*3, 0.0f);
    tower.assign((size_t)nt*ntower*3, 0.0f);
    std::vector<int16_t> buf((size_t)(nz*ny + ntower)*3);
    for(int it=0; it<nt; ++it)
    {
        if(!in.read(reinterpret_cast<char*>(buf.data()), (std::streamsize)(buf.size()*sizeof(int16_t))))
        {
            msg = file + ": ends at time level " + std::to_string(it) + " of " + std::to_string(nt);
            return false;
        }
        for(size_t q=0; q<(size_t)nz*ny; ++q)
            for(int c=0; c<3; ++c)
                data[((size_t)it*nz*ny + q)*3 + c] = (float)((buf[q*3 + c] - so[2*c + 1])/so[2*c]);
        for(int k=0; k<ntower; ++k)
            for(int c=0; c<3; ++c)
                tower[((size_t)it*ntower + k)*3 + c] = (float)((buf[((size_t)nz*ny + k)*3 + c] - so[2*c + 1])/so[2*c]);
    }
    return true;
}

bool wind_grid::write_bts(const std::string &file, const std::string &description, std::string &msg,
                          const wind_field *mean, double t_mean) const
{
    // the velocities to write (with the mean wind of a fluctuation grid)
    std::vector<double> v(data.begin(), data.end()), vt(tower.begin(), tower.end());
    if(mean)
        for(int i=0; i<nz; ++i)
            for(int j=0; j<ny; ++j)
            {
                const double U = mean->speed(eta0() + j*dy, z_bottom + i*dz, t_mean);
                for(int it=0; it<nt; ++it)
                    v[(((size_t)it*nz + i)*ny + j)*3] += U;
            }
    std::ofstream out(file, std::ios::binary);
    if(!out)
    {
        msg = "cannot write '" + file + "'";
        return false;
    }
    float so[6];
    for(int c=0; c<3; ++c)
    {
        double lo = 1e300, hi = -1e300;
        for(size_t q=c; q<v.size(); q+=3) {lo = std::min(lo, v[q]); hi = std::max(hi, v[q]);}
        for(size_t q=c; q<vt.size(); q+=3) {lo = std::min(lo, vt[q]); hi = std::max(hi, vt[q]);}
        const double scl = hi > lo ? 65535.0/(hi - lo) : 1.0;
        so[2*c] = (float)scl;
        so[2*c + 1] = (float)(-32768.0 - scl*lo);
    }
    put<int16_t>(out, periodic ? 8 : 7);
    put<int32_t>(out, nz); put<int32_t>(out, ny); put<int32_t>(out, ntower); put<int32_t>(out, nt);
    put<float>(out, (float)dz); put<float>(out, (float)dy); put<float>(out, (float)dt);
    put<float>(out, (float)(mean ? mean->hub_speed(t_mean) : u_hub)); put<float>(out, (float)z_hub); put<float>(out, (float)z_bottom);
    for(int i=0; i<6; ++i) put<float>(out, so[i]);
    put<int32_t>(out, (int32_t)description.size());
    out.write(description.data(), (std::streamsize)description.size());
    auto q16 = [&](double x, int c)
    {
        const double y = std::nearbyint(so[2*c]*x + so[2*c + 1]);
        return (int16_t)std::clamp(y, -32768.0, 32767.0);
    };
    std::vector<int16_t> buf((size_t)(nz*ny + ntower)*3);
    for(int it=0; it<nt; ++it)
    {
        for(size_t q=0; q<(size_t)nz*ny; ++q)
            for(int c=0; c<3; ++c)
                buf[q*3 + c] = q16(v[((size_t)it*nz*ny + q)*3 + c], c);
        for(int k=0; k<ntower; ++k)
            for(int c=0; c<3; ++c)
                buf[((size_t)nz*ny + k)*3 + c] = q16(vt[((size_t)it*ntower + k)*3 + c], c);
        out.write(reinterpret_cast<const char*>(buf.data()), (std::streamsize)(buf.size()*sizeof(int16_t)));
    }
    if(!out)
    {
        msg = "error writing '" + file + "'";
        return false;
    }
    return true;
}

double wind_grid::kaimal_spectrum(const iec_turbine &tc, int k, double f, double V, double sigma1)
{
    static constexpr double sr[3] = {1.0, 0.8, 0.5}, lr[3] = {8.1, 2.7, 0.66};
    const double s = sr[k]*sigma1, L = lr[k]*tc.Lambda1();
    return 4.0*s*s*L/V/std::pow(1.0 + 6.0*f*L/V, 5.0/3.0);
}

double wind_grid::iec_coherence(const iec_turbine &tc, double r, double f, double V)
{
    const double Lc = 8.1*tc.Lambda1();
    const double a = f*r/V, b = 0.12*r/Lc;
    return std::exp(-12.0*std::sqrt(a*a + b*b));
}

wind_grid wind_grid::kaimal(const iec_turbine &tc, double V, double sigma1, const kaimal_settings &ks)
{
    wind_grid g;
    g.ny = ks.ny; g.nz = ks.nz;
    g.dy = ks.ny > 1 ? ks.width/(ks.ny - 1) : 1.0;
    g.dz = ks.nz > 1 ? ks.height/(ks.nz - 1) : 1.0;
    g.z_hub = tc.z_hub;
    g.z_bottom = tc.z_hub - 0.5*(ks.nz - 1)*g.dz;
    g.u_hub = V;
    g.dt = ks.dt;
    g.periodic = true;
    g.has_mean = false;
    int N = 2;
    while(N*ks.dt < ks.T - 1.0e-9*ks.T) N *= 2;
    g.nt = N;
    const int P = g.ny*g.nz;
    std::vector<double> ye(P), ze(P);
    for(int i=0; i<g.nz; ++i)
        for(int j=0; j<g.ny; ++j)
        {
            ye[i*g.ny + j] = g.eta0() + j*g.dy;
            ze[i*g.ny + j] = g.z_bottom + i*g.dz;
        }
    const double rmin = std::min(g.ny > 1 ? g.dy : 1e30, g.nz > 1 ? g.dz : 1e30);
    const double df = 1.0/(N*ks.dt);
    std::mt19937 gen(ks.seed);
    auto phase = [&gen]() {return 2.0*M_PI*(gen()/4294967296.0);};      // portable (mt19937 output is standard)
    g.data.assign((size_t)N*P*3, 0.0f);
    Eigen::FFT<double> fft;
    std::vector<std::complex<double>> spec(N), sig(N);
    MatX C(P, P);
    std::vector<std::vector<std::complex<double>>> X(P, std::vector<std::complex<double>>(N/2));
    for(int c=0; c<3; ++c)
    {
        const bool coh = c==0 || ks.coherence_vw;
        for(int k=1; k<N/2; ++k)
        {
            const double f = k*df;
            const double S = kaimal_spectrum(tc, c, f, V, sigma1);
            Eigen::VectorXcd xi(P);
            for(int j=0; j<P; ++j)
                xi(j) = std::polar(1.0, phase());
            if(coh && P > 1 && iec_coherence(tc, rmin, f, V) > 1.0e-8)
            {
                for(int a=0; a<P; ++a)
                    for(int b=0; b<=a; ++b)
                    {
                        const double r = std::hypot(ye[a] - ye[b], ze[a] - ze[b]);
                        C(a, b) = C(b, a) = S*iec_coherence(tc, r, f, V);
                    }
                Eigen::LLT<MatX> llt(C);
                const VecX yr = llt.matrixL()*VecX(xi.real()), yi = llt.matrixL()*VecX(xi.imag());
                for(int j=0; j<P; ++j)
                    X[j][k] = std::sqrt(2.0*df)*std::complex<double>(yr(j), yi(j));
            }
            else
                for(int j=0; j<P; ++j)
                    X[j][k] = std::sqrt(2.0*S*df)*xi(j);
        }
        for(int j=0; j<P; ++j)
        {
            std::fill(spec.begin(), spec.end(), std::complex<double>(0.0, 0.0));
            for(int k=1; k<N/2; ++k)
            {
                spec[k] = 0.5*N*X[j][k];
                spec[N - k] = std::conj(spec[k]);
            }
            fft.inv(sig, spec);
            for(int it=0; it<N; ++it)
                g.data[((size_t)it*P + j)*3 + c] = (float)sig[it].real();
        }
    }
    if(ks.scale_hub)
    {
        int jh = 0;
        double best = 1e300;
        for(int j=0; j<P; ++j)
        {
            const double d = std::hypot(ye[j], ze[j] - tc.z_hub);
            if(d < best) {best = d; jh = j;}
        }
        static constexpr double sr[3] = {1.0, 0.8, 0.5};
        for(int c=0; c<3; ++c)
        {
            double m = 0.0, q = 0.0;
            for(int it=0; it<N; ++it)
            {
                const double x = g.data[((size_t)it*P + jh)*3 + c];
                m += x; q += x*x;
            }
            m /= N;
            const double s = std::sqrt(std::max(q/N - m*m, 1e-300));
            const double f = sr[c]*sigma1/s;
            for(size_t q2=c; q2<g.data.size(); q2+=3)
                g.data[q2] = (float)(g.data[q2]*f);
        }
    }
    return g;
}

}

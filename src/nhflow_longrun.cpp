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

#include"nhflow_longrun.h"
#include"reef3d_mpmd.h"
#include"nhflow_amr.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include"increment.h"
#include"reefckpt.h"
#include<cmath>
#include<fstream>
#include<iostream>
#include<iomanip>
#include<limits>
#include<mpi.h>
#include<sys/stat.h>

namespace
{
// constituents in the order of choice, speeds in degrees per hour
const char *cname[nhflow_longrun::NC] = {"M2","S2","K1","O1","M4","N2","K2","P1","Q1","MS4","MN4","M6","2N2"};
const double cspeed[nhflow_longrun::NC] = {28.9841042, 30.0000000, 15.0410686, 13.9430356, 57.9682084, 28.4397295,
                                           30.0821373, 14.9589314, 13.3986609, 58.9841042, 57.4238337, 86.9523127, 27.8953548};

// P 162: astronomical arguments: the mean longitudes s (moon), h (sun), p (lunar perigee), N (lunar node) and
// T = 15 deg per hour of UT counted from midnight (the convention of Foreman 1977 and t_tide/UTide, in which
// tau = T + h - s and K1 = tau + s + 90, O1 = tau - s - 90) in degrees at the time days (days since 2000-01-01
// 12:00 UTC). V0 + u agree with UTide's ut_FUV within 0.7 deg and f within 2.2 % for all constituents except
// 2N2 (5 deg, 12 %: UTide adds the perigee satellites, here f and u of M2)
void astro(double days, double &s, double &h, double &pp, double &N, double &T)
{
    const double c = days/36525.0;
    s = 218.3165 + 481267.8813*c;
    h = 280.4665 + 36000.7698*c;
    pp = 83.3532 + 4069.0137*c;
    N = 125.0445 - 1934.1363*c;
    double ut = std::fmod(days + 0.5, 1.0);
    if(ut<0.0)
    ut += 1.0;
    T = 360.0*ut;
}

// days since 2000-01-01 12:00 UTC of a calendar date (proleptic Gregorian)
double days2000(int y, int mo, int d, int hh, int mi)
{
    y -= mo<=2;
    const int era = (y>=0 ? y : y-399)/400;
    const int yoe = y - era*400;
    const int doy = (153*(mo + (mo>2 ? -3 : 9)) + 2)/5 + d - 1;
    const int doe = yoe*365 + yoe/4 - yoe/100 + doy;
    const long dn = (long)era*146097 + doe - 719468;     // days since 1970-01-01
    return double(dn - 10957) - 0.5 + (hh + mi/60.0)/24.0;
}

// nodal factor f, nodal angle u (deg) and the equilibrium argument V (deg) of constituent k (order of cname)
void nodal(int k, double s, double h, double pp, double N, double T, double &f, double &u, double &V)
{
    const double r = N*M_PI/180.0;
    const double fM2 = 1.0004 - 0.0373*cos(r) + 0.0002*cos(2*r), uM2 = -2.14*sin(r);
    const double fK1 = 1.0060 + 0.1150*cos(r) - 0.0088*cos(2*r) + 0.0006*cos(3*r);
    const double uK1 = -8.86*sin(r) + 0.68*sin(2*r) - 0.07*sin(3*r);
    const double fO1 = 1.0089 + 0.1871*cos(r) - 0.0147*cos(2*r) + 0.0014*cos(3*r);
    const double uO1 = 10.80*sin(r) - 1.34*sin(2*r) + 0.19*sin(3*r);
    const double fK2 = 1.0241 + 0.2863*cos(r) + 0.0083*cos(2*r) - 0.0015*cos(3*r);
    const double uK2 = -17.74*sin(r) + 0.68*sin(2*r) - 0.04*sin(3*r);
    const double tau = T + h - s;

    switch(k)
    {
        case 0:  f = fM2;         u = uM2;     V = 2*tau;                  break;   // M2
        case 1:  f = 1.0;         u = 0.0;     V = 2*T;                    break;   // S2
        case 2:  f = fK1;         u = uK1;     V = tau + s + 90.0;         break;   // K1
        case 3:  f = fO1;         u = uO1;     V = tau - s - 90.0;         break;   // O1
        case 4:  f = fM2*fM2;     u = 2*uM2;   V = 4*tau;                  break;   // M4
        case 5:  f = fM2;         u = uM2;     V = 2*tau - s + pp;         break;   // N2
        case 6:  f = fK2;         u = uK2;     V = 2*T + 2*h;              break;   // K2
        case 7:  f = 1.0;         u = 0.0;     V = T - h + 270.0;          break;   // P1
        case 8:  f = fO1;         u = uO1;     V = tau - 2*s + pp - 90.0;  break;   // Q1
        case 9:  f = fM2;         u = uM2;     V = 2*tau + 2*T;            break;   // MS4
        case 10: f = fM2*fM2;     u = 2*uM2;   V = 4*tau - s + pp;         break;   // MN4
        case 11: f = fM2*fM2*fM2; u = 3*uM2;   V = 6*tau;                  break;   // M6
        default: f = fM2;         u = uM2;     V = 2*tau - 2*s + 2*pp;     break;   // 2N2
    }
}

inline int lij(const lexer *q, int ii, int jj) { return (ii-q->imin)*q->jmax + (jj-q->jmin); }
inline int lijk(const lexer *q, int ii, int jj, int kk) { return (ii-q->imin)*q->jmax*q->kmax + (jj-q->jmin)*q->kmax + kk-q->kmin; }

// Cholesky solve of the symmetric positive definite n x n system A x = b (A row-major, overwritten)
bool chol(vector<double> &A, int n)
{
    for(int j=0; j<n; ++j)
    {
        double s = A[j*n+j];
        for(int k=0; k<j; ++k)
        s -= A[j*n+k]*A[j*n+k];
        if(s<=0.0)
        return false;
        A[j*n+j] = sqrt(s);
        for(int i=j+1; i<n; ++i)
        {
            double t = A[i*n+j];
            for(int k=0; k<j; ++k)
            t -= A[i*n+k]*A[j*n+k];
            A[i*n+j] = t/A[j*n+j];
        }
    }
    return true;
}

void chol_solve(const vector<double> &L, int n, double *x)
{
    for(int i=0; i<n; ++i)
    {
        double t = x[i];
        for(int k=0; k<i; ++k)
        t -= L[i*n+k]*x[k];
        x[i] = t/L[i*n+i];
    }
    for(int i=n-1; i>=0; --i)
    {
        double t = x[i];
        for(int k=i+1; k<n; ++k)
        t -= L[k*n+i]*x[k];
        x[i] = t/L[i*n+i];
    }
}
}

nhflow_longrun::nhflow_longrun(lexer *p) : t0(p->P160_ts), dtout(p->P160_dt), rayleigh(p->P161)
{
    NB = 2*NC + 1;
    NV = 3*NB + 6;
    M.assign(NB*NB,0.0);

    for(int k=0; k<NC; ++k)
    {
    omega.push_back(cspeed[k]*M_PI/180.0/3600.0);
    name.push_back(cname[k]);
    }

    tnext = dtout>0.0 ? t0 + dtout : 1.0e30;

    for(int k=0; k<3*NB+6; ++k)
    {
    char nm[16];
    snprintf(nm,sizeof(nm),"lr_%02d",k);
    fname.push_back(nm);
    }

    if(p->mpirank==0)
    mkdir("./REEF3D_NHFLOW_Tidal",0777);
}

nhflow_longrun::~nhflow_longrun()
{
}

void nhflow_longrun::add_grid(lexer *q, fdm_nhf *d, int i0, int i1, int j0, int j1)
{
    grid g;
    g.q = q;
    g.d = d;
    g.i0 = i0;
    g.i1 = i1;
    g.j0 = j0;
    g.j1 = j1;
    g.s.assign((size_t)q->imax*q->jmax*NV,0.0);

    const int b = 3*NB;
    for(int c=0; c<q->imax*q->jmax; ++c)
    {
    g.s[(size_t)c*NV + b+3] = -1.0e30;      // eta max
    g.s[(size_t)c*NV + b+4] =  1.0e30;      // eta min
    }

    G.push_back(g);
}

void nhflow_longrun::step(lexer *p, fdm_nhf *d, ghostcell *pgc, nhflow_amr *amr)
{
    if(p->simtime<t0)
    return;

    if(!ini)
    {
        add_grid(p,d,0,p->knox,0,p->knoy);

        if(amr!=nullptr)
        amr->longrun_grids(this);

        ini = true;

        if(p->mpirank==0)
        cout<<"NHFLOW P 160: long-run statistics from t = "<<p->simtime<<" s, maps every "<<dtout<<" s"<<endl;
    }

    // the state at the end of the step, weight dt
    const double w = p->dt;
    const double t = p->simtime;
    vector<double> phi(NB);

    phi[0] = 1.0;
    for(int k=0; k<NC; ++k)
    {
    phi[1+2*k] = cos(omega[k]*t);
    phi[2+2*k] = sin(omega[k]*t);
    }

    for(int a=0; a<NB; ++a)
    for(int b=0; b<NB; ++b)
    M[a*NB+b] += w*phi[a]*phi[b];

    Twin += w;

    for(grid &g : G)
    accumulate(g,&phi[0],w);

    const bool last = (p->simtime>=p->N41-1.0e-9) || (p->count>=p->N45);

    if(p->simtime>=tnext || last)
    {
        write(p,pgc);

        while(tnext<=p->simtime)
        tnext += dtout;
    }
}

void nhflow_longrun::accumulate(grid &g, const double *phi, double w)
{
    lexer *q = g.q;
    fdm_nhf *d = g.d;
    const int m = increment::marge;
    const int b = 3*NB;

    for(int ii=0; ii<q->knox; ++ii)
    for(int jj=0; jj<q->knoy; ++jj)
    {
        const int c = lij(q,ii,jj);

        if(q->flagslice4[c]<0)
        continue;

        double *s = &g.s[(size_t)c*NV];

        const bool wet = q->wet[c]==1 && d->WL(ii,jj)>q->A544;
        const double eta = d->eta(ii,jj);
        double u=0.0, v=0.0, qx=0.0, qy=0.0;

        if(wet)
        for(int kk=0; kk<q->knoz; ++kk)
        {
            const int n = lijk(q,ii,jj,kk);
            const double dz = q->DZN[kk+m];
            u += dz*d->U[n];
            v += dz*d->V[n];
            qx += dz*d->UH[n];
            qy += dz*d->VH[n];
        }

        for(int a=0; a<NB; ++a)
        {
        const double f = w*phi[a];
        s[a] += f*eta;
        s[NB+a] += f*u;
        s[2*NB+a] += f*v;
        }

        s[b] += w*qx;
        s[b+1] += w*qy;

        if(wet)
        {
        s[b+2] += w;
        s[b+3] = max(s[b+3],eta);
        s[b+4] = min(s[b+4],eta);
        s[b+5] = max(s[b+5],sqrt(u*u+v*v));
        }
    }
}

// the constituents the record of length T (s) resolves, in the order of choice
int nhflow_longrun::select(double T, vector<int> &sel) const
{
    sel.clear();
    const double Th = T/3600.0;

    for(int k=0; k<NC; ++k)
    {
        const double fk = cspeed[k]/360.0;         // cycles per hour

        if(fk*Th<rayleigh)
        continue;

        // separated from every constituent before it in the order of choice, chosen or not: a
        // constituent the record cannot split from a larger one is left out instead of taking its energy
        bool ok = true;
        for(int l=0; l<k; ++l)
        if(fabs(fk - cspeed[l]/360.0)*Th<rayleigh)
        ok = false;

        if(ok)
        sel.push_back(k);
    }

    return (int)sel.size();
}

void nhflow_longrun::write(lexer *p, ghostcell *pgc)
{
    vector<int> sel;
    const int ns = select(Twin,sel);
    const int n = 2*ns + 1;

    // the normal matrix of the chosen basis
    vector<int> idx(1,0);
    for(int k : sel)
    {
    idx.push_back(1+2*k);
    idx.push_back(2+2*k);
    }

    vector<double> L(n*n);
    for(int a=0; a<n; ++a)
    for(int b=0; b<n; ++b)
    L[a*n+b] = M[idx[a]*NB+idx[b]];

    const bool okc = chol(L,n);

    if(!okc && p->mpirank==0)
    cout<<"NHFLOW P 160: the normal matrix of the harmonic analysis is singular, maps without constituents"<<endl;

    const int nsel = okc ? ns : 0;
    const double nan = std::numeric_limits<double>::quiet_NaN();

    // file names
    int ng = (int)G.size();
    vector<int> all(p->mpi_size,0);
    MPI_Allgather(&ng,1,MPI_INT,&all[0],1,MPI_INT,reef3d_comm());

    vector<double> xe(n), xu(n), xv(n);

    // P 162: nodal factor and Greenwich phase offset (V0 + u) of the chosen constituents
    vector<double> nf(nsel,1.0), nvu(nsel,0.0);
    if(p->P162==1)
    {
        const double ds = days2000(p->P162_y,p->P162_mo,p->P162_d,p->P162_h,p->P162_mi);
        double s0,h0,p0,N0,T0, sm,hm,pm,Nm,Tm;
        astro(ds,s0,h0,p0,N0,T0);
        astro(ds + (t0 + 0.5*Twin)/86400.0,sm,hm,pm,Nm,Tm);
        for(int k=0; k<nsel; ++k)
        {
            double f,u,V,f0,u0,V0;
            nodal(sel[k],sm,hm,pm,Nm,Tm,f,u,V);
            nodal(sel[k],s0,h0,p0,N0,T0,f0,u0,V0);
            nf[k] = f;
            nvu[k] = std::fmod(V0 + u, 360.0);
        }
    }

    for(int gi=0; gi<ng; ++gi)
    {
        grid &g = G[gi];
        lexer *q = g.q;
        fdm_nhf *d = g.d;
        const int nx = g.i1-g.i0, ny = g.j1-g.j0;
        const int b = 3*NB;

        // the fields of the cells, in the order of the file
        vector<string> fname;
        vector<vector<double>> fval;
        fval.reserve(10 + 5*nsel);      // field() hands out references into fval

        auto field = [&](const string &nm) -> vector<double>&
        {
            fname.push_back(nm);
            fval.push_back(vector<double>((size_t)nx*ny,nan));
            return fval.back();
        };

        vector<double> &f_bed = field("bed");
        vector<double> &f_em = field("eta_mean");
        vector<double> &f_um = field("u_mean");
        vector<double> &f_vm = field("v_mean");
        vector<double> &f_qx = field("qx_mean");
        vector<double> &f_qy = field("qy_mean");
        vector<double> &f_emax = field("eta_max");
        vector<double> &f_emin = field("eta_min");
        vector<double> &f_smax = field("speed_max");
        vector<double> &f_wet = field("wet_fraction");

        vector<vector<double>*> fc;
        for(int k=0; k<nsel; ++k)
        {
        const string c = name[sel[k]];
        fc.push_back(&field(c+"_eta_amp"));
        fc.push_back(&field(c+"_eta_pha"));
        fc.push_back(&field(c+"_major"));
        fc.push_back(&field(c+"_minor"));
        fc.push_back(&field(c+"_incl"));
        }

        for(int jj=g.j0; jj<g.j1; ++jj)
        for(int ii=g.i0; ii<g.i1; ++ii)
        {
            const int c = lij(q,ii,jj);
            const size_t o = (size_t)(jj-g.j0)*nx + (ii-g.i0);

            if(q->flagslice4[c]<0)
            continue;

            const double *s = &g.s[(size_t)c*NV];

            f_bed[o] = d->bed(ii,jj);
            f_em[o] = s[0]/Twin;
            f_um[o] = s[NB]/Twin;
            f_vm[o] = s[2*NB]/Twin;
            f_qx[o] = s[b]/Twin;
            f_qy[o] = s[b+1]/Twin;
            f_wet[o] = s[b+2]/Twin;
            f_smax[o] = s[b+2]>0.0 ? s[b+5] : nan;
            f_emax[o] = s[b+2]>0.0 ? s[b+3] : nan;
            f_emin[o] = s[b+2]>0.0 ? s[b+4] : nan;

            if(nsel==0)
            continue;

            for(int a=0; a<n; ++a)
            {
            xe[a] = s[idx[a]];
            xu[a] = s[NB+idx[a]];
            xv[a] = s[2*NB+idx[a]];
            }

            chol_solve(L,n,&xe[0]);
            chol_solve(L,n,&xu[0]);
            chol_solve(L,n,&xv[0]);

            for(int k=0; k<nsel; ++k)
            {
                const double ae = xe[1+2*k], be = xe[2+2*k];
                const double au = xu[1+2*k], bu = xu[2+2*k];
                const double av = xv[1+2*k], bv = xv[2+2*k];

                double pha = atan2(be,ae)*180.0/M_PI;
                if(pha<0.0)
                pha += 360.0;

                // rotary components of (u + i v): anticlockwise wp, clockwise wm
                const double wpr = 0.5*(au+bv), wpi = 0.5*(av-bu);
                const double wmr = 0.5*(au-bv), wmi = 0.5*(av+bu);
                const double ap = sqrt(wpr*wpr+wpi*wpi), am = sqrt(wmr*wmr+wmi*wmi);
                double inc = 0.5*(atan2(wpi,wpr) + atan2(wmi,wmr))*180.0/M_PI;
                while(inc<0.0) inc += 180.0;
                while(inc>=180.0) inc -= 180.0;

                // P 162: H = A/f, G = g + V0 + u
                pha = std::fmod(pha + nvu[k], 360.0);
                if(pha<0.0)
                pha += 360.0;

                (*fc[5*k])[o] = sqrt(ae*ae+be*be)/nf[k];
                (*fc[5*k+1])[o] = pha;
                (*fc[5*k+2])[o] = (ap + am)/nf[k];
                (*fc[5*k+3])[o] = (ap - am)/nf[k];
                (*fc[5*k+4])[o] = inc;
            }
        }

        char fn[256];
        snprintf(fn,sizeof(fn),"./REEF3D_NHFLOW_Tidal/REEF3D-NHFLOW-Tidal-%08i-%04i-%04i.vts",count,p->mpirank+1,gi+1);
        ofstream out(fn);
        out<<"<?xml version=\"1.0\"?>\n";
        out<<"<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
        out<<"<StructuredGrid WholeExtent=\"0 "<<nx<<" 0 "<<ny<<" 0 0\">\n";
        out<<"<FieldData><DataArray type=\"Float64\" Name=\"TimeValue\" NumberOfTuples=\"1\">"<<setprecision(12)<<p->simtime<<"</DataArray>";
        out<<"<DataArray type=\"Float64\" Name=\"record_start\" NumberOfTuples=\"1\">"<<t0<<"</DataArray>";
        out<<"<DataArray type=\"Float64\" Name=\"record_length\" NumberOfTuples=\"1\">"<<Twin<<"</DataArray></FieldData>\n";
        out<<"<Piece Extent=\"0 "<<nx<<" 0 "<<ny<<" 0 0\">\n<CellData Scalars=\"eta_mean\">\n";
        for(size_t f=0; f<fname.size(); ++f)
        {
            out<<"<DataArray type=\"Float64\" Name=\""<<fname[f]<<"\" format=\"ascii\">\n";
            for(double v : fval[f])
            {
                if(std::isnan(v))
                out<<"nan\n";
                else
                out<<setprecision(8)<<v<<"\n";
            }
            out<<"</DataArray>\n";
        }
        out<<"</CellData>\n<Points>\n<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
        for(int jj=g.j0; jj<=g.j1; ++jj)
        for(int ii=g.i0; ii<=g.i1; ++ii)
        out<<setprecision(12)<<q->XN2D[q->nij(ii,jj)]<<" "<<q->YN2D[q->nij(ii,jj)]<<" 0\n";
        out<<"</DataArray>\n</Points>\n</Piece>\n</StructuredGrid>\n</VTKFile>\n";
        out.close();
    }

    if(p->mpirank==0)
    {
        char fn[256];
        snprintf(fn,sizeof(fn),"./REEF3D_NHFLOW_Tidal/REEF3D-NHFLOW-Tidal-%08i.vtm",count);
        ofstream out(fn);
        out<<"<?xml version=\"1.0\"?>\n<VTKFile type=\"vtkMultiBlockDataSet\" version=\"1.0\">\n<vtkMultiBlockDataSet>\n";
        int idx2=0;
        for(int r=0; r<p->mpi_size; ++r)
        for(int gi=0; gi<all[r]; ++gi)
        {
        out<<"<DataSet index=\""<<idx2<<"\" file=\"REEF3D-NHFLOW-Tidal-"<<setw(8)<<setfill('0')<<count<<"-"<<setw(4)<<r+1<<"-"<<setw(4)<<gi+1<<".vts\"/>\n";
        out<<setfill(' ');
        ++idx2;
        }
        out<<"</vtkMultiBlockDataSet>\n</VTKFile>\n";
        out.close();

        snprintf(fn,sizeof(fn),"./REEF3D_NHFLOW_Tidal/REEF3D-NHFLOW-Tidal-%08i.txt",count);
        ofstream txt(fn);
        txt<<"# NHFLOW long-run statistics (P 160)\n";
        txt<<"time "<<setprecision(12)<<p->simtime<<"\nrecord_start "<<t0<<"\nrecord_length_s "<<Twin<<"\nrecord_length_days "<<Twin/86400.0<<"\nrayleigh "<<rayleigh<<"\n";
        if(p->P162==1)
        {
        txt<<"# P 162: the run starts at "<<p->P162_y<<"-"<<setw(2)<<setfill('0')<<p->P162_mo<<"-"<<setw(2)<<p->P162_d<<" "
           <<setw(2)<<p->P162_h<<":"<<setw(2)<<p->P162_mi<<" UTC"<<setfill(' ')<<"; amplitudes H = A/f (nodal factor at the middle of the record),\n";
        txt<<"# phases Greenwich phase lags G = g + V0 + u (V0 at the start of the run, u at the middle of the record)\n";
        txt<<"# nodal f and V0+u (deg):";
        for(int k=0; k<nsel; ++k)
        txt<<" "<<name[sel[k]]<<" "<<setprecision(6)<<nf[k]<<" "<<nvu[k];
        txt<<"\n";
        txt<<"# constituent  speed_deg_per_h  period_h\n";
        }
        else
        txt<<"# constituent  speed_deg_per_h  period_h  (phases relative to t = 0 of the run, no nodal corrections)\n";
        for(int k=0; k<nsel; ++k)
        txt<<name[sel[k]]<<" "<<setprecision(10)<<cspeed[sel[k]]<<" "<<360.0/cspeed[sel[k]]<<"\n";
        txt.close();

        cout<<"NHFLOW P 160: maps of "<<setprecision(4)<<Twin/86400.0<<" days written, constituents:";
        for(int k=0; k<nsel; ++k)
        cout<<" "<<name[sel[k]];
        cout<<endl;
    }

    ++count;
}

// ---------------------------------------------------------------------------
// checkpoints (P 48, I 45)
// ---------------------------------------------------------------------------

// the sums of every grid component-major (a slice of its lexer per component), copied from the sums or zero
void nhflow_longrun::pack_all(bool copy)
{
    pack.assign(G.size(),vector<double>());
    for(size_t gi=0; gi<G.size(); ++gi)
    {
        const grid &g = G[gi];
        const size_t nc = (size_t)g.q->imax*g.q->jmax;
        pack[gi].assign(nc*NV,0.0);
        if(copy)
        for(size_t c=0; c<nc; ++c)
        for(int k=0; k<NV; ++k)
        pack[gi][(size_t)k*nc + c] = g.s[c*NV + k];
    }
}

void nhflow_longrun::unpack_all()
{
    for(size_t gi=0; gi<G.size() && gi<pack.size(); ++gi)
    {
        grid &g = G[gi];
        const size_t nc = (size_t)g.q->imax*g.q->jmax;
        for(size_t c=0; c<nc; ++c)
        for(int k=0; k<NV; ++k)
        g.s[c*NV + k] = pack[gi][(size_t)k*nc + c];
    }
    pack.clear();
}

void nhflow_longrun::ckpt_fields(lexer *q, const std::function<void(const char*, char, double*)> &v)
{
    for(size_t gi=0; gi<G.size() && gi<pack.size(); ++gi)
    if(G[gi].q==q)
    {
        const size_t nc = (size_t)q->imax*q->jmax;
        for(int k=0; k<NV; ++k)
        v(fname[k].c_str(),'S',&pack[gi][(size_t)k*nc]);
        return;
    }
}

// the record: version, started, t0, dt_out, next output, its length, the number of outputs, the basis, then
// the matrix of the normal equations (the same on all ranks); the sums of the grids when it has started
void nhflow_longrun::ckpt_write(lexer *p, reefckpt &ck, nhflow_amr *amr)
{
    vector<double> b = {1.0, ini ? 1.0 : 0.0, t0, dtout, tnext, Twin, double(count), double(NB), double(NV), rayleigh};
    b.insert(b.end(),M.begin(),M.end());
    ck.blob("nhflow_longrun",b);

    if(!ini)
    return;

    pack_all(true);
    ckpt_fields(p,[&](const char *n, char, double *V){ ck.put_S(n,V); });
    if(amr!=nullptr)
    amr->longrun_ckpt_write(ck,this);
    pack.clear();
}

bool nhflow_longrun::ckpt_read(lexer *p, fdm_nhf *d, reefckpt &ck, nhflow_amr *amr, string &msg)
{
    // a checkpoint without the statistics (written without P 160): the record starts at the restart
    if(!ck.has("nhflow_longrun"))
    {
        if(p->simtime>t0)
        {
            t0 = p->simtime;
            tnext = dtout>0.0 ? t0 + dtout : 1.0e30;
        }
        if(p->mpirank==0)
        cout<<"NHFLOW P 160: the checkpoint has no long-run statistics, the record starts at t = "<<t0<<" s"<<endl;
        return true;
    }

    const vector<double> b = ck.get_blob("nhflow_longrun");
    if(b.size()<10 || int(b[7])!=NB || int(b[8])!=NV || b.size()!=size_t(10+NB*NB))
    {
        msg = "the long-run statistics (P 160) of the checkpoint do not fit this version";
        return false;
    }

    // not started at the checkpoint: P 160 of this run
    if(b[1]<0.5)
    return true;

    // the record continues from its start; a new dt_out: the next output on its multiples from t0
    t0 = b[2];
    Twin = b[5];
    count = int(b[6]);
    M.assign(b.begin()+10,b.end());
    if(fabs(b[3]-dtout)<1.0e-9*MAX(1.0,fabs(dtout)))
    tnext = b[4];
    else
    tnext = dtout>0.0 ? t0 + (floor((p->simtime-t0)/dtout + 1.0e-9) + 1.0)*dtout : 1.0e30;

    // the grids (the patches of the blocks are there), then their sums
    if(!ini)
    {
        add_grid(p,d,0,p->knox,0,p->knoy);
        if(amr!=nullptr)
        amr->longrun_grids(this);
        ini = true;
    }

    pack_all(false);
    bool ok = true;
    ckpt_fields(p,[&](const char *n, char, double *V){ ok = ck.has(n) && ck.get_S(n,V) && ok; });
    if(ok && amr!=nullptr)
    ok = amr->longrun_ckpt_read(ck,this,msg);
    if(!ok)
    {
        if(msg.empty())
        msg = "the long-run statistics (P 160) of the checkpoint are incomplete";
        pack.clear();
        return false;
    }
    unpack_all();

    if(p->mpirank==0)
    cout<<"NHFLOW P 160: the record from t = "<<t0<<" s continues ("<<Twin/86400.0<<" days, "<<count<<" outputs so far)"<<endl;
    return true;
}

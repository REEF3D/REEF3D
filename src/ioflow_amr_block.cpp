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

#include"ioflow_amr_block.h"
#include"geo_mesh.h"
#include"wave_interface.h"
#include"lexer.h"
#include"fdm2D.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include"slice.h"
#include<fstream>
#include<sstream>

ioflow_amr_block::ioflow_amr_block(lexer *p, ghostcell *pgc, patchBC_interface *ppBC) : ioflow_v(p,pgc,ppBC)
{
}

ioflow_amr_block::~ioflow_amr_block()
{
    delete pwave;
}

// the open faces of the G 54 lines nin (-i face) and nout (+i face), -1 none
void ioflow_amr_block::open_faces(lexer *p, int nin, int nout)
{
    if(nin>=0)
    {
        in = 1;
        Q = p->G54_v[nin];
        if(geo_mesh::g54_kind(p->G54_t[nin])==3)
        {
            series_read(series_file(3,int(p->G54_v[nin])),sq_t,sq_v);
            Q = series_value(sq_t,sq_v,p->simtime);
        }
    }

    if(nout>=0)
    {
        out = 1;
        wsf = p->G54_v[nout];
        if(geo_mesh::g54_kind(p->G54_t[nout])==4)
        {
            series_read(series_file(4,int(p->G54_v[nout])),sw_t,sw_v);
            wsf = series_value(sw_t,sw_v,p->simtime);
        }
        p->F62 = wsf;
    }
}

// G 54 type 3: discharge-<id>.dat, type 4: waterlevel-<id>.dat
std::string ioflow_amr_block::series_file(int t, int id)
{
    return std::string(t==3 ? "discharge-" : "waterlevel-")+std::to_string(id)+".dat";
}

// lines "t value" ('#' comments, other lines without two numbers are skipped), times increasing; false: no
// file, no value or the times not increasing
bool ioflow_amr_block::series_read(const std::string &name, std::vector<double> &t, std::vector<double> &v)
{
    t.clear();
    v.clear();
    std::ifstream f(name.c_str());
    if(!f)
    return false;

    std::string line;
    while(std::getline(f,line))
    {
        const size_t c = line.find('#');
        if(c!=std::string::npos)
        line.erase(c);
        std::istringstream ss(line);
        double a,b;
        if(!(ss>>a>>b))
        continue;
        if(!t.empty() && a<=t.back())
        return false;
        t.push_back(a);
        v.push_back(b);
    }
    return !t.empty();
}

// linear in time, the first and the last value held beyond the ends
double ioflow_amr_block::series_value(const std::vector<double> &t, const std::vector<double> &v, double time)
{
    if(time<=t.front())
    return v.front();
    if(time>=t.back())
    return v.back();

    size_t a=0, b=t.size()-1;
    while(b-a>1)
    {
        const size_t m = (a+b)/2;
        if(t[m]<=time)
        a = m;
        else
        b = m;
    }
    const double w = (time-t[a])/(t[b]-t[a]);
    return v[a] + w*(v[b]-v[a]);
}

// the values of the time series at the time of the step (the start of the step)
void ioflow_amr_block::series_update(lexer *p)
{
    if(!sq_t.empty())
    Q = series_value(sq_t,sq_v,p->simtime);

    if(!sw_t.empty())
    wsf = p->F62 = series_value(sw_t,sw_v,p->simtime);
}

// water level at the +i face: F 64 ramps it from F 60 over the first steps (as ioflow_f with F 62)
double ioflow_amr_block::level(lexer *p) const
{
    if(p->F64>0)
    {
        double f=1.0;
        if(p->count<p->F64)
        f = 0.5*cos(PI + PI*double(p->count)/double(p->F64)) + 0.5;

        return f*wsf + (1.0-f)*p->F60;
    }

    return wsf;
}

// ------------------------------------------------------------------------------------------- SFLOW
// as ioflow_f::Qin2D / Qout2D (B 60 1) over the cells of the open faces of the block
void ioflow_amr_block::discharge2D(lexer *p, fdm2D *b, ghostcell *pgc)
{
    series_update(p);

    // the backgrounds of the edge zones at the time of the step, the waves (iowave::sflow_bg_update)
    if(ez[0]>=0 || ez[1]>=0)
    bgs.update(p,p->simtime);
    if(pwave!=nullptr)
    {
        wave_run l0;
        wave_in(p,l0);
        p->wavetime = p->simtime;
        pwave->wave_prestep(p,pgc);
        wave_get(p,wown);
        wave_out(p,l0);
    }

    if(in==0)
    return;

    Ai=0.0;
    for(n=0;n<p->gcslin_count;n++)
    {
    i=p->gcslin[n][0];
    j=p->gcslin[n][1];

        if(p->wet[IJ]==1 && p->geo_curv==1)
        {
        const int qf = (i-1-p->imin)*p->jmax + (j-p->jmin);
        Ai += p->LEN1[qf]*b->hp(i-1,j);
        continue;
        }

        if(p->wet[IJ]==1)
        Ai += p->DYN[JP]*b->hp(i-1,j);
    }

    Ui = Q/(Ai>1.0e-20?Ai:1.0e20);
    p->Ui = Ui;

}

void ioflow_amr_block::inflow2D(lexer *p, fdm2D *b, ghostcell *pgc, slice &P, slice &Q, slice &bed, slice &eta)
{
    // the edge zones of the faces (B 522)
    edges2D(p,b,pgc,P,Q);

    if(in==1)
    for(n=0;n<p->gcslin_count;n++)
    {
    i=p->gcslin[n][0];
    j=p->gcslin[n][1];

        // curvilinear grid: inflow along the normal of the inflow face
        if(p->wet[IJ]==1 && p->geo_curv==1)
        {
        const int qf = (i-1-p->imin)*p->jmax + (j-p->jmin);
        const double nx = p->SX1[qf]/p->LEN1[qf], ny = p->SY1[qf]/p->LEN1[qf];

        P(i-1,j)=P(i-2,j)=P(i-3,j)=Ui*nx;
        Q(i-1,j)=Q(i-2,j)=Q(i-3,j)=Ui*ny;
        continue;
        }

        P(i-1,j)=P(i-2,j)=P(i-3,j) = (p->wet[IJ]==1) ? Ui : 0.0;
        Q(i-1,j)=Q(i-2,j)=Q(i-3,j) = 0.0;
    }
}

// as ioflow_f::waterlevel2D with F 62 (F 50 2/3)
void ioflow_amr_block::waterlevel2D(lexer *p, fdm2D *b, ghostcell *pgc, slice &eta)
{
    // the edge zones of the faces (B 522)
    edges2D_eta(p,b,eta);

    if(out==0)
    return;

    const double wsfout = level(p);

    for(n=0;n<p->gcslout_count;n++)
    {
    i=p->gcslout[n][0];
    j=p->gcslout[n][1];

        if(p->wet[IJ]==1)
        {
        eta(i,j)   = wsfout-p->wd;
        eta(i+1,j) = wsfout-p->wd;
        eta(i+2,j) = wsfout-p->wd;
        eta(i+3,j) = wsfout-p->wd;

        b->hp(i,j)   = MAX(eta(i+1,j) + p->wd - b->bed(i,j),0.0);
        b->hp(i+1,j) = MAX(eta(i+1,j) + p->wd - b->bed(i,j),0.0);
        b->hp(i+2,j) = MAX(eta(i+2,j) + p->wd - b->bed(i,j),0.0);
        b->hp(i+3,j) = MAX(eta(i+3,j) + p->wd - b->bed(i,j),0.0);
        }
    }
}

// ------------------------------------------------------------------------------------------- NHFLOW
// as ioflow_f::Qin_nhf: wet area of the -i face and the discharge through it
void ioflow_amr_block::qin_nhflow(lexer *p, fdm_nhf *d, double &A, double &Qi)
{
    A=0.0;
    Qi=0.0;
    for(n=0;n<p->gcin_count;n++)
    if(p->gcin[n][3]>0)
    {
        i=p->gcin[n][0];
        j=p->gcin[n][1];
        k=p->gcin[n][2];

        if(p->geo_curv==1)
        {
        const int qf = (i-1-p->imin)*p->jmax + (j-p->jmin);
        A+=p->LEN1[qf]*p->DZN[KP]*d->WL(i-1,j);
        Qi+=p->DZN[KP]*d->WL(i-1,j)*(d->U[Im1JK]*p->SX1[qf] + d->V[Im1JK]*p->SY1[qf]);
        continue;
        }

        const double area=p->DYN[JP]*p->DZN[KP]*d->WL(i-1,j);
        A+=area;
        Qi+=area*d->U[Im1JK];
    }
}

// as ioflow_f::discharge_nhflow (B 60 1): Ui from the inflow area, Uo the mean outflow velocity (start4V outflow 2)
void ioflow_amr_block::discharge_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    series_update(p);

    // the backgrounds of the edge, relaxation and beach zones at the time of the step (iowave::nhflow_bg_update)
    if(ez[0]>=0 || ez[1]>=0 || rz[0]>=0 || rz[1]>=0)
    bgs.update(p,p->simtime);

    // the waves at the start of the step (iowave::wave_prestep), B 95: the targets of the step until the stages set them
    tstage = -1.0e20;
    if(pwave!=nullptr)
    {
        wave_run l0;
        wave_in(p,l0);
        p->wavetime = p->simtime;
        pwave->wave_prestep(p,pgc);
        wave_get(p,wown);
        wave_out(p,l0);
    }

    if(in==1)
    {
    double Qi;
    qin_nhflow(p,d,Ai,Qi);
    Ui = Q/(Ai>1.0e-20?Ai:1.0e20);
    p->Ui = Ui;
    p->Qi = Qi;

    }

    if(out==1)
    {
    double Ao=0.0;
    p->Qo=0.0;
    for(n=0;n<p->gcout_count;n++)
    if(p->gcout[n][3]>0)
    {
        i=p->gcout[n][0];
        j=p->gcout[n][1];
        k=p->gcout[n][2];

        if(p->wet[IJ]==1 && p->DF[IJK]>0 && p->geo_curv==1)
        {
        Ao+=p->LEN1[IJ]*p->DZN[KP]*d->WL(i,j);
        p->Qo+=p->DZN[KP]*d->WL(i,j)*(d->U[IJK]*p->SX1[IJ] + d->V[IJK]*p->SY1[IJ]);
        continue;
        }

        if(p->wet[IJ]==1 && p->DF[IJK]>0)
        {
        const double area=p->DYN[JP]*p->DZN[KP]*d->WL(i,j);
        Ao+=area;
        p->Qo+=area*d->U[IJK];
        }
    }
    p->Uo=p->Qo/(Ao>1.0e-20?Ao:1.0e20);
    }
}

// curvilinear grid: the inflow speed in U, UH of the ghost cells turned into the normal of the face (ioflow_f::curv_inflow_dir)
void ioflow_amr_block::curv_dir(lexer *p, double *U, double *V, double *UH, double *VH)
{
    const int qf = (i-1-p->imin)*p->jmax + (j-p->jmin);
    const double len = p->LEN1[qf];
    const double nx = len>0.0 ? p->SX1[qf]/len : 1.0;
    const double ny = len>0.0 ? p->SY1[qf]/len : 0.0;

    const int c[3] = {Im1JK,Im2JK,Im3JK};
    for(int q=0; q<3; ++q)
    {
        V[c[q]]  = U[c[q]]*ny;
        U[c[q]]  = U[c[q]]*nx;
        VH[c[q]] = UH[c[q]]*ny;
        UH[c[q]] = UH[c[q]]*nx;
    }
}

// the discharge through the -i face and the edge zones of the faces
void ioflow_amr_block::inflow_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc, double *U, double *V, double *W, double *UH, double *VH, double *WH, slice &WL)
{
    if(in==1)
    inflow_q(p,d,U,V,W,UH,VH,WH);

    edges_nhflow(p,d,pgc,U,V,W,UH,VH,WH);
}

// as ioflow_f::inflow_plain_nhflow (B 61 1) and inflow_log_nhflow (B 61 2, 4, 5)
void ioflow_amr_block::inflow_q(lexer *p, fdm_nhf *d, double *U, double *V, double *W, double *UH, double *VH, double *WH)
{
    if(p->B61==1)
    {
    for(n=0;n<p->gcin_count;n++)
    {
    i=p->gcin[n][0];
    j=p->gcin[n][1];
    k=p->gcin[n][2];

        const bool wet = p->wet[IJ]==1 && p->DF[IJK]>0;
        U[Im1JK]=U[Im2JK]=U[Im3JK] = wet ? Ui : 0.0;
        V[Im1JK]=V[Im2JK]=V[Im3JK] = 0.0;
        W[Im1JK]=W[Im2JK]=W[Im3JK] = 0.0;
        UH[Im1JK]=UH[Im2JK]=UH[Im3JK] = wet ? (d->eta(i,j)+d->depth(i,j))*Ui : 0.0;
        VH[Im1JK]=VH[Im2JK]=VH[Im3JK] = 0.0;
        WH[Im1JK]=WH[Im2JK]=WH[Im3JK] = 0.0;

        if(p->geo_curv==1)
        curv_dir(p,U,V,UH,VH);
    }
    return;
    }

    if(p->B61!=2 && p->B61!=4 && p->B61!=5)
    return;

    // log profile from the bed shear velocity of the mean speed over the deepest water column
    double hmax=-1.0e20;
    for(n=0;n<p->gcin_count;++n)
    {
        i=p->gcin[n][0];
        j=p->gcin[n][1];
        hmax=MAX(hmax,d->WL(i,j));
    }

    const double ks = (p->S10>0) ? p->S20*p->S21 : p->B50;
    const double H=hmax;
    const double shearvel = Ui/(2.5*log((11.0*H/ks)));

    for(n=0;n<p->gcin_count;n++)
    {
    i=p->gcin[n][0];
    j=p->gcin[n][1];
    k=p->gcin[n][2];

        const double zcoor = p->ZSP[IJK]-d->bed(i,j);

        if(p->wet[IJ]==1 && p->DF[IJK]>0)
        U[Im1JK]=U[Im2JK]=U[Im3JK] = shearvel*2.5*log(MAX(30.0*MIN(zcoor,hmax)/ks,1.0));

        if(p->wet[IJ]==0 || p->DF[IJK]<0)
        U[Im1JK]=U[Im2JK]=U[Im3JK] = 0.0;

        V[Im1JK]=V[Im2JK]=V[Im3JK] = 0.0;

        if(p->geo_curv==1)
        curv_dir(p,U,V,UH,VH);
    }

    // discharge corrected to Q
    for(int q=0; q<5; ++q)
    {
    double A,Qi;
    qin_nhflow(p,d,A,Qi);
    p->Qi=Qi;
    const double ratio = (fabs(Qi)<1.0e-20) ? 1.0 : Q/Qi;

        for(n=0;n<p->gcin_count;++n)
        {
        i=p->gcin[n][0];
        j=p->gcin[n][1];
        k=p->gcin[n][2];

        U[Im1JK]*=ratio;
        U[Im2JK]*=ratio;
        U[Im3JK]*=ratio;

        if(p->geo_curv==1)
        {
        V[Im1JK]*=ratio;
        V[Im2JK]*=ratio;
        V[Im3JK]*=ratio;
        }
        }
    }

    for(n=0;n<p->gcin_count;n++)
    {
    i=p->gcin[n][0];
    j=p->gcin[n][1];
    k=p->gcin[n][2];

        const bool wet = p->wet[IJ]==1 && p->DF[IJK]>0;

        if(p->geo_curv==0)
        V[Im1JK]=V[Im2JK]=V[Im3JK] = 0.0;

        W[Im1JK]=W[Im2JK]=W[Im3JK] = 0.0;

        UH[Im1JK] = wet ? U[Im1JK]*d->WL(i,j) : 0.0;
        UH[Im2JK] = wet ? U[Im2JK]*d->WL(i,j) : 0.0;
        UH[Im3JK] = wet ? U[Im3JK]*d->WL(i,j) : 0.0;

        VH[Im1JK]=VH[Im2JK]=VH[Im3JK] = 0.0;

        if(p->geo_curv==1 && wet)
        {
        VH[Im1JK] = V[Im1JK]*d->WL(i,j);
        VH[Im2JK] = V[Im2JK]*d->WL(i,j);
        VH[Im3JK] = V[Im3JK]*d->WL(i,j);
        }

        WH[Im1JK]=WH[Im2JK]=WH[Im3JK] = 0.0;
    }
}

void ioflow_amr_block::rkinflow_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc, double *U, double *V, double *W, double *UH, double *VH, double *WH, slice &WL)
{
    for(n=0;n<p->gcin_count;n++)
    {
    i=p->gcin[n][0];
    j=p->gcin[n][1];
    k=p->gcin[n][2];

        U[Im3JK]=U[Im2JK]=U[Im1JK]=d->U[Im1JK];
        V[Im3JK]=V[Im2JK]=V[Im1JK]=d->V[Im1JK];
        W[Im3JK]=W[Im2JK]=W[Im1JK]=d->W[Im1JK];

        UH[Im3JK]=UH[Im2JK]=UH[Im1JK]=d->UH[Im1JK];
        VH[Im3JK]=VH[Im2JK]=VH[Im1JK]=d->VH[Im1JK];
        WH[Im3JK]=WH[Im2JK]=WH[Im1JK]=d->WH[Im1JK];
    }

    edges_rk(p,d,U,V,W,UH,VH,WH);
}

void ioflow_amr_block::rkinflow_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc, double *F, double *G)
{
    for(n=0;n<p->gcin_count;n++)
    {
    i=p->gcin[n][0];
    j=p->gcin[n][1];
    k=p->gcin[n][2];

        F[Im3JK]=F[Im2JK]=F[Im1JK]=G[Im1JK];
    }
}

// as ioflow_f::fsfinflow_nhflow: the inflow ghost level from the mean surface of the wet inflow cells (F 50 2, 4),
// the outflow level (G 54 2) in the last cells and their ghost cells
void ioflow_amr_block::fsfinflow_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc, slice &WL)
{
    if(in==1 && (p->F50==2 || p->F50==4))
    {
    int cnt=0;
    double zval=0.0;
    for(n=0;n<p->gcslin_count;n++)
    {
        i=p->gcslin[n][0];
        j=p->gcslin[n][1];

        if(p->wet[IJ]==1)
        {
        zval+=d->eta(i,j);
        ++cnt;
        }
    }

    p->Hi = p->wd + (cnt>0 ? zval/double(cnt) : 0.0);

    for(n=0;n<p->gcslin_count;n++)
    {
        i=p->gcslin[n][0];
        j=p->gcslin[n][1];

        if(p->wet[IJ]==1)
        {
        WL(i-1,j) = p->Hi - d->bed(i,j);
        WL(i-2,j) = p->Hi - d->bed(i,j);
        WL(i-3,j) = p->Hi - d->bed(i,j);

        d->eta(i-1,j) = WL(i-1,j) - d->depth(i,j);
        d->eta(i-2,j) = WL(i-2,j) - d->depth(i,j);
        d->eta(i-3,j) = WL(i-3,j) - d->depth(i,j);
        }
    }
    }

    if(out==1)
    {
    const double wsfout = level(p);

    for(n=0;n<p->gcslout_count;n++)
    {
        i=p->gcslout[n][0];
        j=p->gcslout[n][1];

        if(p->wet[IJ]==1)
        {
        WL(i,j)   = wsfout - d->bed(i,j);
        WL(i+1,j) = wsfout - d->bed(i,j);
        WL(i+2,j) = wsfout - d->bed(i,j);
        WL(i+3,j) = wsfout - d->bed(i,j);

        d->eta(i,j)   = WL(i,j)   - d->depth(i,j);
        d->eta(i+1,j) = WL(i+1,j) - d->depth(i,j);
        d->eta(i+2,j) = WL(i+2,j) - d->depth(i,j);
        d->eta(i+3,j) = WL(i+3,j) - d->depth(i,j);
        }
    }
    }

    edges_wl(p,d,WL);
}

// ------------------------------------------------------------------------------------------- NHFLOW edge zones
// the edge zones of the -i (zm) and +i (zp) faces (B 520 methods 3-6 with B 522 b), as iowave::nhflow_open_edges
// on the native grid without waves (the zones and the backgrounds are checked by reefamr::blocks_open)
void ioflow_amr_block::open_edges(lexer *p, const bc_zone *zm, const bc_zone *zp, const background_state &bg)
{
    ezone.clear();
    ez[0] = ez[1] = -1;
    if(zm!=nullptr)
    {
        ez[0] = int(ezone.size());
        ezone.push_back(*zm);
    }
    if(zp!=nullptr)
    {
        ez[1] = int(ezone.size());
        ezone.push_back(*zp);
    }
    if(!ezone.empty())
    {
        bgs = bg;
        bgs.update(p,p->simtime);
    }
}

// depth average of a layer field in the column (i,j) (iowave::nhflow_col_ubar)
double ioflow_amr_block::col_ubar(lexer *p, double *F)
{
    const int kk = k;
    double s = 0.0;

    for(k=0; k<p->knoz; ++k)
    s += F[IJK]*p->DZN[KP];

    k = kk;
    return s;
}

// B 513: profile factor of layer kk (iowave::bg_prof)
double ioflow_amr_block::bg_prof(lexer *p, int b, int kk, double h)
{
    double zc=0.0, s=0.0, fk=1.0;

    for(int m=0; m<p->knoz; ++m)
    {
        const double dz = p->DZN[m+marge];
        const double f = bgs.shape(b,zc+0.5*dz,h);
        s += dz*f;

        if(m==kk)
        fk = f;

        zc += dz;
    }

    return s>0.0 ? fk/s : 1.0;
}

namespace
{
// face f (0: -i, 1: +i) of the edge cell (i,j): the face index of the face lists, the inward unit normal, the
// length and the physical centre (iowave curv_edge_normal, curv_edge_len, curv_edge_centre)
void edge_face(const lexer *p, int f, int i, int j, double &cx, double &cy, double &len, double &xf, double &yf)
{
    const int fi = (f==0) ? i-1 : i;
    const int q = (fi-p->imin)*p->jmax + (j-p->jmin);
    const double sg = (f==0) ? 1.0 : -1.0;
    len = p->LEN1[q];
    cx = len>0.0 ? sg*p->SX1[q]/len : sg;
    cy = len>0.0 ? sg*p->SY1[q]/len : 0.0;
    const int a = p->nij(fi+1,j), c = p->nij(fi+1,j+1);
    xf = 0.5*(p->XN2D[a] + p->XN2D[c]);
    yf = 0.5*(p->YN2D[a] + p->YN2D[c]);
}
}

// the ghost cells of the edges (iowave::nhflow_open_edges, x edges on the curvilinear grid): the characteristics
// of the background and the interior (Riemann), Flather, the level of the background (clamped level) or its
// normal current / the B 525 discharge (clamped discharge)
void ioflow_amr_block::edges_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc_w, double *U, double *V, double *W, double *UH, double *VH, double *WH)
{
    if(ez[0]<0 && ez[1]<0)
    return;

    if(edge_h.empty())
    edge_h.assign(size_t(p->imax)*size_t(p->jmax),0.0);

    const double g = fabs(p->W22);

    for(int f=0; f<2; ++f)
    {
        if(ez[f]<0)
        continue;

        const bc_zone &z = ezone[ez[f]];
        const int sc = (f==0) ? 1 : 4;
        const int di = (f==0) ? -1 : 1;
        const double nn = -double(di);
        const int b = z.bg>0 ? bgs.index(z.bg) : -1;

        // waves of the edge's sources (B 524, Riemann only, iowave::nhflow_open_edges): per column eta and the column
        // means of u, v at the centre of the edge cell, at the -i edge the non-hydrostatic pressure of the waves
        const bool waves = z.method==bc_method::riemann && !z.sources.empty() && wave_on(p,pgc_w);
        wave_run l0;
        if(waves)
        {
            wave_in(p,l0);
            pwave->select_sources(&z.sources);
            p->wavetime = p->simtime;
            if(ew_eta.empty())
            {
                ew_eta.assign(size_t(p->imax)*size_t(p->jmax),0.0);
                ew_u.assign(size_t(p->imax)*size_t(p->jmax),0.0);
                ew_v.assign(size_t(p->imax)*size_t(p->jmax),0.0);
            }
            for(int list=0; list<2; ++list)
            {
                const int cs = list==0 ? p->gcslin_count : p->gcslout_count;
                int **gs = list==0 ? p->gcslin : p->gcslout;
                for(n=0; n<cs; ++n)
                {
                    i = gs[n][0];
                    j = gs[n][1];
                    if(gs[n][3]!=sc)
                    continue;
                    double xg,yg;
                    wave_frame(p,p->XC2D[IJ],p->YC2D[IJ],xg,yg);
                    ew_eta[IJ] = pwave->wave_eta(p,pgc_w,xg,yg);
                    double su=0.0, sv=0.0;
                    for(k=0; k<p->knoz; ++k)
                    {
                        su += pwave->wave_u(p,pgc_w,xg,yg,p->ZSP[IJK]-p->phimean)*p->DZN[KP];
                        sv += pwave->wave_v(p,pgc_w,xg,yg,p->ZSP[IJK]-p->phimean)*p->DZN[KP];
                    }
                    ew_u[IJ] = su;
                    ew_v[IJ] = sv;
                    if(f==0)
                    edge_pressure(p,d,pgc_w,xg,yg);
                }
            }
        }

        // clamped discharge edge with B 525: normal velocity Q / (sum of h ds along the face), ramped
        double uq = 0.0;
        if(z.method==bc_method::clamp_q && z.has_Q)
        {
            double area = 0.0;
            for(int list=0; list<2; ++list)
            {
                const int cs = list==0 ? p->gcslin_count : p->gcslout_count;
                int **gs = list==0 ? p->gcslin : p->gcslout;
                for(n=0; n<cs; ++n)
                if(gs[n][3]==sc)
                {
                    i = gs[n][0];
                    j = gs[n][1];
                    if(p->wet[IJ]==1)
                    {
                        double cx,cy,len,xf,yf;
                        edge_face(p,f,i,j,cx,cy,len,xf,yf);
                        area += d->WL(i,j)*len;
                    }
                }
            }
            const double r = (z.Q_tramp>0.0 && p->simtime<z.Q_tramp) ? 0.5*(1.0-cos(PI*fmax(p->simtime,0.0)/z.Q_tramp)) : 1.0;
            uq = area>1.0e-12 ? r*z.Q/area : 0.0;
        }

        for(int list=0; list<2; ++list)
        {
            const int count = list==0 ? p->gcin_count : p->gcout_count;
            int **gc = list==0 ? p->gcin : p->gcout;

            for(n=0; n<count; ++n)
            {
                i = gc[n][0];
                j = gc[n][1];
                k = gc[n][2];

                if(gc[n][3]!=sc)
                continue;

                const double h0 = d->depth(i,j);
                const double hi = d->WL(i,j);

                if(p->wet[IJ]==0 || hi<=1.0e-6 || h0<=1.0e-6)
                continue;

                // background at the centre of the boundary face, the inward normal (cx, cy)
                double cx,cy,len,xf,yf;
                edge_face(p,f,i,j,cx,cy,len,xf,yf);
                const double eb = b>=0 ? bgs.eta(b,xf,yf) : 0.0;
                double ub=0.0, vb=0.0;
                if(b>=0)
                bgs.vel(b,h0,xf,yf,ub,vb);

                const double ui = nn*(col_ubar(p,d->U)*cx + col_ubar(p,d->V)*cy);

                // the waves of the layer and of the column (ramped), added to the background
                double uwk=0.0, vwk=0.0, wwk=0.0, ewc=0.0, uwc=0.0, vwc=0.0;
                if(waves)
                {
                    double xg,yg;
                    wave_frame(p,p->XC2D[IJ],p->YC2D[IJ],xg,yg);
                    const double zk = p->ZSP[IJK]-p->phimean;
                    const double rw = wramp(p);
                    uwk = rw*pwave->wave_u(p,pgc_w,xg,yg,zk);
                    vwk = rw*pwave->wave_v(p,pgc_w,xg,yg,zk);
                    wwk = rw*pwave->wave_w(p,pgc_w,xg,yg,zk);
                    ewc = rw*ew_eta[IJ];
                    uwc = rw*ew_u[IJ];
                    vwc = rw*ew_v[IJ];
                }

                const double unb = (ub + uwc)*cx + (vb + vwc)*cy;      // target, normal
                const double uni = nn*ui;              // interior, normal
                double ung, hg;

                if(z.method==bc_method::riemann)
                {
                    const double hb = fmax(h0+eb+ewc,1.0e-6);
                    const double Rin  = unb + 2.0*sqrt(g*hb);
                    const double Rout = uni - 2.0*sqrt(g*hi);

                    hg  = pow(Rin-Rout,2.0)/(16.0*g);
                    ung = 0.5*(Rin+Rout);
                }
                else
                if(z.method==bc_method::flather)
                {
                    hg = hi;
                    ung = (ub*cx + vb*cy) - sqrt(g/hi)*((hi-h0) - eb);
                }
                else
                if(z.method==bc_method::clamp_level)
                {
                    hg = fmax(h0+eb,1.0e-6);
                    ung = uni;
                }
                else
                {
                    hg = hi;
                    ung = z.has_Q ? uq : (ub*cx + vb*cy);
                }

                edge_h[IJ] = hg;

                // inflow: tangential velocity of the background, outflow: of the interior; the normal part
                // from the edge, with the vertical profile of the background current (B 513) at inflow
                // (the waves: the normal part their deviation from the column mean, the tangential part with the
                // background at inflow)
                const bool inflow = ung>0.0;
                const double tx = -cy, ty = cx;
                const double fp = (inflow && bgs.profiled(b)) ? bg_prof(p,b,k,fmax(h0+eb,1.0e-6)) : 1.0;
                const double dn = (uwk - uwc)*cx + (vwk - vwc)*cy;
                const double ut = inflow ? (ub*fp + uwk)*tx + (vb*fp + vwk)*ty : U[IJK]*tx + V[IJK]*ty;
                const double ug = (ung*fp + dn)*cx + ut*tx;
                const double vg = (ung*fp + dn)*cy + ut*ty;
                const double wg = inflow ? wwk : W[IJK];

                const int ii=i;
                for(int q=1; q<=3; ++q)
                {
                    i = ii + q*di;
                    U[IJK] = ug;
                    V[IJK] = vg;
                    W[IJK] = wg;
                    UH[IJK] = hg*ug;
                    VH[IJK] = hg*vg;
                    WH[IJK] = hg*wg;
                }
                i = ii;
            }
        }

        if(waves)
        {
            pwave->select_sources(nullptr);
            wave_out(p,l0);
        }
    }
}

// the ghost cells of the RK arrays from the step's ghost cells (iowave::nhflow_open_edges_rk)
void ioflow_amr_block::edges_rk(lexer *p, fdm_nhf *d, double *U, double *V, double *W, double *UH, double *VH, double *WH)
{
    for(int f=0; f<2; ++f)
    {
        if(ez[f]<0)
        continue;

        const int sc = (f==0) ? 1 : 4;
        const int di = (f==0) ? -1 : 1;

        for(int list=0; list<2; ++list)
        {
            const int count = list==0 ? p->gcin_count : p->gcout_count;
            int **gc = list==0 ? p->gcin : p->gcout;

            for(n=0; n<count; ++n)
            {
                i = gc[n][0];
                j = gc[n][1];
                k = gc[n][2];

                if(gc[n][3]!=sc)
                continue;

                const int ii=i;
                for(int q=1; q<=3; ++q)
                {
                    i = ii + q*di;
                    U[IJK]=d->U[IJK];
                    V[IJK]=d->V[IJK];
                    W[IJK]=d->W[IJK];
                    UH[IJK]=d->UH[IJK];
                    VH[IJK]=d->VH[IJK];
                    WH[IJK]=d->WH[IJK];
                }
                i = ii;
            }
        }
    }
}

// the ghost water level of the edges (iowave::nhflow_open_edges_wl): Riemann, clamped level the depth of the
// step, Flather, clamped discharge zero gradient
void ioflow_amr_block::edges_wl(lexer *p, fdm_nhf *d, slice &WL)
{
    for(int f=0; f<2; ++f)
    {
        if(ez[f]<0)
        continue;

        const bc_zone &z = ezone[ez[f]];
        const int sc = (f==0) ? 1 : 4;
        const int di = (f==0) ? -1 : 1;
        const bool own_h = z.method==bc_method::riemann || z.method==bc_method::clamp_level;

        for(int list=0; list<2; ++list)
        {
            const int count = list==0 ? p->gcslin_count : p->gcslout_count;
            int **gc = list==0 ? p->gcslin : p->gcslout;

            for(n=0; n<count; ++n)
            {
                i = gc[n][0];
                j = gc[n][1];

                if(gc[n][3]!=sc)
                continue;

                const double hg = (own_h && !edge_h.empty() && edge_h[IJ]>0.0) ? edge_h[IJ] : WL(i,j);

                for(int q=1; q<=3; ++q)
                {
                    WL(i+q*di,j) = hg;
                    d->eta(i+q*di,j) = hg - d->depth(i,j);
                }
            }
        }
    }
}

// ------------------------------------------------------------------------------------------- SFLOW edge zones
// the ghost cells of the edges of the -i and +i faces (iowave::sflow_open_edges on the native grid, here along the
// normal of the block face): the characteristics of the background (+ the waves of the B 524 sources, averaged
// over B 160 + 1 levels) and the interior (Riemann), Flather, the level of the background (clamped level) or its
// normal current / the B 525 discharge (clamped discharge); U, V, W, UA, VA and eta in the three ghost cells
void ioflow_amr_block::edges2D(lexer *p, fdm2D *b, ghostcell *pgc, slice &U, slice &V)
{
    if(ez[0]<0 && ez[1]<0)
    return;

    if(edge_h.empty())
    edge_h.assign(size_t(p->imax)*size_t(p->jmax),0.0);

    const double g = fabs(p->W22);

    for(int f=0; f<2; ++f)
    {
        if(ez[f]<0)
        continue;

        const bc_zone &z = ezone[ez[f]];
        const int sc = (f==0) ? 1 : 4;
        const int di = (f==0) ? -1 : 1;
        const int bi = z.bg>0 ? bgs.index(z.bg) : -1;
        const int **gsl[2] = {(const int**)p->gcslin,(const int**)p->gcslout};
        const int gsc[2] = {p->gcslin_count,p->gcslout_count};

        // clamped discharge edge with B 525: normal velocity Q / (sum of h ds along the face), ramped
        double uq = 0.0;
        if(z.method==bc_method::clamp_q && z.has_Q)
        {
            double area = 0.0;
            for(int list=0; list<2; ++list)
            for(n=0; n<gsc[list]; ++n)
            if(gsl[list][n][3]==sc)
            {
                i = gsl[list][n][0];
                j = gsl[list][n][1];
                if(p->wet[IJ]==1)
                {
                    double cx,cy,len,xf,yf;
                    edge_face(p,f,i,j,cx,cy,len,xf,yf);
                    area += b->WL(i,j)*len;
                }
            }
            const double r = (z.Q_tramp>0.0 && p->simtime<z.Q_tramp) ? 0.5*(1.0-cos(PI*fmax(p->simtime,0.0)/z.Q_tramp)) : 1.0;
            uq = area>1.0e-12 ? r*z.Q/area : 0.0;
        }

        const bool waves = z.method==bc_method::riemann && !z.sources.empty() && wave_on(p,pgc);
        wave_run l0;
        if(waves)
        {
            wave_in(p,l0);
            pwave->select_sources(&z.sources);
            p->wavetime = p->simtime;
        }

        for(int list=0; list<2; ++list)
        for(n=0; n<gsc[list]; ++n)
        {
            i = gsl[list][n][0];
            j = gsl[list][n][1];

            if(gsl[list][n][3]!=sc)
            continue;

            const double h0 = b->depth(i,j);
            const double hi = b->eta(i,j) + b->depth(i,j);

            if(p->wet[IJ]==0 || hi<=1.0e-6 || h0<=1.0e-6)
            continue;

            // background at the centre of the boundary face, the inward normal (cx, cy)
            double cx,cy,len,xf,yf;
            edge_face(p,f,i,j,cx,cy,len,xf,yf);
            const double eb = bi>=0 ? bgs.eta(bi,xf,yf) : 0.0;
            double ub=0.0, vb=0.0;
            if(bi>=0)
            bgs.vel(bi,h0,xf,yf,ub,vb);

            // waves (ramped): eta and the velocities averaged over B 160 + 1 levels
            double ewc=0.0, uwc=0.0, vwc=0.0, wwc=0.0;
            if(waves)
            {
                double xg,yg;
                wave_frame(p,p->XC2D[IJ],p->YC2D[IJ],xg,yg);
                const double rw = wramp(p);
                const double ew = pwave->wave_eta(p,pgc,xg,yg);
                const double dz = (eb + ew + p->wd - b->bed(i,j))/double(p->B160);

                double su=0.0, sv=0.0, sw=0.0;
                double zl = -p->wd;
                for(int qn=0; qn<=p->B160; ++qn)
                {
                    su += pwave->wave_u(p,pgc,xg,yg,zl);
                    sv += pwave->wave_v(p,pgc,xg,yg,zl);
                    sw += pwave->wave_w(p,pgc,xg,yg,zl);
                    zl += dz;
                }

                ewc = rw*ew;
                uwc = rw*su/double(p->B160+1);
                vwc = rw*sv/double(p->B160+1);
                wwc = rw*sw/double(p->B160+1);
            }

            const double unb = (ub + uwc)*cx + (vb + vwc)*cy;      // target, normal
            const double uni = U(i,j)*cx + V(i,j)*cy;              // interior, normal
            double ung, hg;

            if(z.method==bc_method::riemann)
            {
                const double hb = fmax(h0+eb+ewc,1.0e-6);
                const double Rin  = unb + 2.0*sqrt(g*hb);
                const double Rout = uni - 2.0*sqrt(g*hi);

                hg  = pow(Rin-Rout,2.0)/(16.0*g);
                ung = 0.5*(Rin+Rout);
            }
            else
            if(z.method==bc_method::flather)
            {
                hg = hi;
                ung = (ub*cx + vb*cy) - sqrt(g/hi)*((hi-h0) - eb);
            }
            else
            if(z.method==bc_method::clamp_level)
            {
                hg = fmax(h0+eb,1.0e-6);
                ung = uni;
            }
            else
            {
                hg = hi;
                ung = z.has_Q ? uq : (ub*cx + vb*cy);
            }

            edge_h[IJ] = hg;

            // inflow: tangential velocity of the background (+ waves), outflow: of the interior
            const bool inflow = ung>0.0;
            const double tx = -cy, ty = cx;
            const double ut = inflow ? (ub + uwc)*tx + (vb + vwc)*ty : U(i,j)*tx + V(i,j)*ty;
            const double ug = ung*cx + ut*tx;
            const double vg = ung*cy + ut*ty;
            const double wg = inflow ? wwc : b->W(i,j);
            const double eg = (z.method==bc_method::riemann || z.method==bc_method::clamp_level) ? hg - h0 : b->eta(i,j);

            for(int q=1; q<=3; ++q)
            {
                const int ii = i + q*di;
                U(ii,j) = ug;
                V(ii,j) = vg;
                b->W(ii,j) = wg;
                b->UA(ii,j) = ug;
                b->VA(ii,j) = vg;
                b->eta(ii,j) = eg;
            }
        }

        if(waves)
        {
            pwave->select_sources(nullptr);
            wave_out(p,l0);
        }
    }
}

// the ghost level of the edges after a water level update (iowave::sflow_open_edges_eta): Riemann, clamped level
// h_g of the step, Flather, clamped discharge zero gradient
void ioflow_amr_block::edges2D_eta(lexer *p, fdm2D *b, slice &eta)
{
    for(int f=0; f<2; ++f)
    {
        if(ez[f]<0)
        continue;

        const bc_zone &z = ezone[ez[f]];
        const int sc = (f==0) ? 1 : 4;
        const int di = (f==0) ? -1 : 1;
        const bool own_h = z.method==bc_method::riemann || z.method==bc_method::clamp_level;

        for(int list=0; list<2; ++list)
        {
            const int cs = list==0 ? p->gcslin_count : p->gcslout_count;
            int **gs = list==0 ? p->gcslin : p->gcslout;

            for(n=0; n<cs; ++n)
            {
                i = gs[n][0];
                j = gs[n][1];

                if(gs[n][3]!=sc)
                continue;

                const double eg = (own_h && !edge_h.empty() && edge_h[IJ]>0.0) ? edge_h[IJ] - b->depth(i,j) : eta(i,j);

                for(int q=1; q<=3; ++q)
                eta(i+q*di,j) = eg;
            }
        }
    }
}

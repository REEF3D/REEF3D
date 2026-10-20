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
}

// the open faces of the G 54 lines nin (-i face) and nout (+i face), -1 none
void ioflow_amr_block::open_faces(lexer *p, int nin, int nout)
{
    if(nin>=0)
    {
        in = 1;
        Q = p->G54_v[nin];
        if(p->G54_t[nin]==3)
        {
            series_read(series_file(3,int(p->G54_v[nin])),sq_t,sq_v);
            Q = series_value(sq_t,sq_v,p->simtime);
        }
    }

    if(nout>=0)
    {
        out = 1;
        wsf = p->G54_v[nout];
        if(p->G54_t[nout]==4)
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

// as ioflow_f::inflow_plain_nhflow (B 61 1) and inflow_log_nhflow (B 61 2, 4, 5)
void ioflow_amr_block::inflow_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc, double *U, double *V, double *W, double *UH, double *VH, double *WH, slice &WL)
{
    if(in==0)
    return;

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
}

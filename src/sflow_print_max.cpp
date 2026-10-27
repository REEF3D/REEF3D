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

#include"sflow_print_max.h"
#include"lexer.h"
#include"fdm2D.h"
#include"ghostcell.h"
#include"slice4.h"
#include"met_forcing.h"
#include"runlog.h"
#include<cmath>
#include<limits>
#include<iomanip>
#include<sys/stat.h>
#include<sys/types.h>

namespace
{
    // the state: index of the fields
    enum {ETA=0, TETA, DEPTH, SPEED, TSPEED, DUR, FIRST, CLASS, WIND, PA};

    // a value that does not exist (NaN in the output; not 'none': an enumerator of vtk3D)
    const double unset_value = -1.0e30;

    // P 115: reduced over the ranks (the rank of the point has the value)
    const double absent_value = -1.0e20;

    double out(double v)
    {
        return v<=0.5*unset_value ? std::numeric_limits<double>::quiet_NaN() : v;
    }
}

sflow_print_max::sflow_print_max(lexer *p, fdm2D *b, ghostcell *pgc)
{
    fname[ETA]="max_eta";            fname[TETA]="time_max_eta";
    fname[DEPTH]="max_depth";        fname[SPEED]="max_speed";
    fname[TSPEED]="time_max_speed";  fname[DUR]="flood_duration";
    fname[FIRST]="first_flood_time"; fname[CLASS]="inundation";
    fname[WIND]="max_wind10";        fname[PA]="min_air_pressure";

    // all fields, also without wind and pressure: the same checkpoint for every W 200
    for(int k=0; k<nfield; ++k)
    mx[k] = new slice4(p);

    if(p->mpirank==0)
    mkdir("./REEF3D_SFLOW_MAX",0777);

    // P 115: the cell of each point (a fluid cell of this rank)
    p->Iarray(ploc,MAX(p->P115,1));

    for(int n=0; n<p->P115; ++n)
    {
    ploc[n] = 0;
    i = p->posc_i(p->P115_x[n]);
    j = p->j_dir==1 ? p->posc_j(p->P115_y[n]) : 0;

        if(i>=0 && i<p->knox && j>=0 && j<p->knoy)
        if(p->flagslice4[IJ]>0)
        ploc[n] = 1;
    }

    if(p->mpirank==0)
    cout<<"SFLOW maxima (P 112): from t = "<<p->P113<<" s, flooding depth P 114 "<<p->P114<<" m, "<<p->P115<<" points (P 115)"<<endl;
}

sflow_print_max::~sflow_print_max()
{
    for(int k=0; k<nfield; ++k)
    delete mx[k];
}

bool sflow_print_max::wet(lexer *p, fdm2D *b)
{
    return p->wet[IJ]==1 && b->hp(i,j)>p->P114;
}

void sflow_print_max::start_window(lexer *p, fdm2D *b, ghostcell *pgc)
{
    SLICEBASELOOP
    {
    (*mx[ETA])(i,j) = unset_value;
    (*mx[TETA])(i,j) = unset_value;
    (*mx[DEPTH])(i,j) = 0.0;
    (*mx[SPEED])(i,j) = 0.0;
    (*mx[TSPEED])(i,j) = unset_value;
    (*mx[DUR])(i,j) = 0.0;
    (*mx[FIRST])(i,j) = unset_value;
    (*mx[CLASS])(i,j) = (p->flagslice4[IJ]>0 && wet(p,b)) ? 2.0 : 1.0;
    (*mx[WIND])(i,j) = 0.0;
    (*mx[PA])(i,j) = -unset_value;
    }

    started = true;
}

void sflow_print_max::update(lexer *p, fdm2D *b, ghostcell *pgc)
{
    // a continued run (I 45): the window may have started before (the fields of the checkpoint)
    if(!checked)
    {
        int s=0;
        SLICEBASELOOP
        if((*mx[CLASS])(i,j)>0.5)
        s=1;

        started = pgc->globalimax(s)>0;
        checked = true;
    }

    if(p->simtime<p->P113)
    return;

    if(!started)
    start_window(p,b,pgc);

    const double t = p->simtime;
    const double dt = p->dt_old;    // the step just made (the time step classes set dt_old = dt, then the next dt)
    const double dtw = MIN(dt, t - p->P113);    // its part after the start of the window

    met_forcing *m = b->met;
    has_wind = m!=nullptr && m->wind;
    has_pressure = m!=nullptr && m->pressure;

    SLICELOOP4
    {
        if(has_wind)
        (*mx[WIND])(i,j) = MAX((*mx[WIND])(i,j), sqrt(m->U10(i,j)*m->U10(i,j) + m->V10(i,j)*m->V10(i,j)));

        if(has_pressure)
        (*mx[PA])(i,j) = MIN((*mx[PA])(i,j), m->Pa(i,j));

        if(!wet(p,b))
        continue;

        const double e = b->eta(i,j);
        const double h = b->hp(i,j);
        const double u = sqrt(b->U(i,j)*b->U(i,j) + p->y_dir*b->V(i,j)*b->V(i,j));

        if(e>(*mx[ETA])(i,j))
        {
        (*mx[ETA])(i,j) = e;
        (*mx[TETA])(i,j) = t;
        }

        (*mx[DEPTH])(i,j) = MAX((*mx[DEPTH])(i,j), h);

        if(u>(*mx[SPEED])(i,j) || (*mx[TSPEED])(i,j)<=0.5*unset_value)
        {
        (*mx[SPEED])(i,j) = MAX(u,(*mx[SPEED])(i,j));
        (*mx[TSPEED])(i,j) = t;
        }

        (*mx[DUR])(i,j) += dtw;

        if((*mx[CLASS])(i,j)<1.5)
        {
        (*mx[CLASS])(i,j) = 3.0;
        (*mx[FIRST])(i,j) = t;
        }
    }
}

void sflow_print_max::cellvalues(lexer *p, fdm2D *b, int k, double *v)
{
    // the output value of field k in the cell (i,j); v[0] the value
    if(p->flagslice4[IJ]<0 || !started)
    {
    v[0] = std::numeric_limits<double>::quiet_NaN();
    return;
    }

    const double c = (*mx[CLASS])(i,j);

    if(k==CLASS)
    v[0] = c<1.5 ? 0.0 : (c<2.5 ? 1.0 : 2.0);

    else
    if(k==PA)
    v[0] = (*mx[PA])(i,j)>=-0.5*unset_value ? std::numeric_limits<double>::quiet_NaN() : (*mx[PA])(i,j);

    else
    if((k==DEPTH || k==SPEED) && c<1.5)
    v[0] = std::numeric_limits<double>::quiet_NaN();

    else
    v[0] = out((*mx[k])(i,j));
}

void sflow_print_max::write(lexer *p, fdm2D *b, ghostcell *pgc)
{
    met_forcing *m = b->met;
    has_wind = m!=nullptr && m->wind;
    has_pressure = m!=nullptr && m->pressure;
    nout = 8 + (has_wind?1:0) + (has_pressure?1:0);

    if(p->P112>0)
    {
    write_vtp(p,b,pgc);

        if(p->mpirank==0)
        write_pvtp(p);
    }

    if(p->P115>0)
    write_points(p,b,pgc);

    // the highest water level of the domain
    double emax=unset_value, xm=0.0, ym=0.0, tm=0.0;
    if(started)
    SLICELOOP4
    if((*mx[ETA])(i,j)>emax)
    {
    emax = (*mx[ETA])(i,j);
    xm = p->XP[IP];
    ym = p->YP[JP];
    tm = (*mx[TETA])(i,j);
    }

    const double gmax = pgc->globalmax(emax);
    if(emax<gmax)
    xm=ym=tm=-1.0e20;
    xm = pgc->globalmax(xm);
    ym = pgc->globalmax(ym);
    tm = pgc->globalmax(tm);

    if(p->mpirank==0 && gmax>0.5*unset_value)
    cout<<"SFLOW maxima: highest water level "<<setprecision(4)<<gmax<<" m at x = "<<setprecision(9)<<xm<<" y = "<<ym<<", t = "<<tm<<" s"<<endl;
}

void sflow_print_max::write_vtp(lexer *p, fdm2D *b, ghostcell *pgc)
{
    int offset[50], n=0, iin;
    float ffn;
    double v[1];

    // offsets: points, the fields (cell data), connectivity, offsets
    offset[n]=0;
    ++n;
    offset[n]=offset[n-1]+sizeof(float)*p->pointnum2D*3+sizeof(int);
    ++n;
    for(int k=0; k<nfield; ++k)
    if(k<WIND || (k==WIND && has_wind) || (k==PA && has_pressure))
    {
    offset[n]=offset[n-1]+sizeof(float)*p->polygon_sum+sizeof(int);
    ++n;
    }
    offset[n]=offset[n-1]+sizeof(int)*p->polygon_sum*3+sizeof(int);
    ++n;
    offset[n]=offset[n-1]+sizeof(int)*p->polygon_sum+sizeof(int);
    ++n;

    sprintf(name,"./REEF3D_SFLOW_MAX/REEF3D-SFLOW-MAX-%06i.vtp",p->mpirank+1);
    ofstream result;
    result.open(name, ios::binary);

    vtp3D::beginning(p,result,p->pointnum2D,0,0,0,p->polygon_sum);

    n=0;
    vtp3D::points(result,offset,n);

    result<<"<CellData>\n";
    for(int k=0; k<nfield; ++k)
    if(k<WIND || (k==WIND && has_wind) || (k==PA && has_pressure))
    {
    result<<"<DataArray type=\"Float32\" Name=\""<<fname[k]<<"\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
    }
    result<<"</CellData>\n";

    vtp3D::polys(result,offset,n);

    vtp3D::ending(result);

    // points at the bed
    iin=sizeof(float)*p->pointnum2D*3;
    result.write((char*)&iin, sizeof(int));
    TPSLICELOOP
    {
    ffn=(p->geo_curv==1) ? p->XN2D[p->nij(i+1,j+1)] : p->XN[IP1];
    result.write((char*)&ffn, sizeof(float));

    ffn=(p->geo_curv==1) ? p->YN2D[p->nij(i+1,j+1)] : p->YN[JP1];
    result.write((char*)&ffn, sizeof(float));

    ffn=b->bednode(i,j);
    result.write((char*)&ffn, sizeof(float));
    }

    // the fields: the cell value on both triangles of the cell (the order of the connectivity)
    for(int k=0; k<nfield; ++k)
    if(k<WIND || (k==WIND && has_wind) || (k==PA && has_pressure))
    {
    iin=sizeof(float)*p->polygon_sum;
    result.write((char*)&iin, sizeof(int));
        SLICEBASELOOP
        {
        cellvalues(p,b,k,v);
        ffn=float(v[0]);
        result.write((char*)&ffn, sizeof(float));
        result.write((char*)&ffn, sizeof(float));
        }
    }

    // connectivity (as sflow_vtp_fsf)
    iin=sizeof(int)*p->polygon_sum*3;
    result.write((char*)&iin, sizeof(int));
    SLICEBASELOOP
    {
    iin=int(b->nodeval(i-1,j-1))-1;
    result.write((char*)&iin, sizeof(int));
    iin=int(b->nodeval(i,j-1))-1;
    result.write((char*)&iin, sizeof(int));
    iin=int(b->nodeval(i,j))-1;
    result.write((char*)&iin, sizeof(int));

    iin=int(b->nodeval(i-1,j-1))-1;
    result.write((char*)&iin, sizeof(int));
    iin=int(b->nodeval(i,j))-1;
    result.write((char*)&iin, sizeof(int));
    iin=int(b->nodeval(i-1,j))-1;
    result.write((char*)&iin, sizeof(int));
    }

    iin=sizeof(int)*p->polygon_sum;
    result.write((char*)&iin, sizeof(int));
    for(int q=0; q<p->polygon_sum; ++q)
    {
    iin=(q+1)*3;
    result.write((char*)&iin, sizeof(int));
    }

    vtp3D::footer(result);

    result.close();
}

void sflow_print_max::write_pvtp(lexer *p)
{
    sprintf(name,"./REEF3D_SFLOW_MAX/REEF3D-SFLOW-MAX.pvtp");

    ofstream result;
    result.open(name);

    vtp3D::beginningParallel(p,result);

    vtp3D::pointsParallel(result);

    result<<"<PCellData>\n";
    for(int k=0; k<nfield; ++k)
    if(k<WIND || (k==WIND && has_wind) || (k==PA && has_pressure))
    result<<"<PDataArray type=\"Float32\" Name=\""<<fname[k]<<"\"/>\n";
    result<<"</PCellData>\n";

    char pname[200];
    for(int n=0; n<p->M10; ++n)
    {
    sprintf(pname,"REEF3D-SFLOW-MAX-%06i.vtp",n+1);
    result<<"<Piece Source=\""<<pname<<"\"/>\n";
    }

    vtp3D::endingParallel(result);

    result.close();

    // a file overwritten during the run: registered once, as a static file
    if(p->plog)
    p->plog->table_file(p,"max","maxima",name);
}

void sflow_print_max::write_points(lexer *p, fdm2D *b, ghostcell *pgc)
{
    ofstream result;

    if(p->mpirank==0)
    {
    result.open("./REEF3D_SFLOW_MAX/REEF3D-SFLOW-MAX-points.dat");
    result<<"# REEF3D::SFLOW maxima at the points P 115 (the cell holding the point), t = "<<setprecision(9)<<p->simtime
          <<" s, from t = "<<p->P113<<" s (P 113), flooding depth "<<p->P114<<" m (P 114); nan: no value"<<endl;
    result<<"# inundation: 0 never wet, 1 wet at P 113, 2 flooded later; eta relative to the still water level F 60"<<endl;
    result<<"n x y bed max_eta time_max_eta max_depth max_speed time_max_speed flood_duration first_flood_time inundation";
    if(has_wind)
    result<<" max_wind10";
    if(has_pressure)
    result<<" min_air_pressure";
    result<<endl;
    }

    double v[1];

    for(int n=0; n<p->P115; ++n)
    {
    double val[nfield+1];
    for(int k=0; k<=nfield; ++k)
    val[k] = absent_value;

        if(ploc[n]==1)
        {
        i = p->posc_i(p->P115_x[n]);
        j = p->j_dir==1 ? p->posc_j(p->P115_y[n]) : 0;

        val[0] = b->bed(i,j) - p->wd;
            for(int k=0; k<nfield; ++k)
            {
            cellvalues(p,b,k,v);
            // NaN as a marker below the absent_value value (reduced by the maximum)
            val[k+1] = std::isnan(v[0]) ? -0.5e20 : v[0];
            }
        }

        for(int k=0; k<=nfield; ++k)
        val[k] = pgc->globalmax(val[k]);

        if(p->mpirank==0)
        {
        result<<n+1<<" "<<setprecision(9)<<p->P115_x[n]<<" "<<p->P115_y[n];

            for(int k=0; k<=nfield; ++k)
            {
            if(k==WIND+1 && !has_wind)
            continue;
            if(k==PA+1 && !has_pressure)
            continue;

            if(val[k]<=0.5*absent_value)
            result<<" nan";      // outside the domain, a solid cell, or no value
            else
            result<<" "<<setprecision(7)<<val[k];
            }
        result<<endl;
        }
    }

    if(p->mpirank==0)
    result.close();

    if(p->plog)
    p->plog->table_file(p,"max_points","maxima_points","./REEF3D_SFLOW_MAX/REEF3D-SFLOW-MAX-points.dat");
}

void sflow_print_max::ckpt_fields(const std::function<void(const char*, char, double*)> &fd, const std::function<void(const char*, char, int*)> &fi)
{
    fd("max_eta",'S',mx[ETA]->V);          fd("max_teta",'S',mx[TETA]->V);
    fd("max_depth",'S',mx[DEPTH]->V);      fd("max_speed",'S',mx[SPEED]->V);
    fd("max_tspeed",'S',mx[TSPEED]->V);    fd("max_duration",'S',mx[DUR]->V);
    fd("max_first",'S',mx[FIRST]->V);      fd("max_class",'S',mx[CLASS]->V);
    fd("max_wind",'S',mx[WIND]->V);        fd("max_pa",'S',mx[PA]->V);
}

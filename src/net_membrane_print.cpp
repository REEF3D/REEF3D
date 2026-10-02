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

#include"net_membrane.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include<fstream>
#include<iomanip>
#include<cstdio>

void net_membrane::print_timeseries(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    double ein=0.0, ain=0.0, eout=0.0, aout=0.0, umax=0.0, vol=0.0;

    SLICELOOP4
    if(p->wet[IJ]==1)
    {
        const double A = p->DXN[IP]*p->DYN[JP];

        vol += d->WL(i,j)*A;

        if(inside_footprint(p->XP[IP],p->YP[JP],delta))
        {
            ein += d->eta(i,j)*A;
            ain += A;
        }

        if(outside_footprint(p->XP[IP],p->YP[JP],delta))
        {
            eout += d->eta(i,j)*A;
            aout += A;
        }
    }

    LOOP
    if(p->wet[IJ]==1 && p->DF[IJK]>0)
    umax = MAX(umax, sqrt(d->U[IJK]*d->U[IJK] + d->V[IJK]*d->V[IJK] + d->W[IJK]*d->W[IJK]));

    ein  = pgc->globalsum(ein);
    ain  = pgc->globalsum(ain);
    eout = pgc->globalsum(eout);
    aout = pgc->globalsum(aout);
    umax = pgc->globalmax(umax);
    vol  = pgc->globalsum(vol);

    if(p->mpirank==0)
    {
        const double etain  = ain>0.0  ? ein/ain  : 0.0;
        const double etaout = aout>0.0 ? eout/aout : 0.0;
        const double dhl = etain - etaout;
        
        double zm=0.0, zmin=1.0e20;
        for(int q : flo_)
        {
            zm += x_[q](2);
            zmin = MIN(zmin, x_[q](2));
        }
        zm = flo_.empty() ? 0.0 : zm/double(flo_.size());

        ofstream ts((outdir+"/REEF3D_NHFLOW_Membrane_"+to_string(nMem)+".dat").c_str(), ios::app);
        ts<<setprecision(10)<<p->simtime<<" "<<etain<<" "<<etaout<<" "<<dhl<<" "<<Qleak<<" "
          <<Fx<<" "<<Fy<<" "<<Fz<<" "<<Fzfloor<<" "<<-p->W1*fabs(p->W22)*dhl*Afloor<<" "<<urelmax<<" "<<umax<<" "<<vol<<" "
          <<Fb_(0)+Ffl_(0)<<" "<<Fb_(1)+Ffl_(1)<<" "<<Fb_(2)+Ffl_(2)<<" "<<zm<<" "<<zmin<<" "<<vmax_<<" "<<Tmax_;
        
        // strong coupling: iterations of the time step (all stages), largest relative residual at the end of a stage
        if(iterated())
        ts<<" "<<citstep_<<" "<<cres_;
        
        ts<<"\n";
    }
}

bool net_membrane::print_now(lexer *p)
{
    // membrane.dat 'print dt': own interval, 'print 0': off; default: NHFLOW print control (P 30 or P 20)
    if(prm.printdt==0.0)
    return false;
    
    bool flag = printcount==0;
    double dtp = -1.0;
    
    if(prm.printdt>0.0)
    dtp = prm.printdt;
    
    else if(p->P30>0.0)
    dtp = p->P30;
    
    if(dtp>0.0)
    {
        if(p->simtime>=printtime)
        flag=true;
    }
    else if(p->P20>0 && p->count%p->P20==0)
    flag=true;
    
    if(flag && dtp>0.0)
    while(printtime<=p->simtime)
    printtime += dtp;
    
    return flag;
}

void net_membrane::print_vtp(lexer *p)
{
    // binary vtp of the membrane: node velocity, displacement and load per area, panel pressure jump, load,
    // membrane tension and tag; a pvd collection lists the files with their times
    if(p->mpirank!=0)
    {
        ++printcount;
        return;
    }
    
    const int np = x_.size();
    const int nt = tri_.size();
    
    char name[400];
    snprintf(name,sizeof(name),"%s/REEF3D-NHFLOW-Membrane-%i-%06i.vtp",vtpdir.c_str(),nMem,printcount);
    
    ofstream result(name, ios::binary);
    
    // tension per triangle: mean E t strain of its edges (flexible), compression as 0
    vector<float> tens(nt,0.0f);
    if(prm.structure==2)
    for(int t=0; t<nt; ++t)
    {
        double T=0.0;
        for(int q=0; q<3; ++q)
        {
            const int e = tedge_[t][q];
            const double L = (x_[edge_[e][1]]-x_[edge_[e][0]]).norm();
            T += MAX(0.0, prm.EA*(L - L0_[e])/L0_[e]);
        }
        tens[t] = float(T/3.0);
    }
    
    int offset[20];
    int n=0;
    offset[n]=0;
    ++n;
    
    auto add = [&](int bytes) {offset[n] = offset[n-1] + bytes + sizeof(int); ++n;};
    
    add(sizeof(float)*np*3);    // points
    add(sizeof(float)*np*3);    // velocity
    add(sizeof(float)*np*3);    // displacement
    add(sizeof(float)*np*3);    // load per area
    add(sizeof(float)*nt);      // dp
    add(sizeof(float)*nt*3);    // force
    add(sizeof(float)*nt);      // tension
    add(sizeof(int)*nt);        // tag
    add(sizeof(int)*nt*3);      // connectivity
    add(sizeof(int)*nt);        // offsets
    
    vtp3D::beginning(p, result, np, 0, 0, 0, nt);
    
    n=0;
    vtp3D::points(result, offset, n);
    
    result<<"<PointData>\n";
    result<<"<DataArray type=\"Float32\" Name=\"velocity\" NumberOfComponents=\"3\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
    result<<"<DataArray type=\"Float32\" Name=\"displacement\" NumberOfComponents=\"3\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
    result<<"<DataArray type=\"Float32\" Name=\"load\" NumberOfComponents=\"3\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
    result<<"</PointData>\n";
    
    result<<"<CellData>\n";
    result<<"<DataArray type=\"Float32\" Name=\"dp\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
    result<<"<DataArray type=\"Float32\" Name=\"force\" NumberOfComponents=\"3\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
    result<<"<DataArray type=\"Float32\" Name=\"tension\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
    result<<"<DataArray type=\"Int32\" Name=\"tag\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
    result<<"</CellData>\n";
    
    vtp3D::polys(result, offset, n);
    vtp3D::ending(result);
    
    int iin;
    float ffn;
    
    auto vec3 = [&](const vector<Eigen::Vector3d> &v)
    {
        iin = sizeof(float)*np*3;
        result.write((char*)&iin, sizeof(int));
        for(int q=0; q<np; ++q)
        for(int c=0; c<3; ++c)
        {
            ffn = float(v[q](c));
            result.write((char*)&ffn, sizeof(float));
        }
    };
    
    vec3(x_);
    vec3(xdot_);
    
    {
        vector<Eigen::Vector3d> u(np), f(np);
        for(int q=0; q<np; ++q)
        {
            u[q] = x_[q] - x0_[q];
            f[q] = Eigen::Vector3d(nf_[3*q+0],nf_[3*q+1],nf_[3*q+2])/MAX(an_[q],1.0e-20);
        }
        vec3(u);
        vec3(f);
    }
    
    // dp: pressure jump across the panel, positive when the inside pressure is higher
    iin = sizeof(float)*nt;
    result.write((char*)&iin, sizeof(int));
    for(int t=0; t<nt; ++t)
    {
        ffn = float(Eigen::Vector3d(tf_[3*t+0],tf_[3*t+1],tf_[3*t+2]).dot(tn_[t])/ta_[t]);
        result.write((char*)&ffn, sizeof(float));
    }
    
    iin = sizeof(float)*nt*3;
    result.write((char*)&iin, sizeof(int));
    for(int t=0; t<nt; ++t)
    for(int c=0; c<3; ++c)
    {
        ffn = float(tf_[3*t+c]);
        result.write((char*)&ffn, sizeof(float));
    }
    
    iin = sizeof(float)*nt;
    result.write((char*)&iin, sizeof(int));
    result.write((char*)tens.data(), sizeof(float)*nt);
    
    iin = sizeof(int)*nt;
    result.write((char*)&iin, sizeof(int));
    for(int t=0; t<nt; ++t)
    {
        iin = ttag_[t];
        result.write((char*)&iin, sizeof(int));
    }
    
    iin = sizeof(int)*nt*3;
    result.write((char*)&iin, sizeof(int));
    for(int t=0; t<nt; ++t)
    for(int q=0; q<3; ++q)
    {
        iin = tri_[t][q];
        result.write((char*)&iin, sizeof(int));
    }
    
    iin = sizeof(int)*nt;
    result.write((char*)&iin, sizeof(int));
    for(int t=0; t<nt; ++t)
    {
        iin = 3*(t+1);
        result.write((char*)&iin, sizeof(int));
    }
    
    vtp3D::footer(result);
    result.close();
    
    // collection file for the time series in ParaView
    pvdtime_.push_back(p->simtime);
    
    ofstream pvd((vtpdir+"/REEF3D-NHFLOW-Membrane-"+to_string(nMem)+".pvd").c_str());
    pvd<<"<?xml version=\"1.0\"?>\n<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">\n<Collection>\n";
    for(size_t q=0; q<pvdtime_.size(); ++q)
    {
        char fn[100];
        snprintf(fn,sizeof(fn),"REEF3D-NHFLOW-Membrane-%i-%06i.vtp",nMem,(int)q);
        pvd<<"<DataSet timestep=\""<<setprecision(10)<<pvdtime_[q]<<"\" group=\"\" part=\"0\" file=\""<<fn<<"\"/>\n";
    }
    pvd<<"</Collection>\n</VTKFile>\n";
    
    ++printcount;
}

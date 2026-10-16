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

#include"nhflow_rans_io.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include"nhflow_wall.h"
#include<cstring>

nhflow_rans_io::nhflow_rans_io(lexer *p, fdm_nhf *d) : nhflow_strain(p,d),
									 ke_c_1e(1.44), ke_c_2e(1.92),ke_sigma_k(1.0),ke_sigma_e(1.3),
									 kw_alpha(5.0/9.0), kw_beta(3.0/40.0),kw_sigma_k(2.0),kw_sigma_w(2.0)
{
    p->Darray(KIN,p->imax*p->jmax*(p->kmax+2));
    p->Darray(EPS,p->imax*p->jmax*(p->kmax+2));
    p->Iarray(WALLF,p->imax*p->jmax*(p->kmax+2));
    p->Iarray(WETN,p->imax*p->jmax);
    
    for(int n=0; n<p->imax*p->jmax; ++n)
    WETN[n] = -1;   // set from p->wet in the first wetting_ini
}

nhflow_rans_io::~nhflow_rans_io()
{
}

// 2D output (bed mode 0, free surface mode 1): eddyv, kin, eps/omega of the bottom or top sigma layer at the
// slice vertices, 3 arrays of pointnum2D as offset_ParaView_2D and name_vtp declare (was: a first eddyv array of
// pointnum 3D values, the second eddyv array not written and the layer index left at 0/1). Not called at present.
void nhflow_rans_io::print_2D(lexer* p, fdm_nhf *d, ghostcell *pgc, ofstream &result, int mode)
{
    const int kk = (mode==1) ? p->knoz-1 : 0;
    double *field[3] = {d->EV, KIN, EPS};
    
    for(int q=0; q<3; ++q)
    {
    double *F = field[q];
    
    iin=4*(p->pointnum2D);
    result.write((char*)&iin, sizeof (int));
    
        TPSLICELOOP
        {
        k = kk;
        
        if(p->j_dir==0)
        {
        jj=j;
        j=0;
        ffn=float(0.5*(F[IJK]+F[Ip1JK]));
        j=jj;
        }
        
        if(p->j_dir==1)
        ffn=float(0.25*(F[IJK]+F[Ip1JK]+F[IJp1K]+F[Ip1Jp1K]));
        
        result.write((char*)&ffn, sizeof (float));
        }
    }
}

void nhflow_rans_io::print_3D(lexer* p, fdm_nhf *d, ghostcell *pgc, std::vector<char> &buffer, size_t &m)
{
    // eddyv
    iin=4*(p->pointnum);
    std::memcpy(&buffer[m],&iin,sizeof(int));
    m+=sizeof(int);

    TPLOOP
    {
        if(p->j_dir==0)
        {
            jj=j;
            j=0;
            ffn=float(0.25*(d->EV[IJK]+d->EV[Ip1JK]+d->EV[IJKp1]+d->EV[Ip1JKp1]));
            j=jj;
        }
        else if(p->j_dir==1)
            ffn=float(0.125*(d->EV[IJK]+d->EV[Ip1JK]+d->EV[IJp1K]+d->EV[Ip1Jp1K]
                        +d->EV[IJKp1]+d->EV[Ip1JKp1]+d->EV[IJp1Kp1]+d->EV[Ip1Jp1Kp1])); 
            
        std::memcpy(&buffer[m],&ffn,sizeof(float));
        m+=sizeof(float);
    }

    // kin
    iin=4*(p->pointnum);
    std::memcpy(&buffer[m],&iin,sizeof(int));
    m+=sizeof(int);

    TPLOOP
    {
        if(p->j_dir==0)
        {
            jj=j;
            j=0;
            ffn=float(0.25*(KIN[IJK]+KIN[Ip1JK]+KIN[IJKp1]+KIN[Ip1JKp1]));
            j=jj;
        }
        else if(p->j_dir==1)
            ffn=float(0.125*(KIN[IJK]+KIN[Ip1JK]+KIN[IJp1K]+KIN[Ip1Jp1K]
                        +KIN[IJKp1]+KIN[Ip1JKp1]+KIN[IJp1Kp1]+KIN[Ip1Jp1Kp1]));

        std::memcpy(&buffer[m],&ffn,sizeof(float));
        m+=sizeof(float);
    }

    // eps
    iin=4*(p->pointnum);
    std::memcpy(&buffer[m],&iin,sizeof(int));
    m+=sizeof(int);

    TPLOOP
    {
        if(p->j_dir==0)
        {
            jj=j;
            j=0;
            ffn=float(0.25*(EPS[IJK]+EPS[Ip1JK]+EPS[IJKp1]+EPS[Ip1JKp1]));
            j=jj;
        }
        else if(p->j_dir==1)
            ffn=float(0.125*(EPS[IJK]+EPS[Ip1JK]+EPS[IJp1K]+EPS[Ip1Jp1K]
                        +EPS[IJKp1]+EPS[Ip1JKp1]+EPS[IJp1Kp1]+EPS[Ip1Jp1Kp1]));

        std::memcpy(&buffer[m],&ffn,sizeof(float));
        m+=sizeof(float);
    }
}

double nhflow_rans_io::ccipol_kinval(lexer *p, ghostcell *pgc, double xp, double yp, double zp)
{
    double val=0.0;

    //val=p->ccipol4( kin, xp, yp, zp);

    return val;
}

double nhflow_rans_io::ccipol_epsval(lexer *p, ghostcell *pgc, double xp, double yp, double zp)
{
    double val=0.0;

    //val=p->ccipol4( eps, xp, yp, zp);

    return val;
}

double nhflow_rans_io::ccipol_a_kinval(lexer *p, ghostcell *pgc, double xp, double yp, double zp)
{
    double val=0.0;

    //val=p->ccipol4a( kin, xp, yp, zp);

    return val;
}

double nhflow_rans_io::ccipol_a_epsval(lexer *p, ghostcell *pgc, double xp, double yp, double zp)
{
    double val=0.0;

    //val=p->ccipol4a( eps, xp, yp, zp);

    return val;
}

double nhflow_rans_io::kinval(int ii, int jj, int kk)
{
    double val=0.0;

    //val=kin(ii,jj,kk);

    return val;
}

double nhflow_rans_io::epsval(int ii, int jj, int kk)
{
    double val=0.0;

    //val=eps(ii,jj,kk);

    return val;
}

void nhflow_rans_io::kinget(int ii, int jj, int kk,double val)
{
    i=ii;
    j=jj;
    k=kk;
    
    KIN[IJK]=val;
}

void nhflow_rans_io::epsget(int ii, int jj, int kk,double val)
{
    i=ii;
    j=jj;
    k=kk;
    
    EPS[IJK]=val;
}

void nhflow_rans_io::gcupdate(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    pgc->start4V(p,KIN,20);
    pgc->start4V(p,EPS,30);
}

void nhflow_rans_io::name_pvtp(lexer *p, fdm_nhf *d, ghostcell *pgc, ofstream &result)
{
    result<<"<PDataArray type=\"Float32\" Name=\"eddyv\"/>\n";
    
    result<<"<PDataArray type=\"Float32\" Name=\"kin\"/>\n";
	
	if(p->A560==1 || p->A560==21)
	result<<"<PDataArray type=\"Float32\" Name=\"epsilon\"/>\n";
	if(p->A560==2 || p->A560==22)
    result<<"<PDataArray type=\"Float32\" Name=\"omega\"/>\n";
}

void nhflow_rans_io::name_vtp(lexer *p, fdm_nhf *d, ghostcell *pgc, ofstream &result, int *offset, int &n)
{
    result<<"<DataArray type=\"Float32\" Name=\"eddyv\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
    result<<"<DataArray type=\"Float32\" Name=\"kin\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
	if(p->A560==1 || p->A560==21)
	result<<"<DataArray type=\"Float32\" Name=\"epsilon\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
	if(p->A560==2 || p->A560==22)
    result<<"<DataArray type=\"Float32\" Name=\"omega\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
}

void nhflow_rans_io::offset_ParaView_2D(lexer *p, int *offset, int &n)
{
    offset[n]=offset[n-1]+4*(p->pointnum2D)+4;
	++n;
    offset[n]=offset[n-1]+4*(p->pointnum2D)+4;
	++n;
    offset[n]=offset[n-1]+4*(p->pointnum2D)+4;
	++n;
}

void nhflow_rans_io::name_ParaView_parallel(lexer *p, ofstream &result)
{
    result<<"<PDataArray type=\"Float32\" Name=\"eddyv\"/>\n";
    
    result<<"<PDataArray type=\"Float32\" Name=\"kin\"/>\n";
	
	if(p->A560==1 || p->A560==21)
	result<<"<PDataArray type=\"Float32\" Name=\"epsilon\"/>\n";
	if(p->A560==2 || p->A560==22)
    result<<"<PDataArray type=\"Float32\" Name=\"omega\"/>\n";
}

void nhflow_rans_io::name_ParaView(lexer *p, stringstream &result, int *offset, int &n)
{
    result<<"<DataArray type=\"Float32\" Name=\"eddyv\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
    result<<"<DataArray type=\"Float32\" Name=\"kin\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
	if(p->A560==1 || p->A560==21)
	result<<"<DataArray type=\"Float32\" Name=\"epsilon\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
	if(p->A560==2 || p->A560==22)
    result<<"<DataArray type=\"Float32\" Name=\"omega\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
}

void nhflow_rans_io::offset_ParaView(lexer *p, int *offset, int &n)
{
    offset[n]=offset[n-1]+4*(p->pointnum)+4;
	++n;
	offset[n]=offset[n-1]+4*(p->pointnum)+4;
	++n;
    offset[n]=offset[n-1]+4*(p->pointnum)+4;
	++n;
}


// newly wetted columns (wet now, dry at the last turbulence step): k and eps/omega, and their old time levels
// KN/EN, start from the mean of the wet neighbour columns at the same sigma level (neighbours with k > 0 and
// eps/omega > 0, i.e. not wetted in this step themselves). Before, the column started from k = eps = 0, and k
// diffused in from the neighbours faster than eps/omega, which gave large nu_t = cmu k^2/eps in shallow cells.
// komega: eps holds omega; the mean of k and omega is used the same way.
void nhflow_rans_io::wetting_ini(lexer *p, fdm_nhf *d, double *KN, double *EN, bool komega)
{
    int nn;
    double ksum, esum, cnt;
    
    SLICELOOP4
    {
        if(WETN[IJ]==0 && p->wet[IJ]==1)
        {
            KLOOP
            if(p->DF[IJK]>0)
            {
            ksum=esum=cnt=0.0;
            
            const int nb[4] = {Im1JK, Ip1JK, IJm1K, IJp1K};
            const int nbs[4] = {Im1J, Ip1J, IJm1, IJp1};
            
                for(nn=0; nn<(p->j_dir==1?4:2); ++nn)
                if(p->wet[nbs[nn]]==1 && p->flag4[nb[nn]]>0 && p->DF[nb[nn]]>0 && KIN[nb[nn]]>0.0 && EPS[nb[nn]]>0.0)
                {
                ksum += KIN[nb[nn]];
                esum += EPS[nb[nn]];
                cnt  += 1.0;
                }
                
                if(cnt>0.0)
                {
                KIN[IJK] = KN[IJK] = ksum/cnt;
                EPS[IJK] = EN[IJK] = esum/cnt;
                }
            }
        }
    }
}

// length-scale limit (A 564 >= 1): the turbulence length scale l = cmu^0.75 k^1.5/eps does not exceed kappa h
// (h the local water depth): eps >= cmu^0.75 k^1.5/(kappa h), omega >= k^0.5/(cmu^0.25 kappa h), i.e.
// nu_t <= cmu^0.25 k^0.5 kappa h. In open-channel equilibrium l stays below 0.16 h, so the limit only acts
// where eps/omega -> 0+ with k > 0 (newly wetted and very shallow cells, still water with residual k).
// Also updates WETN for wetting_ini.
void nhflow_rans_io::length_limit(lexer *p, fdm_nhf *d, bool komega)
{
    const double kappa = 0.4;
    
    if(p->A564>=1)
    LOOP
    if(p->wet[IJ]==1 && p->DF[IJK]>0 && KIN[IJK]>0.0 && d->WL(i,j)>1.0e-6)
    {
        const double lmax = kappa*d->WL(i,j);
        
        if(!komega)
        EPS[IJK] = MAX(EPS[IJK], pow(p->cmu,0.75)*pow(KIN[IJK],1.5)/lmax);
        
        if(komega)
        EPS[IJK] = MAX(EPS[IJK], sqrt(KIN[IJK])/(pow(p->cmu,0.25)*lmax));
    }
    
    SLICELOOP4
    WETN[IJ] = p->wet[IJ];
}

// horizontal wall distance for the SST F2 (A 564 2): distance of the cell centre, along its sigma layer, to the nearest
// x or y wall face with a wall function, i.e. side walls with A 519 2 (solid cells and wall boundaries as in
// nhflow_wall_faces). Seeded with half a cell at the wall cells and spread by d = min(d, d_nb + dx_nb) over the wet
// fluid cells, with the MPI exchange in every sweep. Computed once at the first call, so it holds the walls of the
// initial geometry (moving 6DOF bodies and thin bodies are covered locally in sst_walldist). 1e20 without side walls.
// On curvilinear grids (G 50) the distance is taken along the computational spacing DXP/DYP.
void nhflow_rans_io::walldist_h(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    if(WDH!=nullptr)
    return;
    
    const double big = 1.0e20;
    const int size = p->imax*p->jmax*(p->kmax+2);
    int w[6];
    
    p->Darray(WDH,size);
    
    for(int n=0; n<size; ++n)
    WDH[n] = big;
    
    if(p->A519!=2)
    return;
    
    double seeded = 0.0;
    
    LOOP
    if(p->flag4[IJK]>0 && p->DF[IJK]>0)
    {
        nhflow_wall_faces(p,i,j,k,w);
        
        if(w[0] || w[1])
        WDH[IJK] = MIN(WDH[IJK], 0.5*p->DXN[IP]);
        
        if(w[2] || w[3])
        WDH[IJK] = MIN(WDH[IJK], 0.5*p->DYN[JP]);
        
        if(WDH[IJK]<big)
        seeded = 1.0;
    }
    
    if(pgc->globalmax(seeded)<0.5)
    return;
    
    const int maxit = p->gknox + p->gknoy + 10;
    
    for(int it=0; it<maxit; ++it)
    {
        pgc->start4V_par(p,WDH,24);
        
        double change = 0.0;
        
        LOOP
        if(p->flag4[IJK]>0 && p->DF[IJK]>0)
        {
            double dn = WDH[IJK];
            
            if(p->flag4[Im1JK]>0 && p->DF[Im1JK]>0)
            dn = MIN(dn, WDH[Im1JK] + p->DXP[IM1]);
            
            if(p->flag4[Ip1JK]>0 && p->DF[Ip1JK]>0)
            dn = MIN(dn, WDH[Ip1JK] + p->DXP[IP]);
            
            if(p->j_dir==1 && p->flag4[IJm1K]>0 && p->DF[IJm1K]>0)
            dn = MIN(dn, WDH[IJm1K] + p->DYP[JM1]);
            
            if(p->j_dir==1 && p->flag4[IJp1K]>0 && p->DF[IJp1K]>0)
            dn = MIN(dn, WDH[IJp1K] + p->DYP[JP]);
            
            if(dn<WDH[IJK])
            {
            change = MAX(change, WDH[IJK]<0.5*big ? WDH[IJK]-dn : 1.0);   // first value: 1
            WDH[IJK] = dn;
            }
        }
        
        if(pgc->globalmax(change)<=1.0e-12)
        break;
    }
    
    pgc->start4V_par(p,WDH,24);
}

// wall distance for F2: bed, solids (SOLID), side walls (WDH), top walls in the same column (A 519 2) and the
// walls of the cell itself with a wall function (nhflow_wall.h, incl. thin bodies X 330 'mobility sharp')
double nhflow_rans_io::sst_walldist(lexer *p, fdm_nhf *d)
{
    double y = p->ZSP[IJK] - d->bed(i,j);

    if(d->SOLID[IJK] > 0.0)
    y = MIN(y, d->SOLID[IJK]);
    
    if(WDH!=nullptr)
    y = MIN(y, WDH[IJK]);
    
    // top wall: first solid cell above in the column, distance to its lower face
    if(p->A519==2)
    {
        const int ksave = k;
        
        for(k=ksave; k<p->knoz-1; ++k)
        if(p->flag4[IJKp1]<0 || p->DF[IJKp1]<0)
        {
            const double ztop = p->ZN[KP1]*d->WL(i,j) + d->bed(i,j);
            k = ksave;
            y = MIN(y, ztop - p->ZSP[IJK]);
            break;
        }
        
        k = ksave;
    }
    
    double dist,ks,ut;
    
    if(nhflow_turb_wall(p,d,i,j,k,dist,ks,ut)==1)
    y = MIN(y, dist);

    return MAX(y, 1.0e-10);
}

double nhflow_rans_io::sst_F2(lexer *p, fdm_nhf *d, double kval, double wval)
{
    double y = sst_walldist(p,d);

    kval = MAX(kval, 0.0);
    wval = MAX(wval, 1.0e-20);

    double arg2 = MAX( 2.0*sqrt(kval)/(p->cmu*wval*y),
                       500.0*d->VISC[IJK]/(y*y*wval) );

    arg2 = MIN(arg2, 25.0);

    return tanh(arg2*arg2);
}

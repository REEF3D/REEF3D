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
#include<cstring>

nhflow_rans_io::nhflow_rans_io(lexer *p, fdm_nhf *d) : nhflow_strain(p,d),
									 ke_c_1e(1.44), ke_c_2e(1.92),ke_sigma_k(1.0),ke_sigma_e(1.3),
									 kw_alpha(5.0/9.0), kw_beta(3.0/40.0),kw_sigma_k(2.0),kw_sigma_w(2.0),
									 sst_alpha1(5.0/9.0), sst_alpha2(0.44), sst_beta1(3.0/40.0), sst_beta2(0.0828), 
									 sst_sigma_k1(0.85), sst_sigma_k2(1.0), sst_sigma_w1(0.5), sst_sigma_w2(0.856)
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

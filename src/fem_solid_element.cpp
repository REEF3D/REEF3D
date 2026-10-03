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

#include"fem_solid.h"
#include<cmath>
#include<algorithm>

// natural coordinates of the hex8 nodes (VTK order)
static const double xin[8] = {-1, 1, 1,-1,-1, 1, 1,-1};
static const double eta[8] = {-1,-1, 1, 1,-1,-1, 1, 1};
static const double zet[8] = {-1,-1,-1,-1, 1, 1, 1, 1};

void fem_solid::shape_derivatives()
{
    // all elements are identical voxels: reference derivatives are computed once

    // one-point (centre) integration
    for(int a=0; a<8; ++a)
    {
        dN0[a][0] = xin[a]/(4.0*hx);
        dN0[a][1] = eta[a]/(4.0*hy);
        dN0[a][2] = zet[a]/(4.0*hz);
    }

    // 2x2x2 Gauss points
    const double g = 1.0/std::sqrt(3.0);
    for(int q=0; q<8; ++q)
    {
        const double xg = g*xin[q], eg = g*eta[q], zg = g*zet[q];
        for(int a=0; a<8; ++a)
        {
            dNg[q][a][0] = 2.0/hx * xin[a]*(1.0+eta[a]*eg)*(1.0+zet[a]*zg)/8.0;
            dNg[q][a][1] = 2.0/hy * eta[a]*(1.0+xin[a]*xg)*(1.0+zet[a]*zg)/8.0;
            dNg[q][a][2] = 2.0/hz * zet[a]*(1.0+xin[a]*xg)*(1.0+eta[a]*eg)/8.0;
        }
    }

    // hourglass base vectors (for a rectangular box they are orthogonal to
    // linear fields, so the Flanagan-Belytschko gamma vectors equal them)
    for(int a=0; a<8; ++a)
    {
        gam[0][a] = xin[a]*eta[a];
        gam[1][a] = eta[a]*zet[a];
        gam[2][a] = zet[a]*xin[a];
        gam[3][a] = xin[a]*eta[a]*zet[a];
    }
}

double fem_solid::damage_exp(double kappa,double e0,double ef) const
{
    // exponential softening: sigma = ft exp(-(eps-e0)/(ef-e0)) beyond e0
    if(kappa<=e0)
    return 0.0;
    return 1.0 - e0/kappa*std::exp(-(kappa-e0)/(ef-e0));
}

void fem_solid::stress(const material& mt,gpstate& st,const Mat3& F,const Mat3& Fdot,double h,double w,Mat3& P,double& svm,bool& failed)
{
    const Mat3 I = Mat3::Identity();
    const double J = F.determinant();

    failed = false;

    if(J<erode_J || J>1.0/erode_J)
    {
        failed = true;
        P.setZero();
        svm = 0.0;
        return;
    }

    Mat3 E = 0.5*(F.transpose()*F - I);
    Mat3 S;

    if(mt.type==MAT_ELASTIC)
    {
        S = mt.lambda*E.trace()*I + 2.0*mt.mu*E;
    }
    else if(mt.type==MAT_J2)
    {
        Mat3 Ep;
        Ep << st.Ep[0], st.Ep[3], st.Ep[5],
              st.Ep[3], st.Ep[1], st.Ep[4],
              st.Ep[5], st.Ep[4], st.Ep[2];

        const Mat3 Ee = E - Ep;
        const Mat3 Str = mt.lambda*Ee.trace()*I + 2.0*mt.mu*Ee;
        const double pm = Str.trace()/3.0;
        Mat3 s = Str - pm*I;
        const double seq = std::sqrt(1.5*(s.array()*s.array()).sum());
        const double fy = mt.sigy + mt.H*st.ep;

        if(seq>fy && seq>0.0)
        {
            const double dp = (seq-fy)/(3.0*mt.mu + mt.H);
            const Mat3 nrm = 1.5*s/seq;
            Ep += dp*nrm;
            s *= 1.0 - 3.0*mt.mu*dp/seq;
            st.ep += dp;
            wdiss += w*(fy + 0.5*mt.H*dp)*dp;

            st.Ep[0]=Ep(0,0); st.Ep[1]=Ep(1,1); st.Ep[2]=Ep(2,2);
            st.Ep[3]=Ep(0,1); st.Ep[4]=Ep(1,2); st.Ep[5]=Ep(0,2);
        }

        S = s + pm*I;

        if(mt.epsfail>0.0 && st.ep>=mt.epsfail)
        failed = true;
    }
    else
    {
        // concrete: isotropic damage, Rankine tension and compression crushing
        // on the principal effective stresses, crack-band regularised
        const Mat3 S0 = mt.lambda*E.trace()*I + 2.0*mt.mu*E;

        Eigen::SelfAdjointEigenSolver<Mat3> es;
        es.computeDirect(S0,Eigen::EigenvaluesOnly);
        const Eigen::Vector3d sp = es.eigenvalues();   // ascending

        const double et = std::max(sp(2),0.0)/mt.E;
        const double ec = std::max(-sp(0),0.0)/mt.E;
        st.kt = std::max(st.kt,et);
        st.kc = std::max(st.kc,ec);

        double dt = 0.0, dc = 0.0;
        if(mt.ft>0.0)
        {
            const double e0 = mt.ft/mt.E;
            const double ef = 0.5*e0 + mt.Gf/(mt.ft*h);
            dt = damage_exp(st.kt,e0,ef);
        }
        if(mt.fc>0.0)
        {
            const double e0 = mt.fc/mt.E;
            const double ef = 0.5*e0 + mt.Gc/(mt.fc*h);
            dc = damage_exp(st.kc,e0,ef);
        }

        const double dnew = std::min(1.0, 1.0 - (1.0-dt)*(1.0-dc));
        if(dnew>st.d)
        {
            const double psi0 = 0.5*(S0.array()*E.array()).sum();
            wdiss += w*psi0*(dnew-st.d);
            st.d = dnew;
        }

        S = (1.0-st.d)*S0;

        if(st.d>=mt.derode)
        failed = true;
    }

    // first Piola-Kirchhoff stress
    P = F*S;

    // bulk viscosity (compression only), damps the ringing of impacts
    const Mat3 Finv = F.inverse();
    const double trD = (Fdot*Finv).trace();
    if(trD<0.0)
    {
        const double rho = mt.rho/J;
        const double q = rho*h*(bulkq1*mt.cp*(-trD) + bulkq2*h*trD*trD);
        P -= q*J*Finv.transpose();
    }

    // von Mises of the Cauchy stress
    const Mat3 sig = F*S*F.transpose()/J;
    const Mat3 dev = sig - sig.trace()/3.0*I;
    svm = std::sqrt(1.5*(dev.array()*dev.array()).sum());
}

void fem_solid::internal_forces(double dts)
{
    (void)dts;
    const double w = Vel/double(ngp);

    for(int e=0; e<nelem(); ++e)
    {
        element& el = elems[e];
        if(!el.alive)
        continue;

        const material& mt = mats[el.mat];

        double xa[8][3], va[8][3];
        for(int a=0; a<8; ++a)
        for(int i=0; i<3; ++i)
        {
            xa[a][i] = x[el.n[a]](i);
            va[a][i] = v[el.n[a]](i);
        }

        double fe[8][3] = {};
        double svm_sum = 0.0;
        int nfail = 0;
        double dmean = 0.0;
        double Jmin = 1.0e30;

        for(int g=0; g<ngp; ++g)
        {
            const double (*dN)[3] = (ngp==1) ? dN0 : dNg[g];

            Mat3 F = Mat3::Zero(), Fd = Mat3::Zero();
            for(int a=0; a<8; ++a)
            for(int i=0; i<3; ++i)
            for(int J=0; J<3; ++J)
            {
                F(i,J)  += xa[a][i]*dN[a][J];
                Fd(i,J) += va[a][i]*dN[a][J];
            }

            Mat3 P;
            double svm;
            bool failed;
            gpstate& st = gps[e*ngp+g];
            stress(mt,st,F,Fd,helem,w,P,svm,failed);

            if(failed) ++nfail;
            svm_sum += svm;
            dmean += st.d;
            Jmin = std::min(Jmin,F.determinant());

            for(int a=0; a<8; ++a)
            for(int i=0; i<3; ++i)
            fe[a][i] += w*(P(i,0)*dN[a][0] + P(i,1)*dN[a][1] + P(i,2)*dN[a][2]);
        }

        el.svm = svm_sum/double(ngp);
        el.J = Jmin;
        dmean /= double(ngp);

        // erosion: failure in at least half of the integration points
        if(2*nfail>=ngp && nfail>0)
        {
            el.alive = false;
            surf_dirty = true;
            for(int a=0; a<8; ++a)
            --nalive[el.n[a]];
            continue;
        }

        // hourglass control (one-point integration): stiffness on the
        // hourglass modes of the current positions, rotation invariant
        if(ngp==1)
        {
            double bb = 0.0;
            for(int a=0; a<8; ++a)
            bb += dN0[a][0]*dN0[a][0] + dN0[a][1]*dN0[a][1] + dN0[a][2]*dN0[a][2];

            const double k = hg_coef*(1.0-dmean)*(mt.lambda+2.0*mt.mu)*Vel*bb/8.0;

            for(int al=0; al<4; ++al)
            {
                double q[3] = {0.0,0.0,0.0};
                for(int a=0; a<8; ++a)
                for(int i=0; i<3; ++i)
                q[i] += gam[al][a]*xa[a][i];

                for(int a=0; a<8; ++a)
                for(int i=0; i<3; ++i)
                fe[a][i] += k*gam[al][a]*q[i];
            }
        }

        for(int a=0; a<8; ++a)
        for(int i=0; i<3; ++i)
        fint[el.n[a]](i) += fe[a][i];
    }
}

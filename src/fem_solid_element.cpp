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

#include"fem_solid.h"
#include<cmath>
#include<algorithm>

// natural coordinates of the hex8 nodes (VTK order)
static const double xin[8] = {-1, 1, 1,-1,-1, 1, 1,-1};
static const double eta[8] = {-1,-1, 1, 1,-1,-1, 1, 1};
static const double zet[8] = {-1,-1,-1,-1, 1, 1, 1, 1};

bool fem_solid::element_geometry(const Vec3* Xa,egeom& G) const
{
    // isoparametric hex8: Gauss-point gradients, volume, uniform gradient
    // (one-point integration), hourglass vectors orthogonal to linear fields
    const double g = 1.0/std::sqrt(3.0);
    bool ok = true;
    G.V = 0.0;
    for(int a=0; a<8; ++a)
    {
        G.mw[a] = 0.0;
        for(int i=0; i<3; ++i)
        G.dN0[a][i] = 0.0;
    }

    for(int q=0; q<8; ++q)
    {
        const double xg = g*xin[q], eg = g*eta[q], zg = g*zet[q];
        double dxi[8][3], N[8];
        for(int a=0; a<8; ++a)
        {
            dxi[a][0] = xin[a]*(1.0+eta[a]*eg)*(1.0+zet[a]*zg)/8.0;
            dxi[a][1] = eta[a]*(1.0+xin[a]*xg)*(1.0+zet[a]*zg)/8.0;
            dxi[a][2] = zet[a]*(1.0+xin[a]*xg)*(1.0+eta[a]*eg)/8.0;
            N[a] = (1.0+xin[a]*xg)*(1.0+eta[a]*eg)*(1.0+zet[a]*zg)/8.0;
        }

        Mat3 J = Mat3::Zero();
        for(int a=0; a<8; ++a)
        for(int i=0; i<3; ++i)
        for(int j=0; j<3; ++j)
        J(i,j) += Xa[a](i)*dxi[a][j];

        const double detJ = J.determinant();
        if(detJ<=0.0)
        ok = false;

        const Mat3 Jinv = J.inverse();
        for(int a=0; a<8; ++a)
        for(int i=0; i<3; ++i)
        G.dNg[q][a][i] = Jinv(0,i)*dxi[a][0] + Jinv(1,i)*dxi[a][1] + Jinv(2,i)*dxi[a][2];

        G.wg[q] = detJ;
        G.V += detJ;
        for(int a=0; a<8; ++a)
        {
            G.mw[a] += detJ*N[a];
            for(int i=0; i<3; ++i)
            G.dN0[a][i] += detJ*G.dNg[q][a][i];
        }
    }

    if(G.V<=0.0)
    return false;

    G.bb = 0.0;
    for(int a=0; a<8; ++a)
    for(int i=0; i<3; ++i)
    {
        G.dN0[a][i] /= G.V;
        G.bb += G.dN0[a][i]*G.dN0[a][i];
    }

    // hourglass base vectors, made orthogonal to the linear fields
    for(int a=0; a<8; ++a)
    {
        G.gam[0][a] = xin[a]*eta[a];
        G.gam[1][a] = eta[a]*zet[a];
        G.gam[2][a] = zet[a]*xin[a];
        G.gam[3][a] = xin[a]*eta[a]*zet[a];
    }
    for(int al=0; al<4; ++al)
    {
        double gx[3] = {0.0,0.0,0.0};
        for(int b=0; b<8; ++b)
        for(int i=0; i<3; ++i)
        gx[i] += G.gam[al][b]*Xa[b](i);

        double tmp[8];
        for(int a=0; a<8; ++a)
        tmp[a] = G.gam[al][a] - (gx[0]*G.dN0[a][0] + gx[1]*G.dN0[a][1] + gx[2]*G.dN0[a][2]);
        for(int a=0; a<8; ++a)
        G.gam[al][a] = tmp[a];
    }

    // characteristic length: volume / largest face area
    static const int lf[6][4] = {{0,4,7,3},{1,2,6,5},{0,1,5,4},{3,7,6,2},{0,3,2,1},{4,5,6,7}};
    double amax = 0.0;
    for(int f=0; f<6; ++f)
    {
        const Vec3 c = (Xa[lf[f][2]]-Xa[lf[f][0]]).cross(Xa[lf[f][3]]-Xa[lf[f][1]]);
        amax = std::max(amax,0.5*c.norm());
    }
    G.L = amax>0.0 ? G.V/amax : 0.0;
    G.h = std::cbrt(G.V);

    return ok;
}

void fem_solid::shape_derivatives()
{
    // regular voxel (shared by all elements that are not distorted)
    Vec3 Xv[8];
    for(int a=0; a<8; ++a)
    Xv[a] = Vec3(0.5*(1.0+xin[a])*hx, 0.5*(1.0+eta[a])*hy, 0.5*(1.0+zet[a])*hz);
    element_geometry(Xv,vgeo);

    // distorted elements: their own geometry
    geos.clear();
    for(element& el : elems)
    {
        el.geo = -1;
        Vec3 Xa[8];
        bool regular = true;
        const Vec3 base(ox+el.ix*hx, oy+el.iy*hy, oz+el.iz*hz);
        for(int a=0; a<8; ++a)
        {
            Xa[a] = X[el.n[a]];
            if((Xa[a]-base-Xv[a]).squaredNorm() > 1.0e-20*hx*hx)
            regular = false;
        }
        if(regular)
        continue;

        egeom G;
        element_geometry(Xa,G);
        el.geo = (int)geos.size();
        geos.push_back(G);
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

    for(int e=0; e<nelem(); ++e)
    {
        element& el = elems[e];
        if(!el.alive || el.rigid)
        continue;

        const material& mt = mats[el.mat];
        const egeom& G = geom(el);

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
            const double (*dN)[3] = (ngp==1) ? G.dN0 : G.dNg[g];
            const double w = (ngp==1) ? G.V : G.wg[g];

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
            stress(mt,st,F,Fd,G.h,w,P,svm,failed);

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
            const double k = hg_coef*(1.0-dmean)*(mt.lambda+2.0*mt.mu)*G.V*G.bb/8.0;

            for(int al=0; al<4; ++al)
            {
                double q[3] = {0.0,0.0,0.0};
                for(int a=0; a<8; ++a)
                for(int i=0; i<3; ++i)
                q[i] += G.gam[al][a]*xa[a][i];

                for(int a=0; a<8; ++a)
                for(int i=0; i<3; ++i)
                fe[a][i] += k*G.gam[al][a]*q[i];
            }
        }

        for(int a=0; a<8; ++a)
        for(int i=0; i<3; ++i)
        fint[el.n[a]](i) += fe[a][i];
    }
}

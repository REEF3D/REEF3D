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

#include"sflow_boussinesq.h"
#include"lexer.h"
#include"fdm2D.h"
#include"ghostcell.h"
#include"solver2D.h"
#include"sflow_curv.h"

#define HWL (WL(i,j)>p->A244?WL(i,j):1.0e20)

sflow_boussinesq::sflow_boussinesq(lexer *p, fdm2D *b, ghostcell *pgc) : zeta(-0.53), breakratio(0.8),
                                    A(p),B(p),Ax(p),Ay(p),Bx(p),By(p),
                                    U4(p),V4(p),U1p(p),V1p(p),Cx(p),Cy(p),
                                    T1(p),T2(p),etat(p),
                                    vy(p),hvy(p),ux(p),hux(p),
                                    mask(p),io(p),ua_n(p),va_n(p),f(p),
                                    Xp(p),Xn(p),Xs(p),Xr(p),Yp(p),Ye(p),Yw(p),Yr(p),
                                    Ra(p),Rb(p),ucart(p),vcart(p),uy(p),vx(p),huy(p),hvx(p)
{
    // (not for the patches of the mesh refinement, whose lexer has no domain origin)
    if(p->mpirank==0 && p->origin_i>-1000000)
    cout<<"SFLOW Boussinesq equations (FUNWAVE-TVD formulation)"<<endl;
}

sflow_boussinesq::~sflow_boussinesq()
{
    for(auto v : kv)
    delete v;
}

// ---------------------------------------------------------------------------
// helpers
// ---------------------------------------------------------------------------
double sflow_boussinesq::za(lexer *p, fdm2D *b, int ii, int jj)
{
    return zeta*b->depth(ii,jj) + (1.0+zeta)*b->eta(ii,jj);
}

double sflow_boussinesq::ddx(lexer *p, slice &g)
{
    if(p->geo_curv==1)
    return sflow_curv::ddx(p,g,i,j);
    
    return (g(i+1,j) - g(i-1,j))/(p->DXP[IP] + p->DXP[IM1]);
}

double sflow_boussinesq::ddy(lexer *p, slice &g)
{
    if(p->j_dir==0)
    return 0.0;
    
    if(p->geo_curv==1)
    return sflow_curv::ddy(p,g,i,j);
    
    return (g(i,j+1) - g(i,j-1))/(p->DYP[JP] + p->DYP[JM1]);
}

// ---------------------------------------------------------------------------
// Dispersive operators. The switch m enters as A -> m A, B -> m B, and the
// gradients of m A, m B are taken compactly from face values:
//   d/dx(m B)_i = (m_e B_e - m_w B_w)/dx,  B_e = (u_i+1 - u_i)/dx + (v_y)_e
// u2 and V1' then share the same symmetric operator, which keeps the
// semi-discrete system neutrally stable for any spatial variation of m.
// ---------------------------------------------------------------------------

// implicit (u_a in x) part of V1'_x at cell i
void sflow_boussinesq::coef_x(lexer *p, fdm2D *b, double &cE, double &cW, double &cC)
{
    const double z  = za(p,b,i,j);
    const double hC = b->depth(i,j);
    const double hE = b->depth(i+1,j);
    const double hW = b->depth(i-1,j);
    const double ep = 0.5*(b->eta(i,j) + b->eta(i+1,j));
    const double em = 0.5*(b->eta(i,j) + b->eta(i-1,j));
    const double c1p = 0.25*(b->eta(i,j)*b->eta(i,j) + b->eta(i+1,j)*b->eta(i+1,j));
    const double c1m = 0.25*(b->eta(i,j)*b->eta(i,j) + b->eta(i-1,j)*b->eta(i-1,j));
    const double mp = 0.5*(mask(i,j) + mask(i+1,j));
    const double mm = 0.5*(mask(i,j) + mask(i-1,j));
    const double cp = mp/(p->DXP[IP]*p->DXN[IP]);
    const double cm = mm/(p->DXP[IM1]*p->DXN[IP]);
    
    cE =  cp*(0.5*z*z + z*hE - c1p - ep*hE);
    cW =  cm*(0.5*z*z + z*hW - c1m - em*hW);
    cC = -cp*(0.5*z*z + z*hC - c1p - ep*hC)
         -cm*(0.5*z*z + z*hC - c1m - em*hC);
}

void sflow_boussinesq::coef_y(lexer *p, fdm2D *b, double &cN, double &cS, double &cC)
{
    const double z  = za(p,b,i,j);
    const double hC = b->depth(i,j);
    const double hN = b->depth(i,j+1);
    const double hS = b->depth(i,j-1);
    const double ep = 0.5*(b->eta(i,j) + b->eta(i,j+1));
    const double em = 0.5*(b->eta(i,j) + b->eta(i,j-1));
    const double c1p = 0.25*(b->eta(i,j)*b->eta(i,j) + b->eta(i,j+1)*b->eta(i,j+1));
    const double c1m = 0.25*(b->eta(i,j)*b->eta(i,j) + b->eta(i,j-1)*b->eta(i,j-1));
    const double mp = 0.5*(mask(i,j) + mask(i,j+1));
    const double mm = 0.5*(mask(i,j) + mask(i,j-1));
    const double cp = mp/(p->DYP[JP]*p->DYN[JP]);
    const double cm = mm/(p->DYP[JM1]*p->DYN[JP]);
    
    cN =  cp*(0.5*z*z + z*hN - c1p - ep*hN);
    cS =  cm*(0.5*z*z + z*hS - c1m - em*hS);
    cC = -cp*(0.5*z*z + z*hC - c1p - ep*hC)
         -cm*(0.5*z*z + z*hC - c1m - em*hC);
}

double sflow_boussinesq::lx(lexer *p, fdm2D *b, slice &u)
{
    double cE,cW,cC;
    
    coef_x(p,b,cE,cW,cC);
    
    return cE*u(i+1,j) + cW*u(i-1,j) + cC*u(i,j);
}

double sflow_boussinesq::ly(lexer *p, fdm2D *b, slice &v)
{
    double cN,cS,cC;
    
    if(p->j_dir==0)
    return 0.0;
    
    coef_y(p,b,cN,cS,cC);
    
    return cN*v(i,j+1) + cS*v(i,j-1) + cC*v(i,j);
}

// cell-centred A, B and the tangential derivatives needed at the faces
void sflow_boussinesq::cellterms(lexer *p, fdm2D *b, ghostcell *pgc)
{
    // curvilinear grid: full gradients of u_a, v_a, h u_a, h v_a (the faces need the tangential parts)
    if(p->geo_curv==1)
    {
        SLICELOOP4
        {
        sflow_curv::grad(p,b->UA,i,j,ux(i,j),uy(i,j));
        sflow_curv::grad(p,b->VA,i,j,vx(i,j),vy(i,j));
        sflow_curv::grad5(p,i,j,b->depth(i+1,j)*b->UA(i+1,j),b->depth(i-1,j)*b->UA(i-1,j),
                                b->depth(i,j+1)*b->UA(i,j+1),b->depth(i,j-1)*b->UA(i,j-1),hux(i,j),huy(i,j));
        sflow_curv::grad5(p,i,j,b->depth(i+1,j)*b->VA(i+1,j),b->depth(i-1,j)*b->VA(i-1,j),
                                b->depth(i,j+1)*b->VA(i,j+1),b->depth(i,j-1)*b->VA(i,j-1),hvx(i,j),hvy(i,j));
        
        B(i,j) = mask(i,j)*(ux(i,j) + vy(i,j));
        A(i,j) = mask(i,j)*(hux(i,j) + hvy(i,j));
        
            if(p->wet[IJ]==0)
            {
            A(i,j) = 0.0;
            B(i,j) = 0.0;
            }
        }
        
        slice4 *g[10] = {&ux,&uy,&vx,&vy,&hux,&huy,&hvx,&hvy,&A,&B};
        for(int k=0; k<10; ++k)
        pgc->gcsl_start4(p,*g[k],1);
        return;
    }
    
    SLICELOOP4
    {
    ux(i,j)  = ddx(p,b->UA);
    hux(i,j) = (b->depth(i+1,j)*b->UA(i+1,j) - b->depth(i-1,j)*b->UA(i-1,j))/(p->DXP[IP] + p->DXP[IM1]);
    vy(i,j)  = 0.0;
    hvy(i,j) = 0.0;
    
        if(p->j_dir==1)
        {
        vy(i,j)  = ddy(p,b->VA);
        hvy(i,j) = (b->depth(i,j+1)*b->VA(i,j+1) - b->depth(i,j-1)*b->VA(i,j-1))/(p->DYP[JP] + p->DYP[JM1]);
        }
        
    B(i,j) = mask(i,j)*(ux(i,j) + vy(i,j));
    A(i,j) = mask(i,j)*(hux(i,j) + hvy(i,j));
    
        if(p->wet[IJ]==0)
        {
        A(i,j) = 0.0;
        B(i,j) = 0.0;
        }
    }
    
    pgc->gcsl_start4(p,ux,1);
    pgc->gcsl_start4(p,hux,1);
    pgc->gcsl_start4(p,vy,1);
    pgc->gcsl_start4(p,hvy,1);
    pgc->gcsl_start4(p,A,1);
    pgc->gcsl_start4(p,B,1);
}

// gradients of mA, mB (Ax,Ay,Bx,By), u2 = (U4,V4) and V1' = (U1p,V1p)
// requires the ghost cells of UA,VA
void sflow_boussinesq::operators(lexer *p, fdm2D *b, ghostcell *pgc)
{
    double z,h,e,cA,cB;
    double BE,BW,AE,AW,mE,mW,eE,eW,qE,qW,FE,FW;
    
    if(p->geo_curv==1)
    {
    operators_curv(p,b,pgc);
    return;
    }
    
    cellterms(p,b,pgc);
    
    SLICELOOP4
    {
    Ax(i,j)  = Bx(i,j)  = Ay(i,j)  = By(i,j) = 0.0;
    U4(i,j)  = V4(i,j)  = 0.0;
    U1p(i,j) = V1p(i,j) = 0.0;
    
        if(p->wet[IJ]==1)
        {
        z = za(p,b,i,j);
        h = b->depth(i,j);
        e = b->eta(i,j);
        
        cB = 0.5*z*z - (h*h - h*e + e*e)/6.0;
        cA = z + 0.5*(h - e);
        
        // x
        mE = 0.5*(mask(i,j) + mask(i+1,j));
        mW = 0.5*(mask(i,j) + mask(i-1,j));
        eE = 0.5*(b->eta(i,j) + b->eta(i+1,j));
        eW = 0.5*(b->eta(i,j) + b->eta(i-1,j));
        qE = 0.25*(b->eta(i,j)*b->eta(i,j) + b->eta(i+1,j)*b->eta(i+1,j));
        qW = 0.25*(b->eta(i,j)*b->eta(i,j) + b->eta(i-1,j)*b->eta(i-1,j));
        
        BE = (b->UA(i+1,j) - b->UA(i,j))/p->DXP[IP]  + 0.5*(vy(i,j) + vy(i+1,j));
        BW = (b->UA(i,j) - b->UA(i-1,j))/p->DXP[IM1] + 0.5*(vy(i,j) + vy(i-1,j));
        AE = (b->depth(i+1,j)*b->UA(i+1,j) - b->depth(i,j)*b->UA(i,j))/p->DXP[IP]  + 0.5*(hvy(i,j) + hvy(i+1,j));
        AW = (b->depth(i,j)*b->UA(i,j) - b->depth(i-1,j)*b->UA(i-1,j))/p->DXP[IM1] + 0.5*(hvy(i,j) + hvy(i-1,j));
        
        Bx(i,j) = (mE*BE - mW*BW)/p->DXN[IP];
        Ax(i,j) = (mE*AE - mW*AW)/p->DXN[IP];
        
        FE = mE*(qE*BE + eE*AE);
        FW = mW*(qW*BW + eW*AW);
        
        U4(i,j)  = cB*Bx(i,j) + cA*Ax(i,j);
        U1p(i,j) = 0.5*z*z*Bx(i,j) + z*Ax(i,j) - (FE - FW)/p->DXN[IP];
        
            // y
            if(p->j_dir==1)
            {
            mE = 0.5*(mask(i,j) + mask(i,j+1));
            mW = 0.5*(mask(i,j) + mask(i,j-1));
            eE = 0.5*(b->eta(i,j) + b->eta(i,j+1));
            eW = 0.5*(b->eta(i,j) + b->eta(i,j-1));
            qE = 0.25*(b->eta(i,j)*b->eta(i,j) + b->eta(i,j+1)*b->eta(i,j+1));
            qW = 0.25*(b->eta(i,j)*b->eta(i,j) + b->eta(i,j-1)*b->eta(i,j-1));
            
            BE = (b->VA(i,j+1) - b->VA(i,j))/p->DYP[JP]  + 0.5*(ux(i,j) + ux(i,j+1));
            BW = (b->VA(i,j) - b->VA(i,j-1))/p->DYP[JM1] + 0.5*(ux(i,j) + ux(i,j-1));
            AE = (b->depth(i,j+1)*b->VA(i,j+1) - b->depth(i,j)*b->VA(i,j))/p->DYP[JP]  + 0.5*(hux(i,j) + hux(i,j+1));
            AW = (b->depth(i,j)*b->VA(i,j) - b->depth(i,j-1)*b->VA(i,j-1))/p->DYP[JM1] + 0.5*(hux(i,j) + hux(i,j-1));
            
            By(i,j) = (mE*BE - mW*BW)/p->DYN[JP];
            Ay(i,j) = (mE*AE - mW*AW)/p->DYN[JP];
            
            FE = mE*(qE*BE + eE*AE);
            FW = mW*(qW*BW + eW*AW);
            
            V4(i,j)  = cB*By(i,j) + cA*Ay(i,j);
            V1p(i,j) = 0.5*z*z*By(i,j) + z*Ay(i,j) - (FE - FW)/p->DYN[JP];
            }
        }
    }
    
    pgc->gcsl_start4(p,U4,1);
    pgc->gcsl_start4(p,V4,1);
}

// ---------------------------------------------------------------------------
// dispersion switch: breaking (A 246), |eta|/h >= 0.8, Fr >= 1, shoreline, in- and outflow
// ---------------------------------------------------------------------------
void sflow_boussinesq::mask_update(lexer *p, fdm2D *b, ghostcell *pgc, slice &WL)
{
    // cells next to in- and outflow boundaries: 1 inflow (west), 2 outflow (east)
    SLICELOOP4
    io(i,j) = 0.0;
    
    for(n=0;n<p->gcslin_count;n++)
    {
    i=p->gcslin[n][0];
    j=p->gcslin[n][1];
    io(i,j) = 1.0;
    }
    
    for(n=0;n<p->gcslout_count;n++)
    {
    i=p->gcslout[n][0];
    j=p->gcslout[n][1];
    io(i,j) += 2.0;
    }
    
    // raw switch: 1 dispersive, 0 nonlinear shallow water equations
    SLICELOOP4
    {
    f(i,j) = 1.0;
    
    if(p->wet[IJ]==0 || p->deep[IJ]==0 || b->breaking(i,j)>0)
    f(i,j) = 0.0;
    
    if(p->wet[Im1J]==0 || p->wet[Ip1J]==0)
    f(i,j) = 0.0;
    
    if(p->j_dir==1 && (p->wet[IJm1]==0 || p->wet[IJp1]==0))
    f(i,j) = 0.0;
    
    if(b->depth(i,j)<=10.0*p->A244 || WL(i,j)<=10.0*p->A244)
    f(i,j) = 0.0;
    
    // breaking (Tonelli & Petti 2009, FUNWAVE-TVD): |eta|/h >= 0.8
    // (crests and deep troughs / rundown, where the total depth H -> 0)
    if(b->depth(i,j)>0.0 && fabs(b->eta(i,j))/b->depth(i,j) >= breakratio)
    f(i,j) = 0.0;
    
    // supercritical flow (backwash, bores): Fr = |u|/sqrt(g H) >= 1
    if(WL(i,j)>p->A244 && b->U(i,j)*b->U(i,j) + b->V(i,j)*b->V(i,j)*p->y_dir >= fabs(p->W22)*WL(i,j))
    f(i,j) = 0.0;
    
    // in- and outflow boundaries (ghost values from ioflow)
    if(io(i,j)>0.5)
    f(i,j) = 0.0;
    }
    
    pgc->gcsl_start4(p,f,1);
    
    // smooth transition: widen the NSWE region by one cell, then two smoothing passes;
    // the dispersive terms are weighted with mask in [0,1]
    SLICELOOP4
    {
    T1(i,j) = MIN(f(i,j),MIN(f(i-1,j),f(i+1,j)));
    
    if(p->j_dir==1)
    T1(i,j) = MIN(T1(i,j),MIN(f(i,j-1),f(i,j+1)));
    }
    
    pgc->gcsl_start4(p,T1,1);
    
    for(int qn=0;qn<2;++qn)
    {
        SLICELOOP4
        {
        if(p->j_dir==0)
        T2(i,j) = 0.5*T1(i,j) + 0.25*(T1(i-1,j) + T1(i+1,j));
        
        if(p->j_dir==1)
        T2(i,j) = 0.5*T1(i,j) + 0.125*(T1(i-1,j) + T1(i+1,j) + T1(i,j-1) + T1(i,j+1));
        }
        
        pgc->gcsl_start4(p,T2,1);
        
        SLICELOOP4
        T1(i,j) = T2(i,j);
        
        pgc->gcsl_start4(p,T1,1);
    }
    
    SLICELOOP4
    {
    mask(i,j) = f(i,j)*T1(i,j);
    
    if(mask(i,j)<1.0e-6)
    mask(i,j) = 0.0;
    }
    
    pgc->gcsl_start4(p,mask,1);
    
    // FUNWAVE-TVD practice: V is kept where the switch changes (momentum
    // conserving), u_a follows from the next inversion
}

// ---------------------------------------------------------------------------
// wide-stencil (central) gradients of mA, mB, u2 and V1' for the explicit
// terms psi (no response to grid-scale modes in the explicit terms)
// ---------------------------------------------------------------------------
void sflow_boussinesq::explicitterms(lexer *p, fdm2D *b, ghostcell *pgc)
{
    double z,h,e;
    
    cellterms(p,b,pgc);
    
    SLICELOOP4
    {
    T1(i,j) = 0.5*b->eta(i,j)*b->eta(i,j)*B(i,j) + b->eta(i,j)*A(i,j);
    
    Ax(i,j) = ddx(p,A);
    Ay(i,j) = ddy(p,A);
    Bx(i,j) = ddx(p,B);
    By(i,j) = ddy(p,B);
    }
    
    pgc->gcsl_start4(p,T1,1);
    
    SLICELOOP4
    {
    U4(i,j)  = V4(i,j)  = 0.0;
    U1p(i,j) = V1p(i,j) = 0.0;
    
        if(p->wet[IJ]==1)
        {
        z = za(p,b,i,j);
        h = b->depth(i,j);
        e = b->eta(i,j);
        
        U4(i,j) = (0.5*z*z - (h*h - h*e + e*e)/6.0)*Bx(i,j) + (z + 0.5*(h - e))*Ax(i,j);
        V4(i,j) = (0.5*z*z - (h*h - h*e + e*e)/6.0)*By(i,j) + (z + 0.5*(h - e))*Ay(i,j);
        
        U1p(i,j) = 0.5*z*z*Bx(i,j) + z*Ax(i,j) - ddx(p,T1);
        V1p(i,j) = 0.5*z*z*By(i,j) + z*Ay(i,j) - ddy(p,T1);
        }
    }
    
    pgc->gcsl_start4(p,U4,1);
    pgc->gcsl_start4(p,V4,1);
}

// ---------------------------------------------------------------------------
// dispersive source psi (RK stage, stage state in UA,VA,eta; WL stage depth)
// ---------------------------------------------------------------------------
void sflow_boussinesq::source(lexer *p, fdm2D *b, ghostcell *pgc, slice &WL)
{
    double z,e,zx,zy,om0,om2,psix,psiy;
    
    explicitterms(p,b,pgc);
    
    SLICELOOP4
    {
    etat(i,j) = 0.0;
    T1(i,j) = 0.0;
    T2(i,j) = 0.0;
    
        if(p->wet[IJ]==1)
        {
        if(p->geo_curv==1)
        etat(i,j) = -(b->FEx(i,j)*p->LEN1[IJ] - b->FEx(i-1,j)*p->LEN1[Im1J]
                    + b->FEy(i,j)*p->LEN2[IJ] - b->FEy(i,j-1)*p->LEN2[IJm1])/p->AREA2D[IJ];
        
        else
        etat(i,j) = -(b->FEx(i,j) - b->FEx(i-1,j))/p->DXN[IP] 
                    -(b->FEy(i,j) - b->FEy(i,j-1))/p->DYN[JP]*p->y_dir;
        
        z = za(p,b,i,j);
        e = b->eta(i,j);
        
        T1(i,j) = etat(i,j)*(A(i,j) + e*B(i,j));
        
        T2(i,j) = (z - e)*(b->UA(i,j)*Ax(i,j) + b->VA(i,j)*Ay(i,j))
                + 0.5*(z*z - e*e)*(b->UA(i,j)*Bx(i,j) + b->VA(i,j)*By(i,j))
                + 0.5*(A(i,j) + e*B(i,j))*(A(i,j) + e*B(i,j));
        }
    }
    
    pgc->gcsl_start4(p,T1,1);
    pgc->gcsl_start4(p,T2,1);
    
    SLICELOOP4
    if(p->wet[IJ]==1 && mask(i,j)>0.0)
    {
    z  = za(p,b,i,j);
    zx = zeta*ddx(p,b->depth) + (1.0+zeta)*ddx(p,b->eta);
    zy = zeta*ddy(p,b->depth) + (1.0+zeta)*ddy(p,b->eta);
    
    om0 = ddx(p,b->VA) - ddy(p,b->UA);
    om2 = zx*(Ay(i,j) + z*By(i,j)) - zy*(Ax(i,j) + z*Bx(i,j));
    
    psix = etat(i,j)*(U1p(i,j) - U4(i,j))
         + WL(i,j)*( b->UA(i,j)*ddx(p,U4) + b->VA(i,j)*ddy(p,U4)
                   + U4(i,j)*ddx(p,b->UA) + V4(i,j)*ddy(p,b->UA)
                   - ddx(p,T1) - ddx(p,T2)
                   + om0*V4(i,j) + om2*b->VA(i,j));
    
    b->F(i,j) += mask(i,j)*psix;
    
        if(p->j_dir==1)
        {
        psiy = etat(i,j)*(V1p(i,j) - V4(i,j))
             + WL(i,j)*( b->UA(i,j)*ddx(p,V4) + b->VA(i,j)*ddy(p,V4)
                       + U4(i,j)*ddx(p,b->VA) + V4(i,j)*ddy(p,b->VA)
                       - ddy(p,T1) - ddy(p,T2)
                       - om0*U4(i,j) - om2*b->UA(i,j));
        
        b->G(i,j) += mask(i,j)*psiy;
        }
    }
}

// ---------------------------------------------------------------------------
// u_a from V = H (u_a + V1'(u_a)), line-implicit per component, 
// tangential derivatives lagged: V1' = L(u_a) + C, C from the current u_a
// ---------------------------------------------------------------------------
void sflow_boussinesq::invert(lexer *p, fdm2D *b, ghostcell *pgc, solver2D *psolv, slice &Vx, slice &Vy, slice &WL)
{
    invert_prepare(p,b,pgc,Vx,Vy,WL);
    invert_solve(p,b,psolv,pgc);
}

// rows of the two line systems for u_a and v_a: x-x and y-y derivatives implicit, the cross
// derivatives (Cx, Cy) from the u_a at hand
void sflow_boussinesq::invert_prepare(lexer *p, fdm2D *b, ghostcell *pgc, slice &Vx, slice &Vy, slice &WL)
{
    double c1,c2,c3;
    
    if(p->geo_curv==1)
    {
    prepare_curv(p,b,pgc,Vx,Vy,WL);
    return;
    }
    
    operators(p,b,pgc);
    
    SLICELOOP4
    {
    Cx(i,j) = 0.0;
    Cy(i,j) = 0.0;
    
        if(p->wet[IJ]==1)
        {
        Cx(i,j) = U1p(i,j) - lx(p,b,b->UA);
        
        if(p->j_dir==1)
        Cy(i,j) = V1p(i,j) - ly(p,b,b->VA);
        }
    }
    
    // x-component
    SLICELOOP4
    {
    Xp(i,j) = 1.0;
    Xn(i,j) = 0.0;
    Xs(i,j) = 0.0;
    Xr(i,j) = 0.0;
    
        if(p->wet[IJ]==1)
        {
        coef_x(p,b,c1,c2,c3);
        
        Xp(i,j) = 1.0 + c3;
        Xn(i,j) = c1;
        Xs(i,j) = c2;
        Xr(i,j) = Vx(i,j)/HWL - Cx(i,j);
        }
    
        // walls: u antisymmetric (implicit), in- and outflow: ghost values
        if(p->flagslice4[Im1J]<0)
        {
        if(int(io(i,j)+0.5)%2==0)
        Xp(i,j) -= Xs(i,j);
        
        else
        Xr(i,j) -= Xs(i,j)*b->UA(i-1,j);
        
        Xs(i,j) = 0.0;
        }
        
        if(p->flagslice4[Ip1J]<0)
        {
        if(int(io(i,j)+0.5)<2)
        Xp(i,j) -= Xn(i,j);
        
        else
        Xr(i,j) -= Xn(i,j)*b->UA(i+1,j);
        
        Xn(i,j) = 0.0;
        }
    }
    
    // y-component
    if(p->j_dir==1)
    SLICELOOP4
    {
    Yp(i,j) = 1.0;
    Yw(i,j) = 0.0;
    Ye(i,j) = 0.0;
    Yr(i,j) = 0.0;
    
        if(p->wet[IJ]==1)
        {
        coef_y(p,b,c1,c2,c3);
        
        Yp(i,j) = 1.0 + c3;
        Yw(i,j) = c1;
        Ye(i,j) = c2;
        Yr(i,j) = Vy(i,j)/HWL - Cy(i,j);
        }
    
        // side walls: v antisymmetric (implicit)
        if(p->flagslice4[IJm1]<0)
        {
        Yp(i,j) -= Ye(i,j);
        Ye(i,j) = 0.0;
        }
        
        if(p->flagslice4[IJp1]<0)
        {
        Yp(i,j) -= Yw(i,j);
        Yw(i,j) = 0.0;
        }
    }
}

// the line systems of invert_prepare; outside the computed cells (the partition halo, the filled
// cells of a refined patch) f carries u_a, which the patch line solver takes as boundary value,
// so that sflow_amr can repeat the solve with updated values there
void sflow_boussinesq::invert_solve(lexer *p, fdm2D *b, solver2D *psolv, ghostcell *pgc)
{
    if(p->geo_curv==1)
    {
    solve_curv(p,b,pgc);
    return;
    }
    
    for(int ii=p->imin; ii<p->imin+p->imax; ++ii)
    for(int jj=p->jmin; jj<p->jmin+p->jmax; ++jj)
    f(ii,jj) = b->UA(ii,jj);
    
    n=0;
    SLICELOOP4
    {
    b->M.p[n] = Xp(i,j);
    b->M.n[n] = Xn(i,j);
    b->M.s[n] = Xs(i,j);
    b->M.e[n] = 0.0;
    b->M.w[n] = 0.0;
    b->rhsvec.V[n] = Xr(i,j);
    ++n;
    }
    
    psolv->start(p,pgc,f,b->M,b->xvec,b->rhsvec,4);
    
    SLICELOOP4
    b->UA(i,j) = (p->wet[IJ]==1)?f(i,j):0.0;
    
    if(p->j_dir==1)
    {
        for(int ii=p->imin; ii<p->imin+p->imax; ++ii)
        for(int jj=p->jmin; jj<p->jmin+p->jmax; ++jj)
        f(ii,jj) = b->VA(ii,jj);
        
        n=0;
        SLICELOOP4
        {
        b->M.p[n] = Yp(i,j);
        b->M.n[n] = 0.0;
        b->M.s[n] = 0.0;
        b->M.e[n] = Ye(i,j);
        b->M.w[n] = Yw(i,j);
        b->rhsvec.V[n] = Yr(i,j);
        ++n;
        }
        
        psolv->start(p,pgc,f,b->M,b->xvec,b->rhsvec,4);
        
        SLICELOOP4
        b->VA(i,j) = (p->wet[IJ]==1)?f(i,j):0.0;
    }
}

// volume flux M = H (u_a + u2), requires the ghost cells of UA,VA
void sflow_boussinesq::flux(lexer *p, fdm2D *b, ghostcell *pgc, slice &WL)
{
    operators(p,b,pgc);
    
    SLICELOOP4
    {
    b->MX(i,j) = 0.0;
    b->MY(i,j) = 0.0;
    
        if(p->wet[IJ]==1)
        {
        b->MX(i,j) = WL(i,j)*(b->UA(i,j) + U4(i,j));
        b->MY(i,j) = WL(i,j)*(b->VA(i,j) + V4(i,j))*p->y_dir;
        }
    }
}

// V from u_a:  mode 0: cells whose V1' stencil contains a u_a changed since save(), 1: all cells
void sflow_boussinesq::forward(lexer *p, fdm2D *b, ghostcell *pgc, slice &Vx, slice &Vy, slice &WL, int mode)
{
    int sel;
    
    operators(p,b,pgc);
    
    if(mode==0)
    {
        SLICELOOP4
        f(i,j) = (fabs(b->UA(i,j)-ua_n(i,j))>0.0 || fabs(b->VA(i,j)-va_n(i,j))>0.0)?1.0:0.0;
        
        pgc->gcsl_start4(p,f,1);
    }
    
    SLICELOOP4
    {
    sel = 0;
    
    if(mode==1)
    sel = 1;
    
    if(mode==0)
    for(int qi=-1;qi<=1;++qi)
    for(int qj=-p->j_dir;qj<=p->j_dir;++qj)
    if(f(i+qi,j+qj)>0.5)
    sel = 1;
    
        if(sel==1)
        {
        Vx(i,j) = WL(i,j)*(b->UA(i,j) + U1p(i,j));
        Vy(i,j) = WL(i,j)*(b->VA(i,j) + V1p(i,j))*p->y_dir;
        
            if(p->wet[IJ]==0)
            {
            Vx(i,j) = 0.0;
            Vy(i,j) = 0.0;
            }
        }
    }
}

void sflow_boussinesq::save(lexer *p, fdm2D *b)
{
    SLICELOOP4
    {
    ua_n(i,j) = b->UA(i,j);
    va_n(i,j) = b->VA(i,j);
    }
}

// ---------------------------------------------------------------------------
// curvilinear grid (G 50): the compact face operators of the Cartesian kernels with the face
// normals. Face f of cell P to neighbour N, outward face vector S (length L, unit normal n,
// tangent t), centre distance dn along n:
//   grad g at f = n (g_N - g_P)/dn + t (t . mean of the cell gradients of P and N)
//   B_f = (grad u)_x + (grad v)_y,  A_f the same for h u, h v
//   grad(m B)_x = sum_f m_f B_f S_x / A_P,  ... ;  the Cartesian grid gives the Cartesian kernels
// ---------------------------------------------------------------------------
void sflow_boussinesq::operators_curv(lexer *p, fdm2D *b, ghostcell *pgc)
{
    cellterms(p,b,pgc);
    
    const int di[4]={1,-1,0,0}, dj[4]={0,0,1,-1};
    
    SLICELOOP4
    {
    Ax(i,j)  = Bx(i,j)  = Ay(i,j)  = By(i,j) = 0.0;
    U4(i,j)  = V4(i,j)  = 0.0;
    U1p(i,j) = V1p(i,j) = 0.0;
    
        if(p->wet[IJ]==1)
        {
        const double z = za(p,b,i,j);
        const double h = b->depth(i,j);
        const double e = b->eta(i,j);
        const double cB = 0.5*z*z - (h*h - h*e + e*e)/6.0;
        const double cA = z + 0.5*(h - e);
        double bx=0.0, by=0.0, ax=0.0, ay=0.0, fx=0.0, fy=0.0;
        
            for(int d=0; d<4; ++d)
            {
            const int a = i+di[d], c = j+dj[d];
            const int dir = (d<2) ? 1 : 2;
            const int fi = (d==1) ? i-1 : i, fj = (d==3) ? j-1 : j;
            const double sg = (d==0 || d==2) ? 1.0 : -1.0;
            double nx,ny,L;
            sflow_curv::normal(p,dir,fi,fj,nx,ny,L);
            nx*=sg; ny*=sg;
            const double tx = -ny, ty = nx;
            const double sx = nx*L, sy = ny*L;
            const double dn = sflow_curv::dnorm(p,dir,fi,fj);
            
            // face gradients: normal part across the face, tangential part from the cells
            auto fg = [&](slice &g, slice &gx, slice &gy, double &rx, double &ry)
            {
                const double dnv = (g(a,c) - g(i,j))/dn;
                const double tg = tx*0.5*(gx(i,j)+gx(a,c)) + ty*0.5*(gy(i,j)+gy(a,c));
                rx = nx*dnv + tx*tg;
                ry = ny*dnv + ty*tg;
            };
            
            double uxf,uyf,vxf,vyf,huxf,huyf,hvxf,hvyf;
            fg(b->UA,ux,uy,uxf,uyf);
            fg(b->VA,vx,vy,vxf,vyf);
            
            // h u_a, h v_a across the face
            {
            const double dnu = (b->depth(a,c)*b->UA(a,c) - h*b->UA(i,j))/dn;
            const double dnv = (b->depth(a,c)*b->VA(a,c) - h*b->VA(i,j))/dn;
            const double tgu = tx*0.5*(hux(i,j)+hux(a,c)) + ty*0.5*(huy(i,j)+huy(a,c));
            const double tgv = tx*0.5*(hvx(i,j)+hvx(a,c)) + ty*0.5*(hvy(i,j)+hvy(a,c));
            huxf = nx*dnu + tx*tgu; huyf = ny*dnu + ty*tgu;
            hvxf = nx*dnv + tx*tgv; hvyf = ny*dnv + ty*tgv;
            }
            (void)uyf; (void)vxf; (void)huyf; (void)hvxf;
            
            const double Bf = uxf + vyf;
            const double Af = huxf + hvyf;
            const double mf = 0.5*(mask(i,j) + mask(a,c));
            const double ef = 0.5*(e + b->eta(a,c));
            const double qf = 0.25*(e*e + b->eta(a,c)*b->eta(a,c));
            const double Ff = mf*(qf*Bf + ef*Af);
            
            bx += mf*Bf*sx;  by += mf*Bf*sy;
            ax += mf*Af*sx;  ay += mf*Af*sy;
            fx += Ff*sx;     fy += Ff*sy;
            }
        
        const double ar = p->AREA2D[IJ];
        Bx(i,j) = bx/ar; By(i,j) = by/ar;
        Ax(i,j) = ax/ar; Ay(i,j) = ay/ar;
        
        U4(i,j)  = cB*Bx(i,j) + cA*Ax(i,j);
        V4(i,j)  = cB*By(i,j) + cA*Ay(i,j);
        U1p(i,j) = 0.5*z*z*Bx(i,j) + z*Ax(i,j) - fx/ar;
        V1p(i,j) = 0.5*z*z*By(i,j) + z*Ay(i,j) - fy/ar;
        }
    }
    
    pgc->gcsl_start4(p,U4,1);
    pgc->gcsl_start4(p,V4,1);
}

// grid-aligned basis of cell (i,j): e1 along i (through the neighbour centres, one-sided in the
// outermost ghost layers), e2 = e1 turned by 90 deg; on a Cartesian grid (1,0) and (0,1)
void sflow_boussinesq::basis(lexer *p, int ii, int jj, double &ex, double &ey)
{
    const int q = (ii-p->imin)*p->jmax + (jj-p->jmin);
    const int qe = (ii+1<p->imin+p->imax) ? q+p->jmax : q;
    const int qw = (ii-1>=p->imin) ? q-p->jmax : q;
    const double dx = p->XC2D[qe] - p->XC2D[qw];
    const double dy = p->YC2D[qe] - p->YC2D[qw];
    const double l = sqrt(dx*dx + dy*dy);
    ex = l>0.0 ? dx/l : 1.0;
    ey = l>0.0 ? dy/l : 0.0;
}

// implicit part of V1' along e (dir 0: e1, 1: e2) in the component along e: the normal derivatives
// across the four faces, weight m_f (S.e)^2/(A L dn); c: i+1, i-1, j+1, j-1.  With the grid-aligned
// components the own component enters only through the normal derivatives (Cartesian grid: the
// x-x and y-y derivatives of the line systems), the other one through the lagged part
void sflow_boussinesq::coef_curv(lexer *p, fdm2D *b, int dir, double *c, double &cC)
{
    const int di[4]={1,-1,0,0}, dj[4]={0,0,1,-1};
    const double z = za(p,b,i,j);
    const double hP = b->depth(i,j);
    const double e = b->eta(i,j);
    const double ar = p->AREA2D[IJ];
    double ex,ey;
    basis(p,i,j,ex,ey);
    if(dir==1)
    {
    const double t = ex;
    ex = -ey;
    ey = t;
    }
    
    cC = 0.0;
    
    for(int d=0; d<4; ++d)
    {
    const int a = i+di[d], cc = j+dj[d];
    const int fdir = (d<2) ? 1 : 2;
    const int fi = (d==1) ? i-1 : i, fj = (d==3) ? j-1 : j;
    double nx,ny,L;
    sflow_curv::normal(p,fdir,fi,fj,nx,ny,L);
    const double sc = (nx*ex + ny*ey)*L;
    const double dn = sflow_curv::dnorm(p,fdir,fi,fj);
    const double mf = 0.5*(mask(i,j) + mask(a,cc));
    const double ef = 0.5*(e + b->eta(a,cc));
    const double qf = 0.25*(e*e + b->eta(a,cc)*b->eta(a,cc));
    const double w = mf*sc*sc/(ar*L*dn);
    
    c[d] = w*(0.5*z*z - qf + (z - ef)*b->depth(a,cc));
    cC  -= w*(0.5*z*z - qf + (z - ef)*hP);
    }
}

// the inversion on a curvilinear grid (solve_curv, sflow_amr::bous_window_curv): the right-hand
// sides V.e1/H, V.e2/H, the start values a = u_a.e1, b = u_a.e2 in UA, VA (u_a, v_a kept in
// ucart, vcart), and the rows of the preconditioner: the normal derivatives of the own component
// across the faces (coef_curv), walls with the mirrored state (Cartesian grid: u antisymmetric);
// in- and outflow neighbours left out
void sflow_boussinesq::prepare_curv(lexer *p, fdm2D *b, ghostcell *pgc, slice &Vx, slice &Vy, slice &WL)
{
    const int di[4]={1,-1,0,0}, dj[4]={0,0,1,-1};
    double c[4],cC;
    
    for(int n=0; n<p->imax*p->jmax; ++n)
    {
    ucart.V[n] = b->UA.V[n];
    vcart.V[n] = b->VA.V[n];
    }
    
    for(int dir=0; dir<2; ++dir)
    {
    slice4 &Pp = (dir==0) ? Xp : Yp;
    slice4 &Pu = (dir==0) ? Xn : Yw;    // i+1 / j+1
    slice4 &Pd = (dir==0) ? Xs : Ye;    // i-1 / j-1
    slice4 &Rh = (dir==0) ? Ra : Rb;
    
        SLICELOOP4
        {
        Pp(i,j) = 1.0;
        Pu(i,j) = Pd(i,j) = 0.0;
        Rh(i,j) = 0.0;
        
            if(p->wet[IJ]==1)
            {
            double ex,ey;
            basis(p,i,j,ex,ey);
            if(dir==1)
            {
            const double t = ex;
            ex = -ey;
            ey = t;
            }
            
            coef_curv(p,b,dir,c,cC);
            
            Pp(i,j) = 1.0 + cC;
            Rh(i,j) = (Vx(i,j)*ex + Vy(i,j)*ey)/HWL;
            
                for(int d=0; d<4; ++d)
                {
                const int a = i+di[d], cc = j+dj[d];
                const bool line = (dir==0) ? (d<2) : (d>=2);
                
                if(p->flagslice4[(a-p->imin)*p->jmax+(cc-p->jmin)]>0)
                {
                if(line)
                (d==0 || d==2 ? Pu : Pd)(i,j) = c[d];
                continue;
                }
                
                const int iof = int(io(i,j)+0.5);
                if((d==1 && iof%2==1) || (d==0 && iof>=2))
                continue;
                
                double nx,ny,L;
                sflow_curv::normal(p,(d<2)?1:2,(d==1)?i-1:i,(d==3)?j-1:j,nx,ny,L);
                const double ne = nx*ex + ny*ey;
                Pp(i,j) += c[d]*(1.0 - 2.0*ne*ne);
                }
            }
        }
    }
    
    // a, b of every cell with neighbours (own basis): the start values
    for(int ii=p->imin+1; ii<p->imin+p->imax-1; ++ii)
    for(int jj=p->jmin+1; jj<p->jmin+p->jmax-1; ++jj)
    {
    double ex,ey;
    basis(p,ii,jj,ex,ey);
    const double u = b->UA(ii,jj), v = b->VA(ii,jj);
    b->UA(ii,jj) =  u*ex + v*ey;
    b->VA(ii,jj) = -u*ey + v*ex;
    }
}

// ---------------------------------------------------------------------------
// curvilinear grid: the inversion for (a, b) by BiCGStab (right preconditioned), matrix free
// ---------------------------------------------------------------------------
void sflow_boussinesq::set_curv(lexer *p, fdm2D *b, slice &xa, slice &xb, bool all)
{
    for(int n=0; n<p->imax*p->jmax; ++n)
    {
    b->UA.V[n] = ucart.V[n];
    b->VA.V[n] = vcart.V[n];
    }
    
    auto put = [&](int ii, int jj)
    {
        double ex,ey;
        basis(p,ii,jj,ex,ey);
        const double a = xa(ii,jj), c = xb(ii,jj);
        b->UA(ii,jj) = a*ex - c*ey;
        b->VA(ii,jj) = a*ey + c*ex;
    };
    
    if(all)
    {
        for(int ii=p->imin; ii<p->imin+p->imax; ++ii)
        for(int jj=p->jmin; jj<p->jmin+p->jmax; ++jj)
        if(p->flagslice4[(ii-p->imin)*p->jmax+(jj-p->jmin)]>0)
        put(ii,jj);
    }
    else
    SLICELOOP4
    put(i,j);
}

void sflow_boussinesq::apply_curv(lexer *p, fdm2D *b, ghostcell *pgc, slice &xa, slice &xb, slice &ya, slice &yb, bool all)
{
    set_curv(p,b,xa,xb,all);
    
    if(velbc)
    velbc(p,b,pgc);
    
    operators_curv(p,b,pgc);
    
    SLICELOOP4
    {
        if(p->wet[IJ]==1)
        {
        double ex,ey;
        basis(p,i,j,ex,ey);
        const double ux_ = b->UA(i,j) + U1p(i,j), vy_ = b->VA(i,j) + V1p(i,j);
        ya(i,j) =  ux_*ex + vy_*ey;
        yb(i,j) = -ux_*ey + vy_*ex;
        }
        else
        {
        ya(i,j) = xa(i,j);
        yb(i,j) = xb(i,j);
        }
    }
}

void sflow_boussinesq::prec_curv(lexer *p, int dir, slice &r, slice &z)
{
    std::vector<double> cp, dp;
    
    // a along i (Xs: i-1, Xn: i+1), b along j (Ye: j-1, Yw: j+1); runs of computed cells
    slice4 &P = (dir==0) ? Xp : Yp;
    slice4 &U = (dir==0) ? Xn : Yw;
    slice4 &D = (dir==0) ? Xs : Ye;
    const int na = (dir==0) ? p->knoy : p->knox;
    const int nl = (dir==0) ? p->knox : p->knoy;
    
    auto cell = [&](int t, int k, int &ii, int &jj) { ii = (dir==0) ? k : t; jj = (dir==0) ? t : k; };
    auto in = [&](int ii, int jj) { return p->flagslice4[(ii-p->imin)*p->jmax+(jj-p->jmin)]>0; };
    
    for(int t=0; t<na; ++t)
    {
        int k=0;
        while(k<nl)
        {
            int ii,jj;
            cell(t,k,ii,jj);
            if(!in(ii,jj)) { z(ii,jj) = 0.0; ++k; continue; }
            int k0=k;
            while(k<nl) { cell(t,k,ii,jj); if(!in(ii,jj)) break; ++k; }
            const int m = k-k0;
            cp.resize(m); dp.resize(m);
            
            for(int s=0; s<m; ++s)
            {
                cell(t,k0+s,ii,jj);
                const double lo = (s>0) ? D(ii,jj) : 0.0;
                const double den = P(ii,jj) - ((s>0) ? lo*cp[s-1] : 0.0);
                cp[s] = U(ii,jj)/den;
                dp[s] = (r(ii,jj) - ((s>0) ? lo*dp[s-1] : 0.0))/den;
            }
            for(int s=m-2; s>=0; --s)
            dp[s] -= cp[s]*dp[s+1];
            
            for(int s=0; s<m; ++s)
            {
                cell(t,k0+s,ii,jj);
                z(ii,jj) = dp[s];
            }
        }
    }
}

void sflow_boussinesq::solve_curv(lexer *p, fdm2D *b, ghostcell *pgc)
{
    enum { X=0, R, RH, PV, VV, S, T, PH, SH, C0, ZA, ZB, OA, OB, NA, NB, NK };
    while((int)kv.size()<NK)
    kv.push_back(new slice4(p));
    
    auto v = [&](int k) -> slice& { return *kv[k]; };
    auto dot = [&](int k, int l)
    {
        double s=0.0;
        SLICELOOP4
        s += v(k)(i,j)*v(l)(i,j);
        return pgc->globalsum(s);
    };
    
    // a, b of prepare_curv (UA, VA hold them): start and the lagged other component
    for(int n=0; n<p->imax*p->jmax; ++n)
    {
    v(OA).V[n] = b->UA.V[n];
    v(OB).V[n] = b->VA.V[n];
    }
    
    int itsum=0;
    
    for(int dir=0; dir<2; ++dir)
    {
        slice &own = v(dir==0 ? OA : OB);
        slice &Rh = (dir==0) ? Ra : Rb;
        
        // y = own component of the operator for the own component x and the other one o
        auto opr = [&](slice &x, slice &o, slice &y)
        {
            if(dir==0)
            apply_curv(p,b,pgc,x,o,y,v(ZB),false);
            else
            apply_curv(p,b,pgc,o,x,v(ZA),y,false);
        };
        
        for(int n=0; n<p->imax*p->jmax; ++n)
        {
        v(X).V[n] = own.V[n];
        v(T).V[n] = 0.0;
        }
        
        // constant part (in- and outflow ghost values), A x0 with the lagged other component
        opr(v(T),v(T),v(C0));
        opr(v(X),v(dir==0 ? OB : OA),v(VV));
        
        double bn=0.0;
        SLICELOOP4
        {
        v(R)(i,j) = Rh(i,j) - v(VV)(i,j);
        v(RH)(i,j) = v(R)(i,j);
        v(PV)(i,j) = v(VV)(i,j) = 0.0;
        bn += Rh(i,j)*Rh(i,j);
        }
        bn = sqrt(pgc->globalsum(bn));
        double rn = sqrt(dot(R,R));
        
        // the linear part: own component alone, without the constant part
        for(int n=0; n<p->imax*p->jmax; ++n)
        v(T).V[n] = 0.0;
        auto apply = [&](int kx, int ky)
        {
            for(int n=0; n<p->imax*p->jmax; ++n)
            v(ZA).V[n] = v(ZB).V[n] = 0.0;
            if(dir==0)
            apply_curv(p,b,pgc,v(kx),v(ZA),v(ky),v(ZB),false);
            else
            apply_curv(p,b,pgc,v(ZA),v(kx),v(ZB),v(ky),false);
            SLICELOOP4
            v(ky)(i,j) -= v(C0)(i,j);
        };
        
        double rho=1.0, alp=1.0, om=1.0;
        int it=0;
        
        if(bn<=0.0)
        {
            SLICELOOP4
            v(X)(i,j) = 0.0;
        }
        else
        while(rn/bn>p->N44 && it<p->N46)
        {
            double rhon = dot(RH,R);
            if(it>0 && fabs(rhon)<=1.0e-12*rn*sqrt(dot(RH,RH)))
            {
                SLICELOOP4
                {
                v(RH)(i,j) = v(R)(i,j);
                v(PV)(i,j) = v(VV)(i,j) = 0.0;
                }
                rho = alp = om = 1.0;
                rhon = rn*rn;
            }
            const double beta = (rhon/rho)*(alp/om);
            rho = rhon;
            
            SLICELOOP4
            v(PV)(i,j) = v(R)(i,j) + beta*(v(PV)(i,j) - om*v(VV)(i,j));
            prec_curv(p,dir,v(PV),v(PH));
            apply(PH,VV);
            alp = rho/dot(RH,VV);
            
            SLICELOOP4
            v(S)(i,j) = v(R)(i,j) - alp*v(VV)(i,j);
            prec_curv(p,dir,v(S),v(SH));
            apply(SH,T);
            const double tt = dot(T,T);
            om = tt>0.0 ? dot(T,S)/tt : 0.0;
            
            SLICELOOP4
            {
            v(X)(i,j) += alp*v(PH)(i,j) + om*v(SH)(i,j);
            v(R)(i,j) = v(S)(i,j) - om*v(T)(i,j);
            }
            rn = sqrt(dot(R,R));
            ++it;
        }
        itsum += it;
        
        // the solution of this component (the other solve keeps the old one: lagged)
        slice &res = v(dir==0 ? NA : NB);
        for(int n=0; n<p->imax*p->jmax; ++n)
        res.V[n] = own.V[n];
        SLICELOOP4
        res(i,j) = (p->wet[IJ]==1) ? v(X)(i,j) : 0.0;
    }
    
    // u_a, v_a: the solutions in the computed cells, the rest as before the solve
    set_curv(p,b,v(NA),v(NB),false);
    
    last_it = itsum;
}

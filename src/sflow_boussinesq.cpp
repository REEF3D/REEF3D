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

#define HWL (WL(i,j)>p->A244?WL(i,j):1.0e20)

sflow_boussinesq::sflow_boussinesq(lexer *p, fdm2D *b, ghostcell *pgc) : zeta(-0.53), breakratio(0.8),
                                    A(p),B(p),Ax(p),Ay(p),Bx(p),By(p),
                                    U4(p),V4(p),U1p(p),V1p(p),Cx(p),Cy(p),
                                    T1(p),T2(p),etat(p),
                                    vy(p),hvy(p),ux(p),hux(p),
                                    mask(p),io(p),ua_n(p),va_n(p),f(p)
{
    if(p->mpirank==0)
    cout<<"SFLOW Boussinesq equations (FUNWAVE-TVD formulation)"<<endl;
}

sflow_boussinesq::~sflow_boussinesq()
{
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
    return (g(i+1,j) - g(i-1,j))/(p->DXP[IP] + p->DXP[IM1]);
}

double sflow_boussinesq::ddy(lexer *p, slice &g)
{
    if(p->j_dir==0)
    return 0.0;
    
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
// dispersion switch: breaking (A 246), eta/h >= 0.8, eta/h <= -0.5, shoreline, in- and outflow
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
    
    // breaking (Tonelli & Petti 2009, FUNWAVE-TVD): eta/h >= 0.8
    if(b->depth(i,j)>0.0 && b->eta(i,j)/b->depth(i,j) >= breakratio)
    f(i,j) = 0.0;
    
    // strong drawdown: eta/h <= -0.5
    if(b->depth(i,j)>0.0 && b->eta(i,j)/b->depth(i,j) <= -0.5)
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
    double c1,c2,c3;
    
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
    n=0;
    SLICELOOP4
    {
    b->M.p[n] = 1.0;
    b->M.n[n] = 0.0;
    b->M.s[n] = 0.0;
    b->M.e[n] = 0.0;
    b->M.w[n] = 0.0;
    b->rhsvec.V[n] = 0.0;
    
        if(p->wet[IJ]==1)
        {
        coef_x(p,b,c1,c2,c3);
        
        b->M.p[n] = 1.0 + c3;
        b->M.n[n] = c1;
        b->M.s[n] = c2;
        b->rhsvec.V[n] = Vx(i,j)/HWL - Cx(i,j);
        }
    
    f(i,j) = b->UA(i,j);
    ++n;
    }
    
    n=0;
    SLICELOOP4
    {
        // walls: u antisymmetric (implicit), in- and outflow: ghost values
        if(p->flagslice4[Im1J]<0)
        {
        if(int(io(i,j)+0.5)%2==0)
        b->M.p[n] -= b->M.s[n];
        
        else
        b->rhsvec.V[n] -= b->M.s[n]*b->UA(i-1,j);
        
        b->M.s[n] = 0.0;
        }
        
        if(p->flagslice4[Ip1J]<0)
        {
        if(int(io(i,j)+0.5)<2)
        b->M.p[n] -= b->M.n[n];
        
        else
        b->rhsvec.V[n] -= b->M.n[n]*b->UA(i+1,j);
        
        b->M.n[n] = 0.0;
        }
    ++n;
    }
    
    psolv->start(p,pgc,f,b->M,b->xvec,b->rhsvec,4);
    
    SLICELOOP4
    b->UA(i,j) = (p->wet[IJ]==1)?f(i,j):0.0;
    
    // y-component
    if(p->j_dir==1)
    {
        n=0;
        SLICELOOP4
        {
        b->M.p[n] = 1.0;
        b->M.n[n] = 0.0;
        b->M.s[n] = 0.0;
        b->M.e[n] = 0.0;
        b->M.w[n] = 0.0;
        b->rhsvec.V[n] = 0.0;
        
            if(p->wet[IJ]==1)
            {
            coef_y(p,b,c1,c2,c3);
            
            b->M.p[n] = 1.0 + c3;
            b->M.w[n] = c1;
            b->M.e[n] = c2;
            b->rhsvec.V[n] = Vy(i,j)/HWL - Cy(i,j);
            }
        
        f(i,j) = b->VA(i,j);
        ++n;
        }
        
        n=0;
        SLICELOOP4
        {
            // side walls: v antisymmetric (implicit)
            if(p->flagslice4[IJm1]<0)
            {
            b->M.p[n] -= b->M.e[n];
            b->M.e[n] = 0.0;
            }
            
            if(p->flagslice4[IJp1]<0)
            {
            b->M.p[n] -= b->M.w[n];
            b->M.w[n] = 0.0;
            }
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

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

#include"6DOF_obj.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include"vrans_definitions.h"

void sixdof_obj::hydrodynamic_forces_nhflow_volume
(
    lexer *p, fdm_nhf *d, ghostcell *pgc,
    double *FX, double *FY, double *FZ,
    slice &WL, int iter, bool finalize
)
{
    double dV,H;
    double rx,ry,rz;
    double fx,fy,fz;
    double Sfx,Sfy,Sfz,SKx,SKy,SKz;
    double dts,fsf_z;
    double Vsub,xB,yB,zB,Fb;

    Xe=Ye=Ze=Ke=Me=Ne=0.0;
    Sfx=Sfy=Sfz=SKx=SKy=SKz=0.0;

    LOOP
    {
        H = Hsolidface_nhflow(p,d,0,0,0);   // this body's own Heaviside from d->FB

        if(H<1.0e-12)
        continue;

        dV = p->DXN[IP]*p->DYN[JP]*p->DZN[KP]*WL(i,j);

        rx = p->pos_x() - c_(0);
        ry = p->pos_y() - c_(1);
        rz = p->pos_z() - c_(2);          // ZSP[IJK] under G2==1

        // the forcing exactly as the momentum update applied it
        fx = CPORNH*FX[IJK];
        fy = CPORNH*FY[IJK];
        fz = CPORNH*FZ[IJK];

        // reaction on the body
        Xe -= p->W1*fx*dV;
        Ye -= p->W1*fy*dV;
        Ze -= p->W1*fz*dV;

        Ke -= p->W1*(ry*fz - rz*fy)*dV;
        Me -= p->W1*(rz*fx - rx*fz)*dV;
        Ne -= p->W1*(rx*fy - ry*fx)*dV;

        // momentum of the fictitious fluid inside the body region.
        // SAME cells, SAME weights, SAME dV as the sum above - this
        // consistency is what the Kempe & Froehlich correction relies on.
        Sfx += p->W1*H*d->U[IJK]*dV;
        Sfy += p->W1*H*d->V[IJK]*dV;
        Sfz += p->W1*H*d->W[IJK]*dV;

        SKx += p->W1*H*(ry*d->W[IJK] - rz*d->V[IJK])*dV;
        SKy += p->W1*H*(rz*d->U[IJK] - rx*d->W[IJK])*dV;
        SKz += p->W1*H*(rx*d->V[IJK] - ry*d->U[IJK])*dV;
    }

    Xe = pgc->globalsum(Xe);
    Ye = pgc->globalsum(Ye);
    Ze = pgc->globalsum(Ze);
    Ke = pgc->globalsum(Ke);
    Me = pgc->globalsum(Me);
    Ne = pgc->globalsum(Ne);

    Sfx = pgc->globalsum(Sfx);
    Sfy = pgc->globalsum(Sfy);
    Sfz = pgc->globalsum(Sfz);
    SKx = pgc->globalsum(SKx);
    SKy = pgc->globalsum(SKy);
    SKz = pgc->globalsum(SKz);

    // ---- fictitious-mass correction: + rho_f d/dt int_Omega_b u dV
    dts = alpha[iter]*p->dt;

    if(fictmass_ini==1)
    {
        Xe += (Sfx - Sfx_n)/dts;
        Ye += (Sfy - Sfy_n)/dts;
        Ze += (Sfz - Sfz_n)/dts;

        Ke += (SKx - SKx_n)/dts;
        Me += (SKy - SKy_n)/dts;
        Ne += (SKz - SKz_n)/dts;
    }

    Sfx_n=Sfx; Sfy_n=Sfy; Sfz_n=Sfz;
    SKx_n=SKx; SKy_n=SKy; SKz_n=SKz;
    fictmass_ini=1;

    // ---- buoyancy
    // NHFLOW carries no gravity term in the vertical momentum equation
    // (wpgrad is empty, no W22 in the W path), so the hydrostatic pressure
    // is never formed and the forcing reaction contains NO buoyancy.
    // The horizontal hydrostatic (Froude-Krylov) part IS in the sum above,
    // because f_x has to fight -g deta/dx. Only the vertical piece is missing.
    fsf_z = p->wd;      // single plane for the whole body - see note below

    buoyancy_nhflow(p,d,pgc,fsf_z,Vsub,xB,yB,zB);

    Fb = p->W1*fabs(p->W22)*Vsub;

    Ze += Fb;
    Ke += (yB - c_(1))*Fb;
    Me -= (xB - c_(0))*Fb;

    // ---- weight
    Xe += p->W20*Mass_fb;
    Ye += p->W21*Mass_fb;
    Ze += p->W22*Mass_fb;

    if(p->mpirank==0 && finalize==1)
    {
    cout<<"Vsub: "<<Vsub<<" (m/rho_f = "<<Mass_fb/p->W1<<")"
        <<" xB: "<<xB<<" zB: "<<zB<<endl;
    cout<<"Xe: "<<Xe<<" Ze: "<<Ze<<" Me: "<<Me<<endl;
    }
}


void sixdof_obj::buoyancy_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc, double fsf_z,
                                 double &Vsub, double &xB, double &yB, double &zB)
{
    double x0,y0,z0,x1,y1,z1,x2,y2,z2;
    double xc,yc,zc,nx,ny,nz,norm;
    double wpx[4],wpy[4],wpz[4];
    double sx[3],sy[3],sz[3];
    double mpx[3],mpy[3],mpz[3];
    double A_sub,w_q;
    double V,Mx,My,Mz;
    int nwet,q,m;

    V=Mx=My=Mz=0.0;

    for(int n=0; n<tricount; ++n)
    {
        x0=tri_x[n][0]; y0=tri_y[n][0]; z0=tri_z[n][0];
        x1=tri_x[n][1]; y1=tri_y[n][1]; z1=tri_z[n][1];
        x2=tri_x[n][2]; y2=tri_y[n][2]; z2=tri_z[n][2];

        xc=(x0+x1+x2)/3.0;
        yc=(y0+y1+y2)/3.0;
        zc=(z0+z1+z2)/3.0;

        if(xc >= p->originx && xc < p->endx &&
          (p->j_dir==0 || (yc >= p->originy && yc < p->endy)))
        {
            nx = (y1-y0)*(z2-z0) - (y2-y0)*(z1-z0);
            ny = (x2-x0)*(z1-z0) - (x1-x0)*(z2-z0);
            nz = (x1-x0)*(y2-y0) - (x2-x0)*(y1-y0);

            norm = sqrt(nx*nx+ny*ny+nz*nz);
            nz /= norm>1.0e-20?norm:1.0e20;

            nwet = clip_facet_poly(p,x0,y0,z0,x1,y1,z1,x2,y2,z2,fsf_z,wpx,wpy,wpz);

            if(nwet==0)
            continue;

            for(q=0; q<nwet-2; ++q)
            {
                sx[0]=wpx[0];   sy[0]=wpy[0];   sz[0]=wpz[0];
                sx[1]=wpx[q+1]; sy[1]=wpy[q+1]; sz[1]=wpz[q+1];
                sx[2]=wpx[q+2]; sy[2]=wpy[q+2]; sz[2]=wpz[q+2];

                A_sub = triangle_area(p,sx[0],sy[0],sz[0],sx[1],sy[1],sz[1],sx[2],sy[2],sz[2]);

                if(A_sub<1.0e-20)
                continue;

                mpx[0]=0.5*(sx[0]+sx[1]); mpy[0]=0.5*(sy[0]+sy[1]); mpz[0]=0.5*(sz[0]+sz[1]);
                mpx[1]=0.5*(sx[1]+sx[2]); mpy[1]=0.5*(sy[1]+sy[2]); mpz[1]=0.5*(sz[1]+sz[2]);
                mpx[2]=0.5*(sx[2]+sx[0]); mpy[2]=0.5*(sy[2]+sy[0]); mpz[2]=0.5*(sz[2]+sz[0]);

                w_q = A_sub/3.0;

                for(m=0; m<3; ++m)
                {
                    V  +=                (mpz[m]-fsf_z)*nz*w_q;
                    Mx += mpx[m]        *(mpz[m]-fsf_z)*nz*w_q;
                    My += mpy[m]        *(mpz[m]-fsf_z)*nz*w_q;
                    Mz += 0.5*(mpz[m]*mpz[m] - fsf_z*fsf_z)*nz*w_q;
                }
            }
        }
    }

    V  = pgc->globalsum(V);
    Mx = pgc->globalsum(Mx);
    My = pgc->globalsum(My);
    Mz = pgc->globalsum(Mz);

    Vsub = V;

    xB = fabs(V)>1.0e-20 ? Mx/V : c_(0);
    yB = fabs(V)>1.0e-20 ? My/V : c_(1);
    zB = fabs(V)>1.0e-20 ? Mz/V : c_(2);
}

// Clip a facet against the horizontal plane z = fsf_z and return the wetted
// polygon itself rather than just its area and centroid. Winding is preserved,
// so a fan from vertex 0 triangulates it without overlap.
// Returns the number of wetted vertices: 3, 4, or 0 when the facet is dry.
int sixdof_obj::clip_facet_poly(lexer *p,
                                double x0,double y0,double z0,
                                double x1,double y1,double z1,
                                double x2,double y2,double z2,
                                double fsf_z,
                                double *wpx, double *wpy, double *wpz)
{
    double xA,yA,zA,xB,yB,zB,xC,yC,zC;
    double xSb,ySb,zSb,xSc,ySc,zSc;
    double f;
    int nabove;

    nabove = 0;

    if(z0 > fsf_z)
    ++nabove;

    if(z1 > fsf_z)
    ++nabove;

    if(z2 > fsf_z)
    ++nabove;

    // fully dry
    if(nabove==3)
    return 0;

    // fully wet
    if(nabove==0)
    {
    wpx[0]=x0; wpy[0]=y0; wpz[0]=z0;
    wpx[1]=x1; wpy[1]=y1; wpz[1]=z1;
    wpx[2]=x2; wpy[2]=y2; wpz[2]=z2;

    return 3;
    }

    // Relabel so that A is the odd vertex out: the single dry one when
    // nabove==1, the single wet one when nabove==2. The cut points then
    // always lie on edges A->B and A->C. The cyclic order 0,1,2 -> 1,2,0
    // -> 2,0,1 preserves the facet winding.
    if(nabove==1)
    {
        if(z0 > fsf_z)
        {xA=x0;yA=y0;zA=z0; xB=x1;yB=y1;zB=z1; xC=x2;yC=y2;zC=z2;}

        else if(z1 > fsf_z)
        {xA=x1;yA=y1;zA=z1; xB=x2;yB=y2;zB=z2; xC=x0;yC=y0;zC=z0;}

        else
        {xA=x2;yA=y2;zA=z2; xB=x0;yB=y0;zB=z0; xC=x1;yC=y1;zC=z1;}
    }

    else
    {
        if(z0 <= fsf_z)
        {xA=x0;yA=y0;zA=z0; xB=x1;yB=y1;zB=z1; xC=x2;yC=y2;zC=z2;}

        else if(z1 <= fsf_z)
        {xA=x1;yA=y1;zA=z1; xB=x2;yB=y2;zB=z2; xC=x0;yC=y0;zC=z0;}

        else
        {xA=x2;yA=y2;zA=z2; xB=x0;yB=y0;zB=z0; xC=x1;yC=y1;zC=z1;}
    }

    // cut points on A->B and A->C
    f = clip_edge_vol(zA,zB,fsf_z);

    xSb = xA + f*(xB-xA);
    ySb = yA + f*(yB-yA);
    zSb = zA + f*(zB-zA);

    f = clip_edge_vol(zA,zC,fsf_z);

    xSc = xA + f*(xC-xA);
    ySc = yA + f*(yC-yA);
    zSc = zA + f*(zC-zA);

    // two vertices dry: the wet part is the triangle A-Sb-Sc
    if(nabove==2)
    {
    wpx[0]=xA;  wpy[0]=yA;  wpz[0]=zA;
    wpx[1]=xSb; wpy[1]=ySb; wpz[1]=zSb;
    wpx[2]=xSc; wpy[2]=ySc; wpz[2]=zSc;

    return 3;
    }

    // one vertex dry: the wet part is the quad B-C-Sc-Sb
    wpx[0]=xB;  wpy[0]=yB;  wpz[0]=zB;
    wpx[1]=xC;  wpy[1]=yC;  wpz[1]=zC;
    wpx[2]=xSc; wpy[2]=ySc; wpz[2]=zSc;
    wpx[3]=xSb; wpy[3]=ySb; wpz[3]=zSb;

    return 4;
}


double sixdof_obj::clip_edge_vol(double za, double zb, double fsf_z)
{
    double dz,f;

    dz = zb - za;

    if(fabs(dz) < 1.0e-12)
    return 0.0;

    f = (fsf_z - za)/dz;

    return MIN(MAX(f,0.0),1.0);
}
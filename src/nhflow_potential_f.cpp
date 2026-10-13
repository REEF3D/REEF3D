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

#include"nhflow_potential_f.h"
#include"solver.h"
#include"ghostcell.h"
#include"fdm_nhf.h"
#include"lexer.h"
#include"nhflow_curv.h"
#include<iomanip>
#include<cmath>

// Curvilinear grid (G 50): the horizontal Laplacian compact (face lengths over the area and the centre
// distances along the face normals, nhflow_curv::lapcoef), the inflow flux Ui along the normal of the
// inflow faces, the velocities by the chain rule (nhflow_curv::grad5), the ghost cells of the in- and
// outflow faces with the centre distance along the face normal.  Multi-block grid (G 50 1 with the
// blocks of DIVEMesh R 20 1): one solve over all blocks, nhflow_amr::potential (amr_assemble,
// amr_velocities here for every block).



nhflow_potential_f::nhflow_potential_f(lexer* p) 
{
    gcval_pot=49;
}

nhflow_potential_f::~nhflow_potential_f()
{
}

void nhflow_potential_f::start(lexer*p, fdm_nhf *d, solver* psolv, ghostcell* pgc)
{
    int itermem;

    // multi-block grid: one solve over the blocks, after they are set up (nhflow_amr::potential)
    if(p->geo_mb>0)
    {
    if(p->mpirank==0)
    cout<<"potential flow solver: on the blocks of the multi-block grid (mesh refinement ini)"<<endl<<endl;
    return;
    }
    
    p->Darray(PSI,p->imax*p->jmax*(p->kmax+2));
    p->Iarray(BC,p->imax*p->jmax*(p->kmax+3));
    
    if(p->mpirank==0 )
	cout<<"potential flow solver..."<<endl<<endl;
    
    if(fabs(p->Ui)<1.0e-9)
    {
    if(p->mpirank==0)
    cout<<"no inflow...potential flow stop"<<endl;
    goto finalize;
    }
    
    if(fabs(p->Uo)<1.0e-9)
    {
    if(p->mpirank==0)
    cout<<"no outflow...potential flow stop"<<endl;
    goto finalize;
    }
    
    ini_bc(p,d,pgc);

    starttime=pgc->timer();
	
	LOOP
	PSI[IJK] = 0.0;
    
    pgc->start49V(p,PSI,gcval_pot);
    curv_ghost(p);
	
    itermem=p->N46;
    p->N46=5000;
	
    p->final_res=0.0;

    for(int qn=0; qn<1;++qn)
    {
    laplace(p,d,pgc);
    psolv->startV(p,pgc,PSI,d->rhsvec,d->M,44);
    
    pgc->start49V(p,PSI,gcval_pot);
    curv_ghost(p);
    }

    ucalc(p,d);
	vcalc(p,d);
	wcalc(p,d);

    // The potential problem is pure Neumann apart from the in/outflow fluxes. If the wet
    // domain falls apart into regions that are not connected (e.g. a weir or dam whose crest
    // sticks out of the initial water level), each region only sees inflow or outflow, the
    // system is inconsistent and the solve diverges. Start from rest instead of passing
    // non-finite or unconverged velocities on to the first time step.
    {
    int fail = (p->solveriter>=p->N46) ? 1 : 0;

    // a solve that diverged within N 46 iterations (final residual far above the tolerance,
    // REEFMG reports it in final_res, other solvers leave it at 0) or velocities far
    // beyond the in/outflow velocities would otherwise be passed on as the initial state
    if(!std::isfinite(p->final_res) || p->final_res>1.0e-2)
    fail=1;

    const double ulim = 100.0*MAX(1.0,MAX(fabs(p->Ui),fabs(p->Uo)));

    LOOP
    if(!std::isfinite(d->U[IJK]) || !std::isfinite(d->V[IJK]) || !std::isfinite(d->W[IJK])
    || fabs(d->U[IJK])>ulim || fabs(d->V[IJK])>ulim || fabs(d->W[IJK])>ulim)
    fail=1;

    fail = pgc->globalimax(fail);

        if(fail==1)
        {
            if(p->mpirank==0)
            cout<<"potential flow solver: no converged solution (wet domain disconnected or "
                <<"in/outflow unbalanced?) - starting from rest instead"<<endl<<endl;

            LOOP
            {
            d->U[IJK]=d->V[IJK]=d->W[IJK]=0.0;
            d->UH[IJK]=d->VH[IJK]=d->WH[IJK]=0.0;
            }
        }
    }

	pgc->start4V(p,d->U,10);
    pgc->start4V(p,d->V,11);
    pgc->start4V(p,d->W,12);
    
    pgc->start4V(p,d->UH,14);
    pgc->start4V(p,d->VH,15);
    pgc->start4V(p,d->WH,16);
    
    pgc->start4V(p,d->U,10);
    pgc->start4V(p,d->V,11);
    pgc->start4V(p,d->W,12);
    
    pgc->start4V(p,d->UH,14);
    pgc->start4V(p,d->VH,15);
    pgc->start4V(p,d->WH,16);

    endtime=pgc->timer();
    p->laplaceiter=p->solveriter;
	p->laplacetime=endtime-starttime;
	if(p->mpirank==0  && (p->count%p->P12==0))
	cout<<"lapltime: "<<p->laplacetime<<"  lapiter: "<<p->laplaceiter<<endl<<endl;

    p->N46=itermem;

    finalize:
    
    p->del_Iarray(BC,p->imax*p->jmax*(p->kmax+3));
    p->del_Darray(PSI,p->imax*p->jmax*(p->kmax+2));
}

void nhflow_potential_f::ucalc(lexer *p, fdm_nhf *d)
{	
    // curvilinear grid: both horizontal components by the chain rule (vcalc does nothing)
    if(p->geo_curv==1)
    {
    curv_velocities(p,d);
    return;
    }
    
	LOOP
    {
	d->U[IJK] =  (PSI[Ip1JK]-PSI[Im1JK])/(p->DXP[IP]+p->DXP[IM1]) 
    
              + 0.5*(p->sigx[FIJK]+p->sigx[FIJKp1])*((PSI[IJKp1]-PSI[IJKm1])/(p->DZP[KP]+p->DZP[KM1]));
              
    d->UH[IJK] = d->WL(i,j)*d->U[IJK];
    }
    
    
    LOOP
    if(p->DF[IJK]<0 || p->DF[Im1JK]<0 || p->DF[Ip1JK]<0 
    ||p->wet[IJ]==0||p->wet[Im1J]==0||p->wet[Ip1J]==0)
    {
	d->U[IJK] = 0.0;
    d->UH[IJK] = 0.0;
    }
}

void nhflow_potential_f::vcalc(lexer *p, fdm_nhf *d)
{	
    if(p->geo_curv==1)
    return;
    
    LOOP
    {
	d->V[IJK] =  (PSI[IJp1K]-PSI[IJm1K])/(p->DYP[JP]+p->DYP[JP1]) 
    
              + 0.5*(p->sigy[FIJK]+p->sigy[FIJKp1])*((PSI[IJKp1]-PSI[IJKm1])/(p->DZP[KP]+p->DZP[KM1]));
              
    d->VH[IJK] = d->WL(i,j)*d->V[IJK];
    }


    LOOP
    if(p->DF[IJK]<0 || p->DF[IJm1K]<0 || p->DF[IJp1K]<0 
    ||p->wet[IJ]==0||p->wet[IJm1]==0||p->wet[IJp1]==0)
    {
	d->V[IJK] = 0.0;
    d->VH[IJK] = 0.0;
    }
}

void nhflow_potential_f::wcalc(lexer *p, fdm_nhf *d)
{
    LOOP
    {
	d->W[IJK] =  p->sigz[IJ]*(PSI[IJKp1]-PSI[IJKm1])/(p->DZP[KP]+p->DZP[KP1]); 
              
    d->WH[IJK] = d->WL(i,j)*d->W[IJK];
    }
    
    
	LOOP
    if(p->DF[IJK]<0 || p->DF[IJm1K]<0 || p->DF[IJp1K]<0||p->wet[IJ]==0)
    {
	d->W[IJK] = 0.0;
    d->WH[IJK] = 0.0;
    }
}

void nhflow_potential_f::rhs(lexer *p, fdm_nhf *d)
{
    count=0;
    LOOP
    {
    d->rhsvec.V[count] = 0.0;
    count++;
    }
}

void nhflow_potential_f::laplace(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    n=0;
    BASELOOP
    {
    d->M.p[n]  =  1.0;

        d->M.n[n] = 0.0;
        d->M.s[n] = 0.0;

        d->M.w[n] = 0.0;
        d->M.e[n] = 0.0;
        
        d->M.t[n] = 0.0;
        d->M.b[n] = 0.0;
        
        d->rhsvec.V[n] =  0.0;
        
    ++n;
    }

    
    n=0;
    LOOP
	{
        if(p->wet[IJ]==1 && p->DF[IJK]>0 && p->geo_curv==1)
        {
            sigxyz2 = pow(p->sigx[FIJK],2.0) + pow(p->sigy[FIJK],2.0) + pow(p->sigz[IJ],2.0);
            
            double cIp,cIm,cJp,cJm;
            nhflow_curv::lapcoef(p,i,j,cIp,cIm,cJp,cJm);
            
            d->M.p[n]  =  cIp + cIm + cJp + cJm
                        + sigxyz2/(p->DZP[KM1]*p->DZN[KP])
                        + sigxyz2/(p->DZP[KM1]*p->DZN[KM1]);
            
            d->M.n[n] = -cIp;
            d->M.s[n] = -cIm;
            d->M.w[n] = -cJp;
            d->M.e[n] = -cJm;
            
            d->M.t[n] = -sigxyz2/(p->DZP[KM1]*p->DZN[KP])     
                        - p->sigxx[FIJK]/((p->DZN[KP]+p->DZN[KM1]));
                        
            d->M.b[n] = -sigxyz2/(p->DZP[KM1]*p->DZN[KM1]) 
                        + p->sigxx[FIJK]/((p->DZN[KP]+p->DZN[KM1]));
            
            // the cross terms 2 (sigx d/dx + sigy d/dy) dPSI/dsigma with the current PSI
            const double ds = p->DZN[KP]+p->DZN[KM1];
            auto q = [&](int a, int b)
            {
                const int c = (a-p->imin)*p->jmax*p->kmax + (b-p->jmin)*p->kmax + k-p->kmin;
                return (PSI[c+1] - PSI[c-1])/ds;
            };
            double qx,qy;
            nhflow_curv::grad5(p,i,j,q(i+1,j),q(i-1,j),q(i,j+1),q(i,j-1),qx,qy);
            d->rhsvec.V[n] = (p->sigx[FIJK]+p->sigx[FIJKp1])*qx + (p->sigy[FIJK]+p->sigy[FIJKp1])*qy;
        }
        
        if(p->wet[IJ]==1 && p->DF[IJK]>0 && p->geo_curv==0)
        {
            sigxyz2 = pow(p->sigx[FIJK],2.0) + pow(p->sigy[FIJK],2.0) + pow(p->sigz[IJ],2.0);
            
            d->M.p[n]  =  1.0/(p->DXP[IP]*p->DXN[IP])
                        + 1.0/(p->DXP[IM1]*p->DXN[IP])
                        
                        + 1.0/(p->DYP[JP]*p->DYN[JP])*p->y_dir
                        + 1.0/(p->DYP[JM1]*p->DYN[JP])*p->y_dir
                        
                        + sigxyz2/(p->DZP[KM1]*p->DZN[KP])
                        + sigxyz2/(p->DZP[KM1]*p->DZN[KM1]);


            d->M.n[n] = -1.0/(p->DXP[IP]*p->DXN[IP]);
            d->M.s[n] = -1.0/(p->DXP[IM1]*p->DXN[IP]);

            d->M.w[n] = -1.0/(p->DYP[JP]*p->DYN[JP])*p->y_dir;
            d->M.e[n] = -1.0/(p->DYP[JM1]*p->DYN[JP])*p->y_dir;

            d->M.t[n] = -sigxyz2/(p->DZP[KM1]*p->DZN[KP])     
                        - p->sigxx[FIJK]/((p->DZN[KP]+p->DZN[KM1]));
                        
            d->M.b[n] = -sigxyz2/(p->DZP[KM1]*p->DZN[KM1]) 
                        + p->sigxx[FIJK]/((p->DZN[KP]+p->DZN[KM1]));
            
            
            d->rhsvec.V[n] =  (p->sigx[FIJK]+p->sigx[FIJKp1])*(PSI[Ip1JKp1] - PSI[Im1JKp1] - PSI[Ip1JKm1] + PSI[Im1JKm1])
                            /((p->DXP[IP]+p->DXP[IM1])*(p->DZN[KP]+p->DZN[KM1]))
                        
                            + (p->sigx[FIJK]+p->sigx[FIJKp1])*(PSI[IJp1Kp1] - PSI[IJm1Kp1] - PSI[IJp1Km1] + PSI[IJm1Km1])
                            /((p->DYP[JP]+p->DYP[JM1])*(p->DZN[KP]+p->DZN[KM1]))*p->y_dir;

        }
        
        if(p->wet[IJ]==0 || p->DF[IJK]<0)
        {
        d->M.p[n]  =  1.0;


        d->M.n[n] = 0.0;
        d->M.s[n] = 0.0;

        d->M.w[n] = 0.0;
        d->M.e[n] = 0.0;

        d->M.t[n] = 0.0;
        d->M.b[n] = 0.0;
        
        d->rhsvec.V[n] =  0.0;
        }
        
        
	++n;
	}
    
    
    double denom,ab;  
    
    n=0;
	LOOP
	{
        if(p->DF[IJK]>0 && p->flag4[IJK]>0 && p->wet[IJ]==1)
        {
           
		if((p->flag4[Im1JK]<0 || p->DF[Im1JK]<0 ||  p->wet[Im1J]==0) && BC[Im1JK]==0)
		{
		d->M.p[n] += d->M.s[n];
		d->M.s[n] = 0.0;
		}
        
        if((p->flag4[Im1JK]<0 || p->DF[Im1JK]<0) && BC[Im1JK]==1)
        {
        // curvilinear grid: the flux Ui through the face (its length over the area of the cell)
        if(p->geo_curv==1)
        d->rhsvec.V[n] += -p->Ui*p->LEN1[Im1J]/p->AREA2D[IJ];
        else
        d->rhsvec.V[n] += d->M.s[n]*p->Ui*p->DXP[IM1];
        d->M.p[n] += d->M.s[n];
        d->M.s[n] = 0.0;
        }
		
		if((p->flag4[Ip1JK]<0 || p->DF[Ip1JK]<0 || p->wet[Ip1J]==0) && BC[Ip1JK]==0)
		{
		d->M.p[n] += d->M.n[n];
		d->M.n[n] = 0.0;
		}
        
        // outflow: psi = 0 on the outflow face (ghost = -psi), not the flux Uo.  With the
        // prescribed flux the problem was pure Neumann, and the sigma operator (bed slope
        // terms sigxx, KBEDBC) is not in conservation form, so inflow = outflow did not make
        // it consistent over real bathymetry: the solve stalled or diverged.  The Dirichlet
        // outflow removes the null space; the outflow profile follows from the potential.
        if((p->flag4[Ip1JK]<0 || p->DF[Ip1JK]<0) && BC[Ip1JK]==2)
        {
        d->M.p[n] -= d->M.n[n];
        d->M.n[n] = 0.0;
        }
		
		if(p->flag4[IJm1K]<0 || p->DF[IJm1K]<0 || p->wet[IJm1]==0)
		{
		d->M.p[n] += d->M.e[n];
		d->M.e[n] = 0.0;
		}
		
		if(p->flag4[IJp1K]<0 || p->DF[IJp1K]<0 || p->wet[IJp1]==0)
		{
		d->M.p[n] += d->M.w[n];
		d->M.w[n] = 0.0;
		}
		
		if(p->flag4[IJKm1]>0 && p->DF[IJKm1]<0)
		{
		d->M.p[n] += d->M.b[n];
		d->M.b[n] = 0.0;
		}
		
		if(p->flag4[IJKp1]<0 || p->DF[IJKp1]<0)
		{
		d->M.p[n] += d->M.t[n];
		d->M.t[n] = 0.0;
		}
        
        
         // KBEDBC
            if(p->flag4[IJKm1]<0)
            {
            sigxyz2 = pow(p->sigx[FIJK],2.0) + pow(p->sigy[FIJK],2.0) + pow(p->sigz[IJ],2.0);
            
            ab = -(sigxyz2/(p->DZP[KM1]*p->DZN[KM1]) - p->sigxx[FIJK]/(p->DZN[KP]+p->DZN[KM1]));
            
            denom = p->sigz[IJ] + d->Bx(i,j)*p->sigx[FIJK] + d->By(i,j)*p->sigy[FIJK];

                // curvilinear grid: dPSI/dx, dPSI/dy at the bed by the chain rule (sflow_curv::grad5)
                if(p->geo_curv==1)
                {
                    if(p->wet[Ip1J]==1 && p->wet[Im1J]==1 && p->wet[IJp1]==1 && p->wet[IJm1]==1
                    && p->flag4[Ip1JK]>0 && p->flag4[Im1JK]>0 && p->flag4[IJp1K]>0 && p->flag4[IJm1K]>0)
                    {
                    auto xc = [&](int a, int b) { return p->XC2D[(a-p->imin)*p->jmax + (b-p->jmin)]; };
                    auto yc = [&](int a, int b) { return p->YC2D[(a-p->imin)*p->jmax + (b-p->jmin)]; };
                    const double xe = xc(i+1,j) - xc(i-1,j), ye = yc(i+1,j) - yc(i-1,j);
                    const double xn = xc(i,j+1) - xc(i,j-1), yn = yc(i,j+1) - yc(i,j-1);
                    const double det = xe*yn - ye*xn;
                    if(fabs(det)>0.0)
                    {
                    const double c = ab*2.0*p->DZN[KP]/(denom*det);
                    const double bx = d->Bx(i,j), by = d->By(i,j);
                    d->M.n[n] += c*( bx*yn - by*xn);
                    d->M.s[n] += c*(-bx*yn + by*xn);
                    d->M.w[n] += c*(-bx*ye + by*xe);
                    d->M.e[n] += c*( bx*ye - by*xe);
                    }
                    }
                }
                else
                {
                    if(p->wet[Ip1J]==1 && p->wet[Im1J]==1)
                    {
                    d->M.n[n] +=  ab*2.0*p->DZN[KP]*(d->Bx(i,j))/(denom*(p->DXP[IP] + p->DXP[IM1]));
                    d->M.s[n] += -ab*2.0*p->DZN[KP]*(d->Bx(i,j))/(denom*(p->DXP[IP] + p->DXP[IM1]));
                    }
                    
                    if(p->wet[IJp1]==1 && p->wet[IJm1]==1)
                    {
                    d->M.w[n] +=  ab*2.0*p->DZN[KP]*d->By(i,j)/(denom*(p->DYP[JP] + p->DYP[JM1]));
                    d->M.e[n] += -ab*2.0*p->DZN[KP]*d->By(i,j)/(denom*(p->DYP[JP] + p->DYP[JM1]));
                    }
                }

                d->M.t[n] += ab;
                d->M.b[n] = 0.0;
            }
        }

	++n;
	}
    
    
    
    /*
    count=0;
    if(p->mpirank==2)
	for(int q=0;q<p->gcpara4_count;++q)
	{
    i=p->gcpara4[q][0];
    j=p->gcpara4[q][1];
    k=p->gcpara4[q][2];
        
        cout<<"DF2: "<<p->DF[Ip1JK]<<" "<<p->DF[Ip2JK]<<endl;
	}
    
    
    count=0;
    if(p->mpirank==3)
	for(int q=0;q<p->gcpara1_count;++q)
	{
    i=p->gcpara1[q][0];
    j=p->gcpara1[q][1];
    k=p->gcpara1[q][2];
        
        cout<<"DF3: "<<p->DF[IJK]<<" "<<p->DF[Ip1JK]<<endl;
	}*/
}

void nhflow_potential_f::ini_bc(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    BASELOOP
    BC[IJK]=0;
    
    LOOP
    {
        if(p->flag4[Im1JK]<0)
		BC[Im1JK]=0;
		
		if(p->flag4[Ip1JK]<0)
		BC[Ip1JK]=0;
		
		if(p->flag4[IJm1K]<0)
		BC[IJm1K]=0;
		
		if(p->flag4[IJp1K]<0)
		BC[IJp1K]=0;
		
		if(p->flag4[IJKm1]<0)
		BC[IJKm1]=0;
		
		if(p->flag4[IJKp1]<0)
		BC[IJKp1]=0;
    }
    
    // -----
    
    GC4LOOP
    {
        if(p->gcb4[n][4]==1 || p->gcb4[n][4]==6)
        {
            i=p->gcb4[n][0]; 
            j=p->gcb4[n][1];
            k=p->gcb4[n][2];  
            
            WETDRY
            {
            if(p->gcb4[n][3]==1)
            BC[Im1JK]=1;
            
            if(p->gcb4[n][3]==3)
            BC[IJm1K]=1;
            
            if(p->gcb4[n][3]==2)
            BC[IJp1K]=1;
            
            if(p->gcb4[n][3]==4)
            BC[Ip1JK]=1;
            }
 
        }
        
        if(p->gcb4[n][4]==2 || p->gcb4[n][4]==7 || p->gcb4[n][4]==8)
        {
            i=p->gcb4[n][0]; 
            j=p->gcb4[n][1];
            k=p->gcb4[n][2];  
            
            WETDRY
            {
            if(p->gcb4[n][3]==1)
            BC[Im1JK]=2;
            
            if(p->gcb4[n][3]==3)
            BC[IJm1K]=2;
            
            if(p->gcb4[n][3]==2)
            BC[IJp1K]=2;
            
            if(p->gcb4[n][3]==4)
            BC[Ip1JK]=2;
            }
 
        }
    }
}

// curvilinear grid: the ghost cells of the in- (BC 1) and outflow faces (BC 2) of start49V with the
// centre distance along the face normal; the other i-walls zero gradient as the j-walls
void nhflow_potential_f::curv_ghost(lexer *p)
{
    if(p->geo_curv!=1)
    return;

    auto xc = [&](int a, int b) { return p->XC2D[(a-p->imin)*p->jmax + (b-p->jmin)]; };
    auto yc = [&](int a, int b) { return p->YC2D[(a-p->imin)*p->jmax + (b-p->jmin)]; };
    auto dn = [&](int a, int b)     // cell (a,b) to (a+1,b) along the normal of face a+1/2
    {
        double nx,ny,len;
        nhflow_curv::normal(p,1,a,b,nx,ny,len);
        return (xc(a+1,b)-xc(a,b))*nx + (yc(a+1,b)-yc(a,b))*ny;
    };

    LOOP
    {
        if(p->flag4[Im1JK]<0)
        PSI[Im1JK] = (BC[Im1JK]==1) ? PSI[IJK] - p->Ui*dn(i-1,j) : PSI[IJK];

        if(p->flag4[Ip1JK]<0)
        PSI[Ip1JK] = (BC[Ip1JK]==2) ? PSI[IJK] + p->Uo*dn(i,j) : PSI[IJK];
    }
}

// curvilinear grid: U, V from the gradient along the sigma layer (chain rule) and the sigma terms
void nhflow_potential_f::curv_velocities(lexer *p, fdm_nhf *d)
{
    LOOP
    {
        double px,py;
        nhflow_curv::grad5(p,i,j,PSI[Ip1JK],PSI[Im1JK],PSI[IJp1K],PSI[IJm1K],px,py);
        const double dz = (PSI[IJKp1]-PSI[IJKm1])/(p->DZP[KP]+p->DZP[KM1]);
        
        d->U[IJK] = px + 0.5*(p->sigx[FIJK]+p->sigx[FIJKp1])*dz;
        d->V[IJK] = py + 0.5*(p->sigy[FIJK]+p->sigy[FIJKp1])*dz;
        d->UH[IJK] = d->WL(i,j)*d->U[IJK];
        d->VH[IJK] = d->WL(i,j)*d->V[IJK];
    }
    
    LOOP
    if(p->DF[IJK]<0 || p->DF[Im1JK]<0 || p->DF[Ip1JK]<0 || p->DF[IJm1K]<0 || p->DF[IJp1K]<0
    || p->wet[IJ]==0 || p->wet[Im1J]==0 || p->wet[Ip1J]==0 || p->wet[IJm1]==0 || p->wet[IJp1]==0)
    {
    d->U[IJK] = d->V[IJK] = 0.0;
    d->UH[IJK] = d->VH[IJK] = 0.0;
    }
}

// multi-block grid (nhflow_amr::potential): the rows of this grid (its fdm M and rhsvec, the LOOP order
// of the Poisson rows), the potential 0 in the F-layout array of the composite solve
void nhflow_potential_f::amr_assemble(lexer *p, fdm_nhf *d, ghostcell *pgc, double *PSIF)
{
    p->Darray(PSI,p->imax*p->jmax*(p->kmax+2));
    p->Iarray(BC,p->imax*p->jmax*(p->kmax+3));
    
    ini_bc(p,d,pgc);
    
    LOOP
    PSI[IJK] = 0.0;
    
    // the ghost cells as start() has them for the rows (the cross terms)
    pgc->start49V(p,PSI,gcval_pot);
    curv_ghost(p);
    
    laplace(p,d,pgc);
    
    for(int q=0; q<p->imax*p->jmax*(p->kmax+2); ++q)
    PSIF[q] = 0.0;
    
    p->del_Iarray(BC,p->imax*p->jmax*(p->kmax+3));
    p->del_Darray(PSI,p->imax*p->jmax*(p->kmax+2));
}

// multi-block grid: the velocities of this grid from the solution (F layout, the cells around the grid
// filled), ghost cells as start(); returns 1 where the solution is not usable
int nhflow_potential_f::amr_velocities(lexer *p, fdm_nhf *d, ghostcell *pgc, const double *PSIF)
{
    p->Darray(PSI,p->imax*p->jmax*(p->kmax+2));
    p->Iarray(BC,p->imax*p->jmax*(p->kmax+3));
    
    ini_bc(p,d,pgc);
    
    for(i=p->imin; i<p->imin+p->imax; ++i)
    for(j=p->jmin; j<p->jmin+p->jmax; ++j)
    for(k=p->kmin; k<p->kmin+p->kmax; ++k)
    PSI[IJK] = PSIF[FIJK];
    
    pgc->start49V(p,PSI,gcval_pot);
    curv_ghost(p);
    
    ucalc(p,d);
	vcalc(p,d);
	wcalc(p,d);
    
    int fail=0;
    const double ulim = 100.0*MAX(1.0,MAX(fabs(p->Ui),fabs(p->Uo)));
    LOOP
    if(!std::isfinite(d->U[IJK]) || !std::isfinite(d->V[IJK]) || !std::isfinite(d->W[IJK])
    || fabs(d->U[IJK])>ulim || fabs(d->V[IJK])>ulim || fabs(d->W[IJK])>ulim)
    fail=1;
    
    p->del_Iarray(BC,p->imax*p->jmax*(p->kmax+3));
    p->del_Darray(PSI,p->imax*p->jmax*(p->kmax+2));
    
    return fail;
}

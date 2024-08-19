/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTIBILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
Author: Alexander Hanke
--------------------------------------------------------------------*/

#include "sedpart_movement.h"
#include "particles_obj.h"

#include "lexer.h"
#include "fdm.h"
#include "ghostcell.h"
#include "field4a.h"
#include "sediment_fdm.h"
#include "turbulence.h"

namespace sediment_particle::movement
{
    Tavouktsoglou::Tavouktsoglou(lexer *p) : drho(p->W1/p->S22) ,invKinVis(p->W1/p->W2), 
                                            Ps(p->Q14),beta(p->Q15),epsilon(p->Q16),theta_crit(p->Q17)
    {
        p->Darray(stressTensor,p->imax*p->jmax*p->kmax);
        p->Darray(cellSum,p->imax*p->jmax*p->kmax);
        p->Darray(cellSumTopo,p->imax*p->jmax*p->kmax);
        p->Darray(columnSum,p->imax*p->jmax);
    }

    Tavouktsoglou::~Tavouktsoglou()
    {
        delete[] stressTensor;
        stressTensor = nullptr;
        delete[] cellSum;
        cellSum = nullptr;
        delete[] cellSumTopo;
        cellSumTopo = nullptr;
        delete[] columnSum;
        columnSum = nullptr;
    }

    /**
     * @brief Sets up the Tavouktsoglou class.
     *
     * This function is responsible for setting up the Tavouktsoglou class by calculating
     * the cellSumTopo and columnSum values based on the given parameters.
     *
     * @param p A pointer to the lexer object.
     * @param a A reference to the fdm object.
     * @param diameter The diameter value.
     */
    void Tavouktsoglou::setup(lexer *p, fdm &a, double &diameter)
    {
        PLAINLOOP
        {
            cellSumTopo[IJK] = maxParticlesPerCell(p,a,diameter);
            columnSum[IJ] += cellSumTopo[IJK];
        }
    }

    void Tavouktsoglou::setupState(lexer *p, fdm &a, ghostcell &pgc, particles_obj &PP)
    {
        particlePerCell(p,pgc,PP);
        PLAINLOOP
        columnSum[IJ] += cellSumTopo[IJK]+cellSum[IJK];
        particleStressTensor(p,a,pgc,PP);
    }

    /**
     * @brief Determines if a particle should be seeded in the current cell.
     *
     * This function is responsible for determining if a particle should be seeded in the current cell.
     * It does so by checking if the cellSumTopo value is greater than 0 and if the cellSum value is less than the maximum value.
     *
     * @param p A pointer to the lexer object.
     * @param PP A reference to the particles_obj object.
     * @param index The index of the particle.
     * @param max The maximum value.
     * @return True if the particle should be seeded, false otherwise.
     */
    seedReturn Tavouktsoglou::seeding(lexer *p, particles_obj &PP, size_t &index, double max, bool free)
    {
        if(free)
        {
            cellSum[IJK] += PP.PackingFactor[index];
        }
        else
        {
            if(cellSumTopo[IJK]>=PP.PackingFactor[index])
                cellSumTopo[IJK] -= PP.PackingFactor[index];
            else if (cellSumTopo[IJK]>0)
            {
                PP.PackingFactor[index] = cellSumTopo[IJK];
                cellSumTopo[IJK] -= PP.PackingFactor[index];
            }
            else
            {
                return seedReturn::REMOVE;
            }
            cellSum[IJK] += PP.PackingFactor[index];
            
        }

        if(cellSum[IJK]>=max)
            return seedReturn::STOP;
        return seedReturn::CONTINUE;
    }

    /**
     * @brief Transfers the particle to the current cell.
     *
     * This function is responsible for transferring the particle to the current cell.
     * It does so by adding the PackingFactor value to the cellSum value.
     *
     * @param p A pointer to the lexer object.
     * @param PP A reference to the particles_obj object.
     * @param index The index of the particle.
     */
    void Tavouktsoglou::transfer(lexer *p, particles_obj &PP, size_t &index)
    {
        cellSum[IJK] += PP.PackingFactor[index];
    }

    /**
     * @brief Removes the particle from the current cell.
     *
     * This function is responsible for removing the particle from the current cell.
     * It does so by subtracting the PackingFactor value from the cellSum value.
     *
     * @param p A pointer to the lexer object.
     * @param PP A reference to the particles_obj object.
     * @param index The index of the particle.
     */
    void Tavouktsoglou::remove(lexer *p, particles_obj &PP, size_t &index)
    {
        cellSum[IJK] -= PP.PackingFactor[index];
    }

    /**
     * @brief Moves the particles with the flow field.
     *
     * This function is responsible for moving the particle with the flow field.
     * It does so by calculating the experienced accelerations and the resulting new velocity and position of the particle.
     *
     * @param p A pointer to the lexer object.
     * @param a A reference to the fdm object.
     * @param pgc A reference to the ghostcell object.
     * @param PP A reference to the particles_obj object.
     */
    void Tavouktsoglou::move(lexer *p, fdm &a, ghostcell &pgc, particles_obj &PP, sediment_fdm &s, turbulence &pturb)
    {
        double RKu,RKv,RKw;
        double u,v,w;
        double du1, du2, du3, dv1, dv2, dv3, dw1, dw2, dw3;
        /// @brief Difference between flowfield and particle velocity
        double du, dv, dw;
        double Dp, thetas;
        double pressureDivX=0, pressureDivY=0, pressureDivZ=0;
        double stressDivX=0, stressDivY=0, stressDivZ=0;
        double netBuoyX=(1.0-drho)*p->W20, netBuoyY=(1.0-drho)*p->W21, netBuoyZ=(1.0-drho)*p->W22;

        bool limited = false;
        bool debugPrint = false;
        bool bedLoad = false;
        bool shearVel = true;

        particlePerCell(p,pgc,PP);
        particleStressTensor(p,a,pgc,PP);

        for(size_t n=0;n<PP.loopindex;n++)
        {
            if(PP.Flag[n]>0) // INT32_MIN
            {
                // Prep
                if(p->global_xmin+p->Q73>PP.X[n])
                limited = true;
                else
                limited = false;

                i=p->posc_i(PP.X[n]);
                j=p->posc_j(PP.Y[n]);
                k=p->posc_k(PP.Z[n]);

                thetas=theta_s(p,a,PP,i,j,k);

                // Non interpolation leads to blockyness
                // stressDivX = (stressTensor[Ip1JK] - stressTensor[IJK])/(p->DXN[IP]);
                // stressDivX = (stressTensor[IJK] - stressTensor[Im1JK])/(p->DXN[IM1]);
                // stressDivY = (0.5*(stressTensor[IJp1K]+stressTensor[Ip1Jp1K]) - 0.5*(stressTensor[IJm1K]+stressTensor[Ip1Jm1K]))/(p->DYN[JM1]+p->DYN[JP]);
                // stressDivZ = (0.5*(stressTensor[IJKp1]+stressTensor[Ip1JKp1]) - 0.5*(stressTensor[IJKm1]+stressTensor[Ip1JKm1]))/(p->DYN[KM1]+p->DYN[KP]);

                stressDivX = (p->ccipol4c(stressTensor,PP.X[n]+0.5*p->DXN[IP],PP.Y[n],PP.Z[n]) - p->ccipol4c(stressTensor,PP.X[n]-0.5*p->DXN[IP],PP.Y[n],PP.Z[n]))/p->DXN[IP];
                stressDivY = (p->ccipol4c(stressTensor,PP.X[n],PP.Y[n]+0.5*p->DYN[JP],PP.Z[n]) - p->ccipol4c(stressTensor,PP.X[n],PP.Y[n]-0.5*p->DYN[JP],PP.Z[n]))/p->DYN[JP];
                stressDivZ = (p->ccipol4c(stressTensor,PP.X[n],PP.Y[n],PP.Z[n]+0.5*p->DZN[KP]) - p->ccipol4c(stressTensor,PP.X[n],PP.Y[n],PP.Z[n]-0.5*p->DZN[KP]))/p->DXN[KP];

                // stressDivX = (stressTensor[Ip1JK] - stressTensor[Im1JK])/(p->DXN[IP]+p->DXN[IM1]);
                // stressDivY = (stressTensor[IJp1K] - stressTensor[IJm1K])/(p->DYN[JP]+p->DYN[JM1]);
                // stressDivZ = (stressTensor[IJKp1] - stressTensor[IJKm1])/(p->DZN[KP]+p->DZN[KM1]);

                // if(u<0)
                    // stressDivX = (stressTensor[Ip1JK] - p->ccipol4c(stressTensor,PP.X[n],PP.Y[n],PP.Z[n]))/(0.5*p->DXN[IP1]+p->XN[IP1]-PP.X[n]);
                // else
                    // stressDivX = (p->ccipol4c(stressTensor,PP.X[n],PP.Y[n],PP.Z[n])-stressTensor[Im1JK])/(PP.X[n]-(0.5*p->DXN[IM1]+p->XN[IM1]));
                // if(v<0)
                //     stressDivY = (stressTensor[IJp1K] - p->ccipol4c(stressTensor,PP.X[n],PP.Y[n],PP.Z[n]))/(0.5*p->DYN[IP1]+p->YN[IP1]-PP.Y[n]);
                // else
                //     stressDivY = (p->ccipol4c(stressTensor,PP.X[n],PP.Y[n],PP.Z[n])-stressTensor[IJm1K])/(PP.Y[n]-(0.5*p->DYN[IM1]+p->YN[IM1]));
                // if(w<0)
                //     stressDivZ = (stressTensor[IJKp1] - p->ccipol4c(stressTensor,PP.X[n],PP.Y[n],PP.Z[n]))/(0.5*p->DZN[IP1]+p->ZN[IP1]-PP.Z[n]);
                // else
                //     stressDivZ = (p->ccipol4c(stressTensor,PP.X[n],PP.Y[n],PP.Z[n])-stressTensor[IJm1K])/(PP.Z[n]-(0.5*p->DZN[IM1]+p->ZN[IM1]));

                // pressureDivX = (a.press(i+1,j,k) - a.press(i,j,k))/(p->DXN[IP]);
                pressureDivX = (p->ccipol4_c(a.press,PP.X[n]+0.5*p->DXN[IP],PP.Y[n],PP.Z[n]) - p->ccipol4_c(a.press,PP.X[n]-0.5*p->DXN[IP],PP.Y[n],PP.Z[n]))/p->DXN[IP];
                pressureDivY = (0.5*(a.press(i,j+1,k)+a.press(i+1,j+1,k)) - 0.5*(a.press(i,j-1,k)+a.press(i+1,j-1,k)))/(p->DYN[JM1]+p->DYN[JP]);
                pressureDivZ = (0.5*(a.press(i,j,k+1)+a.press(i+1,j,k+1)) - 0.5*(a.press(i,j,k-1)+a.press(i+1,j,k-1)))/(p->DYN[KM1]+p->DYN[KP]);

                if(p->ccipol4(a.topo,PP.X[n],PP.Y[n],PP.Z[n])<PP.d50*10)
                    bedLoad=true;

                if(!bedLoad)
                {
                    u=p->ccipol1c(a.u,PP.X[n],PP.Y[n],PP.Z[n]);
                    v=p->ccipol2c(a.v,PP.X[n],PP.Y[n],PP.Z[n]);
                    w=p->ccipol3c(a.w,PP.X[n],PP.Y[n],PP.Z[n]);
                }
                else
                {
                    double uvel,vvel,u_abs;
                    double signx,signy;
                    
                    double ux1,vx1,ux2,vx2,uy1,vy1,uy2,vy2;
                    double sgx1,sgx2,sgy1,sgy2;
                    double ux1_abs,ux2_abs,uy1_abs,uy2_abs;
                    
                    uvel=p->ccslipol4(s.P,PP.X[n],PP.Y[n]);
                    vvel=p->ccslipol4(s.Q,PP.X[n],PP.Y[n]);
                    
                    u_abs = sqrt(uvel*uvel + vvel*vvel);
                    signx=fabs(u_abs)>1.0e-10?uvel/fabs(u_abs):0.0;
                    signy=fabs(u_abs)>1.0e-10?vvel/fabs(u_abs):0.0;
                    
                    if(shearVel)
                    {
                        // u=sqrt(tau/rho)
                        u = sqrt(fabs(p->ccslipol4(s.tau_eff,PP.X[n],PP.Y[n]))/p->W1)*(uvel>=0?1:-1)*signx;
                        v = sqrt(fabs(p->ccslipol4(s.tau_eff,PP.X[n],PP.Y[n]))/p->W1)*(vvel>=0?1:-1)*signy;
                        w = 0.0;
                    }
                    else
                    {
                        // alternative
                        du1 = du2 = du3 = p->ccslipol4(s.tau_eff,PP.X[n],PP.Y[n])*3/(PP.d50*PP.PackingFactor[n]*2*PP.density)*signx;
                        dv1 = dv2 = dv3 = p->ccslipol4(s.tau_eff,PP.X[n],PP.Y[n])*3/(PP.d50*PP.PackingFactor[n]*2*PP.density)*signy;
                        dw1 = dw2 = dw3 = 0.0;
                    }
                }
                
                // RK3 step 1
                if(shearVel)
                {
                    du=u-PP.U[n];
                    dv=v-PP.V[n];
                    dw=w-PP.W[n];

                    Dp=drag_model(p,PP.d50*PP.PackingFactor[n],du,dv,dw,thetas);

                    du1=Dp*du;
                    dv1=Dp*dv;
                    dw1=Dp*dw;
                }

                du1+=netBuoyX-pressureDivX/p->S22-stressDivX/(thetas*p->S22)+p->ccipol1c(a.fbh1,PP.X[n],PP.Y[n],PP.Z[n])*(0.0-PP.U[n])/p->dt;
                dv1+=netBuoyY-pressureDivY/p->S22-stressDivY/(thetas*p->S22)+p->ccipol2c(a.fbh2,PP.X[n],PP.Y[n],PP.Z[n])*(0.0-PP.V[n])/p->dt;
                dw1+=netBuoyZ-pressureDivZ/p->S22-stressDivZ/(thetas*p->S22)+p->ccipol3c(a.fbh3,PP.X[n],PP.Y[n],PP.Z[n])*(0.0-PP.W[n])/p->dt;

                if(debugPrint)
                {
                    cout<<"Z-dir1:drag: "<<Dp*dw<<" buoy: "<<netBuoyZ<<" press: "<<-pressureDivZ/p->S22<<" stress: "<<-stressDivZ/((1-thetas)*p->S22)<<"\n";
                }

                RKu=PP.U[n]+du1*p->dt;
                RKv=PP.V[n]+dv1*p->dt;
                RKw=PP.W[n]+dw1*p->dt;
                
                // RK step 2
                if(shearVel)
                {
                    du=u-RKu;
                    dv=v-RKv;
                    dw=w-RKw;

                    Dp=drag_model(p,PP.d50*PP.PackingFactor[n],du,dv,dw,thetas);

                    du2=Dp*du;
                    dv2=Dp*dv;
                    dw2=Dp*dw;
                }

                du2+=netBuoyX-pressureDivX/p->S22-stressDivX/(thetas*p->S22)+p->ccipol1c(a.fbh1,PP.X[n],PP.Y[n],PP.Z[n])*(0.0-RKu)/p->dt;
                dv2+=netBuoyY-pressureDivY/p->S22-stressDivY/(thetas*p->S22)+p->ccipol2c(a.fbh2,PP.X[n],PP.Y[n],PP.Z[n])*(0.0-RKv)/p->dt;
                dw2+=netBuoyZ-pressureDivZ/p->S22-stressDivZ/(thetas*p->S22)+p->ccipol3c(a.fbh3,PP.X[n],PP.Y[n],PP.Z[n])*(0.0-RKw)/p->dt;

                if(debugPrint)
                {
                    cout<<"Z-dir2:drag: "<<Dp*dw<<" buoy: "<<netBuoyZ<<" press: "<<-pressureDivZ/p->S22<<" stress: "<<-stressDivZ/((1-thetas)*p->S22)<<"\n";
                }

                du2=0.25*du2+0.25*du1;
                dv2=0.25*dv2+0.25*dv1;
                dw2=0.25*dw2+0.25*dw1;

                RKu=PP.U[n]+du2*p->dt;
                RKv=PP.V[n]+dv2*p->dt;
                RKw=PP.W[n]+dw2*p->dt;
                
                // RK step 3
                if(shearVel)
                {
                    du=u-RKu;
                    dv=v-RKv;
                    dw=w-RKw;

                    Dp=drag_model(p,PP.d50*PP.PackingFactor[n],du,dv,dw,thetas);

                    du3=Dp*du;
                    dv3=Dp*dv;
                    dw3=Dp*dw;
                }

                du3+=netBuoyX-pressureDivX/p->S22-stressDivX/(thetas*p->S22)+p->ccipol1c(a.fbh1,PP.X[n],PP.Y[n],PP.Z[n])*(0.0-RKu)/p->dt;
                dv3+=netBuoyY-pressureDivY/p->S22-stressDivY/(thetas*p->S22)+p->ccipol2c(a.fbh2,PP.X[n],PP.Y[n],PP.Z[n])*(0.0-RKv)/p->dt;
                dw3+=netBuoyZ-pressureDivZ/p->S22-stressDivZ/(thetas*p->S22)+p->ccipol3c(a.fbh3,PP.X[n],PP.Y[n],PP.Z[n])*(0.0-RKw)/p->dt;

                if(debugPrint)
                {
                    cout<<"Z-dir3:drag: "<<Dp*dw<<" buoy: "<<netBuoyZ<<" press: "<<-pressureDivZ/p->S22<<" stress: "<<-stressDivZ/((1-thetas)*p->S22)<<endl;
                    debugPrint=false;
                }

                if(du2!=du2||du3!=du3)
                {
                    cerr<<"Particle velocity component u resulted in NaN.\n"
                    <<du2<<","<<du3<<"|"<<Dp<<","<<netBuoyX<<","<<pressureDivX<<","<<stressDivX
                    <<endl;
                    exit(1);
                }
                else
                    PP.U[n] += ((2.0/3.0)*du2 + (2.0/3.0)*du3)*p->dt;
                if(!limited)
                {
                    if(dv2!=dv2||dv3!=dv3)
                    {
                        cerr<<"Particle velocity component v resulted in NaN.\n"
                        <<dv2<<","<<dv3<<"|"<<Dp<<","<<netBuoyY<<","<<pressureDivY<<","<<stressDivY
                        <<endl;
                        exit(1);
                    }
                    else
                        PP.V[n] += ((2.0/3.0)*dv2 + (2.0/3.0)*dv3)*p->dt;
                    if(dw2!=dw2||dw3!=dw3)
                    {
                        cerr<<"Particle velocity component w resulted in NaN.\n"
                        <<dw2<<","<<dw3<<"|"<<Dp<<","<<netBuoyZ<<","<<pressureDivZ<<","<<stressDivZ
                        <<endl;
                        exit(1);
                    }
                    else
                        PP.W[n] += ((2.0/3.0)*dw2 + (2.0/3.0)*dw3)*p->dt;
                }
                
                // Pos update
                PP.X[n] += PP.U[n]*p->dt;
                PP.Y[n] += PP.V[n]*p->dt;
                PP.Z[n] += PP.W[n]*p->dt;

                // Sum update
                cellSum[IJK]-=PP.PackingFactor[n];
                if(cellSum[IJK]<0)
                {
                    clog<<"cellSum is below zero in cell ("<<p->XN[IP]<<"-"<<p->XN[IP1]<<","<<p->YN[JP]<<"-"<<p->YN[JP1]<<","<<p->ZN[KP]<<"-"<<p->ZN[KP1]<<") of partition "<<p->mpirank<<" for particle "<<n<<" at ("<<PP.X[n]<<","<<PP.Y[n]<<","<<PP.Z[n]<<") "<<"."<<endl;
                    cellSum[IJK]=0;
                }
                particleStressTensorUpdateIJK(p,a,PP);
                i=p->posc_i(PP.X[n]);
                j=p->posc_j(PP.Y[n]);
                k=p->posc_k(PP.Z[n]);
                cellSum[IJK]+=PP.PackingFactor[n];
                particleStressTensorUpdateIJK(p,a,PP);
            }
        }
    }

    /**
     * @brief Updates the topo field.
     *
     * This function is responsible for updating the topo field.
     * It does so by calculating the difference between the columnSum and the count value.
     *
     * @param p A pointer to the lexer object.
     * @param pgc A reference to the ghostcell object.
     * @param topo A reference to the field4a object.
     * @param d50 The diameter value.
     */
    void Tavouktsoglou::update(lexer *p, ghostcell &pgc, field4a &topo, double &d50)
    {
        double count;
        ILOOP
        JLOOP
        {
            count = 0.0;
            KLOOP
            {
                count += cellSum[IJK] + cellSumTopo[IJK];
                if(k>0 && cellSumTopo[IJKm1]==0)
                break;
            }
            if(count != columnSum[IJ])
            {
                KLOOP
                topo(i,j,k) -= (count-columnSum[IJ])*4.0/3.0*PI*pow(d50/2.0,3)/(p->DXN[IP]*p->DYN[JP]);
            }
            columnSum[IJ] = count;
        }
        pgc.start4a(p,topo,150);
    }

    void Tavouktsoglou::debug(lexer *p, fdm &a, ghostcell &pgc, particles_obj &PP, sediment_fdm &s)
    {
        // double sumCell = 0;
        // double sumTopo = 0;
        PLAINLOOP
        {
            {
                double uvel,vvel,u_abs;
                double signx,signy;
                
                double ux1,vx1,ux2,vx2,uy1,vy1,uy2,vy2;
                double sgx1,sgx2,sgy1,sgy2;
                double ux1_abs,ux2_abs,uy1_abs,uy2_abs;
                
                uvel=0.5*(s.P(i,j)+s.P(i-1,j));
                vvel=0.5*(s.Q(i,j)+s.Q(i,j-1));
                
                u_abs = sqrt(uvel*uvel + vvel*vvel);
                signx=fabs(u_abs)>1.0e-10?uvel/fabs(u_abs):0.0;
                signy=fabs(u_abs)>1.0e-10?vvel/fabs(u_abs):0.0;

                ux1=s.P(i-1,j);
                vx1=0.25*(s.Q(i,j)+s.Q(i-1,j)+s.Q(i,j-1)+s.Q(i-1,j-1)); 
                
                ux2=s.P(i,j);
                vx2=0.25*(s.Q(i,j)+s.Q(i+1,j)+s.Q(i,j-1)+s.Q(i+1,j-1)); 
                
                
                uy1=0.25*(s.P(i,j-1)+s.P(i,j)+s.P(i-1,j-1)+s.P(i-1,j));
                vy1=s.Q(i,j-1); 
                
                uy2=0.25*(s.P(i,j)+s.P(i,j+1)+s.P(i-1,j)+s.P(i-1,j+1));
                vy2=s.Q(i,j); 
                
                
                ux1_abs = sqrt(ux1*ux1 + vx1*vx1);
                ux2_abs = sqrt(ux2*ux2 + vx2*vx2);
                
                uy1_abs = sqrt(uy1*uy1 + vy1*vy1);
                uy2_abs = sqrt(uy2*uy2 + vy2*vy2);
                    
                sgx1=fabs(ux1_abs)>1.0e-10?ux1/fabs(ux1_abs):0.0;
                sgx2=fabs(ux2_abs)>1.0e-10?ux2/fabs(ux2_abs):0.0;
                
                sgy1=fabs(uy1_abs)>1.0e-10?vy1/fabs(uy1_abs):0.0;
                sgy2=fabs(uy2_abs)>1.0e-10?vy2/fabs(uy2_abs):0.0;
                
                // tau * A = F, F/m = a, a*dt = v
                // a.test(i,j,k) = (s.tau_eff(i+1,j)*sgx2-s.tau_eff(i-1,j)*sgx1)*pow(PP.d50/2.0,2)/(4.0/3.0*pow(PP.d50/2.0,3)*PP.density)*p->dt;
            }


            // for(int n=0;n<=k;n++)
            // a.test(i,j,k) += cellSum[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + n-p->kmin]+cellSumTopo[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + n-p->kmin];
            // a.fb(i,j,k) = cellSum[IJK];
            // sumCell += cellSum[IJK];
            // a.vof(i,j,k) = cellSumTopo[IJK];
            // sumTopo += cellSumTopo[IJK];
        }
        // std::cout<<p->mpirank<<": Sum of cellSum: "<<sumCell<<" : "<<sumTopo<<std::endl;
    }

    double Tavouktsoglou::volume(lexer *p, fdm &a, particles_obj &PP)
    {
        double sum=0;
        ILOOP
            JLOOP
                sum += columnSum[IJ];

        return PI*pow(PP.d50,3.0)*(sum)/6.0;
    }

    /// @brief Writes the state of the Tavouktsoglou class to file.
    /// @ToDo Write cellSumTopo 
    void Tavouktsoglou::writeState(lexer *p, ofstream &result)
    {
        float ffn;
        PLAINLOOP
        {
            ffn=cellSum[IJK];
            result.write((char*)&ffn, sizeof (float));
        }
        result.write((char*)&ffn, sizeof (float));  
    }

    /// Reads the state of the Tavouktsoglou class from file.
    /// Reconstructs cellSum, columnSum and stressTensor
    void Tavouktsoglou::readState(lexer *p, ifstream &result)
    {
        float ffn;
        PLAINLOOP
        {
            result.read((char*)&ffn, sizeof (float));
            cellSum[IJK]=double(ffn);
        }

    }

    /// @brief Topo volume in cell div. by particle volume
    /// Uses i,j&k from increment to pass cell identifier
    /// @param d50 Sauter diameter of particles
    /// @return Ceil of number of particles in cell IJK
    double Tavouktsoglou::maxParticlesPerCell(lexer *p, fdm &a, double d50, bool topo, bool cell)
    {   
        double DZN=topo?0:p->DZN[KP];

        if(topo)
        {
            if (a.topo(i,j,k)<=-0.5*p->DZN[KP]+1.0e-13)
            DZN=p->DZN[KP];
            else if(a.topo(i,j,k)<0.5*p->DZN[KP] -5.0e-18)
            DZN=(p->DZN[KP]*0.5 + a.topo(i,j,k));
        }

        return 6.0*p->DXN[IP]*p->DYN[JP]*DZN*(1.0+(cell?0:-a.porosity(i,j,k)))/(PI*pow(d50,3.0));
    }

    /// @brief Calculate complete intra-particle stress trensor
    void Tavouktsoglou::particleStressTensor(lexer *p, fdm &a, ghostcell &pgc, particles_obj &PP)
    {
        double theta;
        int i,j,k;

        PLAINLOOP
        {
            updateParticleStressTensor(p,a,PP,i,j,k);
        }
        pgc.start4V_par(p,stressTensor,10);
    }

    /// @brief Calculate intra-particle stress trensor for cells around (`increment::i`,`increment::j`,`increment::k`)
    void Tavouktsoglou::particleStressTensorUpdateIJK(lexer *p, fdm &a, particles_obj &PP)
    {
        double theta;
        int i,j,k;

        for (int n=-2; n<3; n++)
            for (int m=-2; m<3; m++)
                for (int l=-2; l<3; l++)
                {
                    i=increment::i+n;
                    j=increment::j+m;
                    k=increment::k+l;

                    updateParticleStressTensor(p,a,PP,i,j,k);
                }
    }

    /// @brief Calculate intra-particle stress trensor for cell ( \p i , \p j , \p k )
    void Tavouktsoglou::updateParticleStressTensor(lexer *p, fdm &a, particles_obj &PP, int i, int j, int k)
    {
        double theta=theta_s(p,a,PP,i,j,k);
        stressTensor[IJK]=Ps*pow(theta,beta)/max(theta_crit-theta,epsilon*(1.0-theta));
    }

    /// @brief Calculate solid volume fraction for cell ( \p i , \p j , \p k )
    double Tavouktsoglou::theta_s(lexer *p, fdm &a, particles_obj &PP, int i, int j, int k) const
    {   
        double theta = PI*pow(PP.d50,3.0)*(cellSum[IJK]+cellSumTopo[IJK])/(6.0*p->DXN[IP]*p->DYN[JP]*p->DYN[KP]);
        if(theta>1)
        theta=1;
        return theta;
    }    

    /// @brief Calculate drag force parameter
    double Tavouktsoglou::drag_model(lexer* p, double d, double du, double dv, double dw, double thetas) const
    {
        double thetaf = 1.0-thetas;
        if(thetaf>1.0-theta_crit) // Saveguard
        thetaf=1.0-theta_crit;

        const double dU=sqrt(du*du+dv*dv+dw*dw);
        if(dU==0) // Saveguard
        return 0;

        const double Rep=dU*d*invKinVis;

        const double Cd=24.0*(pow(thetaf,-2.65)+pow(Rep,2.0/3.0)*pow(thetaf,-1.78)/6.0)/Rep;
        const double Dp=Cd*3.0*drho*dU/d/4.0;

        if(Dp!=Dp)
        cout<<thetaf<<","<<dU<<","<<Rep<<","<<Cd<<"|"<<(dU==0)<<endl;

        return Dp;
    }

    /// @brief Calculate number of particles in cell ( \p i , \p j , \p k )
    void Tavouktsoglou::particlePerCell(lexer *p, ghostcell &pgc, particles_obj &PP)
    {
        PLAINLOOP
        cellSum[IJK]=0;

        for(size_t n=0;n<PP.loopindex;n++)
            if(PP.Flag[n]>INT32_MIN)
            {
                i=p->posc_i(PP.X[n]);
                j=p->posc_j(PP.Y[n]);
                k=p->posc_k(PP.Z[n]);
                cellSum[IJK] += PP.PackingFactor[n];
            }
        
        pgc.start4V_par(p,cellSum,11);
        pgc.start4V_par(p,cellSumTopo,11);
    }

    void Tavouktsoglou::make_moving(lexer *p, fdm &a, particles_obj &PP)
    {
        for(size_t n=0;n<PP.loopindex;n++)
            if(PP.Flag[n]==0)
                PP.Flag[n]=1;
    }

    void Tavouktsoglou::erode(lexer *p, fdm &a, particles_obj &PP)
    {
        ILOOP
        {
            JLOOP
            {
                KLOOP
                {
                    if(cellSum[IJK]>0 || cellSumTopo[IJK] == 0)
                    {
                        --k;
                        // for(int qn=0;qn<ppcell*1000;++qn)
                        // {
                        //     x = p->XN[IP] + p->DXN[IP]*double(rand() % irand)/drand;
                        //     y = p->YN[JP] + p->DYN[JP]*double(rand() % irand)/drand;
                        //     z = p->ZN[KP] + p->DZN[KP]*double(rand() % irand)/drand;

                        //     ipolTopo = p->ccipol4_b(a->topo,x,y,z);
                        //     ipolSolid = p->ccipol4_b(a->solid,x,y,z);

                        //     if (!(ipolTopo>tolerance||ipolTopo<-p->Q102*p->DZN[KP]||ipolSolid<0))
                        //     {
                        //         index=PP.add(x,y,z,flag,0,0,0,p->Q41);
                        //         if(movement->seeding(p, PP, index, ppcell))
                        //             break;
                        //     }
                        // }


                        break;
                    }
                }
            }
        }
    }
};

int sediment_particle::state::solid_clean(lexer* p, particles_obj &PP, sediment_particle::movement::base &movement)
{
    int removed = 0;
    for(size_t n=0;n<PP.loopindex;n++)
    if(PP.Flag[n]>0)
    {
        i = p->posc_i(PP.X[n]);
        j = p->posc_j(PP.Y[n]);
        k = p->posc_k(PP.Z[n]);
        movement.remove(p,PP,n);
        PP.erase(n);
        removed++;
    }
    return removed;
}
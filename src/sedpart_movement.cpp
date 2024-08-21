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
    particleStressBased_T2021::particleStressBased_T2021(lexer *p) : drho(p->W1/p->S22) ,invKinVis(p->W1/p->W2), 
                                            Ps(p->Q14),beta(p->Q15),epsilon(p->Q16),theta_crit(p->Q17)
    {
        p->Darray(stressTensor,p->imax*p->jmax*p->kmax);
        p->Darray(cellSum,p->imax*p->jmax*p->kmax);
        p->Darray(cellSumTopo,p->imax*p->jmax*p->kmax);
        p->Darray(columnSum,p->imax*p->jmax);
        dx=p->global_xmax-p->global_xmin;
        LOOP
        {
            dx = min(dx,MIN3(p->DXN[IP],p->DYN[JP],p->DZN[KP]));
        }
        time = p->simtime;
    }

    particleStressBased_T2021::~particleStressBased_T2021()
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
     * @brief Sets up the particleStressBased_T2021 class.
     *
     * This function is responsible for setting up the particleStressBased_T2021 class by calculating
     * the cellSumTopo and columnSum values based on the given parameters.
     *
     * @param p A pointer to the lexer object.
     * @param a A reference to the fdm object.
     * @param diameter The diameter value.
     */
    void particleStressBased_T2021::setup(lexer *p, fdm &a, double &diameter)
    {
        BLOOP
        {
            cellSumTopo[IJK] = maxParticlesPerCell(p,a,diameter);
            columnSum[IJ] += cellSumTopo[IJK];
        }
    }

    void particleStressBased_T2021::setupState(lexer *p, fdm &a, ghostcell &pgc, particles_obj &PP)
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
    seedReturn particleStressBased_T2021::seeding(lexer *p, particles_obj &PP, size_t &index, double max, bool free)
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
    void particleStressBased_T2021::transfer(lexer *p, particles_obj &PP, size_t &index)
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
    void particleStressBased_T2021::remove(lexer *p, particles_obj &PP, size_t &index)
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
    void particleStressBased_T2021::move(lexer *p, fdm &a, ghostcell &pgc, particles_obj &PP, sediment_fdm &s, turbulence &pturb)
    {
        double RKu,RKv,RKw;
        double u,v,w;
        double du1, du2, du3, dv1, dv2, dv3, dw1, dw2, dw3;
        /// @brief Difference between flowfield and particle velocity
        
       
        double pressureDivX=0, pressureDivY=0, pressureDivZ=0;
        double stressDivX=0, stressDivY=0, stressDivZ=0;
        double netBuoyX=(1.0-drho)*p->W20, netBuoyY=(1.0-drho)*p->W21, netBuoyZ=(1.0-drho)*p->W22;

        bool limited = false;
        bool debugPrint = false;
        bool bedLoad = false;
        bool shearVel = true;

        particlePerCell(p,pgc,PP);
        particleStressTensor(p,a,pgc,PP);
        double timeStep = timestep(p,pgc,PP);
        double RKtimeStep = 0.5*timeStep;


        double Du,Dv,Dw;
        double thetas;
        double DragCoeff;
        double du,dv,dw;

        for(int m=0;m<2;m++)
        {
            for(size_t n=0;n<PP.loopindex;n++)
            {
                if(PP.Flag[n]>0) // INT32_MIN
                {
                    if(p->global_xmin+p->Q73>PP.X[n])
                    limited = true;
                    else
                    limited = false;

                    i=p->posc_i(PP.X[n]);
                    j=p->posc_j(PP.Y[n]);
                    k=p->posc_k(PP.Z[n]);

                    thetas=theta_s(p,a,PP,i,j,k);

                    stressDivX = (stressTensor[Ip1JK] - stressTensor[Im1JK])/(p->DXP[IM1]+p->DXP[IP]);
                    stressDivY = (stressTensor[IJp1K] - stressTensor[IJm1K])/(p->DYP[JM1]+p->DYP[JP]);
                    stressDivZ = (stressTensor[IJKp1] - stressTensor[IJKm1])/(p->DZP[KM1]+p->DZP[KP]);

                    pressureDivX = ((a.press(i+1,j,k)-a.phi(i+1,j,k)*a.ro(i+1,j,k)*fabs(p->W22)) - ((a.press(i-1,j,k)-a.phi(i-1,j,k)*a.ro(i-1,j,k)*fabs(p->W22))))/(p->DXP[IM1]+p->DXP[IP]);
                    pressureDivY = ((a.press(i,j+1,k)-a.phi(i,j+1,k)*a.ro(i,j+1,k)*fabs(p->W22)) - ((a.press(i,j-1,k)-a.phi(i,j-1,k)*a.ro(i,j-1,k)*fabs(p->W22))))/(p->DYP[JM1]+p->DYP[JP]);
                    pressureDivZ = ((a.press(i,j,k+1)-a.phi(i,j,k+1)*a.ro(i,j,k+1)*fabs(p->W22)) - ((a.press(i,j,k-1)-a.phi(i,j,k-1)*a.ro(i,j,k-1)*fabs(p->W22))))/(p->DZP[KM1]+p->DZP[KP]);

                    // if(p->ccipol4(a.topo,PP.X[n],PP.Y[n],PP.Z[n])<PP.d50*10)
                    //     bedLoad=true;

                    u=p->ccipol1c(a.u,PP.X[n],PP.Y[n],PP.Z[n]);
                    v=p->ccipol2c(a.v,PP.X[n],PP.Y[n],PP.Z[n]);
                    w=p->ccipol3c(a.w,PP.X[n],PP.Y[n],PP.Z[n]);
                    
                    Du=u-PP.U[n];
                    Dv=v-PP.V[n];
                    Dw=w-PP.W[n];

                    DragCoeff=drag_model(p,PP.d50,Du,Dv,Dw,thetas);

                    // Acceleration
                    du=DragCoeff*Du;
                    dv=DragCoeff*Dv;
                    dw=DragCoeff*Dw;

                    du+=netBuoyX-pressureDivX/p->S22-stressDivX/(thetas*p->S22);
                    dv+=netBuoyY-pressureDivY/p->S22-stressDivY/(thetas*p->S22);
                    dw+=netBuoyZ-pressureDivZ/p->S22-stressDivZ/(thetas*p->S22);

                    if(dw!=dw)
                    {
                    cout<<"NaN detected.\nu: "<<w<<" up: "<<PP.W[n]<<"\npos: "<<PP.X[n]<<","<<PP.Y[n]<<","<<PP.Z[n]<<"\n drag: "<<DragCoeff<<endl;
                    exit(1);
                    }

                    // Vel update
                    PP.U[n]=0.5*(PP.U[n]+p->ccipol1c(a.fbh1,PP.X[n],PP.Y[n],PP.Z[n])*(0.0-PP.U[n]))+du*RKtimeStep;
                    PP.V[n]=0.5*(PP.V[n]+p->ccipol2c(a.fbh2,PP.X[n],PP.Y[n],PP.Z[n])*(0.0-PP.V[n]))+dv*RKtimeStep;
                    PP.W[n]=0.5*(PP.W[n]+p->ccipol3c(a.fbh3,PP.X[n],PP.Y[n],PP.Z[n])*(0.0-PP.W[n]))+dw*RKtimeStep;

                    if(PP.X[n]!=PP.X[n])
                    {
                    cout<<"NaN detected.\npos: "<<PP.X[n]<<endl;
                    exit(1);
                    }
                    
                    // Pos update
                    PP.X[n] += PP.U[n]*RKtimeStep;
                    PP.Y[n] += PP.V[n]*RKtimeStep;
                    if(!limited)
                    PP.Z[n] += PP.W[n]*RKtimeStep;

                    if(PP.X[n]!=PP.X[n])
                    {
                    cout<<"NaN detected.\npos: "<<PP.X[n]<<endl;
                    exit(1);
                    }

                    // Sum update
                    cellSum[IJK]-=PP.PackingFactor[n];
                    i=p->posc_i(PP.X[n]);
                    j=p->posc_j(PP.Y[n]);
                    k=p->posc_k(PP.Z[n]);
                    cellSum[IJK]+=PP.PackingFactor[n];
                    // particleStressTensorUpdateIJK(p,a,PP);
                }
            }
            particleStressTensor(p,a,pgc,PP);
        }
        if(p->mpirank==0)
        {
            time += timeStep;
            cout<<"Sediment time: "<<time<<" time step: "<<timeStep<<endl;
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
    void particleStressBased_T2021::update(lexer *p, ghostcell &pgc, field4a &topo, double &d50)
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

    void particleStressBased_T2021::debug(lexer *p, fdm &a, ghostcell &pgc, particles_obj &PP, sediment_fdm &s)
    {
        // double sumCell = 0;
        // double sumTopo = 0;
        double thetas=0;
        double pressureDivX=0, pressureDivY=0, pressureDivZ=0;
        double stressDivX=0, stressDivY=0, stressDivZ=0;
        bool initalNaN=true;
        double netBuoyX=(1.0-drho)*p->W20, netBuoyY=(1.0-drho)*p->W21, netBuoyZ=(1.0-drho)*p->W22;
        PLAINLOOP
        {

            a.test(i,j,k) = 0;
            for(int n=0;n<=k;n++)
            a.test(i,j,k) += cellSum[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + n-p->kmin]+cellSumTopo[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + n-p->kmin];
            a.fb(i,j,k) = cellSum[IJK];
            // sumCell += cellSum[IJK];
            a.vof(i,j,k) = cellSumTopo[IJK];
            // sumTopo += cellSumTopo[IJK];
            a.Fi(i,j,k)=stressTensor[IJK];

            pressureDivX = (a.press(i+1,j,k)-a.phi(i+1,j,k)*a.ro(i+1,j,k)*fabs(p->W22) - (a.press(i,j,k)-a.phi(i,j,k)*a.ro(i,j,k)*fabs(p->W22)))/(p->DXN[IP]);
            pressureDivY = (0.5*((a.press(i,j+1,k)-a.phi(i,j+1,k)*a.ro(i,j+1,k)*fabs(p->W22))+(a.press(i+1,j+1,k)-a.phi(i+1,j+1,k)*a.ro(i+1,j+1,k)*fabs(p->W22))) - 0.5*((a.press(i,j-1,k)-a.phi(i,j-1,k)*a.ro(i,j-1,k)*fabs(p->W22))+(a.press(i+1,j-1,k)-a.phi(i+1,j-1,k)*a.ro(i+1,j-1,k)*fabs(p->W22))))/(p->DYN[JM1]+p->DYN[JP]);
            pressureDivZ = (0.5*((a.press(i,j,k+1)-a.phi(i,j,k+1)*a.ro(i,j,k+1)*fabs(p->W22))+(a.press(i+1,j,k+1)-a.phi(i+1,j,k+1)*a.ro(i+1,j,k+1)*fabs(p->W22))) - 0.5*((a.press(i,j,k-1)-a.phi(i,j,k-1)*a.ro(i,j,k-1)*fabs(p->W22))+(a.press(i+1,j,k-1)-a.phi(i+1,j,k-1)*a.ro(i+1,j,k-1)*fabs(p->W22))))/(p->DZN[KM1]+p->DZN[KP]);
            a.test1(i,j,k)=+pressureDivX/p->S22;
            a.test2(i,j,k)=+pressureDivY/p->S22;
            a.test3(i,j,k)=+pressureDivZ/p->S22;
            thetas=theta_s(p,a,PP,i,j,k);
            // stressDivX = (p->ccipol4c(stressTensor,p->XN[IP]+p->DXN[IP],p->YN[JP]+0.5*p->DYN[JP],p->ZN[KP]+0.5*p->DZN[KP]) - p->ccipol4c(stressTensor,p->XN[IP]-0.5*p->DXN[IP],p->YN[JP]+0.5*p->DYN[JP],p->ZN[KP]+0.5*p->DZN[KP]))/p->DXN[IP];
            // stressDivY = (p->ccipol4c(stressTensor,p->XN[IP]+0.5*p->DXN[IP],p->YN[JP]+p->DYN[JP],p->ZN[KP]+0.5*p->DZN[KP]) - p->ccipol4c(stressTensor,p->XN[IP]+0.5*p->DXN[IP],p->YN[JP]-0.5*p->DYN[JP],p->ZN[KP]+0.5*p->DZN[KP]))/p->DYN[JP];
            stressDivX = (stressTensor[Ip1JK] - stressTensor[Im1JK])/(p->DXN[IM1]+p->DXN[IP]);
            stressDivY = (0.5*(stressTensor[IJp1K]+stressTensor[Ip1Jp1K]) - 0.5*(stressTensor[IJm1K]+stressTensor[Ip1Jm1K]))/(p->DYN[JM1]+p->DYN[JP]);
            stressDivZ = (0.5*(stressTensor[IJKp1]+stressTensor[Ip1JKp1]) - 0.5*(stressTensor[IJKm1]+stressTensor[Ip1JKm1]))/(p->DZN[KM1]+p->DZN[KP]);
            a.test4(i,j,k)=-stressDivX/(thetas*p->S22);
            a.test5(i,j,k)=-stressDivY/(thetas*p->S22);
            a.test6(i,j,k)=-stressDivZ/(thetas*p->S22);
            a.test7(i,j,k)=a.test4(i,j,k)+netBuoyX;
            a.test8(i,j,k)=a.test5(i,j,k)+netBuoyY;
            a.test9(i,j,k)=a.test6(i,j,k)+netBuoyZ;
            // if(a.test6(i,j,k)!=a.test6(i,j,k) && initalNaN)
            // {
            //     cout<<"NaN detected: "<<stressTensor[IJKp1]<<"+"<<stressTensor[Ip1JKp1]<<" - "<<stressTensor[IJKm1]<<"+"<<stressTensor[Ip1JKm1]<<" / "<<p->DZN[KM1]<<"+"<<p->DZN[KP]<<"\n";
            //     cout<<"thetas: "<<thetas<<" density: "<<p->S22<<endl;
            //     initalNaN=false;
            // }
        }
        // double netBuoyX=(1.0-drho)*p->W20, netBuoyY=(1.0-drho)*p->W21, netBuoyZ=(1.0-drho)*p->W22;
        if(p->mpirank==0)
        cout<<"NetBuoy: "<<netBuoyX<<","<<netBuoyY<<","<<netBuoyZ<<endl;
        // std::cout<<p->mpirank<<": Sum of cellSum: "<<sumCell<<" : "<<sumTopo<<std::endl;
    }

    double particleStressBased_T2021::volume(lexer *p, fdm &a, particles_obj &PP)
    {
        double sum=0;
        ILOOP
            JLOOP
                sum += columnSum[IJ];

        return PI*pow(PP.d50,3.0)*(sum)/6.0;
    }

    /// @brief Writes the state of the particleStressBased_T2021 class to file.
    /// @ToDo Write cellSumTopo 
    void particleStressBased_T2021::writeState(lexer *p, ofstream &result)
    {
        float ffn;
        PLAINLOOP
        {
            ffn=cellSum[IJK];
            result.write((char*)&ffn, sizeof (float));
        }
        result.write((char*)&ffn, sizeof (float));  
    }

    /// Reads the state of the particleStressBased_T2021 class from file.
    /// Reconstructs cellSum, columnSum and stressTensor
    void particleStressBased_T2021::readState(lexer *p, ifstream &result)
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
    double particleStressBased_T2021::maxParticlesPerCell(lexer *p, fdm &a, double d50, bool topo, bool cell)
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
    void particleStressBased_T2021::particleStressTensor(lexer *p, fdm &a, ghostcell &pgc, particles_obj &PP)
    {
        double theta;
        int i,j,k;

        BLOOP
        {
            updateParticleStressTensor(p,a,PP,i,j,k);
        }
        pgc.start4V_par(p,stressTensor,10);
    }

    /// @brief Calculate intra-particle stress trensor for cells around (`increment::i`,`increment::j`,`increment::k`)
    void particleStressBased_T2021::particleStressTensorUpdateIJK(lexer *p, fdm &a, particles_obj &PP)
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
    void particleStressBased_T2021::updateParticleStressTensor(lexer *p, fdm &a, particles_obj &PP, int i, int j, int k)
    {
        double theta=theta_s(p,a,PP,i,j,k);
        stressTensor[IJK]=Ps*pow(theta,beta)/max(theta_crit-theta,epsilon*(1.0-theta));
    }

    /// @brief Calculate solid volume fraction for cell ( \p i , \p j , \p k )
    double particleStressBased_T2021::theta_s(lexer *p, fdm &a, particles_obj &PP, int i, int j, int k) const
    {   
        double theta = PI*pow(PP.d50,3.0)*(cellSum[IJK]+cellSumTopo[IJK])/(6.0*p->DXN[IP]*p->DYN[JP]*p->DYN[KP]);
        if(theta>1)
        theta=1;
        if(theta<0)
        theta=0;
        // if(theta!=theta)
        // theta=0;
        return theta;
    }    

    /// @brief Calculate drag force parameter
    double particleStressBased_T2021::drag_model(lexer* p, double d, double du, double dv, double dw, double thetas) const
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

        // if(Dp!=Dp)
        // cout<<thetaf<<","<<dU<<","<<Rep<<","<<Cd<<"|"<<(dU==0)<<endl;

        return Dp;
    }

    /// @brief Calculate number of particles in cell ( \p i , \p j , \p k )
    void particleStressBased_T2021::particlePerCell(lexer *p, ghostcell &pgc, particles_obj &PP)
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

    void particleStressBased_T2021::make_moving(lexer *p, fdm &a, particles_obj &PP)
    {
        for(size_t n=0;n<PP.loopindex;n++)
            if(PP.Flag[n]==0)
                PP.Flag[n]=1;
    }

    void particleStressBased_T2021::erode(lexer *p, fdm &a, particles_obj &PP)
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

    void particleStressBased_T2021::bedReDistribution(lexer *p, fdm &a, ghostcell &pgc, particles_obj &PP)
    {
        /// per i,j,k which particles n are in cell
        std::vector<std::vector<size_t>> particlesInCell(p->imax*p->jmax*p->kmax);
        for(size_t n=0;n<PP.loopindex;n++)
        {
            if(PP.Flag[n]>=0) // INT32_MIN
            {
                i=p->posc_i(PP.X[n]);
                j=p->posc_j(PP.Y[n]);
                k=p->posc_k(PP.Z[n]);
                particlesInCell[IJK].push_back(n);
            }
        }
        double tolerance = 5e-18;
        double x,y,z,ipolTopo,ipolSolid;
        int irand=10000;
        double drand=10000;
        PLAINLOOP
        for (auto index : particlesInCell[IJK])
        {
            x = p->XN[IP] + p->DXN[IP]*double(rand() % irand)/drand;
            y = p->YN[JP] + p->DYN[JP]*double(rand() % irand)/drand;
            z = p->ZN[KP] + p->DZN[KP]*double(rand() % irand)/drand;

            ipolTopo = p->ccipol4_b(a.topo,x,y,z);
            ipolSolid = p->ccipol4_b(a.solid,x,y,z);

            if (!(ipolTopo>tolerance||ipolTopo<-p->Q102*p->DZN[KP]||ipolSolid<0))
            {
                PP.X[index] = x;
                PP.Y[index] = y;
                PP.Z[index] = z;
            }
        }
        
    }

    double particleStressBased_T2021::timestep(lexer *p, ghostcell &pgc, particles_obj &PP)
    {
        double maxVelU=0,maxVelV=0,maxVelW=0;
        for(size_t n=0;n<PP.loopindex;n++)
        {
            if(PP.Flag[n]>0)
            {
                maxVelU=max(maxVelU,fabs(PP.U[n]));
                maxVelV=max(maxVelV,fabs(PP.V[n]));
                maxVelW=max(maxVelW,fabs(PP.W[n]));
            }
        }
        dx = pgc.globalmin(dx);
        // if(p->mpirank==0)
        // cout<<"TimeStep: dx: "<<dx<<" vel: "<<sqrt(maxVelU*maxVelU+maxVelV*maxVelV+maxVelW*maxVelW)<<endl;
        double meanVel=sqrt(maxVelU*maxVelU+maxVelV*maxVelV+maxVelW*maxVelW);
        if(meanVel==0)
        meanVel=dx*p->S14/p->S13;
        double dt = p->S14 * (dx/meanVel);
        dt = min(dt,p->S13);
        dt = pgc.globalmin(dt);
        return dt;
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
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
                                            Ps(p->Q14),beta(p->Q15),epsilon(p->Q16),theta_crit(p->Q17),bedChange(p)
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
        // particlePerCell(p,pgc,PP);
        // PLAINLOOP
        // columnSum[IJ] += cellSumTopo[IJK]+cellSum[IJK];
        // particleStressTensor(p,a,pgc,PP);
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
        double u=0,v=0,w=0;
        double du1, du2, du3, dv1, dv2, dv3, dw1, dw2, dw3;
        double ws=0;
        double topoDist=0;
        /// @brief Difference between flowfield and particle velocity
        
       
        double pressureDivX=0, pressureDivY=0, pressureDivZ=0;
        double stressDivX=0, stressDivY=0, stressDivZ=0;
        double netBuoyX=0, netBuoyY=0, netBuoyZ=0;
        // double netBuoyX=(1.0-drho)*p->W20, netBuoyY=(1.0-drho)*p->W21, netBuoyZ=(1.0-drho)*p->W22;

        bool limited = true;
        bool debugPrint = false;
        bool bedLoad = false;
        bool shearVel = true;

        particlePerCell(p,pgc,PP);
        // particleStressTensor(p,a,pgc,PP);
        timestep(p,pgc,PP);
        double RKtimeStep = 0.5*p->dtsed;


        double Du=0,Dv=0,Dw=0;
        double thetas=0;
        double DragCoeff=0;
        double du=0,dv=0,dw=0;
        // int counter = 0;

        for(int m=0;m<2;m++)
        {
            for(size_t n=0;n<PP.loopindex;n++)
            {
                if(PP.Flag[n]>0) // INT32_MIN
                {
                    // ++counter;
                    // if(p->global_xmin+p->Q73>PP.X[n])
                    // limited = true;
                    // else
                    // limited = false;

                    i=p->posc_i(PP.X[n]);
                    j=p->posc_j(PP.Y[n]);
                    k=p->posc_k(PP.Z[n]);

                    thetas=theta_s(p,a,PP,i,j,k);

                    // stressDivX = (stressTensor[Ip1JK] - stressTensor[Im1JK])/(p->DXP[IM1]+p->DXP[IP]);
                    // stressDivY = (stressTensor[IJp1K] - stressTensor[IJm1K])/(p->DYP[JM1]+p->DYP[JP]);
                    // stressDivZ = (stressTensor[IJKp1] - stressTensor[IJKm1])/(p->DZP[KM1]+p->DZP[KP]);

                    // pressureDivX = ((a.press(i+1,j,k)-a.phi(i+1,j,k)*a.ro(i+1,j,k)*fabs(p->W22)) - ((a.press(i-1,j,k)-a.phi(i-1,j,k)*a.ro(i-1,j,k)*fabs(p->W22))))/(p->DXP[IM1]+p->DXP[IP]);
                    // pressureDivY = ((a.press(i,j+1,k)-a.phi(i,j+1,k)*a.ro(i,j+1,k)*fabs(p->W22)) - ((a.press(i,j-1,k)-a.phi(i,j-1,k)*a.ro(i,j-1,k)*fabs(p->W22))))/(p->DYP[JM1]+p->DYP[JP]);
                    // pressureDivZ = ((a.press(i,j,k+1)-a.phi(i,j,k+1)*a.ro(i,j,k+1)*fabs(p->W22)) - ((a.press(i,j,k-1)-a.phi(i,j,k-1)*a.ro(i,j,k-1)*fabs(p->W22))))/(p->DZP[KM1]+p->DZP[KP]);

                    // if(p->ccipol4(a.topo,PP.X[n],PP.Y[n],PP.Z[n])<PP.d50*10)
                    //     bedLoad=true;

                    topoDist=p->ccipol4(a.topo,PP.X[n],PP.Y[n],PP.Z[n]);

                    if(topoDist<velDist*p->DZP[KP])
                    {
                        u=p->ccipol1c(a.u,PP.X[n],PP.Y[n],PP.Z[n]+velDist*p->DZP[KP]-topoDist);
                        v=p->ccipol2c(a.v,PP.X[n],PP.Y[n],PP.Z[n]+velDist*p->DZP[KP]-topoDist);
                        // w=p->ccipol3c(a.w,PP.X[n],PP.Y[n],PP.Z[n]+velDist*p->DZP[KP]-topoDist);
                        if(debugPrint)
                        {
                            cout<<PP.Z[n]+velDist*p->DZP[KP]-topoDist<<endl;
                            debugPrint=false;
                        }
                    }
                    else
                    {
                        u=p->ccipol1c(a.u,PP.X[n],PP.Y[n],PP.Z[n]);
                        v=p->ccipol2c(a.v,PP.X[n],PP.Y[n],PP.Z[n]);
                        // w=p->ccipol3c(a.w,PP.X[n],PP.Y[n],PP.Z[n]);
                    }

                    // PP.Uf[n]=u;
                    // PP.Vf[n]=v;
                    // PP.Wf[n]=w;
                    
                    Du=u-PP.U[n];
                    Dv=v-PP.V[n];
                    Dw=w-PP.W[n];

                    DragCoeff=drag_model(p,PP.d50,Du,Dv,Dw,thetas);
                    // if(debugPrint)
                    // {
                    //     cout<<DragCoeff<<endl;
                    //     debugPrint=false;
                    // }

                    PP.drag[n]=DragCoeff;

                    // sedimentation
                    // if(topoDist>2.5*PP.d50)
                    // {
                    //     ws = sedimentation_velocity(p,PP.d50,Du,Dv,Dw,thetas);
                    //     // if(fabs(topoDist)<p->DZP[KP])
                    //     // {
                    //     //     ws *=topoDist/p->DZP[KP];
                    //     // }
                    //     Dw-=ws;
                    // }

                    // Acceleration
                    du=DragCoeff*Du;
                    dv=DragCoeff*Dv;
                    // dw=DragCoeff*Dw;

                    du+=netBuoyX-pressureDivX/p->S22-stressDivX/(thetas*p->S22);
                    dv+=netBuoyY-pressureDivY/p->S22-stressDivY/(thetas*p->S22);
                    // dw+=netBuoyZ-pressureDivZ/p->S22-stressDivZ/(thetas*p->S22);

                    // if(debugPrint)
                    // {
                    //     cout<<netBuoyZ<<" "<<-stressDivZ/(thetas*p->S22)<<endl;
                    //     debugPrint=false;
                    // }

                    // if(dw!=dw)
                    // {
                    // cout<<"NaN detected.\nu: "<<w<<" up: "<<PP.W[n]<<"\npos: "<<PP.X[n]<<","<<PP.Y[n]<<","<<PP.Z[n]<<"\n drag: "<<DragCoeff<<endl;
                    // exit(1);
                    // }

                    // if(p->mpirank==1&&n==50431&&p->count>=212)
                    // {
                    // cout<<"pos: "<<PP.X[n]<<","<<PP.Y[n]<<","<<PP.Z[n]<<"|"<<PP.Flag[n]
                    // <<"\ndrag: "<<DragCoeff<<" w: "<<w<<" wp: "<<PP.W[n]<<" dw: "<<dw<<" df: "<<0.5*p->ccipol3c(a.fbh3,PP.X[n],PP.Y[n],PP.Z[n])*(0.0-PP.W[n])
                    // <<"\nstress: "<<-stressDivZ/(thetas*p->S22)<<" press: "<<-pressureDivZ/p->S22<<" buoy: "<<netBuoyZ
                    // <<endl;
                    // // StressDiv explodes
                    // PP.Flag[n]=10;}

                    // Vel update
                    PP.U[n]=0.5*(PP.U[n]+p->ccipol1c(a.fbh1,PP.X[n],PP.Y[n],PP.Z[n])*(0.0-PP.U[n]))+du*RKtimeStep;
                    PP.V[n]=0.5*(PP.V[n]+p->ccipol2c(a.fbh2,PP.X[n],PP.Y[n],PP.Z[n])*(0.0-PP.V[n]))+dv*RKtimeStep;
                    // PP.W[n]=0.5*(PP.W[n]+p->ccipol3c(a.fbh3,PP.X[n],PP.Y[n],PP.Z[n])*(0.0-PP.W[n]))+dw*RKtimeStep;

                    if(PP.U[n]!=PP.U[n])
                    {
                    cout<<"NaN detected.\naccel: "<<du<<" df: "<<p->ccipol1c(a.fbh1,PP.X[n],PP.Y[n],PP.Z[n])*(0.0-PP.U[n])<<endl;
                    exit(1);
                    }
                    if(PP.V[n]!=PP.V[n])
                    {
                    cout<<"NaN detected.\naccel: "<<dv<<" df: "<<p->ccipol2c(a.fbh2,PP.X[n],PP.Y[n],PP.Z[n])*(0.0-PP.V[n])<<endl;
                    exit(1);
                    }
                    if(PP.W[n]!=PP.W[n])
                    // ||p->count==213)
                    {
                    cout<<p->mpirank<<"NaN detected.\naccel: "<<dw<<" df: "<<p->ccipol3c(a.fbh3,PP.X[n],PP.Y[n],PP.Z[n])*(0.0-PP.W[n])<<"\n"<<w<<" : "<<DragCoeff<<"pos: "<<PP.X[n]<<","<<PP.Y[n]<<","<<PP.Z[n]<<endl;
                    exit(1);
                    }
                    
                    // Pos update
                    PP.X[n] += PP.U[n]*RKtimeStep;
                    PP.Y[n] += PP.V[n]*RKtimeStep;
                    if(!limited)
                    PP.Z[n] += PP.W[n]*RKtimeStep;

                    // Sum update
                    cellSum[IJK]-=PP.PackingFactor[n];
                    i=p->posc_i(PP.X[n]);
                    j=p->posc_j(PP.Y[n]);
                    k=p->posc_k(PP.Z[n]);
                    cellSum[IJK]+=PP.PackingFactor[n];
                    // particleStressTensorUpdateIJK(p,a,PP);
                }
            }
            // particleStressTensor(p,a,pgc,PP);
        }
        if(p->mpirank==0)
        {
            p->sedtime += p->dtsed;
            cout<<"Sediment time: "<<p->sedtime<<" time step: "<<p->dtsed<<endl;
            // int Flag0=0,Flag1=0;
            // for(size_t n=0;n<PP.loopindex;n++)
            // {
            //     if(PP.Flag[n]==0)
            //     Flag0++;
            //     else if(PP.Flag[n]==1)
            //     Flag1++;
            // }
            // cout<<"Particles: "<<PP.size<<" Flag0: "<<Flag0<<" Flag1: "<<Flag1<<" moved: "<<counter/2<<endl;
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
    void particleStressBased_T2021::update(lexer *p, fdm &a, ghostcell &pgc, particles_obj &PP)
    {
        ILOOP
            JLOOP
            {
                KLOOP
                    a.topo(i,j,k) -= bedChange[IJ]*1.0/6.0*PI*pow(PP.d50,3)/(p->DXN[IP]*p->DYN[JP]);
                columnSum[IJ] += bedChange[IJ];
                activateNew(p,a,PP);
                bedChange[IJ] = 0;
            }
        pgc.start4a(p,a.topo,150);
    }

    void particleStressBased_T2021::activateNew(lexer *p, fdm &a, particles_obj &PP)
    {
        double tolerance = 5e-18;
        double x,y,z,ipolTopo,ipolSolid;
        int flag=0;
        size_t index;

        double counter=0;
        int maxTries=1000;
        int tries=0;

        while(counter<bedChange[IJ] && tries<maxTries)
        {   
            x = p->XN[IP] + p->DXN[IP]*double(rand() % 10000)/10000.0;
            y = p->YN[JP] + p->DYN[JP]*double(rand() % 10000)/10000.0;
            k = 0;
            z = p->ZN[KP]-a.topo(i,j,k) - 1.2*p->DZN[KP]*double(rand() % 10000)/10000.0;
            k = p->posc_k(z);

            ipolTopo = p->ccipol4_b(a.topo,x,y,z);
            ipolSolid = p->ccipol4_b(a.solid,x,y,z);

            if (!(ipolTopo>tolerance||ipolTopo<-p->Q102*p->DZN[KP]||ipolSolid<0) && cellSumTopo[IJK]>=p->Q41)
            {
                index = PP.add(x,y,z,flag,0,0,0,p->Q41);
                counter += PP.PackingFactor[index];
                cellSumTopo[IJK] -= PP.PackingFactor[index];
            }
            ++tries;
        }
    }

    void particleStressBased_T2021::debug(lexer *p, fdm &a, ghostcell &pgc, particles_obj &PP, sediment_fdm &s)
    {
        PLAINLOOP
        {
            a.test(i,j,k) = cellSumTopo[IJK];
        }
        double topoDist=0;
        double u=0,v=0;
        for(size_t n=0;n<PP.loopindex;n++)
            if(PP.Flag[n]>INT32_MIN)
            {
                if(PP.Flag[n]==0)
                {
                    PP.shear_eff[n]=p->ccslipol4(s.tau_eff,PP.X[n],PP.Y[n]);
                    PP.shear_crit[n]=p->ccslipol4(s.tau_crit,PP.X[n],PP.Y[n]);
                }

                topoDist=p->ccipol4(a.topo,PP.X[n],PP.Y[n],PP.Z[n]);

                if(topoDist<velDist*p->DZP[KP])
                {
                    u=p->ccipol1c(a.u,PP.X[n],PP.Y[n],PP.Z[n]+velDist*p->DZP[KP]-topoDist);
                    v=p->ccipol2c(a.v,PP.X[n],PP.Y[n],PP.Z[n]+velDist*p->DZP[KP]-topoDist);
                }
                else
                {
                    u=p->ccipol1c(a.u,PP.X[n],PP.Y[n],PP.Z[n]);
                    v=p->ccipol2c(a.v,PP.X[n],PP.Y[n],PP.Z[n]);
                }

                PP.Uf[n]=u;
                PP.Vf[n]=v;
            }
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
        // cout<<du<<","<<dv<<","<<dw<<endl;
        if(dU==0) // Saveguard
        return 0;

        const double Rep=dU*d*invKinVis;

        const double Cd=24.0*(pow(thetaf,-2.65)+pow(Rep,2.0/3.0)*pow(thetaf,-1.78)/6.0)/Rep;
        const double Dp=Cd*3.0*drho*dU/d/4.0;

        // if(Dp!=Dp)
        // cout<<thetaf<<","<<dU<<","<<Rep<<","<<Cd<<"|"<<(dU==0)<<endl;

        return Dp;
    }

    double particleStressBased_T2021::sedimentation_velocity(lexer *p, double d, double du, double dv, double dw, double thetas) const
    {
        const double dU=sqrt(du*du+dv*dv+dw*dw);
        if(dU==0) // Saveguard
        return 0;

        const double Rep=dU*d*invKinVis;

        const double Cd=24.0/Rep+4.0/sqrt(Rep)+0.4;
        const double ws_single=sqrt(4.0/3.0*(p->S22-p->W1)/p->W1*d*fabs(p->W22)/Cd);
        double n;
        if(Rep<=0.2)
        n=4.65;
        else if(Rep<=1)
        n=4.35*pow(Rep,-0.03);
        else if(Rep<=500)
        n=4.45*pow(Rep,-0.1);
        else
        n=2.39;
        const double ws_swarm = ws_single*pow((1.0-thetas),n);
        return ws_swarm;
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

    void particleStressBased_T2021::erode(lexer *p, fdm &a, particles_obj &PP, sediment_fdm &s)
    {
        double shear_eff;
        double shear_crit;
        for(size_t n=0;n<PP.loopindex;n++)
            if(PP.Flag[n]==0)
            {
                shear_eff=p->ccslipol4(s.tau_eff,PP.X[n],PP.Y[n]);
                shear_crit=p->ccslipol4(s.tau_crit,PP.X[n],PP.Y[n]);
                if(shear_eff>shear_crit)
                if(p->ccipol4_b(a.topo,PP.X[n],PP.Y[n],PP.Z[n])+5.0*PP.d50>0)
                {
                    PP.Flag[n]=1;

                    PP.shear_eff[n]=shear_eff;
                    PP.shear_crit[n]=shear_crit;

                    i=p->posc_i(PP.X[n]);
                    j=p->posc_j(PP.Y[n]);
                    bedChange[IJ] -= PP.PackingFactor[n];
                }
            }
    }

    void particleStressBased_T2021::deposit(lexer *p, fdm &a, particles_obj &PP, sediment_fdm &s)
    {
        double shear_eff;
        double shear_crit;
        for(size_t n=0;n<PP.loopindex;n++)
            if(PP.Flag[n]==1)
            {
                shear_eff=p->ccslipol4(s.tau_eff,PP.X[n],PP.Y[n]);
                shear_crit=p->ccslipol4(s.tau_crit,PP.X[n],PP.Y[n]);
                if(shear_crit>shear_eff)
                if(p->ccipol4_b(a.topo,PP.X[n],PP.Y[n],PP.Z[n])<PP.d50)
                {
                    PP.Flag[n]=0;
                    PP.U[n]=0;
                    PP.V[n]=0;
                    PP.W[n]=0;

                    PP.Uf[n]=0;
                    PP.Vf[n]=0;
                    PP.Wf[n]=0;
                    // PP.shear_eff[n]=shear_eff;
                    // PP.shear_crit[n]=shear_crit;
                    PP.drag[n]=0;

                    i=p->posc_i(PP.X[n]);
                    j=p->posc_j(PP.Y[n]);
                    bedChange[IJ] += PP.PackingFactor[n];
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

    void particleStressBased_T2021::timestep(lexer *p, ghostcell &pgc, particles_obj &PP)
    {
        double maxVelU=0,maxVelV=0,maxVelW=0;
        for(size_t n=0;n<PP.loopindex;n++)
        {
            if(PP.Flag[n]>0)
            {
                maxVelU=max(maxVelU,fabs(PP.U[n]));
                maxVelV=max(maxVelV,fabs(PP.V[n]));
                maxVelW=max(maxVelW,fabs(PP.W[n]));
                // if(PP.W[n]>14000)
                // cout<<p->mpirank<<":"<<n<<":"<<PP.W[n]<<endl;
            }
        }

        // if(p->count==379)
        // cout<<p->mpirank<<":"<<maxVelU<<","<<maxVelV<<","<<maxVelW<<endl;
        dx = pgc.globalmin(dx);
        if(p->mpirank==0)
        cout<<"TimeStep: dx: "<<dx<<" vel: "<<sqrt(maxVelU*maxVelU+maxVelV*maxVelV+maxVelW*maxVelW)<<endl;
        double meanVel=sqrt(maxVelU*maxVelU+maxVelV*maxVelV+maxVelW*maxVelW);
        if(meanVel==0)
        meanVel=dx*p->S14/p->S13;
        p->dtsed = p->S14 * (dx/meanVel);
        p->dtsed = min(p->dtsed,p->S13);
        p->dtsed = pgc.globalmin(p->dtsed);
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
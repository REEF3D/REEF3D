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

partres::partres(lexer *p) : drho(p->W1/p->S22) ,invKinVis(p->W1/p->W2), 
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

partres::~partres()
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
     * @brief Sets up the partres class.
     *
     * This function is responsible for setting up the partres class by calculating
     * the cellSumTopo and columnSum values based on the given parameters.
     *
     * @param p A pointer to the lexer object.
     * @param a A reference to the fdm object.
     * @param diameter The diameter value.
     */
void partres::setup(lexer *p, fdm &a, double &diameter)
{
        BLOOP
        {
            cellSumTopo[IJK] = maxParticlesPerCell(p,a,diameter);
            columnSum[IJ] += cellSumTopo[IJK];
        }
}

void partres::setupState(lexer *p, fdm &a, ghostcell &pgc, particles_obj &PP)
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
seedReturn partres::seeding(lexer *p, particles_obj &PP, size_t &index, double max, bool free)
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
void partres::transfer(lexer *p, particles_obj &PP, size_t &index)
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
    void partres::remove(lexer *p, particles_obj &PP, size_t &index)
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
void partres::move(lexer *p, fdm &a, ghostcell &pgc, particles_obj &PP, sediment_fdm &s, turbulence &pturb)
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
        double dU=0;

        for(int m=0;m<2;m++)
        {
            for(size_t n=0;n<PP.loopindex;n++)
            {
                if(PP.Flag[n]>0)
                {

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
                    // Dw=w-PP.W[n];

                    switch (p->Q202)
                    {
                        case 0:
                        {
                            DragCoeff=drag_model(p,PP.d50,Du,Dv,Dw,thetas);
                            break;
                        }
                        case 1:
                        {
                            relative_velocity(p,a,PP,n,Du,Dv,Dw);
                            dU=sqrt(Du*Du+Dv*Dv+Dw*Dw);
                            const double Re_p = dU*PP.d50/(p->W2/p->W1);
                            DragCoeff=drag_coefficient(Re_p);
                            break;
                        }
                    }
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
                    switch (p->Q202)
                    {
                        case 0:
                        {
                            du=DragCoeff*Du;
                            dv=DragCoeff*Dv;
                            // dw=DragCoeff*Dw;
                            du+=netBuoyX-pressureDivX/p->S22-stressDivX/(thetas*p->S22);
                            dv+=netBuoyY-pressureDivY/p->S22-stressDivY/(thetas*p->S22);
                            // dw+=netBuoyZ-pressureDivZ/p->S22-stressDivZ/(thetas*p->S22);
                            break;
                        }
                        case 1:
                        {
                            const double Fd = DragCoeff * PI/8.0 * pow(PP.d50,2) * p->W1 * pow(dU,2);
                            DragCoeff = Fd /(p->S22*PI*pow(PP.d50,3.0)/6.0);
                            du=DragCoeff;
                            dv=DragCoeff;
                            // dw=DragCoeff;
                            break;
                        }
                    }

                    // Vel update
                    PP.U[n]=0.5*(PP.U[n]+p->ccipol1c(a.fbh1,PP.X[n],PP.Y[n],PP.Z[n])*(0.0-PP.U[n]))+du*RKtimeStep;
                    PP.V[n]=0.5*(PP.V[n]+p->ccipol2c(a.fbh2,PP.X[n],PP.Y[n],PP.Z[n])*(0.0-PP.V[n]))+dv*RKtimeStep;
                    // PP.W[n]=0.5*(PP.W[n]+p->ccipol3c(a.fbh3,PP.X[n],PP.Y[n],PP.Z[n])*(0.0-PP.W[n]))+dw*RKtimeStep;
                    PP.W[n]=0.0;

                    if(PP.U[n]!=PP.U[n] || PP.V[n]!=PP.V[n] || PP.W[n]!=PP.W[n])
                    {
                    cout<<"NaN detected.\nDu: "<<Du<<" Dv: "<<Dv<<" Dw: "<<Dw<<"\nDrag: "<<DragCoeff<<endl;
                    exit(1);
                    }
                    
                    // Pos update
                    PP.X[n] += PP.U[n]*RKtimeStep;
                    PP.Y[n] += PP.V[n]*RKtimeStep;
                    // PP.Z[n] += PP.W[n]*RKtimeStep;

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
    void partres::update(lexer *p, fdm &a, ghostcell &pgc, particles_obj &PP)
    {
        int count=0;
        ILOOP
            JLOOP
            {
                KLOOP
                {
                    a.topo(i,j,k) -= bedChange[IJ]*1.0/6.0*PI*pow(PP.d50,3)/(p->DXN[IP]*p->DYN[JP]);
                    a.fb(i,j,k) = bedChange[IJ];
                }
                columnSum[IJ] += bedChange[IJ];
                if(bedChange[IJ]<0)
                count+=activateNew(p,a,PP);
                bedChange[IJ] = 0;
            }
        pgc.start4a(p,a.topo,150);
        if(count>0)
        cout<<"On partion "<<p->mpirank<<" were "<<count<<" additional particles activated."<<endl;
    }

int partres::activateNew(lexer *p, fdm &a, particles_obj &PP)
{
        double tolerance = 5e-18;
        double x,y,z,ipolTopo,ipolSolid;
        int flag=0;
        size_t index;

        double counter=0;
        int maxTries=1000;
        int tries=0;
        int count=0;

        if(PP.size-bedChange[IJ]>0.9*PP.capacity)
            PP.reserve();

        while(counter<-bedChange[IJ] && tries<maxTries)
        {   
            x = p->XN[IP] + p->DXN[IP]*double(rand() % 10000)/10000.0;
            y = p->YN[JP] + p->DYN[JP]*double(rand() % 10000)/10000.0;
            k = 0;
            z = p->ZN[KP]+0.5*p->DZP[KP]-a.topo(i,j,k) - 5.0*PP.d50*double(rand() % 10000)/10000.0;
            k = p->posc_k(z);

            ipolTopo = p->ccipol4_b(a.topo,x,y,z);
            ipolSolid = p->ccipol4_b(a.solid,x,y,z);

            if (!(ipolTopo>tolerance||ipolTopo<-p->Q102*p->DZN[KP]||ipolSolid<0))
                if(cellSumTopo[IJK]>=p->Q41)
                {
                    index = PP.add(x,y,z,flag,0,0,0,p->Q41);
                    counter += PP.PackingFactor[index];
                    cellSumTopo[IJK] -= PP.PackingFactor[index];
                    ++count;
                }
                else if (cellSumTopo[IJK]+cellSumTopo[IJKm1]>=p->Q41)
                {
                    index = PP.add(x,y,z,flag,0,0,0,p->Q41);
                    counter += PP.PackingFactor[index];
                    cellSumTopo[IJK] -= PP.PackingFactor[index];
                    cellSumTopo[IJKm1] += cellSumTopo[IJK];
                    cellSumTopo[IJK] = 0;

                    ++count;
                    break;
                }
            
            ++tries;
        }
        return count;
}

void partres::debug(lexer *p, fdm &a, ghostcell &pgc, particles_obj &PP, sediment_fdm &s)
{
        PLAINLOOP
        {
            a.test(i,j,k) = cellSumTopo[IJK];
            a.vof(i,j,k) = (s.tau_eff[IJ]>s.tau_crit[IJ]?-1:1);
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

double partres::volume(lexer *p, fdm &a, particles_obj &PP)
{
        double sum=0;
        ILOOP
            JLOOP
                sum += columnSum[IJ];

        return PI*pow(PP.d50,3.0)*(sum)/6.0;
}

/// @brief Writes the state of the partres class to file.
/// @ToDo Write cellSumTopo 
void partres::writeState(lexer *p, ofstream &result)
{
        float ffn;
        PLAINLOOP
        {
            ffn=cellSum[IJK];
            result.write((char*)&ffn, sizeof (float));
        }
        result.write((char*)&ffn, sizeof (float));  
}

/// Reads the state of the partres class from file.
/// Reconstructs cellSum, columnSum and stressTensor
void partres::readState(lexer *p, ifstream &result)
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
    double partres::maxParticlesPerCell(lexer *p, fdm &a, double d50, bool topo, bool cell)
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
    void partres::particleStressTensor(lexer *p, fdm &a, ghostcell &pgc, particles_obj &PP)
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
    void partres::particleStressTensorUpdateIJK(lexer *p, fdm &a, particles_obj &PP)
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
    void partres::updateParticleStressTensor(lexer *p, fdm &a, particles_obj &PP, int i, int j, int k)
    {
        double theta=theta_s(p,a,PP,i,j,k);
        stressTensor[IJK]=Ps*pow(theta,beta)/max(theta_crit-theta,epsilon*(1.0-theta));
    }

    /// @brief Calculate solid volume fraction for cell ( \p i , \p j , \p k )
    double partres::theta_s(lexer *p, fdm &a, particles_obj &PP, int i, int j, int k) const
    {   
        double theta = PI*pow(PP.d50,3.0)*(cellSum[IJK]+cellSumTopo[IJK])/(6.0*p->DXN[IP]*p->DYN[JP]*p->DYN[KP]);
        if(theta>1)
        theta=1;
        if(theta<0)
        theta=0;
        return theta;
    }    

    /// @brief Calculate drag force parameter
    double partres::drag_model(lexer* p, double d, double du, double dv, double dw, double thetas) const
    {
        double thetaf = 1.0-thetas;
        if(thetaf>1.0-theta_crit) // Saveguard
        thetaf=1.0-theta_crit;

        const double dU=sqrt(du*du+dv*dv+dw*dw);
        if(dU==0) // Saveguard
        return 0;

        const double Rep=dU*d*invKinVis;

        // const double Cd=24.0*(pow(thetaf,-2.65)+pow(Rep,2.0/3.0)*pow(thetaf,-1.78)/6.0)/Rep;
        const double Cd=24.0/Rep+4.0/sqrt(Rep)+0.4;
        const double Dp=Cd*3.0*drho*dU/d/4.0;

        return Dp;
    }

    double partres::sedimentation_velocity(lexer *p, double d, double du, double dv, double dw, double thetas) const
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
    void partres::particlePerCell(lexer *p, ghostcell &pgc, particles_obj &PP)
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

    void partres::make_moving(lexer *p, fdm &a, particles_obj &PP)
    {
        for(size_t n=0;n<PP.loopindex;n++)
            if(PP.Flag[n]==0)
                PP.Flag[n]=1;
    }

    void partres::erode(lexer *p, fdm &a, particles_obj &PP, sediment_fdm &s)
    {
        double shear_eff;
        double shear_crit;
        int counter=0;
        double topoDist,u,v,w,du,dv,dw;
        for(size_t n=0;n<PP.loopindex;n++)
            if(PP.Flag[n]==0)
            {
                switch (p->Q200)
                {
                    case 0:
                    {
                        // if(p->ccipol4_b(a.solid,PP.X[n],PP.Y[n],PP.Z[n])<0.6)
                        if(PP.X[n]>=p->S71 && PP.X[n]<=p->S72)
                        if(PP.Y[n]>=p->S77_xs && PP.Y[n]<=p->S77_xe)
                        {
                            shear_eff=p->ccslipol4(s.tau_eff,PP.X[n],PP.Y[n]);
                            shear_crit=p->ccslipol4(s.tau_crit,PP.X[n],PP.Y[n]);
                            if(shear_eff>shear_crit)
                                if(p->ccipol4_b(a.topo,PP.X[n],PP.Y[n],PP.Z[n])+2.0*PP.d50>0)
                                {
                                    PP.Flag[n]=1;
                                    ++counter;

                                    PP.shear_eff[n]=shear_eff;
                                    PP.shear_crit[n]=shear_crit;

                                    i=p->posc_i(PP.X[n]);
                                    j=p->posc_j(PP.Y[n]);
                                    bedChange[IJ] -= PP.PackingFactor[n];
                                }
                        }
                    }
                    break;
                    case 1:
                    {
                        relative_velocity(p,a,PP,n,du,dv,dw);
                        const double dU=sqrt(du*du+dv*dv+dw*dw);
                        const double Re_p = dU*PP.d50/(p->W2/p->W1);
                        const double Cd = drag_coefficient(Re_p);
                        const double Fd = Cd * PI/8.0 * pow(PP.d50,2) * p->W1 * pow(dU,2);
                        const double mu_s = tan(p->S81);
                        const double W = p->W1 * sqrt(p->W20*p->W20+p->W21*p->W21+p->W22*p->W22) * (p->S22/p->W1-1) * PI/6.0 *pow(PP.d50,3);
                        const double Fs = W * mu_s;
                        if(Fd > W * (mu_s*cos(s.teta[IJ])-sin(s.teta[IJ])))
                        {
                            PP.Flag[n]=1;
                            ++counter;
                            bedChange[IJ] -= PP.PackingFactor[n];
                        }
                    }
                    break;
                }
            }
        // if(counter>0)
        // cout<<"On rank "<<p->mpirank<<" were "<<counter<<" particles eroded."<<endl;
    }

    void partres::deposit(lexer *p, fdm &a, particles_obj &PP, sediment_fdm &s)
    {
        double shear_eff;
        double shear_crit;
        double topoDist,u,v,w,du,dv,dw;
        for(size_t n=0;n<PP.loopindex;n++)
            if(PP.Flag[n]==1)
            {
                switch (p->Q201)
                {
                    case 0:
                    {
                        shear_eff=p->ccslipol4(s.tau_eff,PP.X[n],PP.Y[n]);
                        shear_crit=p->ccslipol4(s.tau_crit,PP.X[n],PP.Y[n]);
                        if(shear_crit>shear_eff)
                        if(p->ccipol4_b(a.topo,PP.X[n],PP.Y[n],PP.Z[n])<2.0*PP.d50)
                        {
                            PP.Flag[n]=0;
                            PP.U[n]=0;
                            PP.V[n]=0;
                            PP.W[n]=0;

                            PP.Uf[n]=0;
                            PP.Vf[n]=0;
                            PP.Wf[n]=0;
                            PP.drag[n]=0;

                            i=p->posc_i(PP.X[n]);
                            j=p->posc_j(PP.Y[n]);
                            bedChange[IJ] += PP.PackingFactor[n];
                        }
                    }
                    break;
                    case 1:
                    {
                        relative_velocity(p,a,PP,n,du,dv,dw);
                        const double dU=sqrt(du*du+dv*dv+dw*dw);
                        const double Re_p = dU*PP.d50/(p->W2/p->W1);
                        const double Cd = drag_coefficient(Re_p);
                        const double Fd = Cd * PI/8.0 * pow(PP.d50,2) * p->W1 * pow(dU,2);
                        const double mu_s = tan(p->S81);
                        const double W = p->W1 * sqrt(p->W20*p->W20+p->W21*p->W21+p->W22*p->W22) * (p->S22/p->W1-1) * PI/6.0 *pow(PP.d50,3);
                        const double Fs = W * mu_s;
                        if(Fd < W * (mu_s*cos(s.teta[IJ])-sin(s.teta[IJ])))
                        {
                            PP.Flag[n]=0;
                            bedChange[IJ] += PP.PackingFactor[n];
                            PP.U[n]=0;
                            PP.V[n]=0;
                            PP.W[n]=0;
                        }
                    }
                    break;
                }
            }
    }

    void partres::timestep(lexer *p, ghostcell &pgc, particles_obj &PP)
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
        if(p->mpirank==0)
        cout<<"TimeStep: dx: "<<dx<<" vel: "<<sqrt(maxVelU*maxVelU+maxVelV*maxVelV+maxVelW*maxVelW)<<endl;
        double meanVel=sqrt(maxVelU*maxVelU+maxVelV*maxVelV+maxVelW*maxVelW);
        if(meanVel==0)
        meanVel=dx*p->S14/p->S13;
        p->dtsed = p->S14 * (dx/meanVel);
        p->dtsed = min(p->dtsed,p->S13);
        p->dtsed = min(p->dtsed,p->dt);
        p->dtsed = pgc.globalmin(p->dtsed);
    }

    void partres::relative_velocity(lexer *p, fdm &a, particles_obj &PP, size_t n, double &du, double &dv, double &dw)
    {
        k=p->posc_k(PP.Z[n]);
        double topoDist=p->ccipol4(a.topo,PP.X[n],PP.Y[n],PP.Z[n]);
        double u=p->ccipol1c(a.u,PP.X[n],PP.Y[n],PP.Z[n]+velDist*p->DZP[KP]-topoDist);
        double v=p->ccipol2c(a.v,PP.X[n],PP.Y[n],PP.Z[n]+velDist*p->DZP[KP]-topoDist);
        double w=p->ccipol3c(a.w,PP.X[n],PP.Y[n],PP.Z[n]+velDist*p->DZP[KP]-topoDist);

        du=u-PP.U[n];
        dv=v-PP.V[n];
        dw=w-PP.W[n];
    }

// // 10.1002/wrcr.20303

// #include <Eigen/Dense>
// #include <cmath>

// void sediment_particle::movement::doi10_1002_wrcr_20303::move(lexer *p, fdm *a, particles_obj &PP, sediment_fdm &s)
// {
//     const double kin_vis = p->W2/p->W1;
//     const double d = PP.d50;
//     const double DeltaT = p->dtsed; // time step
//     const Eigen::Vector3d g(p->W20,p->W21,p->W22);

//     bool bedLoad[p->imax*p->jmax];
//     SLICEBASELOOP
//     bedLoad[IJ] = true;


//     // Loop over all particles
//     if(PP.Flag[n] == 0)
//     {
//         i=p->posc_i(PP.X[n]);
//         j=p->posc_j(PP.Y[n]);
//         k=p->posc_k(PP.Z[n]);

//         // Define the normal vector of the plane
//         Eigen::Vector3d normal(0,0,1); // bedslope
//         normal.normalize(); // Ensure the normal vector is a unit vector
//         // Drag
//         const Eigen::Vector3d uf(a->u(i,j,k),a->v(i,j,k),0); // vel from wall law
//         const Eigen::Vector3d up(0,0,0); // bed
//         const Eigen::Vector3d du = uf - up;
//         const double Re_p = du.norm()*d/kin_vis;
//         double Cd = drag_coefficient(Re_p);
        

//         const Eigen::Vector3d Fd = Cd * PI/8.0* pow(d,2)*p->W1*du*du.norm();

//         // Lift
//         // const double G = 1; // Norm of dU/dy
//         // const double Re_shear = G * pow(d,2)/kin_vis; // shear Reynold number
//         const double tau = s.tau_eff[IJ]; // shear stress
//         const double Re_shear = sqrt(tau/p->S22)*d/kin_vis;
//         const double epsilon = sqrt(Re_shear)/Re_p;
//         const double ClSa = 12.92/PI*epsilon;
//         const double J = 0.6765 *(1.0 +tanh(2.5*log(epsilon+0.191)))*(0.667+tanh(6*(epsilon-0.32))); // approx.
//         const double Cl = 0.443*J*ClSa;
//         const Eigen::Vector3d liftDirection = normal; // orto to flow dir = to bed normal?
//         const Eigen::Vector3d Fl = Cl * PI/8.0 * pow(d,2) * p->W1 * pow(du.norm(),2) * liftDirection;

//         // Grav + buoy
        
//         const Eigen::Vector3d Fg = PI/6 * pow(d,3) * (p->S22-p->W1) * g;

//         // total
//         const Eigen::Vector3d Ftot = Fg + Fd + Fl;

//         // Calculate the normal component of the force
//         Eigen::Vector3d F_normal = (Ftot.dot(normal)) * normal;

//         // Calculate the tangential component of the force
//         Eigen::Vector3d F_tangential = Ftot - F_normal;

//         // Calculate the angle in radians
//         const double phi = std::atan(F_tangential.norm() / F_normal.norm());

//         const double moment = PP.d50/2*(F_tangential.norm()*sin(phi)+F_normal.norm()*cos(phi));

//         if(moment > 0)
//         {
//             // initial vel
//             const Eigen::Vector3d u_p_ini = Ftot.normalized() * PP.d50*(phi-phi_old[IJ])/DeltaT;
//             phi_old[IJ] = phi;

//             if(PI/3<=phi && phi <= 2.0*PI/3)
//             {
//                 // activate particles
//                 // claculate new for particles position

//                 --bedParticleNumber[IJ];
//                 PP.Flag[n] = 1;
//                 PP.U[n] = u_p_ini.x();
//                 PP.V[n] = u_p_ini.y();
//                 PP.W[n] = u_p_ini.z();
//             }
//             else if (bedLoad)
//             {
//                 // bed load transport

//                 // bedload
//                 const double psi = 3.67e-4;
//                 const double E = psi * p->S22 * sqrt((1-p->W1/p->S22)*fabs(g.norm())*d*(s.shields_eff[IJ]-s.shields_crit[IJ])); // entrainment rate [kg/m^3*m/s]=[kg/m^2/s]
//                 const double S = p->DXP[IP]+p->DYP[JP]; // pickup area [m^2]
//                 const double VdotE = E * S / p->S22; // instant volume pickup rate - used for transfer with area shift [m^3/s]
//                 const double VE = VdotE * DeltaT; // volume of particles to be transported [m^3]
//                 // Transport in flow direction
//                 const double m_p = p->S22/6*pow(PP.d50,3); // [kg]
//                 const double delta = 0.5 * F_tangential.norm()/(m_p)* pow(DeltaT,2) + u_p_ini.norm() * DeltaT;
//                 double deltaX = delta * F_tangential.normalized().x();
//                 double deltaY = delta * F_tangential.normalized().y();
//                 // Number of particles to be transported per cell and time step
//                 // const double PS = E/m_p*S; // [1/s]

//                 if(deltaX>=0)
//                 {
//                     //IP
//                     p->DXP[IP]-deltaX;
//                 }

//                 if(deltaX>=0 && deltaY>=0)
//                 {
//                     bedParticleNumber[IJ] -= (1-(p->DXP[IP]-deltaX)*(p->DYP[JP]-deltaY)/(p->DXP[IP]*p->DYP[JP]))*VE;
//                     bedParticleNumber[Ip1J] += deltaX*(p->DYP[JP]-deltaY)/(p->DXP[IP]*p->DYP[JP])*VE;
//                     bedParticleNumber[IJp1] += (p->DXP[IP]-deltaX)*deltaY/(p->DXP[IP]*p->DYP[JP])*VE;
//                     bedParticleNumber[Ip1Jp1] += deltaX*deltaY/(p->DXP[IP]*p->DYP[JP])*VE;
//                 }
//                 else if (deltaX>=0 && deltaY<0)
//                 {
//                     deltaY *= -1;
//                     bedParticleNumber[IJ] -= (1-(p->DXP[IP]-deltaX)*(p->DYP[JP]-deltaY)/(p->DXP[IP]*p->DYP[JP]))*VE;
//                     bedParticleNumber[Ip1J] += deltaX*(p->DYP[JP]-deltaY)/(p->DXP[IP]*p->DYP[JP])*VE;
//                     bedParticleNumber[IJm1] += (p->DXP[IP]-deltaX)*deltaY/(p->DXP[IP]*p->DYP[JP])*VE;
//                     bedParticleNumber[Ip1Jm1] += deltaX*deltaY/(p->DXP[IP]*p->DYP[JP])*VE;
//                 }
//                 else if (deltaX<0 && deltaY>=0)
//                 {
//                     deltaX *= -1;
//                     bedParticleNumber[IJ] -= (1-(p->DXP[IP]-deltaX)*(p->DYP[JP]-deltaY)/(p->DXP[IP]*p->DYP[JP]))*VE;
//                     bedParticleNumber[Im1J] += deltaX*(p->DYP[JP]-deltaY)/(p->DXP[IP]*p->DYP[JP])*VE;
//                     bedParticleNumber[IJp1] += (p->DXP[IP]-deltaX)*deltaY/(p->DXP[IP]*p->DYP[JP])*VE;
//                     bedParticleNumber[Im1Jp1] += deltaX*deltaY/(p->DXP[IP]*p->DYP[JP])*VE;
//                 }
//                 else if (deltaX<0 && deltaY<0)
//                 {
//                     deltaX *= -1;
//                     deltaY *= -1;
//                     bedParticleNumber[IJ] -= (1-(p->DXP[IP]-deltaX)*(p->DYP[JP]-deltaY)/(p->DXP[IP]*p->DYP[JP]))*VE;
//                     bedParticleNumber[Im1J] += deltaX*(p->DYP[JP]-deltaY)/(p->DXP[IP]*p->DYP[JP])*VE;
//                     bedParticleNumber[IJm1] += (p->DXP[IP]-deltaX)*deltaY/(p->DXP[IP]*p->DYP[JP])*VE;
//                     bedParticleNumber[Im1Jm1] += deltaX*deltaY/(p->DXP[IP]*p->DYP[JP])*VE;
//                 }
                
//                 // Place particles in cell randomly
//                 bedLoad[IJ] = false;
//             }
//         }            
//     }

//     // Transport?

//     if(PP.Flag[n] = 1)
//     {
//         const Eigen::Vector3d uf(a->u(i,j,k),a->v(i,j,k),0); // vel from wall law

//         const Eigen::Vector3d up(PP.U[n],PP.V[n],PP.Z[n]);
//         const Eigen::Vector3d du = uf - up;
//         const double Re_p = du.norm()*d/kin_vis;
//         double Cd = drag_coefficient(Re_p);

//         // suspended
//         const Eigen::Vector3d DufDt; // delta /delta t + nabla uf
//         const Eigen::Vector3d a_p = (1-p->W1/p->S22)*g + DufDt + 0.5*(DufDt-a_p) + 0.75*Cd/d*du.norm()*du + 0.75*Cl/d*pow(du.norm(),2)*normal;
//     }
// }

double partres::drag_coefficient(double Re_p) const
{
    double Cd;
    if(Re_p<0.1)
    Cd = 24.0/Re_p;
    else if(Re_p<1.0)
    Cd = 22.73/Re_p+0.0903/pow(Re_p,2),+ 3.69;
    else if(Re_p<10.0)
    Cd = 29.1667/Re_p-3.8889/pow(Re_p,2)+ 1.222;
    else if(Re_p<100.0)
    Cd = 46.5/Re_p- 116.67/pow(Re_p,2)+0.6167;
    else if(Re_p<1000.0)
    Cd = 98.33/Re_p- 2778/pow(Re_p,2) +0.3644;
    else if(Re_p<5000.0)
    Cd = 148.62/Re_p-4.75e4/pow(Re_p,2)+0.357;
    else if(Re_p<10000.0)
    Cd = -490.546/Re_p +57.87e4/pow(Re_p,2) +0.46;
    // else if(Re_p<50000.0)
    else
    Cd = -1662.5/Re_p+5.4167e6/pow(Re_p,2)+0.5191;

    return Cd;
}
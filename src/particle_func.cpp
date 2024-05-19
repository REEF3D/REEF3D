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

#include "particle_func.h"

#include "lexer.h"
#include "fdm.h"
#include "ghostcell.h"
#include "boundarycheck.h"

#define PARTICLELOOP for(size_t n=0;n<PP->loopindex;n++) if(PP->Flag[n]>INT32_MIN)

/// @brief Functions to manipulate particle containing objects
particle_func::particle_func(lexer* p) : kinVis(p->W1/p->W2), drho(p->W1/p->S22),
                            Ps(p->Q14),beta(p->Q15),epsilon(p->Q16),theta_crit(p->Q17)
{
    p->Darray(stressTensor,p->imax*p->jmax*p->kmax);
    p->Darray(cellSum,p->imax*p->jmax*p->kmax);
    p->Darray(cellSumTopo,p->imax*p->jmax*p->kmax);
    p->Darray(topoVolumeChange,p->imax*p->jmax);
}
particle_func::~particle_func()
{
    delete[] stressTensor;
    delete[] cellSum;
    delete[] cellSumTopo;
    delete[] topoVolumeChange;
}

/// @brief Removes all tracers with have left the partions boundaries
/// @param p partition object
/// @param PP tracers_obj contains tracer information
/// @return Number of removed tracers
int particle_func::remove(lexer* p, tracers_obj* PP)
{
    bool inBounds=false;
    int removed=0;
    int i,j,k;
    boundarycheck bounderies;

    PARTICLELOOP
        if(PP->Flag[n]>0)
        {
            i = p->posc_i(PP->X[n]);
            j = p->posc_j(PP->Y[n]);
            k = p->posc_k(PP->Z[n]);

            inBounds=bounderies.minboundcheck(p,i,j,k,1);
            if (inBounds)
                inBounds=bounderies.maxboundcheck(p,i,j,k,1);

			// remove out of bounds particles
            if(!inBounds)
            {
                cellSum[IJK]--;
                PP->erase(n);
                removed++;
            }
        }
    
    return removed;
}

/// \copydoc particle_func::remove
int particle_func::remove(lexer* p, particles_obj* PP)
{
    bool inBounds=false;
    int removed=0;
    int i,j,k;
    boundarycheck bounderies;

    PARTICLELOOP
        if(PP->Flag[n]>0)
        {
            i = p->posc_i(PP->X[n]);
            j = p->posc_j(PP->Y[n]);
            k = p->posc_k(PP->Z[n]);

            inBounds=bounderies.minboundcheck(p,i,j,k,1);
            if (inBounds)
                inBounds=bounderies.maxboundcheck(p,i,j,k,1);

			// remove out of bounds particles
            if(!inBounds)
            {
                cellSum[IJK]-=PP->PackingFactor[n];
                PP->erase(n);
                removed++;
            }
        }
    
    return removed;
}

int particle_func::solid_clean(lexer* p, particles_obj* PP)
{
    int removed = 0;
    // PARTICLELOOP
    // if(PP->Flag[n]>0)
    // {
    //     i = p->posc_i(PP->X[n]);
    //     j = p->posc_j(PP->Y[n]);
    //     k = p->posc_k(PP->Z[n]);
    //     cellSum[IJK]-=PP->PackingFactor[n];
    //     PP->erase(n);
    //     removed++;
    // }
    return removed;
}

/// @brief Transfer function
/** Is responsible to transfer tracers to one of the surrounding partitions*/
/// @param p partition object
/// @param pgc ghostcell object
/// @param PP tracers_obj contains tracer information
/// @param maxcount maximum number of tracers which could be transfered
/// @return number of send of tracers
int particle_func::transfer(lexer* p, ghostcell* pgc, tracers_obj* PP, int maxcount)
{
    int xchange=0;

    tracers_obj Send[6]={tracers_obj(maxcount),tracers_obj(maxcount),tracers_obj(maxcount),tracers_obj(maxcount),tracers_obj(maxcount),tracers_obj(maxcount)};
    tracers_obj Recv[6]={tracers_obj(maxcount),tracers_obj(maxcount),tracers_obj(maxcount),tracers_obj(maxcount),tracers_obj(maxcount),tracers_obj(maxcount)};

    int i,j,k;

    PARTICLELOOP
        if(PP->Flag[n]==1)
        {
            i = p->posc_i(PP->X[n]);
            j = p->posc_j(PP->Y[n]);
            k = p->posc_k(PP->Z[n]);


            if(p->flag5[IJK]<0 && p->flag5[IJK]>-10)
            {
                switch (p->flag5[IJK])
                {
                    case -1:
                    {
                        Send[0].add(PP->X[n],PP->Y[n],PP->Z[n],1);
                        break;
                    }

                    case -2:
                    {
                        Send[1].add(PP->X[n],PP->Y[n],PP->Z[n],1);
                        break;
                    }

                    case -3:
                    {
                        Send[2].add(PP->X[n],PP->Y[n],PP->Z[n],1);
                        break;
                    }

                    case -4:
                    {
                        Send[3].add(PP->X[n],PP->Y[n],PP->Z[n],1);
                        break;
                    }

                    case -5:
                    {
                        Send[4].add(PP->X[n],PP->Y[n],PP->Z[n],1);
                        break;
                    }

                    case -6:
                    {
                        Send[5].add(PP->X[n],PP->Y[n],PP->Z[n],1);
                        break;
                    }
                }
                PP->erase(n);
                cellSum[IJK]--;
                ++xchange;
            }
        }

    pgc->para_tracersobj(p,Send,Recv);


    size_t sum=PP->size;
    for(int n=0;n<6;n++)
        sum += Recv[n].size;
    if(sum>PP->capacity)
        PP->reserve(sum);

    for(int n=0;n<6;n++)
    {   
        for(size_t n=0;n<Recv[n].loopindex;n++)
        {
            i = p->posc_i(Recv[n].X[n]);
            j = p->posc_j(Recv[n].Y[n]);
            k = p->posc_k(Recv[n].Z[n]);
            cellSum[IJK]++;
        }
        PP->add_obj(&Recv[n]);
    }

    return xchange;
}


/// @brief Transfer function
/// Is responsible to transfer particles to one of the surrounding partitions
/// @param p partition object
/// @param pgc ghostcell object
/// @param PP particles_obj contains particle information
/// @param maxcount maximum number of particles which could be transfered
/// @return number of send of particles
int particle_func::transfer(lexer* p, ghostcell* pgc, particles_obj* PP, int maxcount)
{
    int xchange=0;

    particles_obj Send[6]={particles_obj(maxcount,PP->d50,PP->density,1),particles_obj(maxcount,PP->d50,PP->density,1),particles_obj(maxcount,PP->d50,PP->density,1),
    particles_obj(maxcount,PP->d50,PP->density,1),particles_obj(maxcount,PP->d50,PP->density,1),particles_obj(maxcount,PP->d50,PP->density,1)};
    particles_obj Recv[6]={particles_obj(maxcount,PP->d50,PP->density,1),particles_obj(maxcount,PP->d50,PP->density,1),particles_obj(maxcount,PP->d50,PP->density,1),
    particles_obj(maxcount,PP->d50,PP->density,1),particles_obj(maxcount,PP->d50,PP->density,1),particles_obj(maxcount,PP->d50,PP->density,1)};;

    int i,j,k;

    PARTICLELOOP
        if(PP->Flag[n]==1)
        {
            i = p->posc_i(PP->X[n]);
            j = p->posc_j(PP->Y[n]);
            k = p->posc_k(PP->Z[n]);


            if(p->flag5[IJK]<0 && p->flag5[IJK]>-10)
            {
                switch (p->flag5[IJK])
                {
                    case -1:
                    {
                        Send[0].add(PP->X[n],PP->Y[n],PP->Z[n],1,PP->U[n],PP->V[n],PP->W[n],PP->PackingFactor[n]);
                        break;
                    }

                    case -2:
                    {
                        Send[1].add(PP->X[n],PP->Y[n],PP->Z[n],1,PP->U[n],PP->V[n],PP->W[n],PP->PackingFactor[n]);
                        break;
                    }

                    case -3:
                    {
                        Send[2].add(PP->X[n],PP->Y[n],PP->Z[n],1,PP->U[n],PP->V[n],PP->W[n],PP->PackingFactor[n]);
                        break;
                    }

                    case -4:
                    {
                        Send[3].add(PP->X[n],PP->Y[n],PP->Z[n],1,PP->U[n],PP->V[n],PP->W[n],PP->PackingFactor[n]);
                        break;
                    }

                    case -5:
                    {
                        Send[4].add(PP->X[n],PP->Y[n],PP->Z[n],1,PP->U[n],PP->V[n],PP->W[n],PP->PackingFactor[n]);
                        break;
                    }

                    case -6:
                    {
                        Send[5].add(PP->X[n],PP->Y[n],PP->Z[n],1,PP->U[n],PP->V[n],PP->W[n],PP->PackingFactor[n]);
                        break;
                    }
                }
                cellSum[IJK]-=PP->PackingFactor[n];
                PP->erase(n);
                ++xchange;
            }
        }

    pgc->para_tracersobj(p,Send,Recv);


    size_t sum=PP->size;
    for(int n=0;n<6;n++)
        sum += Recv[n].size;
    if(sum>PP->capacity)
        PP->reserve(sum);

    for(int n=0;n<6;n++)
    {
        for(size_t n=0;n<Recv[n].loopindex;n++)
        {
            i = p->posc_i(Recv[n].X[n]);
            j = p->posc_j(Recv[n].Y[n]);
            k = p->posc_k(Recv[n].Z[n]);
            cellSum[IJK]+=Recv[n].PackingFactor[n];
        }
        PP->add_obj(&Recv[n]);
    }

    return xchange;
}

/// @brief Particle Reynolds number
/// Calculates particle reynolds number for particle[ \p index ]
/// @return Local particle Reynolds number
double particle_func::reynolds(lexer* p,fdm* a, particles_obj* PP, int index)
{
    const double u=p->ccipol1(a->u,PP->X[index],PP->Y[index],PP->Z[index]);
    const double v=p->ccipol2(a->v,PP->X[index],PP->Y[index],PP->Z[index]);
    const double w=p->ccipol3(a->w,PP->X[index],PP->Y[index],PP->Z[index]);

    const double mean_vel=sqrt(u*u+v*v+w*w);

    const double Re=mean_vel*PP->d50/p->W2;
    // Change to particle diameter once implemented

    return Re;
}

/// @brief Settling velocity
/// Calculates settling velocity using drag coefficent
/// g is assumed to act only in negative z direction
/// @return settling velocity of particle \p index
double particle_func::settling_vel(lexer* p,fdm* a, particles_obj* PP, int index)
{
    return sqrt(4.0/3.0*(PP->density/p->W1-1.0)*fabs(p->W22)*PP->d50/drag_coefficient(p,a,PP,index));
}

/// @brief Drag coefficent Cd
/// Calculates drag coefficent from particle Reynolds number based on emperical 
/// @return Cd
double particle_func::drag_coefficient(lexer* p,fdm* a, particles_obj* PP, int index)
{
    const double Re=reynolds(p,a,PP,index);
    if(Re<0.5)
        return 24.0/Re;
    if(Re<1000.0)
        return 18.5*pow(Re,-0.6);
    if(Re<200000.0)
        return 0.44;
    return NAN;


    /// Andrews and O’Rourke (1996)
    const double Rep=0;
    const double theta_s=0;// vol sol/vol cell
    const double theta_f=1-theta_s;
    return 24.0*(pow(theta_f,-2.65)+pow(Rep,2.0/3.0)*pow(theta_f,-1.78)/6.0)/Rep;
}

/// @brief Set flag of particle to \p minflag aka stationary
void particle_func::make_stationary(lexer* p, fdm* a, tracers_obj* PP, int minflag)
{
    int i,j;
    PARTICLELOOP
        if (p->ccipol4_b(a->topo,PP->X[n],PP->Y[n],PP->Z[n])<0)
            PP->Flag[n]=minflag;
}

/// @brief Set flag of particle to \p minflag aka stationary and remove velocities
void particle_func::make_stationary(lexer* p, fdm* a, particles_obj* PP)
{
    int i,j;
    PARTICLELOOP
    if(PP->Flag[n]>0)
    {
        i=p->posc_i(PP->X[n]);
        j=p->posc_j(PP->Y[n]);
        if (p->ccipol4_b(a->topo,PP->X[n],PP->Y[n],PP->Z[n])-volume(PP,n)/(p->DXN[IP]*p->DYN[JP])<=0||p->ccipol4_b(a->solid,PP->X[n],PP->Y[n],PP->Z[n])<=0)
        {
            PP->Flag[n]=0;
            if(p->count!=0)
            {
                topoVolumeChange[IJ]+=volume(PP,n);
            }
            if(PP->entries>PP->tracers_obj::entries)
            {
                PP->U[n]=0;
                PP->V[n]=0;
                PP->W[n]=0;
            }
            if(p->ccipol4_b(a->solid,PP->X[n],PP->Y[n],PP->Z[n])<=0)
            {
                PP->Flag[n]=-1;
            }
        }
    }
}

/// @brief Calculates volume of particles
/// @return Volume of partice with \p index
double particle_func::volume(particles_obj* PP, int index)
{
    return PI*pow(PP->d50,3.0)*PP->PackingFactor[index]/6.0;
}

/// @brief Cleanup container \note To implement
void particle_func::cleanup(lexer* p, fdm* a, tracers_obj* PP, int max)
{

}

/// @brief Topo volume in cell div. by particle volume
/// Uses i,j&k from increment to pass cell identifier
/// @param d50 Sauter diameter of particles
/// @return Ceil of number of particles in cell IJK
double particle_func::maxParticlesPerCell(lexer* p, fdm* a, double d50, bool topo, bool cell)
{   
    // return 6.0*p->DXN[IP]*p->DYN[JP]*(p->DZN[KP]-(topo?(0.5*p->DZN[KP]+a->topo(i,j,k)):0.0))/((PI*pow(d50,3.0)));
    double DZN=topo?0:p->DZN[KP];

    if(topo)
    {
        if (a->topo(i,j,k)<=-0.5*p->DZN[KP]+1.0e-13)
        DZN=p->DZN[KP];
        else if(a->topo(i,j,k)<0.5*p->DZN[KP] -5.0e-18)
        DZN=(p->DZN[KP]*0.5 + a->topo(i,j,k));
    }

    return 6.0*p->DXN[IP]*p->DYN[JP]*DZN*(1.0+(cell?0:-a->porosity(i,j,k)))/(PI*pow(d50,3.0));
}

/// @brief Max particles in plane
int particle_func::maxParticlesPerXY(lexer* p, fdm* a, double d50)
{
    return ceil(p->DXN[IP]*p->DYN[JP]*(1-p->S24)*6.0/(PI*pow(d50,2.0)));
}

/// @brief Count particles in cell
void particle_func::particlesPerCell(lexer* p, fdm* a, ghostcell* pgc, particles_obj* PP)
{
    PLAINLOOP
        cellSum[IJK]=0;
    
    PARTICLELOOP
    {
        i=p->posc_i(PP->X[n]);
        j=p->posc_j(PP->Y[n]);
        k=p->posc_k(PP->Z[n]);
        cellSum[IJK]+=PP->PackingFactor[n];
    }

    int gcv=11;
    double *f = cellSum;
    double starttime=pgc->timer();
    pgc->gcparaxV(p, f, gcv);
    pgc->gcparacoxV(p, f, gcv);
    pgc->gcparacoxV(p, f, gcv);
    pgc->gcparacoxV(p, f, gcv);
    p->xtime+=pgc->timer()-starttime;
    f = cellSumTopo;
    starttime=pgc->timer();
    pgc->gcparaxV(p, f, gcv);
    pgc->gcparacoxV(p, f, gcv);
    pgc->gcparacoxV(p, f, gcv);
    pgc->gcparacoxV(p, f, gcv);
    p->xtime+=pgc->timer()-starttime;
}

/// @brief Calculate complete intra-particle stress trensor
void particle_func::particleStressTensor(lexer* p, fdm* a, ghostcell* pgc, particles_obj* PP)
{
    double theta;
    int i,j,k;

    PLAINLOOP
    {
        updateParticleStressTensor(p,a,PP,i,j,k);
    }
    pgc->start4V(p,stressTensor,10);
}

/// @brief Calculate intra-particle stress trensor for cells around (`increment::i`,`increment::j`,`increment::k`)
void particle_func::particleStressTensorUpdateIJK(lexer* p, fdm* a, particles_obj* PP)
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
void particle_func::updateParticleStressTensor(lexer* p, fdm* a, particles_obj* PP, int i, int j, int k)
{
    double theta=theta_s(p,a,PP,i,j,k);
    stressTensor[IJK]=Ps*pow(theta,beta)/max(theta_crit-theta,epsilon*(1.0-theta));
}

/// @brief Calculate solid volume fraction for cell ( \p i , \p j , \p k )
double particle_func::theta_s(lexer* p, fdm* a, particles_obj* PP, int i, int j, int k)
{   
    double theta = PI*pow(PP->d50,3.0)*(cellSum[IJK]+cellSumTopo[IJK])/(6.0*p->DXN[IP]*p->DYN[JP]*p->DYN[KP]);
    if(theta>1)
    theta=1;
    return theta;
}    

/// @brief Calculate drag force parameter
double particle_func::drag_model(lexer* p, double d, double du, double dv, double dw, double thetas) const
{
    double thetaf = 1.0-thetas;
    if(thetaf>1.0-theta_crit) // Saveguard
    thetaf=1.0-theta_crit;

    const double dU=sqrt(du*du+dv*dv+dw*dw);
    if(dU==0) // Saveguard
    return 0;

    const double Rep=dU*d*kinVis;

    const double Cd=24.0*(pow(thetaf,-2.65)+pow(Rep,2.0/3.0)*pow(thetaf,-1.78)/6.0)/Rep;
    const double Dp=Cd*3.0*drho*dU/d/4.0;

    if(Dp!=Dp)
    cout<<thetaf<<","<<dU<<","<<Rep<<","<<Cd<<"|"<<(dU==0)<<endl;

    return Dp;
}

/// @brief Determine wether or not a particle would move and set its flag accordingly
void particle_func::make_moving(lexer* p, fdm* a, particles_obj* PP)
{
    double RKu,RKv,RKw;
    double u,v,w;
    double du1, dv1, dw1;
    /// @brief Difference between flowfield and particle velocity
    double du, dv, dw;
    double Dp, thetas;
    double pressureDivX=0, pressureDivY=0, pressureDivZ=0;
    double stressDivX=0,stressDivY=0,stressDivZ=0;
    double netBuoyX=(1.0-drho)*p->W20, netBuoyY=(1.0-drho)*p->W21, netBuoyZ=(1.0-drho)*p->W22;

    PARTICLELOOP
    if(PP->Flag[n]==0)
    {
        i=p->posc_i(PP->X[n]);
        j=p->posc_j(PP->Y[n]);
        k=p->posc_k(PP->Z[n]);

        thetas=theta_s(p,a,PP,i,j,k);

        u=p->ccipol1(a->u,PP->X[n],PP->Y[n],PP->Z[n]);
        v=p->ccipol1(a->v,PP->X[n],PP->Y[n],PP->Z[n]);
        w=p->ccipol1(a->w,PP->X[n],PP->Y[n],PP->Z[n]);

        stressDivX = (stressTensor[Ip1JK] - stressTensor[IJK])/(p->DXN[IP]);
        stressDivY = (0.5*(stressTensor[IJp1K]+stressTensor[Ip1Jp1K]) - 0.5*(stressTensor[IJm1K]+stressTensor[Ip1Jm1K]))/(p->DYN[JM1]+p->DYN[JP]);
        stressDivZ = (0.5*(stressTensor[IJKp1]+stressTensor[Ip1JKp1]) - 0.5*(stressTensor[IJKm1]+stressTensor[Ip1JKm1]))/(p->DYN[KM1]+p->DYN[KP]);

        pressureDivX = (a->press(i+1,j,k) - a->press(i,j,k))/(p->DXN[IP]);
        pressureDivY = (0.5*(a->press(i,j+1,k)+a->press(i+1,j+1,k)) - 0.5*(a->press(i,j-1,k)+a->press(i+1,j-1,k)))/(p->DYN[JM1]+p->DYN[JP]);
        pressureDivZ = (0.5*(a->press(i,j,k+1)+a->press(i+1,j,k+1)) - 0.5*(a->press(i,j,k-1)+a->press(i+1,j,k-1)))/(p->DYN[KM1]+p->DYN[KP]);

        // RK3 step 1
        du=u-PP->U[n];
        dv=v-PP->V[n];
        dw=w-PP->W[n];

        Dp=drag_model(p,PP->d50,du,dv,dw,thetas);

        du1=Dp*du+netBuoyX-pressureDivX/p->S22-stressDivX/((1-thetas)*p->S22);
        dv1=Dp*dv+netBuoyY-pressureDivY/p->S22-stressDivY/((1-thetas)*p->S22);
        dw1=Dp*dw+netBuoyZ-pressureDivZ/p->S22-stressDivZ/((1-thetas)*p->S22);

        double tolerance=0.0;
        if (fabs(du1)>tolerance||fabs(dv1)>tolerance||dw1>tolerance)
        {
            PP->Flag[n]=1;
            if(p->count!=0)
            {
                topoVolumeChange[IJ]-=volume(PP,n);
            }
        }
    }
}

/// @brief All debug code
void particle_func::debug(lexer* p, fdm* a, ghostcell* pgc, tracers_obj* PP)
{
    PLAINLOOP
    a->test(i,j,k)=cellSum[IJK];
}

/// @brief Moves stuck stationary particles out of `fdm::solid`
void particle_func::fixPos(lexer* p, fdm* a, particles_obj* PP)
{
    double solid;
    PARTICLELOOP
    if(PP->Flag[n]==1&&p->ccipol4_b(a->solid,PP->X[n],PP->Y[n],PP->Z[n])<0)
    {
        i=p->posc_i(PP->X[n]);
        j=p->posc_j(PP->Y[n]);
        k=p->posc_k(PP->Z[n]);
        cellSum[IJK]-=PP->PackingFactor[n];

        if(p->flag1[Ip1JK]==SOLID_FLAG&&p->flag1[IJK]==WATER_FLAG)
        {
            while(!(p->flag1[Im1JK]==SOLID_FLAG&&p->flag1[IJK]==WATER_FLAG))
            i++;
            solid=p->ccipol4_b(a->solid,p->XN[IP],PP->Y[n],PP->Z[n]);
            PP->X[n]=p->XN[IP]-solid;
        }
        else
        {
            solid=p->ccipol4_b(a->solid,PP->X[n],PP->Y[n],PP->Z[n]);
            PP->X[n]-=solid;
        }
        i=p->posc_i(PP->X[n]);
        cellSum[IJK]+=PP->PackingFactor[n];
    }
}
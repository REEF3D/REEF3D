#include"VOF_PLIC.h"
#include"gradient.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"convection.h"
#include"solver.h"
#include"ghostcell.h"
#include"ioflow.h"
#include"fluid_update_vof.h"
#include"heat.h"
#include"hires.h"
#include"weno_hj.h"
#include"hric.h"

void VOF_PLIC::advectPhi_Bonn
(
    fdm* a,
    lexer* p,
    int nSweep,
    int sweep
)
{
    double Gp, Gm, uip, uim, dtdxi,
    if(sweep==0)
    {
        uip=a->u(i,j,k);
        uim=a->u(i-1,j,k);
        dtdxi=p->dt/p->DXN[IP];
        if(uip>=0.0)
            Gp=uip*(phistep(i,j,k)+0.5*p->DXP[IP]*(1.0-uip*(p->dt/p->DXP[IP]))*(phistep(i+1,j,k)-phistep(i-1,j,k))/(p->DXP[IP]+p->DXP[IM1]));
        else
            Gp=uip*(phistep(i,j,k)+0.5*p->DXP[IP]*(1.0-uip*(p->dt/p->DXP[IP]))*(phistep(i+2,j,k)-phistep(i,j,k))/(p->DXN[IP]+p->DXP[IP1]));
            
        if(uim<=0.0)
            Gm=uim*(phistep(i,j,k)+0.5*p->DXP[IM1]*(-1.0-uim*(p->dt/p->DXP[IM1]))*(phistep(i+1,j,k)-phistep(i-1,j,k))/(p->DXP[IP]+p->DXP[IM1]));
        else
            Gm=uim*(phistep(i,j,k)+0.5*p->DXP[IM1]*(-1.0-uim*(p->dt/p->DXP[IM1]))*(phistep(i,j,k)-phistep(i-2,j,k))/(p->DXP[IM2]+p->DXP[IM1]));
    }
    else if(sweep==1)
    {
        uip=a->v(i,j,k);
        uim=a->v(i-1,j,k);
        dtdxi=p->dt/p->DYN[JP];
        if(uip>=0.0)
            Gp=uip*(phistep(i,j,k)+0.5*p->DYP[JP]*(1.0-uip*(p->dt/p->DYP[JP]))*(phistep(i,j+1,k)-phistep(i,j-1,k))/(p->DYP[JP]+p->DYP[JM1]));
        else
            Gp=uip*(phistep(i,j,k)+0.5*p->DYP[JP]*(1.0-uip*(p->dt/p->DYP[JP]))*(phistep(i,j+2,k)-phistep(i,j,k))/(p->DYP[JP]+p->DYP[JP1]));
            
        if(uim<=0.0)
            Gm=uim*(phistep(i,j,k)+0.5*p->DYP[JM1]*(-1.0-uim*(p->dt/p->DYP[JM1]))*(phistep(i,j+1,k)-phistep(i,j-1,k))/(p->DYP[JP]+p->DYP[JM1]));
        else
            Gm=uim*(phistep(i,j,k)+0.5*p->DYP[JM1]*(-1.0-uim*(p->dt/p->DYP[JM1]))*(phistep(i,j,k)-phistep(i,j-1,k))/(p->DYP[JM1]+p->DYP[JM2]));
            
    }
    else
    {
        uip=a->w(i,j,k);
        uim=a->w(i,j,k-1);
        dtdxi=p->dt/p->DZN[KP];
        if(uip>=0.0)
            Gp=uip*(phistep(i,j,k)+0.5*p->DZP[KP]*(1.0-uip*(p->dt/p->DZP[IP]))*(phistep(i,j,k+1)-phistep(i,j,k-1))/(p->DZP[IP]+p->DZP[IM1]));
        else
            Gp=uip*(phistep(i,j,k)+0.5*p->DZP[KP]*(1.0-uip*(p->dt/p->DZP[IP]))*(phistep(i,j,k+2)-phistep(i,j,k))/(p->DZP[IP1]+p->DZP[IP]));
            
        if(uim<=0.0)
            Gm=uim*(phistep(i,j,k)+0.5*p->DZP[IM1]*(-1.0-uim*(p->dt/p->DZP[IM1]))*(phistep(i,j,k+1)-phistep(i,j,k-1))/(p->DZP[IP]+p->DZP[IM1]));
        else
            Gm=uim*(phistep(i,j,k)+0.5*p->DZP[IM1]*(-1.0-uim*(p->dt/p->DZP[IM1]))*(phistep(i,j,k)-phistep(i,j,k-2))/(p->DZP[IM2]+p->DZP[IM1]));
    }
    
    if(nSweep==0)
        phiS0(i,j,k)=(a->phi(i,j,k)-dtdxi*(Gp-Gm))/(1.0-dtdxi*(uip-uim));
        
    else if(nSweep==1)
        phiS1(i,j,k)=phiS0(i,j,k)*(1.0+dtdxi*(uip-uim))-dtdxi*(Gp-Gm);
    
    else
        phiS2(i,j,k)=phiS1(i,j,k)+phiS0(i,j,k)*dtdxi*(uip-uim)-dtdxi*(Gp-Gm);
}
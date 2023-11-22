
#include"iowave.h"
#include"lexer.h"
#include"fdm_ptf.h"
#include"ghostcell.h"

void iowave::dirichlet_wavegen_ptf(lexer *p, fdm_ptf* e, ghostcell* pgc, field& u, field& v, field& w)
{
        count=0;
		for(n=0;n<p->gcin_count;n++)
		{
		i=p->gcin[n][0];
		j=p->gcin[n][1];
		k=p->gcin[n][2];		

        uvel=uval[count]*ramp(p);
        vvel=vval[count]*ramp(p);
        wvel=wval[count]*ramp(p);
            
			if(e->phi(i-1,j,k)>=0.0)
			{
			u(i-1,j,k)=uvel+p->Ui;
			u(i-2,j,k)=uvel+p->Ui;
			u(i-3,j,k)=uvel+p->Ui;
            
            v(i-1,j,k)=vvel;
			v(i-2,j,k)=vvel;
			v(i-3,j,k)=vvel;
			
			w(i-1,j,k)=wvel;
			w(i-2,j,k)=wvel;
			w(i-3,j,k)=wvel;
			}

			if(e->phi(i-1,j,k)<0.0 && e->phi(i-1,j,k)>=-p->F45*p->DZP[KP])
			{
			fac= p->B122*(1.0 - fabs(e->phi(i-1,j,k))/(p->F45*p->DZP[KP]));
            
			u(i-1,j,k)=uvel*fac + p->Ui;
			u(i-2,j,k)=uvel*fac + p->Ui;
			u(i-3,j,k)=uvel*fac + p->Ui;
            
            v(i-1,j,k)=vvel*fac;
			v(i-2,j,k)=vvel*fac;
			v(i-3,j,k)=vvel*fac;
            
			w(i-1,j,k)=wvel*fac;
			w(i-2,j,k)=wvel*fac;
			w(i-3,j,k)=wvel*fac;
			}

			if(e->phi(i-1,j,k)<-p->F45*p->DXM)
			{
			u(i-1,j,k)=0.0 + p->Ui;
			u(i-2,j,k)=0.0 + p->Ui;
			u(i-3,j,k)=0.0 + p->Ui;
            
            if(p->W50_air==1)
            {
            u(i-1,j,k)+=p->W50;
            u(i-2,j,k)+=p->W50;
            u(i-3,j,k)+=p->W50;
            }
            
            v(i-1,j,k)=0.0;
			v(i-2,j,k)=0.0;
			v(i-3,j,k)=0.0;

			w(i-1,j,k)=0.0;
			w(i-2,j,k)=0.0;
			w(i-3,j,k)=0.0;
			}
        ++count;
		}
        
        
        if(p->B98==3||p->B98==4||p->B99==3||p->B99==4||p->B99==5)
		{
		for(int q=0;q<4;++q)
		for(n=0;n<p->gcin_count;++n)
		{
		i=p->gcin[n][0]+q;
		j=p->gcin[n][1];
		k=p->gcin[n][2];
        }
        }

        
    // PTF
   if(p->A10==4)
    {
        for(n=0;n<p->gcslin_count;n++)
        {
        i=p->gcslin[n][0];
        j=p->gcslin[n][1];
        k=e->etaloc(i,j);
        {
        e->Fifsf(i-1,j) = e->Fifsf(i,j) - u(i-1,j,k)*1.0*p->DXP[IM1];
        e->Fifsf(i-2,j) = e->Fifsf(i,j) - u(i-1,j,k)*2.0*p->DXP[IM1];
        e->Fifsf(i-3,j) = e->Fifsf(i,j) - u(i-1,j,k)*3.0*p->DXP[IM1];
        }
        }
    }
}
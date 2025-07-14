/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

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

#include"iowave.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

void iowave::u_relax(lexer *p, fdm *a, ghostcell *pgc, field& uvel)
{
    starttime=pgc->timer();
    
    count=0;
    
    ULOOP
    {
        dg = distgen(p);    
        db = distbeach(p);
        
        phival = 0.5*(a->phi(i,j,k)+a->phi(i-1,j,k));

        if(phival>=-psi)
		{
		H=1.0;
		G=1.0;
		}

		if(phival<-epsi)
		{
		H=0.0;
		G=0.0;
		}

		if(phival>=-epsi && phival<-psi)
		{
		H=0.5*(1.0 + phival/epsi + (1.0/PI)*sin((PI*phival)/epsi));
		G=H;
		}
        
        if(phival>=0.0)
        {
        if(p->pos_z()<=p->phimean)
        z=-(fabs(p->phimean-p->pos_z()));
		
		if(p->pos_z()>p->phimean)
        z=(fabs(p->phimean-p->pos_z()));
        }
        
        if(phival<0.0)
        z = 0.5*(eta(i,j)+eta(i+1,j));
		
		// Wave Generation
		if(p->B98==2 && u_switch==1)
        {
            // Zone 1
            if(dg<1.0e20)
            {
            uvel(i,j,k) = (1.0-relax1_wg(i,j))*ramp(p)*uval[count] * H + relax1_wg(i,j)*H*uvel(i,j,k) + (1.0-G)*uvel(i,j,k);
            uvel(i,j,k) += (1.0-relax1_wg(i,j))*ramp(p)*p->W50*(1.0-H) + relax1_wg(i,j)*(1.0-H) *p->W50;
            ++count;
            }
		}
		
		// Numerical Beach
        if(p->B99==1||p->B99==2||beach_relax==1)
		{
            // Zone 2
            if(db<1.0e19)
            uvel(i,j,k) = relax1_nb(i,j)*uvel(i,j,k);
        }
    }
    
    p->wavecalctime+=pgc->timer()-starttime;
}

void iowave::v_relax(lexer *p, fdm *a, ghostcell *pgc, field& vvel)
{
    starttime=pgc->timer();
    
    count=0;
    VLOOP
    {
        dg = distgen(p);    
        db = distbeach(p);
        
        phival = 0.5*(a->phi(i,j,k)+a->phi(i,j-1,k));

        if(phival>=-psi)
		{
		H=1.0;
		G=1.0;
		}

		if(phival<-epsi)
		{
		H=0.0;
		G=0.0;
		}

		if(phival>=-epsi && phival<-psi)
		{
		H=0.5*(1.0 + phival/epsi + (1.0/PI)*sin((PI*phival)/epsi));
		G=H;
		}
		
        if(phival>=0.0)
        {
        if(p->pos_z()<=p->phimean)
        z=-(fabs(p->phimean-p->pos_z()));
		
		if(p->pos_z()>p->phimean)
        z=(fabs(p->phimean-p->pos_z()));
        }
        
        if(phival<0.0)
        z = 0.5*(eta(i,j)+eta(i,j+1));

		// Wave Generation
		if(p->B98==2 && v_switch==1)
        {
            // Zone 1
            if(dg<1.0e20)
            {
            vvel(i,j,k) = (1.0-relax2_wg(i,j))*ramp(p)*vval[count] * H + relax2_wg(i,j)*H*vvel(i,j,k) + (1.0-G)*vvel(i,j,k);
            ++count;
            }
		}
		
		// Numerical Beach
		if(p->B99==1||p->B99==2||beach_relax==1)
		{	
            // Zone 2
            if(db<1.0e20)
            vvel(i,j,k) = relax2_nb(i,j)*vvel(i,j,k);
        }
    }
    
    p->wavecalctime+=pgc->timer()-starttime;
}

void iowave::w_relax(lexer *p, fdm *a, ghostcell *pgc, field& wvel)
{
    starttime=pgc->timer();
    
    count=0;
    WLOOP
    {
        dg = distgen(p);    
        db = distbeach(p);
        
        phival = 0.5*(a->phi(i,j,k)+a->phi(i,j,k-1));

        if(phival>=-psi)
		{
		H=1.0;
		G=1.0;
		}

		if(phival<-epsi)
		{
		H=0.0;
		G=0.0;
		}
		

		if(phival>=-epsi && phival<-psi)
		{
		H=0.5*(1.0 + phival/epsi + (1.0/PI)*sin((PI*phival)/epsi));
		G=H;
		}
		
        if(phival>=0.0)
        {
        if(p->pos_z()<=p->phimean)
        z=-(fabs(p->phimean-p->pos3_z()));
		
		if(p->pos_z()>p->phimean)
        z=(fabs(p->phimean-p->pos3_z()));
        }
        
        if(phival<0.0)
        z = eta(i,j);

		// Wave Generation
		if(p->B98==2 && w_switch==1)
        {
            // Zone 1
            if(dg<1.0e20)
            {
            wvel(i,j,k) = (1.0-relax4_wg(i,j)) * ramp(p)* wval[count] * H + relax4_wg(i,j)*H*wvel(i,j,k) + (1.0-G)*wvel(i,j,k);
            ++count;
            }

		}
		
		// Numerical Beach
        if(p->B99==1||p->B99==2||beach_relax==1)
		{
            // Zone 2
            if(db<1.0e20)
            wvel(i,j,k) = relax4_nb(i,j)*wvel(i,j,k);
        }
    }
    
    p->wavecalctime+=pgc->timer()-starttime;
}



void iowave::p_relax(lexer *p, fdm *a, ghostcell *pgc, field& press)
{
    starttime=pgc->timer();
    
    LOOP
    {
        dg = distgen(p);
        db = distbeach(p);
        
        // Numerical Beach
        if(p->B99==1||p->B99==2||beach_relax==1)
        {            
            // Zone 2
            if(db<1.0e20)
            press(i,j,k) = (1.0-relax4_nb(i,j))*((p->phimean - p->pos_z())*a->ro(i,j,k)*fabs(p->W22)) + relax4_nb(i,j)*press(i,j,k);
        }
    }

p->wavecalctime+=pgc->timer()-starttime;		
}

void iowave::phi_relax(lexer *p, ghostcell *pgc, field& f)
{
    starttime=pgc->timer();
    
    count=0;
    LOOP
    {
        dg = distgen(p);    
        db = distbeach(p);

        if(p->pos_z()<=p->phimean)
        z=-(fabs(p->phimean-p->pos_z()));
            
        if(p->pos_z()>p->phimean)
        z=(fabs(p->phimean-p->pos_z()));
            
        // Wave Generation
        if(p->B98==2 && h_switch==1)
        {
            // Zone 1
            if(dg<1.0e20)
            {
            f(i,j,k) = (1.0-relax4_wg(i,j))*ramp(p) * lsval[count] + relax4_wg(i,j)*f(i,j,k);
            ++count;
            }
        }
            
        // Numerical Beach    
        if(p->B99==2)
        {
            // Zone 2
            if(db<1.0e20)
            f(i,j,k) = (1.0-relax4_nb(i,j)) * (p->phimean-p->pos_z()) + relax4_nb(i,j)*f(i,j,k);
        }
    }
    
    p->wavecalctime+=pgc->timer()-starttime;
}

void iowave::vof_relax(lexer *p, fdm* a, ghostcell *pgc, field& f)
{
    
    starttime=pgc->timer();
    SLICELOOP4
    {
        vofheight(i,j)=0.0;
        genheight(i,j)=0.0;
        KLOOP
        {   
            vofheight(i,j)+=f(i,j,k)*p->DZN[KP];
            
           /* if(f(i,j,k)>p->F94 && f(i,j,k+1)<p->F93)
                vofheight(i,j)=MAX(vofheight(i,j),p->pos_z()+0.5*p->DZN[KP]);
            else if(f(i,j,k)<=p->F94 && f(i,j,k)>=p->F93)
            {
                if((a->nZ(i,j,k)>1E-06 || a->nZ(i,j,k)<-1E-06) && a->Alpha(i,j,k)<1E05)
                {
                    if(fabs(a->Alpha(i,j,k)/a->nZ(i,j,k))<0.5*p->DZN[KP])
                        vofheight(i,j)=MAX(vofheight(i,j),p->pos_z()+a->Alpha(i,j,k)/a->nZ(i,j,k));
                    else
                    {
                        if(a->nZ(i,j,k)>0.0)
                            vofheight(i,j)=MAX(vofheight(i,j),p->pos_z()-0.5*p->DZN[KP]);
                        else
                            vofheight(i,j)=MAX(vofheight(i,j),p->pos_z()+0.5*p->DZN[KP]);
                    }
                }
                else
                {
                    cout<<"surface cell in relax func does not have a plane"<<endl;
                    vofheight(i,j)=MAX(vofheight(i,j),(p->pos_z()-0.5*p->DZN[KP])+f(i,j,k)*p->DZN[KP]);
                }
            }*/
        }
    }
    count=0;
    
    LOOP
    {
    dg = distgen(p);
    db = distbeach(p);
    if(dg<1.0e20)
        genheight(i,j) = (1.0-relax4_wg(i,j))*ramp(p) * (eta(i,j)+p->phimean) + relax4_wg(i,j)*vofheight(i,j);
    else if (db<1.0e20)
        genheight(i,j) = (1.0-relax4_wg(i,j))*ramp(p) * (p->phimean) + relax4_wg(i,j)*vofheight(i,j);
    else
        genheight(i,j)=vofheight(i,j);
    }
    
    pgc->gcsl_start4(p,genheight,1);
    
    LOOP
    {
        dg = distgen(p);
        db = distbeach(p);

            
            
        
        // Wave Generation
        if(p->B98==2 && h_switch==1)
        {
            // Zone 1
            if(dg<1.0e20)
            {
            if(p->pos_z()+0.5*p->DZN[KP]<=genheight(i,j))
                f(i,j,k)=1.0;
            else if(p->pos_z()-0.5*p->DZN[KP]>=genheight(i,j))
                f(i,j,k)=0.0;
            else
            {
               /* double nx_eta,ny_eta,nz_eta,alpha_eta,nsum_eta;
                double nx_vof,ny_vof,nz_vof,alpha_vof,nsum_vof;
                
                nz_vof=1.0;
                nx_vof=-((genheight(i+1,j))-(genheight(i-1,j)))/(p->DXP[IP]+p->DXP[IM1]);
                ny_vof=0.0;
                nsum_vof=sqrt(nx_vof*nx_vof+ny_vof*ny_vof+nz_vof*nz_vof);
                nx_vof=nx_vof/nsum_vof;
                ny_vof=ny_vof/nsum_vof;
                nz_vof=nz_vof/nsum_vof;
                alpha_vof=(genheight(i,j)-p->pos_z())*nz_vof;*/
                
                /*nz_eta=1.0;
                nx_eta=((eta(i+1,j))-(eta(i-1,j)))/(p->DXP[IP]+p->DXP[IM1]);
                ny_eta=((eta(i,j+1))-(eta(i,j-1))/(p->DYP[JP]+p->DYP[JM1]);
                nsum_eta=sqrt(nx_eta*nx_eta+ny_eta*ny_eta+nz_eta*nz_eta);
                nx_eta=nx_eta/nsum_eta;
                ny_eta=ny_eta/nsum_eta;
                nz_eta=nz_eta/nsum_eta;
                if(fabs(eta(i,j)+p->phimean-p->pos_z())<=0.5*p->DZN[KP])
                    alpha_eta=(eta(i,j)+p->phimean-p->pos_z())*nz_eta;
                else if(((eta(i,j)+p->phimean)>p->pos_z()+0.5*p->DZN[KP]) && ((eta(i,j)+p->phimean)<=p->pos_z()+0.5*p->DZN[KP]+p->DZN[KP1]))
                    alpha_eta=(eta(i,j)+p->phimean-(p->pos_z()+p->DZP[KP]))*nz_eta;
                else if(((eta(i,j)+p->phimean)<p->pos_z()-0.5*p->DZN[KP]) && ((eta(i,j)+p->phimean)>=p->pos_z()-0.5*p->DZN[KP]-p->DZN[KM1]))
                    alpha_eta=(eta(i,j)+p->phimean-(p->pos_z()-p->DZP[KM1]))*nz_eta;
                else
                {
                    //cout<<"eta out of relax normal bounds"<<endl;
                    alpha_eta=1E06;
                }
                
                if(a->Alpha(i,j,k)<1E05)
                {
                    nx_vof=a->nX(i,j,k);
                    ny_vof=a->nY(i,j,k);
                    nz_vof=a->nZ(i,j,k);
                    alpha_vof=a->Alpha(i,j,k);
                }
                else if(a->Alpha(i,j,k+1)<1E05)
                {
                    nx_vof=a->nX(i,j,k+1);
                    ny_vof=a->nY(i,j,k+1);
                    nz_vof=a->nZ(i,j,k+1);
                    alpha_vof=a->Alpha(i,j,k+1);
                }
                else if(a->Alpha(i,j,k-1)<1E05)
                {
                    nx_vof=a->nX(i,j,k-1);
                    ny_vof=a->nY(i,j,k-1);
                    nz_vof=a->nZ(i,j,k-1);
                    alpha_vof=a->Alpha(i,j,k-1);
                }
                else
                {
                    cout<<"vofheight out of relax normal bounds"<<endl;
                    alpha_vof=1E06;
                }
    
                if(alpha_vof<1E05 && alpha_eta<1E05)
                {
                    nx_vof = (1.0-relax4_wg(i,j))*ramp(p) * nx_eta + relax4_wg(i,j)*nx_vof;
                    ny_vof = (1.0-relax4_wg(i,j))*ramp(p) * ny_eta + relax4_wg(i,j)*ny_vof;
                    nz_vof = (1.0-relax4_wg(i,j))*ramp(p) * nz_eta + relax4_wg(i,j)*nz_vof;
                    alpha_vof = (1.0-relax4_wg(i,j))*ramp(p) * alpha_eta + relax4_wg(i,j)*alpha_vof;
                    nsum_vof=sqrt(nx_vof*nx_vof+ny_vof*ny_vof+nz_vof*nz_vof);
                    nx_vof=nx_vof/nsum_vof;
                    ny_vof=ny_vof/nsum_vof;
                    nz_vof=nz_vof/nsum_vof;
                }
                else if(alpha_eta<1E05)
                {
                    nx_vof=nx_eta;
                    ny_vof=ny_eta;
                    nz_vof=nz_eta;
                    alpha_vof=alpha_eta;
                }
                */
                
              /*  if(alpha_vof<1E05)
                {
                    f(i,j,k)=V0Calc_PLIC(p,a,nx_vof,ny_vof,nz_vof,alpha_vof);
                }
                else*/
                //{
                    //cout<<"both eta and vof out of relax normal bounds"<<endl;
                    f(i,j,k)=(genheight(i,j)-(p->pos_z()-0.5*p->DZN[KP]))/p->DZN[KP];
               // }
                
            }
            //    f(i,j,k)=(localheight-(p->pos_z()-0.5*p->DZN[KP]))/p->DZN[KP];
           // f(i,j,k) = (1.0-relax4_wg(i,j))*ramp(p) * vofgen(i,j,k) + relax4_wg(i,j)*f(i,j,k);
            ++count;
            }
        }
            
        // Numerical Beach    
        if(p->B99==2)
        {
            // Zone 2
            if(db<1.0e20)
            {
            if(p->pos_z()+0.5*p->DZN[KP]<=genheight(i,j))
                f(i,j,k)=1.0;
            else if(p->pos_z()-0.5*p->DZN[KP]>=genheight(i,j))
                f(i,j,k)=0.0;
            else
                f(i,j,k)=(genheight(i,j)-(p->pos_z()-0.5*p->DZN[KP]))/p->DZN[KP];
            }
        }
    }
    
    
    p->wavecalctime+=pgc->timer()-starttime;

}

void iowave::turb_relax(lexer *p, fdm *a, ghostcell *pgc, field &f)
{
    starttime=pgc->timer();
    
    LOOP
    {
        dg = distgen(p);    
        db = distbeach(p);

        phival = -0.5*(a->phi(i,j,k)+a->phi(i-1,j,k));

        if(phival>=-psi)
		 H=1.0;
         
		if(phival<-epsi)
		H=0.0;

		if(phival>=-epsi && phival<-psi)
		H=0.5*(1.0 + phival/epsi + (1.0/PI)*sin((PI*phival)/epsi));        
        
		// Wave Generation
		if(p->B98==2 && u_switch==1)
        {
            // Zone 1
            if(dg<1.0e20)
            f(i,j,k) = relax4_wg(i,j)*f(i,j,k) + (1.0-H)*f(i,j,k);

		}
    }
    
    p->wavecalctime+=pgc->timer()-starttime;
}

void iowave::fi_relax(lexer *p, ghostcell *pgc, field& f, field& phi)
{
}

void iowave::fivec_relax(lexer *p, ghostcell *pgc, double *f)
{
}


void iowave::visc_relax(lexer *p, ghostcell *pgc, slice& f)
{
    /*count=0;
    SLICELOOP4
    {
        dg = distgen(p);
		db = distbeach(p);
        
		// Numerical Beach
        if(p->B99==1||p->B99==2||beach_relax==1)
		{
            // Zone 2
            if(db<1.0e20)
            f(i,j) = relax4_nb(i,j)*f(i,j) + (1.0-relax4_nb(i,j))*1.86;
        }
    }*/
}


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
        
        //LS version
        //if(p->F80!=4)
        {
            phival = 0.5*(a->phi(i,j,k)+a->phi(i+1,j,k));

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
        }
        
        //PLIC version
        /*else if(p->F80==4)
        {
            if(p->F92==3||p->F92==32)
            {
                if(p->j_dir>0)
                    H=(0.5*p->DXN[IP]*(0.25*a->vof_nte(i,j,k)+0.25*a->vof_ntw(i,j,k)+0.25*a->vof_nbe(i,j,k)+0.25*a->vof_nbw(i,j,k))
                        +0.5*p->DXN[IP1]*(0.25*a->vof_ste(i+1,j,k)+0.25*a->vof_stw(i+1,j,k)+0.25*a->vof_sbe(i+1,j,k)+0.25*a->vof_sbw(i+1,j,k))
                        )/p->DXP[IP];
                else
                    H=(0.25*p->DXN[IP]*a->vof_nt(i,j,k)+0.25*p->DXN[IP]*a->vof_nb(i,j,k)+0.25*p->DXN[IP1]*a->vof_st(i+1,j,k)+0.25*p->DXN[IP1]*a->vof_sb(i+1,j,k))/p->DXP[IP];

            }
            else
                H=(0.5*a->vof(i+1,j,k)*p->DXN[IP1]+0.5*a->vof(i,j,k)*p->DXN[IP])/p->DXP[IP];
        
            G=H;
        
            if(H>=0.5)
            {
                if(p->pos_z()<=p->phimean)
                    z=-(fabs(p->phimean-p->pos_z()));
		
                if(p->pos_z()>p->phimean)
                    z=(fabs(p->phimean-p->pos_z()));
            }
        
            if(H<0.0)
                z = 0.5*(eta(i,j)+eta(i+1,j));
        }*/
		
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
        
       // if(p->F80!=4)
        {
            phival = 0.5*(a->phi(i,j,k)+a->phi(i,j+1,k));

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
        }
        /*else if(p->F80==4)
        {
            if(p->F92==3||p->F92==32)
            {
                if(p->j_dir>0)
                    H=(0.5*p->DXN[IP]*(0.25*a->vof_nte(i,j,k)+0.25*a->vof_ntw(i,j,k)+0.25*a->vof_nbe(i,j,k)+0.25*a->vof_nbw(i,j,k))
                        +0.5*p->DXN[IP1]*(0.25*a->vof_ste(i+1,j,k)+0.25*a->vof_stw(i+1,j,k)+0.25*a->vof_sbe(i+1,j,k)+0.25*a->vof_sbw(i+1,j,k))
                        )/p->DXP[IP];
                else
                    H=(0.25*p->DXN[IP]*a->vof_nt(i,j,k)+0.25*p->DXN[IP]*a->vof_nb(i,j,k)+0.25*p->DXN[IP1]*a->vof_st(i+1,j,k)+0.25*p->DXN[IP1]*a->vof_sb(i+1,j,k))/p->DXP[IP];

            }
            else
                H=(0.5*a->vof(i+1,j,k)*p->DXN[IP1]+0.5*a->vof(i,j,k)*p->DXN[IP])/p->DXP[IP];
        
            G=H;
        
            if(H>=0.5)
            {
                if(p->pos_z()<=p->phimean)
                    z=-(fabs(p->phimean-p->pos_z()));
		
                if(p->pos_z()>p->phimean)
                    z=(fabs(p->phimean-p->pos_z()));
            }
        
            if(H<0.0)
                z = 0.5*(eta(i,j)+eta(i+1,j));
        }*/
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
        
       // if(p->F80!=4)
        {
            phival = 0.5*(a->phi(i,j,k)+a->phi(i,j,k+1));
        
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
        }
     /*   else if(p->F80==4)
        {
            if(p->F92==3||p->F92==32)
            {
                if(p->j_dir>0)
                    H=(0.5*p->DZN[KP]*(0.25*a->vof_nte(i,j,k)+0.25*a->vof_ntw(i,j,k)+0.25*a->vof_ste(i,j,k)+0.25*a->vof_stw(i,j,k))
                        +0.5*p->DZN[KP1]*(0.25*a->vof_nbe(i,j,k+1)+0.25*a->vof_nbw(i,j,k+1)+0.25*a->vof_sbe(i,j,k+1)+0.25*a->vof_sbw(i,j,k+1))
                        )/p->DZP[KP];
                else
                    H=(0.25*p->DZN[KP]*a->vof_nt(i,j,k)+0.25*p->DZN[KP]*a->vof_st(i,j,k)+0.25*p->DZN[KP1]*a->vof_nb(i,j,k+1)+0.25*p->DZN[KP1]*a->vof_sb(i,j,k+1))/p->DZP[KP];

            }
            else
                H=(0.5*a->vof(i,j,k+1)*p->DZN[KP1]+0.5*a->vof(i,j,k)*p->DZN[KP])/p->DZP[KP];
            
            if(H>=0.5)
            {
                if(p->pos_z()<=p->phimean)
                    z=-(fabs(p->phimean-p->pos3_z()));
		
                if(p->pos_z()>p->phimean)
                    z=(fabs(p->phimean-p->pos3_z()));
            }
        
            if(H<0.5)
                z = eta(i,j);
        }*/

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
                    f(i,j,k)=(genheight(i,j)-(p->pos_z()-0.5*p->DZN[KP]))/p->DZN[KP];
            }
            
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
        
      //  if(p->F80!=4)
        {
            phival = 0.5*(a->phi(i,j,k)+a->phi(i-1,j,k));

            if(phival>=-psi)
                H=1.0;
         
            if(phival<-epsi)
                H=0.0;

            if(phival>=-epsi && phival<-psi)
                H=0.5*(1.0 + phival/epsi + (1.0/PI)*sin((PI*phival)/epsi));        
        }
     /*   else if(p->F80==4)
        {
            H=a->vof(i,j,k);
        }*/
		// Wave Generation
		if(p->B98==2 && u_switch==1)
        {
            // Zone 1
            if(dg<1.0e20)
            f(i,j,k) = relax4_wg(i,j)*f(i,j,k);// + (1.0-H)*f(i,j,k);
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


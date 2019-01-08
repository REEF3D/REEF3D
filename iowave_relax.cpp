/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2019 Hans Bihs

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
--------------------------------------------------------------------*/

#include"iowave.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

void iowave::u_relax(lexer *p, fdm *a, ghostcell *pgc, field& uvel)
{
    count=0;
    ULOOP
    {
		xg = xgen1(p);
        yg = ygen1(p);
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
		
		if(p->B121==0)
		G=1.0;
        
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
        if(p->B98==1 && u_switch==1)
        {
            // Zone 1
            if(dg<dist1)
            {
            uvel(i,j,k) = ra1(p,dg) * uval[count] * H + (1.0-ra1(p,dg)*G)*uvel(i,j,k);
            ++count;
            }

            // Zone 2
            if(dg>=dist1 && dg<dist2)
            {
            uvel(i,j,k) = ra2(p,dg) * uval[count] * H + (1.0-ra2(p,dg)*G)*uvel(i,j,k);
            ++count;
            }
		}
		
		if(p->B98==2 && u_switch==1)
        {
            // Zone 1
            if(dg<dist1)
            {
            uvel(i,j,k) = (1.0-rb1(p,dg))*ramp(p)*uval[count] * H + rb1(p,dg)*H*uvel(i,j,k) + (1.0-G)*uvel(i,j,k);
            ++count;
            }
		}
		
		// Numerical Beach
		if(p->B99==1)
		{
            // Zone 3
            if(db<dist3)
            uvel(i,j,k) = (1.0-ra3(p,db))*0.0 + ra3(p,db)*uvel(i,j,k);
        }

        if(p->B99==2||p->B99==4)
		{
            // Zone 3
            if(db<dist3)
            uvel(i,j,k) = (1.0-rb3(p,db))*0.0 + rb3(p,db)*uvel(i,j,k);
        }
    }
}

void iowave::v_relax(lexer *p, fdm *a, ghostcell *pgc, field& vvel)
{
    count=0;
    VLOOP
    {
        xg = xgen2(p);
        yg = ygen2(p);
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
		
		if(p->B121==0)
		G=1.0;
        
        
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
        if(p->B98==1 && v_switch==1)
        {
            // Zone 1
            if(dg<dist1)
            {
            vvel(i,j,k) = ra1(p,dg) * vval[count] * H + (1.0-ra1(p,dg)*G)*vvel(i,j,k);
            ++count;
            }

            // Zone 2
            if(dg>=dist1 && dg<dist2)
            {
            vvel(i,j,k) = ra2(p,dg) * vval[count] * H + (1.0-ra2(p,dg)*G)*vvel(i,j,k);
            ++count;
            }
		}
		
		if(p->B98==2 && v_switch==1)
        {
            // Zone 1
            if(dg<dist1)
            {
            vvel(i,j,k) = (1.0-rb1(p,dg))*ramp(p)*vval[count] * H + rb1(p,dg)*H*vvel(i,j,k) + (1.0-G)*vvel(i,j,k);
            ++count;
            }
		}
		
		// Numerical Beach
		if(p->B99==1)
		{
            // Zone 3
            if(db<dist3)
            vvel(i,j,k) = (1.0-ra3(p,db))*0.0 + ra3(p,db)*vvel(i,j,k);
        }

		if(p->B99==2||p->B99==4)
		{	
            // Zone 3
            if(db<dist3)
            vvel(i,j,k) = (1.0-rb3(p,db))*0.0 + rb3(p,db)*vvel(i,j,k);
        }
    }
}

void iowave::w_relax(lexer *p, fdm *a, ghostcell *pgc, field& wvel)
{
    count=0;
    WLOOP
    {
        xg = xgen(p);
        yg = ygen(p);
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
		
		if(p->B121==0)
		G=1.0;
        
        
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
        if(p->B98==1 && w_switch==1)
        {
            // Zone 1
            if(dg<dist1)
            {
            wvel(i,j,k) = ra1(p,dg) * wval[count] * H + (1.0-ra1(p,dg)*G)*wvel(i,j,k);
            ++count;
            }

            // Zone 2
            if(dg>=dist1 && dg<dist2)
            {
            wvel(i,j,k) = ra2(p,dg) * wval[count] * H + (1.0-ra2(p,dg)*G)*wvel(i,j,k);
            ++count;
            }
		}
		
		if(p->B98==2 && w_switch==1)
        {
            // Zone 1
            if(dg<dist1)
            {
            wvel(i,j,k) = (1.0-rb1(p,dg)) * ramp(p)* wval[count] * H + rb1(p,dg)*H*wvel(i,j,k) + (1.0-G)*wvel(i,j,k);
            ++count;
            }

		}
		
		// Numerical Beach
		if(p->B99==1)
		{
            // Zone 3
            if(db<dist3)
            wvel(i,j,k) = (1.0-ra3(p,db))*0.0 + ra3(p,db)*wvel(i,j,k);
        }

        if(p->B99==2||p->B99==4)
		{
            // Zone 3
            if(db<dist3)
            wvel(i,j,k) = (1.0-rb3(p,db))*0.0 + rb3(p,db)*wvel(i,j,k);
        }
    }		
}

void iowave::p_relax(lexer *p, fdm *a, ghostcell *pgc, field& press)
{
	
    LOOP
    {
		dg = distgen(p);
		db = distbeach(p); 
		
		if(p->B99==1)
        {
			// Zone 3   
            if(db<dist3)
			{
            if(p->D38==0)
            press(i,j,k) = (1.0-ra3(p,db))*((p->phimean - p->pos_z())*a->ro(i,j,k)*fabs(p->W22)) + ra3(p,db)*press(i,j,k);
            
            if(p->D38>0)
            press(i,j,k) = (1.0-ra3(p,db))*0.0 + ra3(p,db)*press(i,j,k);
			}
		}
		

        if(p->B99==2||p->B99==4)
        {            
            // Zone 3
            if(db<dist3)
			{
            if(p->D38==0)
            press(i,j,k) = (1.0-rb3(p,db))*((p->phimean - p->pos_z())*a->ro(i,j,k)*fabs(p->W22)) + rb3(p,db)*press(i,j,k);
            
            if(p->D38>0)
            press(i,j,k) = (1.0-rb3(p,db))*0.0 + rb3(p,db)*press(i,j,k);
			}
        }
    }	 
}

void iowave::phi_relax(lexer *p, ghostcell *pgc, field& f)
{
    count=0;
    FLUIDLOOP
    {
        xg = xgen(p);
        yg = ygen(p);
        dg = distgen(p);
        db = distbeach(p);
        
        if(p->pos_z()<=p->phimean)
        z=-(fabs(p->phimean-p->pos_z()));
            
        if(p->pos_z()>p->phimean)
        z=(fabs(p->phimean-p->pos_z()));
            
        // Wave Generation
        if(p->B98==1 && h_switch==1)
        {
            // Zone 1
            if(dg<dist1)
            {
            f(i,j,k) = ra1(p,dg) * lsval[count] + (1.0-ra1(p,dg))*f(i,j,k);
            ++count;
            }
            
            // Zone 2
            if(dg>=dist1 && dg<dist2)
            {
            f(i,j,k) = ra2(p,dg) * lsval[count] + (1.0-ra2(p,dg))*f(i,j,k);
            ++count;
            }
        }

        if(p->B98==2 && h_switch==1)
        {
            // Zone 1
            if(dg<dist1)
            {
            f(i,j,k) = (1.0-rb1(p,dg))*ramp(p) * lsval[count] + rb1(p,dg) * f(i,j,k);
            ++count;
            }
        }
            
        // Numerical Beach
        if(p->B99==1)
        {
            // Zone 3
            if(db<dist3)
            f(i,j,k) = (1.0-ra3(p,db)) * (p->phimean-p->pos_z()) + ra3(p,db)*f(i,j,k);
        }
            
        if(p->B99==2)
        {
            // Zone 3
            if(db<dist3)
            f(i,j,k) = (1.0-rb3(p,db)) * (p->phimean-p->pos_z()) + rb3(p,db)*f(i,j,k);
        }
    }
}


void iowave::vof_relax(lexer *p, ghostcell *pgc, field& f)
{
    count=0;
    FLUIDLOOP
    {
        xg = xgen(p);
        yg = ygen(p);
		dg = distgen(p);
		db = distbeach(p);
        
		if(p->pos_z()<=p->phimean)
        z=-(fabs(p->phimean-p->pos_z()));
		
		if(p->pos_z()>p->phimean)
        z=(fabs(p->phimean-p->pos_z()));	
  
        double fl = (eta(i-1,j)*p->DXN[IP] + eta(i,j)*p->DXN[IM1])/(p->DXN[IP] + p->DXN[IM1]) + p->phimean;
        double fr = (eta(i,j)*p->DXN[IP] + eta(i+1,j)*p->DXN[IP1])/(p->DXN[IP] + p->DXN[IP1]) + p->phimean;
        double fc = (fl + fr)/2.0;
                
		// Wave Generation
        if(p->B98==1 && h_switch==1)
        {
            // Zone 1
            if(dg<dist1)
            {
                if (MIN(fl,fr) >= p->pos_z() + p->DZN[KP]/2.0)
                {
                    lsval[count] = 1.0;
                }
                else if (MAX(fl,fr) <= p->pos_z() - p->DZN[KP]/2.0)
                {
                    lsval[count] = 0.0;
                }
                else
                {
                    lsval[count] = (fc - p->pos_z() + p->DZN[KP]/2.0)/p->DZN[KP];
                }

                f(i,j,k) = ra1(p,dg) * lsval[count] + (1.0-ra1(p,dg))*f(i,j,k);
                ++count;
            }

            // Zone 2
            if(dg>=dist1 && dg<dist2)
            {
                if (MIN(fl,fr) >= p->pos_z() + p->DZN[KP]/2.0)
                {
                    lsval[count] = 1.0;
                }
                else if (MAX(fl,fr) <= p->pos_z() - p->DZN[KP]/2.0)
                {
                    lsval[count] = 0.0;
                }
                else
                {
                    lsval[count] = (fc - p->pos_z() + p->DZN[KP]/2.0)/p->DZN[KP];
                }
               
                f(i,j,k) = ra2(p,dg) * lsval[count] + (1.0-ra2(p,dg))*f(i,j,k);
                ++count;
            }
		}

        if(p->B98==2 && h_switch==1)
        {
            // Zone 1
            if(dg<dist1)
            {
                if (MIN(fl,fr) >= p->pos_z() + p->DZN[KP]/2.0)
                {
                    lsval[count] = 1.0;
                }
                else if (MAX(fl,fr) <= p->pos_z() - p->DZN[KP]/2.0)
                {
                    lsval[count] = 0.0;
                }
                else
                {
                    lsval[count] = (fc - p->pos_z() + p->DZN[KP]/2.0)/p->DZN[KP];
                }
                
                f(i,j,k) = (1.0-rb1(p,dg))*ramp(p) * lsval[count] + rb1(p,dg) * f(i,j,k);
                ++count;
            }
		}
        
        
        fc = p->phimean;
		double value;
        
		// Numerical Beach
		if(p->B99==1)
		{
		      // Zone 3
            if(db<dist3)
            {
                if (fc >= p->pos_z() + p->DZN[KP]/2.0)
                {
                    value = 1.0;
                }
                else if (fc <= p->pos_z() - p->DZN[KP]/2.0)
                {
                    value = 0.0;
                }
                else 
                {
                    value = (fc - p->pos_z() + p->DZN[KP]/2.0)/p->DZN[KP];
                }
            
                f(i,j,k) = (1.0-ra3(p,db)) * value + ra3(p,db)*f(i,j,k);
            }
        }
		
		if(p->B99==2)
		{
            // Zone 3
            if(db<dist3)
            {
                if (fc >= p->pos_z() + p->DZN[KP]/2.0)
                {
                    value = 1.0;
                }
                else if (fc <= p->pos_z() - p->DZN[KP]/2.0)
                {
                    value = 0.0;
                }
                else 
                {
                    value = (fc - p->pos_z() + p->DZN[KP]/2.0)/p->DZN[KP];
                }

                f(i,j,k) = (1.0-rb3(p,db)) * value + rb3(p,db)*f(i,j,k);
            }
        }
    }
}



void iowave::fi_relax(lexer *p, ghostcell *pgc, field& f, field& phi)
{
    count=0;
    FLUIDLOOP
    {
		xg = xgen(p);
        yg = ygen(p);
		dg = distgen(p);
		db = distbeach(p); 


        phival = phi(i,j,k);
        

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
		H=0.5*(1.0 + phival/fabs(epsi) + (1.0/PI)*sin((PI*phival)/fabs(epsi)));
		G=H;
		}
		
		if(p->B121==0)
		G=1.0;
        
        if(p->A300==1)
        H=G=1.0;
		
		// Wave Generation
        if(p->B98==1 && f_switch==1)
        {
            // Zone 1
            if(dg<dist1)
            {
            f(i,j,k) = ra1(p,dg) * Fival[count] * H + (1.0-ra1(p,dg)*G)*f(i,j,k);
            ++count;
            }

            // Zone 2
            if(dg>=dist1 && dg<dist2)
            {
            f(i,j,k) = ra2(p,dg) * Fival[count] * H + (1.0-ra2(p,dg)*G)*f(i,j,k);
            ++count;
            }
			
		}
		
		if(p->B98==2 && f_switch==1)
        {
            // Zone 1
            if(dg<dist1)
            {
            f(i,j,k) = (1.0-rb1(p,dg))*ramp(p)*Fival[count] * H + rb1(p,dg)*H*f(i,j,k) + (1.0-G)*f(i,j,k);
            ++count;
            }
		}
		
		// Numerical Beach
		if(p->B99==1)
		{
            // Zone 3
            if(db<dist3)
            f(i,j,k) = (1.0-ra3(p,db))*0.0 + ra3(p,db)*f(i,j,k);
        }

        if(p->B99==2||p->B99==4)
		{
            // Zone 3
            if(db<dist3)
            f(i,j,k) = (1.0-rb3(p,db))*0.0 + rb3(p,db)*f(i,j,k);
        }
    }
}

void iowave::fivec_relax(lexer *p, ghostcell *pgc, double *f)
{
    int dbcount=0;
    count=0;
    FLOOP
    {
		dg = distgen(p);
		db = distbeach(p); 
		
		// Wave Generation
        if(p->B98==1 && f_switch==1)
        {
            // Zone 1
            if(dg<dist1)
            {
            f[FIJK] = ra1(p,dg) * Fival[count] + (1.0-ra1(p,dg))*f[FIJK];
            ++count;
            }

            // Zone 2
            if(dg>=dist1 && dg<dist2)
            {
            f[FIJK] = ra2(p,dg) * Fival[count] + (1.0-ra2(p,dg))*f[FIJK];
            ++count;
            }
			
		}
		
		if(p->B98==2 && f_switch==1)
        {
            // Zone 1
            if(dg<dist1)
            {
            f[FIJK] = (1.0-rb1val[count])*ramp(p)*Fival[count] + rb1val[count]*f[FIJK];
            ++count;
            }
		}
		
		// Numerical Beach
		if(p->B99==1)
		{
            // Zone 3
            if(db<dist3)
            f[FIJK] = (1.0-ra3(p,db))*0.0 + ra3(p,db)*f[FIJK];
        }

        if(p->B99==2||p->B99==4)
		{
            // Zone 3
            if(db<dist3)
            {
            f[FIJK] = rb3val[dbcount]*f[FIJK];
            ++dbcount;
            }
        }
    }
}

void iowave::fifsf_relax(lexer *p, ghostcell *pgc, slice& f)
{
    count=0;
    SLICELOOP4
    {
		xg = xgen(p);
        yg = ygen(p);
		dg = distgen(p);
		db = distbeach(p); 
        
        H=G=1.0;
        
        z = eta(i,j);
		
		// Wave Generation
        if(p->B98==1 && f_switch==1)
        {
            // Zone 1
            if(dg<dist1)
            {
            f(i,j) = ra1(p,dg) * Fifsfval[count]  + (1.0-ra1(p,dg))*f(i,j);
            ++count;
            }

            // Zone 2
            if(dg>=dist1 && dg<dist2)
            {
            f(i,j) = ra2(p,dg) * Fifsfval[count]  + (1.0-ra2(p,dg))*f(i,j);
            ++count;
            }
			
		}
		
		if(p->B98==2 && f_switch==1)
        {
            // Zone 1
            if(dg<dist1)
            {
            f(i,j) = (1.0-rb1(p,dg))*ramp(p)*Fifsfval[count]  + rb1(p,dg)*f(i,j);
            ++count;
            }
		}
		
		// Numerical Beach
		if(p->B99==1)
		{
            // Zone 3
            if(db<dist3)
            f(i,j) = ra3(p,db)*f(i,j);
        }

        if(p->B99==2||p->B99==4)
		{
            // Zone 3
            if(db<dist3)
            f(i,j) = rb3(p,db)*f(i,j);
        }
    }
}

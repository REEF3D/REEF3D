/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2023 Hans Bihs

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

#include"lexer.h"
#include"ghostcell.h"

void lexer::lexer_gridspacing(ghostcell *pgc)
{
    Darray(XP,knox+1+4*marge);
    Darray(YP,knoy+1+4*marge);
    Darray(ZP,knoz+1+4*marge);
    
    Darray(RP,knox+1+4*marge);
    Darray(SP,knoy+1+4*marge);
    Darray(TP,knoz+1+4*marge);
    
    Darray(DXN,knox+1+4*marge);
    Darray(DYN,knoy+1+4*marge);
    Darray(DZN,knoz+1+4*marge);
    
    Darray(DXP,knox+1+4*marge);
    Darray(DYP,knoy+1+4*marge);
    Darray(DZP,knoz+1+4*marge);
    
    Darray(DRDXN,knox+1+4*marge);
    Darray(DSDYN,knoy+1+4*marge);
    Darray(DTDZN,knoz+1+4*marge);
    
    Darray(DRDXP,knox+1+4*marge);
    Darray(DSDYP,knoy+1+4*marge);
    Darray(DTDZP,knoz+1+4*marge);
    
    Darray(DDRDDXN,knox+1+4*marge);
    Darray(DDSDDYN,knoy+1+4*marge);
    Darray(DDTDDZN,knoz+1+4*marge);
    
    Darray(DDRDDXP,knox+1+4*marge);
    Darray(DDSDDYP,knoy+1+4*marge);
    Darray(DDTDDZP,knoz+1+4*marge);
    
    Darray(ZSN,imax*jmax*(kmax+1));
    Darray(ZSP,imax*jmax*kmax);

    // XP,YP,ZP
    for(i=-marge;i<knox+marge;++i)
    XP[IP] = 0.5*(XN[IP]+XN[IP1]);
    
    for(j=-marge;j<knoy+marge;++j)
    YP[JP] = 0.5*(YN[JP]+YN[JP1]);
    
    for(k=-marge;k<knoz+marge;++k)
    ZP[KP] = 0.5*(ZN[KP]+ZN[KP1]);
    
    // RP,SP,TP
    for(i=-marge;i<knox+marge;++i)
    RP[IP] = 0.5*(RN[IP]+RN[IP1]);
    
    for(j=-marge;j<knoy+marge;++j)
    SP[JP] = 0.5*(SN[JP]+SN[JP1]);
    
    for(k=-marge;k<knoz+marge;++k)
    TP[KP] = 0.5*(TN[KP]+TN[KP1]);
    
    //dx
    for(i=-marge;i<knox+marge;++i)
    DXN[IP] = XN[IP1]-XN[IP];
    
    for(j=-marge;j<knoy+marge;++j)
    DYN[JP] = YN[JP1]-YN[JP];
    
    for(k=-marge;k<knoz+marge;++k)
    DZN[KP] = ZN[KP1]-ZN[KP];
    
    
    
    // dxn
    
    for(i=-marge;i<knox+marge;++i)
    DXP[IP] = 0.5*(XN[IP2]+XN[IP1]) - 0.5*(XN[IP1]+XN[IP]);

    for(j=-marge;j<knoy+marge;++j)
    DYP[JP] = 0.5*(YN[JP2]+YN[JP1]) - 0.5*(YN[JP1]+YN[JP]);
    
    for(k=-marge;k<knoz+marge;++k)
    DZP[KP] = 0.5*(ZN[KP2]+ZN[KP1]) - 0.5*(ZN[KP1]+ZN[KP]);
    
    DXM = DXD = DYD = 0.0;
    
    int count=0;
    int xcount=0;
    int ycount=0;
    
    for(i=0;i<knox;++i)
    {
    DXM += DXP[IP];
    DXD += DXP[IP];
    ++count;
    ++xcount;
    }
    
    if(j_dir==1)
    for(j=0;j<knoy;++j)
    {
    DXM += DYP[JP];
    DYD += DYP[JP];
    ++count;
    ++ycount;
    }
    
    for(k=0;k<knoz;++k)
    {
    DXM += DZP[KP];
    ++count;
    }
    
    count = pgc->globalisum(count);
    xcount = pgc->globalisum(xcount);
    ycount = pgc->globalisum(ycount);
    
    DXM = pgc->globalsum(DXM);
    DXD = pgc->globalsum(DXD);
    DYD = pgc->globalsum(DYD);
    
    DXM /= double(count); 
    DXD /= double(xcount); 
    DYD /= double(ycount); 
    
    
    DXM = pgc->globalmin(DXM);
    DXD = pgc->globalmin(DXD);
    DYD = pgc->globalmin(DYD);
    
    
    // transformation
    // 1st derivative
    for(i=-1;i<knox+2;++i)
    DRDXN[IP] =  (-RN[IP2] + 8.0*RN[IP1] - 8.0*RN[IM1] + RN[IM2])
                /(-XN[IP2] + 8.0*XN[IP1] - 8.0*XN[IM1] + XN[IM2]);  
                
    for(j=-1;j<knoy+2;++j)
    DSDYN[JP] =  (-SN[JP2] + 8.0*SN[JP1] - 8.0*SN[JM1] + SN[JM2])
                /(-YN[JP2] + 8.0*YN[JP1] - 8.0*YN[JM1] + YN[JM2]);  

    
                
    for(k=-1;k<knoz+2;++k)
    DTDZN[KP] =  (-TN[KP2] + 8.0*TN[KP1] - 8.0*TN[KM1] + TN[KM2])
                /(-ZN[KP2] + 8.0*ZN[KP1] - 8.0*ZN[KM1] + ZN[KM2]);  
                
                
                
    for(i=-1;i<knox+1;++i)
    DRDXP[IP] =  (-RP[IP2] + 8.0*RP[IP1] - 8.0*RP[IM1] + RP[IM2])
                /(-XP[IP2] + 8.0*XP[IP1] - 8.0*XP[IM1] + XP[IM2]);  
                
    for(j=-1;j<knoy+1;++j)
    DSDYP[JP] =  (-SP[JP2] + 8.0*SP[JP1] - 8.0*SP[JM1] + SP[JM2])
                /(-YP[JP2] + 8.0*YP[JP1] - 8.0*YP[JM1] + YP[JM2]);  
                
    for(k=-1;k<knoz+1;++k)
    DTDZP[KP] =  (-TP[KP2] + 8.0*TP[KP1] - 8.0*TP[KM1] + TP[KM2])
                /(-ZP[KP2] + 8.0*ZP[KP1] - 8.0*ZP[KM1] + ZP[KM2]);  
                
    // 2nd derivative
    double den;
    for(i=0;i<knox+2;++i)
    {
    den = (-XN[IP2] + 16.0*XN[IP1] - 30.0*XN[IP] + 16.0*XN[IM1] - XN[IM2]);
    
    DDRDDXN[IP] =  (-RN[IP2] + 16.0*RN[IP1] - 30.0*RN[IP] + 16.0*RN[IM1] - RN[IM2])
                  /(fabs(den)>1.0e-10?den:1.0e20); 
    }
                  
    for(j=-1;j<knoy+2;++j)
    {
    den = (-YN[JP2] + 16.0*YN[JP1] - 30.0*YN[JP] + 16.0*YN[JM1] - YN[JM2]);
    
    DDSDDYN[JP] =  (-SN[JP2] + 16.0*SN[JP1] - 30.0*SN[JP] + 16.0*SN[JM1] - SN[JM2])
                  /(fabs(den)>1.0e-10?den:1.0e20); 
    }
    
    for(k=-1;k<knoz+2;++k)
    {
    den = (-ZN[KP2] + 16.0*ZN[KP1] - 30.0*ZN[KP] + 16.0*ZN[KM1] - ZN[KM2]);

    DDTDDZN[KP] =  (-TN[KP2] + 16.0*TN[KP1] - 30.0*TN[KP] + 16.0*TN[KM1] - TN[KM2])
                  /(fabs(den)>1.0e-10?den:1.0e20); 
    }
            
      
                  
    for(i=-1;i<knox+1;++i)
    {
    den = (-XP[IP2] + 16.0*XP[IP1] - 30.0*XP[IP] + 16.0*XP[IM1] - XP[IM2]); 
    
    DDRDDXP[IP] =  (-RP[IP2] + 16.0*RP[IP1] - 30.0*RP[IP] + 16.0*RP[IM1] - RP[IM2])
                  /(fabs(den)>1.0e-10?den:1.0e20); 
    }
                  
    for(j=-1;j<knoy+1;++j)
    {
    den = (-YP[JP2] + 16.0*YP[JP1] - 30.0*YP[JP] + 16.0*YP[JM1] - YP[JM2]);
    
    DDSDDYP[JP] =  (-SP[JP2] + 16.0*SP[JP1] - 30.0*SP[JP] + 16.0*SP[JM1] - SP[JM2])
                  /(fabs(den)>1.0e-10?den:1.0e20); 
    }
   
    for(k=-1;k<knoz+1;++k)
    {
    den = (-ZP[KP2] + 16.0*ZP[KP1] - 30.0*ZP[KP] + 16.0*ZP[KM1] - ZP[KM2]);
    
    DDTDDZP[KP] =  (-TP[KP2] + 16.0*TP[KP1] - 30.0*TP[KP] + 16.0*TP[KM1] - TP[KM2])
                  /(fabs(den)>1.0e-10?den:1.0e20); 
                  
    }
    
    

    
}

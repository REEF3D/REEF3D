/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

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

#define NLOOP1 for(n=p->sizeM1[0]; n<p->sizeM1[1]; ++n)
#define NLOOP2 for(n=p->sizeM2[0]; n<p->sizeM2[1]; ++n)
#define NLOOP3 for(n=p->sizeM3[0]; n<p->sizeM3[1]; ++n)
#define NLOOP4 for(n=p->sizeM4[0]; n<p->sizeM4[1]; ++n)
#define NLOOP4A for(n=p->sizeM4a[0]; n<p->sizeM4a[1]; ++n)	
#define NLOOP6 for(n=p->sizeM6[0]; n<p->sizeM6[1]; ++n)
#define NLOOP9 for(n=p->sizeM9[0]; n<p->sizeM9[1]; ++n)
#define NLOOP for(n=sizeM[0]; n<sizeM[1]; ++n)
#define VECLOOP for(n=0; n<p->veclength; ++n)


#define PORVAL1 (0.5*(a->porosity(i+1,j,k) + a->porosity(i,j,k)))
#define PORVAL2 (0.5*(a->porosity(i,j+1,k) + a->porosity(i,j,k)))
#define PORVAL3 (0.5*(a->porosity(i,j,k+1) + a->porosity(i,j,k)))
#define PORVAL4 a->porosity(i,j,k)

#define PORVAL4px a->porosity(i+1,j,k)
#define PORVAL4py a->porosity(i,j+1,k)
#define PORVAL4pz a->porosity(i,j,k+1)

#define CPOR4px   (1.0/(1.0+(p->B260*(PORVAL4px<1.0?1.0:0.0))))
#define CPOR4py   (1.0/(1.0+(p->B260*(PORVAL4py<1.0?1.0:0.0))))
#define CPOR4pz   (1.0/(1.0+(p->B260*(PORVAL4pz<1.0?1.0:0.0))))

#define CPOR1   (1.0/(1.0+(p->B260*(PORVAL1<1.0?1.0:0.0))))
#define CPOR2   (1.0/(1.0+(p->B260*(PORVAL2<1.0?1.0:0.0))))
#define CPOR3   (1.0/(1.0+(p->B260*(PORVAL3<1.0?1.0:0.0))))
#define CPOR4   (1.0/(1.0+(p->B260*(PORVAL4<1.0?1.0:0.0))))

#define PORVAL1m (0.5*(a->porosity(i,j,k) + a->porosity(i-1,j,k)))
#define PORVAL2m (0.5*(a->porosity(i,j,k) + a->porosity(i,j-1,k)))
#define PORVAL3m (0.5*(a->porosity(i,j,k) + a->porosity(i,j,k-1)))

#define CPOR1m   (1.0/(1.0+(p->B260*(PORVAL1m<1.0?1.0:0.0))))
#define CPOR2m   (1.0/(1.0+(p->B260*(PORVAL2m<1.0?1.0:0.0))))
#define CPOR3m   (1.0/(1.0+(p->B260*(PORVAL3m<1.0?1.0:0.0))))

#define PORVAL1p (0.5*(a->porosity(i+2,j,k) + a->porosity(i+1,j,k)))
#define PORVAL2p (0.5*(a->porosity(i,j+2,k) + a->porosity(i,j+1,k)))
#define PORVAL3p (0.5*(a->porosity(i,j,k+2) + a->porosity(i,j,k+1)))

#define CPOR1p   (1.0/(1.0+(p->B260*(PORVAL1p<1.0?1.0:0.0))))
#define CPOR2p   (1.0/(1.0+(p->B260*(PORVAL2p<1.0?1.0:0.0))))
#define CPOR3p   (1.0/(1.0+(p->B260*(PORVAL3p<1.0?1.0:0.0))))


#define PORVALNH d->porosity[IJK]
#define PORVALNHm d->porosity[IJK]
#define PORVALNHp d->porosity[IJK]
#define CPORNH  (1.0/(1.0+(p->B260*(PORVALNH<1.0?1.0:0.0))))
#define CPORNHm (1.0/(1.0+(p->B260*(PORVALNHm<1.0?1.0:0.0))))
#define CPORNHp (1.0/(1.0+(p->B260*(PORVALNHp<1.0?1.0:0.0))))


//-----------------------------------------------------------


#define IJK  (i-p->imin)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin

#define Im1JK  (i-p->imin-1)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin
#define Ip1JK  (i-p->imin+1)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin
#define IJm1K  (i-p->imin)*p->jmax*p->kmax + (j-p->jmin-1)*p->kmax + k-p->kmin
#define IJp1K  (i-p->imin)*p->jmax*p->kmax + (j-p->jmin+1)*p->kmax + k-p->kmin
#define IJKm1  (i-p->imin)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin-1
#define IJKp1  (i-p->imin)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin+1

#define Ip1Jp1K  (i-p->imin+1)*p->jmax*p->kmax + (j-p->jmin+1)*p->kmax + k-p->kmin

#define Ip1JKp1  (i-p->imin+1)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin+1
#define Ip1JKm1  (i-p->imin+1)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin-1
#define Im1JKp1  (i-p->imin-1)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin+1
#define Im1JKm1  (i-p->imin-1)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin-1

#define IJp1Kp1  (i-p->imin)*p->jmax*p->kmax + (j-p->jmin+1)*p->kmax + k-p->kmin+1
#define IJp1Km1  (i-p->imin)*p->jmax*p->kmax + (j-p->jmin+1)*p->kmax + k-p->kmin-1
#define IJm1Kp1  (i-p->imin)*p->jmax*p->kmax + (j-p->jmin-1)*p->kmax + k-p->kmin+1
#define IJm1Km1  (i-p->imin)*p->jmax*p->kmax + (j-p->jmin-1)*p->kmax + k-p->kmin-1

#define IJp1Kp1  (i-p->imin)*p->jmax*p->kmax + (j-p->jmin+1)*p->kmax + k-p->kmin+1
#define IJp1Kp1  (i-p->imin)*p->jmax*p->kmax + (j-p->jmin+1)*p->kmax + k-p->kmin+1

#define Im1Jm1Km1 (i-p->imin-1)*p->jmax*p->kmax + (j-p->jmin-1)*p->kmax + k-p->kmin-1
#define Ip1Jm1Km1 (i-p->imin+1)*p->jmax*p->kmax + (j-p->jmin-1)*p->kmax + k-p->kmin-1
#define Ip1Jp1Km1 (i-p->imin+1)*p->jmax*p->kmax + (j-p->jmin+1)*p->kmax + k-p->kmin-1
#define Im1Jp1Km1 (i-p->imin-1)*p->jmax*p->kmax + (j-p->jmin+1)*p->kmax + k-p->kmin-1

#define Im1Jm1Kp1 (i-p->imin-1)*p->jmax*p->kmax + (j-p->jmin-1)*p->kmax + k-p->kmin+1
#define Ip1Jm1Kp1 (i-p->imin+1)*p->jmax*p->kmax + (j-p->jmin-1)*p->kmax + k-p->kmin+1
#define Ip1Jp1Kp1 (i-p->imin+1)*p->jmax*p->kmax + (j-p->jmin+1)*p->kmax + k-p->kmin+1
#define Im1Jp1Kp1 (i-p->imin-1)*p->jmax*p->kmax + (j-p->jmin+1)*p->kmax + k-p->kmin+1

#define Ip1Jp2K  (i-p->imin+1)*p->jmax*p->kmax + (j-p->jmin+2)*p->kmax + k-p->kmin

#define Im1JKp2  (i-p->imin-1)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin+2
#define Ip1JKp2  (i-p->imin+1)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin+2

#define IJm1Kp2  (i-p->imin)*p->jmax*p->kmax + (j-p->jmin-1)*p->kmax + k-p->kmin+2
#define IJp1Kp2  (i-p->imin)*p->jmax*p->kmax + (j-p->jmin+1)*p->kmax + k-p->kmin+2


#define Im2JK  (i-p->imin-2)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin
#define Ip2JK  (i-p->imin+2)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin
#define IJm2K  (i-p->imin)*p->jmax*p->kmax + (j-p->jmin-2)*p->kmax + k-p->kmin
#define IJp2K  (i-p->imin)*p->jmax*p->kmax + (j-p->jmin+2)*p->kmax + k-p->kmin
#define IJKm2  (i-p->imin)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin-2
#define IJKp2  (i-p->imin)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin+2

#define Im3JK  (i-p->imin-3)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin
#define Ip3JK  (i-p->imin+3)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin
#define IJm3K  (i-p->imin)*p->jmax*p->kmax + (j-p->jmin-3)*p->kmax + k-p->kmin
#define IJp3K  (i-p->imin)*p->jmax*p->kmax + (j-p->jmin+3)*p->kmax + k-p->kmin
#define IJKm3  (i-p->imin)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin-3
#define IJKp3  (i-p->imin)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin+3

#define Im4JK  (i-p->imin-4)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin
#define Ip4JK  (i-p->imin+4)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin
#define IJKm4  (i-p->imin)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin-4
#define IJKm5  (i-p->imin)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin-5

#define Im5JK  (i-p->imin-5)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin
#define Ip5JK  (i-p->imin+5)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin

#define Im1Jm1K  (i-p->imin-1)*p->jmax*p->kmax + (j-p->jmin-1)*p->kmax + k-p->kmin
#define Im1Jp1K  (i-p->imin-1)*p->jmax*p->kmax + (j-p->jmin+1)*p->kmax + k-p->kmin
#define Ip1Jm1K  (i-p->imin+1)*p->jmax*p->kmax + (j-p->jmin-1)*p->kmax + k-p->kmin
#define Ip1Jp1K  (i-p->imin+1)*p->jmax*p->kmax + (j-p->jmin+1)*p->kmax + k-p->kmin


#define Im1JKm1  (i-p->imin-1)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin-1
#define IJm1Km1  (i-p->imin)*p->jmax*p->kmax + (j-p->jmin-1)*p->kmax + k-p->kmin-1
#define Im1Jm1Km1  (i-p->imin-1)*p->jmax*p->kmax + (j-p->jmin-1)*p->kmax + k-p->kmin-1

#define IJK6 (i-p->imin)*p->jmax6*p->kmax6 + (j-p->jmin)*p->kmax6 + k-p->kmin

#define Im1JK6 (i-p->imin-1)*p->jmax6*p->kmax6 + (j-p->jmin)*p->kmax6 + k-p->kmin
#define Ip1JK6 (i-p->imin+1)*p->jmax6*p->kmax6 + (j-p->jmin)*p->kmax6 + k-p->kmin
#define IJm1K6 (i-p->imin)*p->jmax6*p->kmax6 + (j-p->jmin-1)*p->kmax6 + k-p->kmin
#define IJp1K6 (i-p->imin)*p->jmax6*p->kmax6 + (j-p->jmin+1)*p->kmax6 + k-p->kmin
#define IJKm16 (i-p->imin)*p->jmax6*p->kmax6 + (j-p->jmin)*p->kmax6 + k-p->kmin-1
#define IJKp16 (i-p->imin)*p->jmax6*p->kmax6 + (j-p->jmin)*p->kmax6 + k-p->kmin+1

#define IJK6cv (i-p->imin6)*p->jmax6*p->kmax6 + (j-p->jmin6)*p->kmax6 + k-p->kmin6

#define Im1JK6cv (i-p->imin6-1)*p->jmax6*p->kmax6 + (j-p->jmin6)*p->kmax6 + k-p->kmin6
#define Ip1JK6cv (i-p->imin6+1)*p->jmax6*p->kmax6 + (j-p->jmin6)*p->kmax6 + k-p->kmin6
#define IJm1K6cv (i-p->imin6)*p->jmax6*p->kmax6 + (j-p->jmin6-1)*p->kmax6 + k-p->kmin6
#define IJp1K6cv (i-p->imin6)*p->jmax6*p->kmax6 + (j-p->jmin6+1)*p->kmax6 + k-p->kmin6
#define IJKm16cv (i-p->imin6)*p->jmax6*p->kmax6 + (j-p->jmin6)*p->kmax6 + k-p->kmin6-1
#define IJKp16cv (i-p->imin6)*p->jmax6*p->kmax6 + (j-p->jmin6)*p->kmax6 + k-p->kmin6+1


#define FIJK  (i-p->imin)*p->jmax*p->kmaxF + (j-p->jmin)*p->kmaxF + k-p->kmin

#define FIm1JK  (i-p->imin-1)*p->jmax*p->kmaxF + (j-p->jmin)*p->kmaxF + k-p->kmin
#define FIp1JK  (i-p->imin+1)*p->jmax*p->kmaxF + (j-p->jmin)*p->kmaxF + k-p->kmin
#define FIJm1K  (i-p->imin)*p->jmax*p->kmaxF + (j-p->jmin-1)*p->kmaxF + k-p->kmin
#define FIJp1K  (i-p->imin)*p->jmax*p->kmaxF + (j-p->jmin+1)*p->kmaxF + k-p->kmin
#define FIJKm1  (i-p->imin)*p->jmax*p->kmaxF + (j-p->jmin)*p->kmaxF + k-p->kmin-1
#define FIJKp1  (i-p->imin)*p->jmax*p->kmaxF + (j-p->jmin)*p->kmaxF + k-p->kmin+1

#define FIm2JK  (i-p->imin-2)*p->jmax*p->kmaxF + (j-p->jmin)*p->kmaxF + k-p->kmin
#define FIp2JK  (i-p->imin+2)*p->jmax*p->kmaxF + (j-p->jmin)*p->kmaxF + k-p->kmin
#define FIJm2K  (i-p->imin)*p->jmax*p->kmaxF + (j-p->jmin-2)*p->kmaxF + k-p->kmin
#define FIJp2K  (i-p->imin)*p->jmax*p->kmaxF + (j-p->jmin+2)*p->kmaxF + k-p->kmin
#define FIJKm2  (i-p->imin)*p->jmax*p->kmaxF + (j-p->jmin)*p->kmaxF + k-p->kmin-2
#define FIJKp2  (i-p->imin)*p->jmax*p->kmaxF + (j-p->jmin)*p->kmaxF + k-p->kmin+2

#define FIm3JK  (i-p->imin-3)*p->jmax*p->kmaxF + (j-p->jmin)*p->kmaxF + k-p->kmin
#define FIp3JK  (i-p->imin+3)*p->jmax*p->kmaxF + (j-p->jmin)*p->kmaxF + k-p->kmin
#define FIJm3K  (i-p->imin)*p->jmax*p->kmaxF + (j-p->jmin-3)*p->kmaxF + k-p->kmin
#define FIJp3K  (i-p->imin)*p->jmax*p->kmaxF + (j-p->jmin+3)*p->kmaxF + k-p->kmin
#define FIJKm3  (i-p->imin)*p->jmax*p->kmaxF + (j-p->jmin)*p->kmaxF + k-p->kmin-3
#define FIJKp3  (i-p->imin)*p->jmax*p->kmaxF + (j-p->jmin)*p->kmaxF + k-p->kmin+3
#define FIJKm4  (i-p->imin)*p->jmax*p->kmaxF + (j-p->jmin)*p->kmaxF + k-p->kmin-4
#define FIJKp4  (i-p->imin)*p->jmax*p->kmaxF + (j-p->jmin)*p->kmaxF + k-p->kmin+4
#define FIJKm5  (i-p->imin)*p->jmax*p->kmaxF + (j-p->jmin)*p->kmaxF + k-p->kmin-5
#define FIJKm6  (i-p->imin)*p->jmax*p->kmaxF + (j-p->jmin)*p->kmaxF + k-p->kmin-6

#define FIp1JKp1  (i-p->imin+1)*p->jmax*p->kmaxF + (j-p->jmin)*p->kmaxF + k-p->kmin+1
#define FIm1JKp1  (i-p->imin-1)*p->jmax*p->kmaxF + (j-p->jmin)*p->kmaxF + k-p->kmin+1
#define FIp1JKm1  (i-p->imin+1)*p->jmax*p->kmaxF + (j-p->jmin)*p->kmaxF + k-p->kmin-1
#define FIm1JKm1  (i-p->imin-1)*p->jmax*p->kmaxF + (j-p->jmin)*p->kmaxF + k-p->kmin-1

#define FIJp1Kp1  (i-p->imin)*p->jmax*p->kmaxF + (j-p->jmin+1)*p->kmaxF + k-p->kmin+1
#define FIJm1Kp1  (i-p->imin)*p->jmax*p->kmaxF + (j-p->jmin-1)*p->kmaxF + k-p->kmin+1
#define FIJp1Km1  (i-p->imin)*p->jmax*p->kmaxF + (j-p->jmin+1)*p->kmaxF + k-p->kmin-1
#define FIJm1Km1  (i-p->imin)*p->jmax*p->kmaxF + (j-p->jmin-1)*p->kmaxF + k-p->kmin-1

#define FIp1Jp1K  (i-p->imin+1)*p->jmax*p->kmaxF + (j-p->jmin+1)*p->kmaxF + k-p->kmin
#define FIm1Jp1K  (i-p->imin-1)*p->jmax*p->kmaxF + (j-p->jmin+1)*p->kmaxF + k-p->kmin
#define FIp1Jm1K  (i-p->imin+1)*p->jmax*p->kmaxF + (j-p->jmin-1)*p->kmaxF + k-p->kmin
#define FIm1Jm1K  (i-p->imin-1)*p->jmax*p->kmaxF + (j-p->jmin-1)*p->kmaxF + k-p->kmin

#define FIp1Jp1Kp1  (i-p->imin+1)*p->jmax*p->kmaxF + (j-p->jmin+1)*p->kmaxF + k-p->kmin+1
#define FIm1Jm1Kp1  (i-p->imin-1)*p->jmax*p->kmaxF + (j-p->jmin-1)*p->kmaxF + k-p->kmin+1

#define FIp1Jp2K  (i-p->imin+1)*p->jmax*p->kmaxF + (j-p->jmin+2)*p->kmaxF + k-p->kmin

#define FIm1JKp2  (i-p->imin-1)*p->jmax*p->kmaxF + (j-p->jmin)*p->kmaxF + k-p->kmin+2
#define FIp1JKp2  (i-p->imin+1)*p->jmax*p->kmaxF + (j-p->jmin)*p->kmaxF + k-p->kmin+2

#define FIJm1Kp2  (i-p->imin)*p->jmax*p->kmaxF + (j-p->jmin-1)*p->kmaxF + k-p->kmin+2
#define FIJp1Kp2  (i-p->imin)*p->jmax*p->kmaxF + (j-p->jmin+1)*p->kmaxF + k-p->kmin+2


/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2023 Hans Bihs

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
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"6DOF_gc.h"
#include"mooring_void.h"
#include"mooring_barQuasiStatic.h"
#include"mooring_Catenary.h"
#include"mooring_Spring.h"
#include"mooring_dynamic.h"
#include"net_void.h"
#include"net_barDyn.h"
#include"net_barQuasiStatic.h"
#include"net_sheet.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include<sys/stat.h>

void sixdof_gc::initialize(lexer *p, fdm *a, ghostcell *pgc, vector<net*>& pnet)
{
    if(p->mpirank==0)
    cout<<"6DOF_gc_ini "<<endl;
    
    if(p->X50==1)
	print_ini_vtp(p,a,pgc);
    
    if(p->X50==2)
	print_ini_stl(p,a,pgc);
    
	ini_parameter(p,a,pgc);
	
    objects(p,a,pgc);
	ray_cast(p,a,pgc);
	reini_AB2(p,a,pgc,a->fb);
	geometry_ini(p,a,pgc);
	
    if (p->X13 == 0)
    {
		position_ini(p,a,pgc);
	}
	else if (p->X13 > 0)
	{
		position_ini_quaternion(p,a,pgc);		
	}
	
	ray_cast(p,a,pgc);
	reini_AB2(p,a,pgc,a->fb);
	pgc->start4a(p,a->fb,50);

	interface(p,true);
	maxvel(p,a,pgc);
	pgc->gcfb_update(p,a);
	
    if(p->X50==1)
    print_vtp(p,a,pgc);
    
    if(p->X50==2)
    print_stl(p,a,pgc);
    
    if(p->X221==1)
    read_motionvec(p,a,pgc);
    
	
	if(p->X310==0)
	{
		pmooring.push_back(new mooring_void());
	}
	else
	{
		MPI_Bcast(&p->mooring_count,1,MPI_DOUBLE,0,pgc->mpi_comm);	

		Xme.resize(p->mooring_count);
		Yme.resize(p->mooring_count);
		Zme.resize(p->mooring_count);
		Kme.resize(p->mooring_count);
		Mme.resize(p->mooring_count);
		Nme.resize(p->mooring_count);

		if(p->mpirank==0 && p->P14==1)
		{
			mkdir("./REEF3D_CFD_6DOF_Mooring",0777);	
		}		

		pmooring.reserve(p->mooring_count);
		X311_xen.resize(p->mooring_count,0.0);
		X311_yen.resize(p->mooring_count,0.0);
		X311_zen.resize(p->mooring_count,0.0);
			
		for (int i=0; i < p->mooring_count; i++)
		{
			if(p->X310==1)
			{
				pmooring.push_back(new mooring_Catenary(i));
			}	
			else if(p->X310==2)
			{
				pmooring.push_back(new mooring_barQuasiStatic(i)); 
			}	
			else if(p->X310==3)
			{
				pmooring.push_back(new mooring_dynamic(i));
			}
			else if(p->X310==4)
			{
				pmooring.push_back(new mooring_Spring(i));
			}
		
			X311_xen[i] = p->X311_xe[i] - xg;
			X311_yen[i] = p->X311_ye[i] - yg;
			X311_zen[i] = p->X311_ze[i] - zg;
		
			pmooring[i]->initialize(p,a,pgc);
		}
	}	

    if (p->X320 == 0)
    {
        pnet.push_back(new net_void());
    }
    else
    {
		MPI_Bcast(&p->net_count,1,MPI_DOUBLE,0,pgc->mpi_comm);
		
        Xne.resize(p->net_count);
		Yne.resize(p->net_count);
		Zne.resize(p->net_count);
		Kne.resize(p->net_count);
		Mne.resize(p->net_count);
		Nne.resize(p->net_count);
    
        if(p->mpirank==0)
        {
            if(p->P14==1)
            {
                mkdir("./REEF3D_CFD_6DOF_Net",0777);	
            }
        }
        else
        {
            p->X320_type = new int[p->net_count];
        }
		
        pnet.reserve(p->net_count);	
  
		for (int ii=0; ii < p->net_count; ii++)
		{
            MPI_Bcast(&p->X320_type[ii],1,MPI_INT,0,pgc->mpi_comm);
			
            if(p->X320_type[ii] > 10)
			{
				pnet.push_back(new net_barDyn(ii,p));
			}
            else if (p->X320_type[ii] < 4)
            {
                pnet.push_back(new net_barQuasiStatic(ii,p));
            }
            else
            {
                 pnet.push_back(new net_sheet(ii,p));
            }
			
            pnet[ii]->initialize(p,a,pgc);
		}
    } 
}

void sixdof_gc::ini_parameter(lexer *p, fdm *a, ghostcell *pgc)
{
    R_.resize(3);
    I_.resize(3);
    e_.resize(13,0.0);
    en_.resize(13,0.0);
    enn_.resize(13,0.0);
    ennn_.resize(13,0.0);
    ennnn_.resize(13,0.0);
    trunc_.resize(13,0.0);
	for (int i=0; i<3; i++)
	{
		R_[i].resize(3,0.0);
        I_[i].resize(3,0.0);
    }
	
	if (p->X102 == 1)
	{
		e_[10] = p->X102_u;
		e_[11] = p->X102_v;
		e_[12] = p->X102_w;
	} 
    
	if (p->X103 == 1)
	{
		e_[4] = p->X103_p;
		e_[5] = p->X103_q;
		e_[6] = p->X103_r;
	}  
	
	p->ufbn=p->vfbn=p->wfbn=0.0;
	p->pfbn=p->qfbn=p->rfbn=0.0;
        
    
	if(p->X210==1)
	{
	p->ufbn = p->X210_u;
	p->vfbn = p->X210_v;
	p->wfbn = p->X210_w;
	}
	
	if(p->X211==1)
	{
	p->pfbn = p->X211_p;
	p->qfbn= p->X211_q;
	p->rfbn = p->X211_r;
	}
    
	
	dxg_sum=dxg_sum=dxg_sum=0.0;
	dphi_sum=dtheta_sum=dpsi_sum=0.0;
	phi=theta=psi=0.0;
	printtime=0.0;
	
	Vfb=Mfb=Rfb=0.0;
	xg=yg=zg=0.0;
	xgn=ygn=zgn=0.0;
	xg_s=yg_s=zg_s=0.0;
	xg_sn=yg_sn=zg_sn=0.0;
	dxg=dyg=dzg=0.0;
	Ix=Iy=Iz=0.0;
	xorig=yorig=zorig=0.0;


	dxg_sum=dyg_sum=dzg_sum=0.0;
	dphi_sum=dtheta_sum=dpsi_sum=0.0;
	phi=theta=psi=0.0;
	dphi=dtheta=dpsi=0.0;
	phi1=theta1=psi1=0.0;
	phi2=theta2=psi2=0.0;
	Xe=Ye=Ze=0.0;
	Ke=Me=Ne=0.0;
	Xs=Ys=Zs=0.0;
	Ks=Ms=Ns=0.0;
	
	Ue=Ve=We=0.0;
	Pe=Qe=Re=0.0;
	dUe=dVe=dWe=0.0;
	dPe=dQe=dRe=0.0;
	
	
	Us=Vs=Ws=0.0;
	Ps=Qs=Rs=0.0;
	
	Usn=Vsn=Wsn=0.0;
	Psn=Qsn=Rsn=0.0;
	Usnn=Vsnn=Wsnn=0.0;
	Psnn=Qsnn=Rsnn=0.0;
	Usnnn=Vsnnn=Wsnnn=0.0;
	Psnnn=Qsnnn=Rsnnn=0.0;
	dUs=dVs=dWs=0.0;
	dPs=dQs=dRs=0.0;
	dUsn=dVsn=dWsn=0.0;
	dPsn=dQsn=dRsn=0.0;
	dUsnn=dVsnn=dWsnn=0.0;
	dPsnn=dQsnn=dRsnn=0.0;
	dUsnnn=dVsnnn=dWsnnn=0.0;
	dPsnnn=dQsnnn=dRsnnn=0.0;
	
	phi_s=theta_s=psi_s=0.0;
	phi_sn=theta_sn=psi_sn=0.0;
	phi_en=theta_en=psi_en=0.0;
	Uen=Ven=Wen=0.0;
	Pen=Qen=Ren=0.0;
	dUen=dVen=dWen=0.0;
	dUenn=dVenn=dWenn=0.0;
	dUennn=dVennn=dWennn=0.0;	
	Uenn=Venn=Wenn=0.0;
	Penn=Qenn=Renn=0.0;
	Uennn=Vennn=Wennn=0.0;
	Pennn=Qennn=Rennn=0.0;
	
	Uext=Vext=Wext=0.0;
	Pext=Qext=Rext=0.0;
}

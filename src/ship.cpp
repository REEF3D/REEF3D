/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2026 Hans Bihs

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
Architect: Hans Bihs
--------------------------------------------------------------------*/

#include"ship.h"
#include"ship_hull.h"
#include"ship_models.h"
#include"6DOF_rigidbody.h"
#include"6DOF_geometry.h"
#include"6DOF_output_dir.h"
#include"lexer.h"
#include<iostream>
#include<sstream>
#include<cstdio>
#include<cmath>

namespace
{
const double DEG = 3.14159265358979323846/180.0;

double wrap_angle(double a)
{
    const double pi = 3.14159265358979323846;
    while(a> pi) a -= 2.0*pi;
    while(a<-pi) a += 2.0*pi;
    return a;
}
}

ship::ship(lexer *p, int number) : id(number), initialized(false),
                                   lpp(0.0), S(0.0), k(0.0), nu(p->W2), B44(0.0), B44q(0.0), Cd(0.0),
                                   thrust(0.0), xthrust(0.0), zthrust(0.0),
                                   friction(1), nstrip(40), lpp_in(false), S_in(false),
                                   xa(0.0), xf(0.0), zw(0.0),
                                   prop(false), xp(0.0), yp(0.0), zp(0.0), Dp(0.0), hub(0.2), nrps(0.0), thick(0.0),
                                   wake(0.0), sample_d(0.0), tded(0.0),
                                   sense(1), inflow_mode(0), psource((p->A10==5 || p->A10==6) ? 1 : 0),
                                   Va(0.0), J(0.0), KT(0.0), KQ(0.0), Tp(0.0), Qp(0.0),
                                   rud(false), rmmg_in(false), rmode(0), rcmd(0.0), zz_d(0.0), zz_psi(0.0), zz_t0(0.0),
                                   ap_psi(0.0), ap_Kp(0.0), ap_Kd(0.0), ap_Ki(0.0),
                                   rrate(15.0*DEG), rmax(35.0*DEG),
                                   delta(0.0), tlast(-1.0), psi_c(0.0), psi_last(0.0), psi0(0.0), eint(0.0),
                                   zz_sign(1), zz_started(false),
                                   alphaR(0.0), UR(0.0), FN(0.0), XR(0.0), YR(0.0), NR(0.0), KR(0.0),
                                   mmg(false), mmg_am(false), pwake_mmg(false), mmg_draft_in(false), mmg_fluid(1),
                                   mmg_mx(0.0), mmg_my(0.0), mmg_Jz(0.0), mmg_d(0.0), xm(0.0), Umin(0.0), rho_am(1000.0),
                                   wC1(0.0), wC2p(1.0), wC2n(1.0), wxP(0.0), XH(0.0), YH(0.0), NH(0.0),
                                   ub(0.0), vb(0.0), wb(0.0), pb(0.0), qb(0.0), rb_(0.0),
                                   Re(0.0), CF(0.0), XF(0.0), Ycf(0.0), Ncf(0.0), Kroll(0.0)
{
    for(int q=0; q<3; ++q)
    kt[q] = kq[q] = 0.0;
    
    Ucur[0] = Ucur[1] = 0.0;
    
    // MMG defaults (KVLCC2, Yasukawa & Yoshimura 2015); xH, lR in units of L until ini
    rp.AR = rp.Lambda = rp.xR = rp.zR = 0.0;
    rp.tR = 0.39;
    rp.aH = 0.3;
    rp.xH = -0.45;
    rp.eps = 1.1;
    rp.kappa = 0.5;
    rp.lR = -0.9;
    rp.gammaR = 0.4;
    rp.gammaRp = -1.0;
    rp.falpha = 0.0;
    
    for(int q=0; q<17; ++q)
    mmgc[q] = 0.0;
    
    // CFD resolves the wall shear (viscous momentum equations): no correlation-line friction by default
    if(p->A10==6)
    friction = 0;
    
    read(p);
    
    // NHFLOW X 38 1: the hull loads already contain a local skin friction
    if(friction==1 && p->X38==1)
    {
        friction = 0;
        
        if(p->mpirank==0)
        cout<<"ship: X 38 1 (local skin friction in the hull loads), ITTC-1957 friction switched off"<<endl;
    }
    
    // the actuator disk is coupled to NHFLOW and CFD
    if(psource==1 && p->A10!=5 && p->A10!=6)
    {
        psource = 0;
        
        if(p->mpirank==0)
        cout<<"ship: propeller_source 1 needs NHFLOW or CFD, the thrust acts on the hull only"<<endl;
    }
}

ship::~ship()
{
    if(out.is_open())
    out.close();
}

void ship::read(lexer *p)
{
    ifstream f("ship.dat");
    
    if(!f.is_open())
    {
        if(p->mpirank==0)
        cout<<"ship: X 350 1 but no ship.dat, defaults are used"<<endl;
        
        return;
    }
    
    string line;
    
    while(getline(f,line))
    {
        const size_t c = line.find('#');
        
        if(c!=string::npos)
        line = line.substr(0,c);
        
        istringstream ls(line);
        string key;
        
        if(!(ls>>key))
        continue;
        
        // hull
        if(key=="lpp")
        {
            ls>>lpp;
            lpp_in = true;
        }
        else if(key=="wetted_surface")
        {
            ls>>S;
            S_in = true;
        }
        else if(key=="form_factor")
        ls>>k;
        
        else if(key=="friction")
        ls>>friction;
        
        else if(key=="viscosity")
        ls>>nu;
        
        else if(key=="roll_damping")
        ls>>B44>>B44q;
        
        else if(key=="crossflow")
        ls>>Cd;
        
        else if(key=="strips")
        ls>>nstrip;
        
        else if(key=="current")
        ls>>Ucur[0]>>Ucur[1];
        
        else if(key=="thrust")
        {
            ls>>thrust;
            
            if(!(ls>>xthrust>>zthrust))
            xthrust = zthrust = 0.0;
        }
        
        // propeller
        else if(key=="propeller")
        {
            ls>>xp>>yp>>zp>>Dp;
            prop = Dp>0.0;
        }
        
        else if(key=="propeller_hub")
        ls>>hub;
        
        else if(key=="propeller_kt")
        ls>>kt[0]>>kt[1]>>kt[2];
        
        else if(key=="propeller_kq")
        ls>>kq[0]>>kq[1]>>kq[2];
        
        else if(key=="propeller_rps")
        ls>>nrps;
        
        else if(key=="propeller_sense")
        ls>>sense;
        
        else if(key=="propeller_thickness")
        ls>>thick;
        
        else if(key=="propeller_inflow")
        {
            string mode;
            ls>>mode;
            
            if(mode=="sample")
            {
                inflow_mode = 1;
                ls>>sample_d;
            }
            else
            {
                inflow_mode = 0;
                ls>>wake;
            }
        }
        
        else if(key=="propeller_source")
        ls>>psource;
        
        else if(key=="thrust_deduction")
        ls>>tded;
        
        // rudder
        else if(key=="rudder")
        {
            ls>>rp.xR>>rp.zR>>rp.AR>>rp.Lambda;
            rud = rp.AR>0.0 && rp.Lambda>0.0;
        }
        
        else if(key=="rudder_mmg")
        {
            ls>>rp.tR>>rp.aH>>rp.xH>>rp.eps>>rp.kappa>>rp.lR>>rp.gammaR;
            rmmg_in = true;
        }
        
        else if(key=="rudder_angle")
        {
            string mode;
            ls>>mode;
            
            if(mode=="zigzag")
            {
                rmode = 1;
                ls>>zz_d>>zz_psi;
                
                if(!(ls>>zz_t0))
                zz_t0 = 0.0;
                
                zz_d *= DEG;
                zz_psi *= DEG;
            }
            else if(mode=="autopilot")
            {
                rmode = 2;
                ls>>ap_psi>>ap_Kp>>ap_Kd;
                
                if(!(ls>>ap_Ki))
                ap_Ki = 0.0;
                
                ap_psi *= DEG;
            }
            else
            {
                rmode = 0;
                ls>>rcmd;
                rcmd *= DEG;
            }
        }
        
        else if(key=="rudder_rate")
        {
            ls>>rrate;
            rrate *= DEG;
        }
        
        else if(key=="rudder_max")
        {
            ls>>rmax;
            rmax *= DEG;
        }
        
        else if(key=="rudder_falpha")
        ls>>rp.falpha;
        
        else if(key=="rudder_gamma")
        ls>>rp.gammaR>>rp.gammaRp;
        
        // MMG manoeuvring model
        else if(key=="mmg_hull")
        {
            for(int q=0; q<17; ++q)
            ls>>mmgc[q];
            
            mmg = true;
        }
        
        else if(key=="mmg_added_mass")
        {
            ls>>mmg_mx>>mmg_my>>mmg_Jz;
            mmg_am = true;
        }
        
        else if(key=="mmg_fluid")
        ls>>mmg_fluid;
        
        else if(key=="mmg_draft")
        {
            ls>>mmg_d;
            mmg_draft_in = true;
        }
        
        else if(key=="propeller_wake_mmg")
        {
            ls>>wC1>>wC2p>>wC2n>>wxP;
            pwake_mmg = true;
        }
        
        else if(p->mpirank==0)
        cout<<"ship: unknown keyword in ship.dat: "<<key<<endl;
    }
}

void ship::ini(lexer *p, const sixdof_rigidbody &b, const sixdof_geometry &g)
{
    // still water level in the body frame of the initial position (level hull)
    zw = p->F60 - b.c(2);
    
    ship_hull::waterline_extent(g.tri_x0,g.tri_y0,g.tri_z0,g.tricount,zw,xa,xf);
    
    if(!lpp_in)
    lpp = xf - xa;
    
    if(!S_in)
    S = ship_hull::wetted_surface(g.tri_x0,g.tri_y0,g.tri_z0,g.tricount,zw,p->j_dir==0);
    
    if(Cd>0.0)
    ship_hull::draft_strips(g.tri_x0,g.tri_y0,g.tri_z0,g.tricount,zw,xa,xf,nstrip,xs,dx,T);
    
    // midship: middle of the waterline, in the CoG frame
    xm = 0.5*(xa + xf);
    
    // MMG positions in units of L; xH from midship
    rp.xH = rp.xH*lpp + xm;
    rp.lR *= lpp;
    
    // MMG draft: deepest hull point below the still water level
    if(!mmg_draft_in)
    {
        double zmin = zw;
        
        for(int n=0; n<g.tricount; ++n)
        for(int q=0; q<3; ++q)
        zmin = g.tri_z0[n][q]<zmin ? g.tri_z0[n][q] : zmin;
        
        mmg_d = zw - zmin;
    }
    
    rho_am = p->W1;
    
    // velocity limit in the denominators of the MMG polynomials
    Umin = 0.01*sqrt(9.81*(lpp>0.0 ? lpp : 1.0));
    
    // actuator disk thickness: resolved by at least 4 cells
    if(prop && thick<=0.0)
    thick = 0.2*Dp>4.0*p->DXM ? 0.2*Dp : 4.0*p->DXM;
    
    // heading
    psi_last = psi_c = b.psi;
    
    if(p->mpirank==0)
    {
        cout<<"ship "<<id<<": L = "<<lpp<<" m, S = "<<S<<" m^2, waterline x = "<<xa<<" .. "<<xf<<" m (from the CoG)"<<endl;
        cout<<"ship "<<id<<": friction "<<friction<<" (1+k) = "<<1.0+k<<", nu = "<<nu
            <<", roll damping "<<B44<<" "<<B44q<<", cross-flow Cd = "<<Cd<<", thrust "<<thrust<<endl;
        
        if(prop)
        cout<<"ship "<<id<<": propeller D = "<<Dp<<" m at ("<<xp<<", "<<yp<<", "<<zp<<"), n = "<<nrps<<" 1/s, disk thickness "<<thick
            <<" m, inflow "<<(inflow_mode==1 ? "sampled" : "wake fraction")<<", actuator disk in the fluid "<<psource<<endl;
        
        if(rud)
        cout<<"ship "<<id<<": rudder AR = "<<rp.AR<<" m^2, Lambda = "<<rp.Lambda<<" at xR = "<<rp.xR<<" m, mode "
            <<(rmode==0 ? "fixed" : (rmode==1 ? "zigzag" : "autopilot"))<<endl;
        
        if(mmg || mmg_am)
        cout<<"ship "<<id<<": MMG hull "<<mmg<<", added mass "<<mmg_am<<", L = "<<lpp<<" m, d = "<<mmg_d
            <<" m, midship at x = "<<xm<<" m from the CoG, fluid loads in surge/sway/yaw "<<(mmg_fluid==1 ? "on" : "off")<<endl;
    }
    
    initialized = true;
}

void ship::steering(lexer *p, double psi, double r)
{
    // rudder machine: once per time step
    if(!rud || p->simtime==tlast)
    return;
    
    const double dt = tlast<0.0 ? 0.0 : p->simtime - tlast;
    tlast = p->simtime;
    
    double cmd = 0.0;
    
    if(rmode==0)
    cmd = rcmd;
    
    if(rmode==1 && p->simtime>=zz_t0)
    {
        // zig-zag: starboard rudder first (heading decreases), counter rudder at -psi / +psi
        if(!zz_started)
        {
            zz_started = true;
            psi0 = psi;
            zz_sign = 1;
        }
        
        const double dpsi = psi - psi0;
        
        if(zz_sign>0 && dpsi<=-zz_psi)
        zz_sign = -1;
        
        else if(zz_sign<0 && dpsi>=zz_psi)
        zz_sign = 1;
        
        cmd = double(zz_sign)*zz_d;
    }
    
    if(rmode==2)
    {
        // heading error > 0: turn to port, i.e. rudder to port (delta < 0)
        const double e = wrap_angle(ap_psi - psi);
        eint += e*dt;
        cmd = -ap_Kp*e + ap_Kd*r - ap_Ki*eint;
    }
    
    cmd = cmd> rmax ?  rmax : cmd;
    cmd = cmd<-rmax ? -rmax : cmd;
    
    // rudder rate
    const double dmax = rrate*dt;
    const double dd = cmd - delta;
    
    delta += dd>dmax ? dmax : (dd<-dmax ? -dmax : dd);
}

double ship::propeller_inflow(const sixdof_rigidbody &b, sixdof_fluid *fluid, const Eigen::Vector3d &centre, const Eigen::Vector3d &axis)
{
    // relative axial inflow on three rings (0.4, 0.7, 0.9 R) of 8 points, sample_d ahead of the disk
    Eigen::Matrix<double,6,1> u6;
    b.velocity(u6);
    const Eigen::Vector3d U(u6(0),u6(1),u6(2)), W(u6(3),u6(4),u6(5));
    
    // two unit vectors normal to the axis
    Eigen::Vector3d e1 = axis.cross(Eigen::Vector3d(0.0,0.0,1.0));
    
    if(e1.norm()<1.0e-6)
    e1 = axis.cross(Eigen::Vector3d(0.0,1.0,0.0));
    
    e1.normalize();
    const Eigen::Vector3d e2 = axis.cross(e1);
    
    const double rr[3] = {0.4,0.7,0.9};
    const int npt = 24;
    double xyz[3*npt], uvw[3*npt];
    
    for(int qr=0; qr<3; ++qr)
    for(int qa=0; qa<8; ++qa)
    {
        const double ang = 2.0*3.14159265358979323846*double(qa)/8.0;
        const Eigen::Vector3d x = centre + sample_d*axis + 0.5*Dp*rr[qr]*(cos(ang)*e1 + sin(ang)*e2);
        const int q = 8*qr + qa;
        
        xyz[3*q] = x(0);
        xyz[3*q+1] = x(1);
        xyz[3*q+2] = x(2);
    }
    
    fluid->velocity(npt,xyz,uvw);
    
    double va=0.0;
    
    for(int q=0; q<npt; ++q)
    {
        const Eigen::Vector3d x(xyz[3*q],xyz[3*q+1],xyz[3*q+2]);
        const Eigen::Vector3d ubody = U + W.cross(x - b.c);
        const Eigen::Vector3d uf(uvw[3*q],uvw[3*q+1],uvw[3*q+2]);
        
        va += (ubody - uf).dot(axis);
    }
    
    return va/double(npt);
}

void ship::add_load(lexer *p, const sixdof_rigidbody &b, const sixdof_geometry &g, sixdof_fluid *fluid, double *F)
{
    if(!initialized)
    ini(p,b,g);
    
    // velocity of the CoG and angular velocity in the ship frame
    Eigen::Matrix<double,6,1> u6;
    b.velocity(u6);
    
    // velocity relative to the water: a uniform current (inertial frame) is subtracted
    const Eigen::Vector3d uI(u6(0) - Ucur[0], u6(1) - Ucur[1], u6(2));
    const Eigen::Vector3d wI(u6(3),u6(4),u6(5));
    const Eigen::Vector3d u = b.R.transpose()*uI;
    const Eigen::Vector3d w = b.R.transpose()*wI;
    
    ub = u(0); vb = u(1); wb = u(2);
    pb = w(0); qb = w(1); rb_ = w(2);
    
    // continuous heading
    psi_c += wrap_angle(b.psi - psi_last);
    psi_last = b.psi;
    
    // loads in the ship frame
    Eigen::Vector3d Fs(0.0,0.0,0.0), Ms(0.0,0.0,0.0);
    
    // hull
    XF = 0.0;
    if(friction==1 && S>0.0 && lpp>0.0)
    XF = ship_models::friction(p->W1,nu,S,lpp,k,ub,Re,CF);
    
    Ycf = Ncf = 0.0;
    if(Cd>0.0)
    ship_models::crossflow(p->W1,Cd,xs,dx,T,vb,rb_,Ycf,Ncf);
    
    Kroll = ship_models::roll_damping(B44,B44q,pb);
    
    Fs(0) = XF + thrust;
    Fs(1) = Ycf;
    
    Ms(0) = Kroll;
    Ms(1) = zthrust*thrust;
    Ms(2) = Ncf;
    
    // MMG hull forces about midship (xG = -xm: CoG ahead of midship), moved to the CoG
    const double xG = -xm;
    const double vm = vb - xG*rb_;
    
    XH = YH = NH = 0.0;
    
    if(mmg)
    {
        ship_models::mmg_hull(p->W1,lpp,mmg_d,mmgc,ub,vm,rb_,Umin,XH,YH,NH);
        
        Fs(0) += XH;
        Fs(1) += YH;
        Ms(2) += NH - xG*YH;
    }
    
    // MMG added mass: velocity terms of the equations of motion (the acceleration terms are
    // solved by the coupling with added_mass()); the Munk moment is part of N'v
    if(mmg_am)
    {
        const double fm = 0.5*p->W1*lpp*lpp*mmg_d;
        const double mx = mmg_mx*fm, my = mmg_my*fm;
        
        Fs(0) += my*vm*rb_ - mx*vb*rb_;
        Fs(1) += (my - mx)*ub*rb_;
        Ms(2) += xG*(mx - my)*ub*rb_;
    }
    
    // propeller
    if(prop)
    {
        const Eigen::Vector3d pos(xp,yp,zp), axb(1.0,0.0,0.0);
        const Eigen::Vector3d centre = b.c + b.R*pos;
        const Eigen::Vector3d axis = b.R*axb;
        
        if(inflow_mode==1 && fluid!=nullptr)
        Va = propeller_inflow(b,fluid,centre,axis);
        
        else if(pwake_mmg)
        {
            const double U = sqrt(ub*ub + vm*vm);
            const double betaP = U>1.0e-10 ? atan2(vm,ub) + wxP*rb_*lpp/U : 0.0;
            
            Va = (1.0 - ship_models::mmg_wake(wake,wC1,wC2p,wC2n,betaP))*ub;
        }
        
        else
        Va = (1.0-wake)*ub;
        
        ship_models::propeller(p->W1,nrps,Dp,kt,kq,Va,J,KT,KQ,Tp,Qp);
        
        // thrust on the hull: with the actuator disk the thrust deduction comes from the flow
        const double Th = psource==1 ? Tp : (1.0-tded)*Tp;
        const Eigen::Vector3d Fp = Th*axb;
        
        Fs += Fp;
        Ms += pos.cross(Fp);
        
        // shaft torque reaction on the hull
        Ms += -double(sense)*Qp*axb;
        
        disk.centre = centre;
        disk.axis = axis;
        disk.R = 0.5*Dp;
        disk.Rh = hub*0.5*Dp;
        disk.thickness = thick;
        disk.T = Tp;
        disk.Q = Qp;
        disk.sense = sense;
    }
    
    // rudder
    if(rud)
    {
        steering(p,psi_c,rb_);
        
        const double uP = prop ? Va : (1.0-wake)*ub;
        
        ship_models::rudder_mmg(p->W1,rp,ub,vm,rb_,delta,prop ? Dp : 0.0,nrps,KT,uP,XR,YR,NR,KR,alphaR,UR,FN);
        
        Fs(0) += XR;
        Fs(1) += YR;
        Ms(0) += KR;
        Ms(2) += NR;
    }
    
    // inertial frame
    const Eigen::Vector3d FI = b.R*Fs;
    const Eigen::Vector3d MI = b.R*Ms;
    
    for(int n=0; n<3; ++n)
    {
        F[n]   += FI(n);
        F[n+3] += MI(n);
    }
}

bool ship::added_mass(const sixdof_rigidbody &b, Eigen::Matrix<double,6,6> &A) const
{
    // MMG added masses about midship moved to the CoG (kinetic energy with vm = v - xG r):
    // surge mx, sway my, yaw Jz + my xG^2, sway-yaw -my xG
    if(!mmg_am || !initialized)
    return false;
    
    const double f = 0.5*rho_am*lpp*lpp*mmg_d;
    const double mx = mmg_mx*f, my = mmg_my*f, Jz = mmg_Jz*f*lpp*lpp;
    const double xG = -xm;
    
    A(0,0) += mx;
    A(1,1) += my;
    A(5,5) += Jz + my*xG*xG;
    A(1,5) += -my*xG;
    A(5,1) += -my*xG;
    
    return true;
}

void ship::fluid_mask(double *w) const
{
    if((mmg || mmg_am) && mmg_fluid==0)
    w[0] = w[1] = w[5] = 0.0;
}

void ship::actuator_disks(vector<sixdof_actuator_disk> &d) const
{
    if(prop && psource==1 && initialized)
    d.push_back(disk);
}

void ship::print(lexer *p)
{
    if(p->mpirank!=0 || p->count%p->X19!=0)
    return;
    
    if(!out.is_open())
    {
        char name[1000];
        snprintf(name,sizeof(name),"%s/REEF3D_ship_%i.dat",sixdof_output_dir(p),id);
        out.open(name);
        out<<"time \t u [m/s] \t v [m/s] \t r [rad/s] \t psi [deg] \t Re \t C_F \t X_F [N] \t Y_cf [N] \t N_cf [Nm] \t K_roll [Nm] \t T_const [N]"
           <<" \t n [1/s] \t Va [m/s] \t J \t KT \t KQ \t T [N] \t Q [Nm]"
           <<" \t delta [deg] \t alpha_R [deg] \t U_R [m/s] \t F_N [N] \t X_R [N] \t Y_R [N] \t N_R [Nm]";
        
        if(mmg)
        out<<" \t X_H [N] \t Y_H [N] \t N_H [Nm] (midship)";
        
        out<<endl;
    }
    
    out<<p->simtime<<" \t "<<ub<<" \t "<<vb<<" \t "<<rb_<<" \t "<<psi_c/DEG<<" \t "<<Re<<" \t "<<CF<<" \t "<<XF<<" \t "<<Ycf<<" \t "<<Ncf<<" \t "<<Kroll<<" \t "<<thrust
       <<" \t "<<nrps<<" \t "<<Va<<" \t "<<J<<" \t "<<KT<<" \t "<<KQ<<" \t "<<Tp<<" \t "<<Qp
       <<" \t "<<delta/DEG<<" \t "<<alphaR/DEG<<" \t "<<UR<<" \t "<<FN<<" \t "<<XR<<" \t "<<YR<<" \t "<<NR;
    
    if(mmg)
    out<<" \t "<<XH<<" \t "<<YH<<" \t "<<NH;
    
    out<<endl;
}

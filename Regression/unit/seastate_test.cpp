// Architect: Hans Bihs
// Standalone verification of the REEF3D::SEASTATE kernels: spectral grid, block-sparse action
// storage, integrated wave parameters, dispersion relation, SWAN spectrum files, source terms (incl. vegetation), surfbeat boundary generator, forcing files,
// structure formulas (d'Angremond, porous), the DIA of a sweep window (Phase 8). No MPI, no REEF3D binary.
// Build:  g++ -O2 -std=c++20 -I../../src seastate_test.cpp ../../src/seastate_grid.cpp ../../src/seastate_store.cpp ../../src/seastate_param.cpp ../../src/seastate_dispersion.cpp ../../src/seastate_swan_spc.cpp ../../src/seastate_source.cpp ../../src/seastate_surfbeat.cpp ../../src/seastate_forcing.cpp ../../src/seastate_bathy.cpp ../../src/seastate_structure.cpp -o seastate_test
// Run:    ./seastate_test
#include"seastate_bathy.h"
#include"seastate_grid.h"
#include"seastate_store.h"
#include"seastate_param.h"
#include"seastate_dispersion.h"
#include"seastate_swan_spc.h"
#include"seastate_source.h"
#include"seastate_surfbeat.h"
#include"seastate_forcing.h"
#include"seastate_structure.h"
#include<cstdio>
#include<fstream>
#include<cmath>
#include<iostream>
#include<string>
#include<vector>

static int nfail = 0;
static void check(bool ok, const std::string &what)
{
    std::cout<<(ok ? "  PASS  " : "  FAIL  ")<<what<<std::endl;
    if(!ok) ++nfail;
}

static const double pi = 3.14159265358979323846;
static const double g  = 9.81;

static bool close(double a, double b, double rtol) {return std::fabs(a-b)<=rtol*std::fabs(b);}

// ---------------------------------------------------------------------------------------------
// spectral grid
// ---------------------------------------------------------------------------------------------
static void test_grid()
{
    std::cout<<"seastate_grid"<<std::endl;

    seastate_grid sg(32,0.04,1.0,36);
    check(sg.valid(),"valid grid for 32 x 36 bins");
    check(sg.nbin==32*36,"nbin = nsig*ndir");
    check(sg.f[0]==0.04 && sg.f[31]==1.0,"end frequencies are fmin and fmax exactly");
    check(close(sg.ratio,std::pow(25.0,1.0/31.0),1e-14),"logarithmic ratio (fmax/fmin)^(1/(nsig-1))");

    bool mono=true;
    for(int l=1;l<sg.nsig;++l)
        mono = mono && close(sg.f[l]/sg.f[l-1],sg.ratio,1e-12);
    check(mono,"constant ratio between neighbouring frequencies");

    double sum=0.0;
    for(int l=0;l<sg.nsig;++l) sum+=sg.dsig[l];
    const double sr=std::sqrt(sg.ratio);
    check(close(sum,sg.sig[31]*sr-sg.sig[0]/sr,1e-12),"bin widths tile [sig_min/sqrt(r), sig_max*sqrt(r)] without gaps");
    check(close(sg.dsig[10],sg.sig[10]*(sr-1.0/sr),1e-12),"dsig_l = sig_l (sqrt(r) - 1/sqrt(r))");

    check(close(sg.dtheta*sg.ndir,2.0*pi,1e-14),"directions cover the full circle");
    check(sg.theta[0]==0.0 && close(sg.theta[9],0.5*pi,1e-14),"theta_0 = 0 (+x), theta_9 = 90 deg (+y) for 36 directions");
    check(close(sg.costh[18],-1.0,1e-14) && std::fabs(sg.sinth[18])<1e-14,"cos/sin tables");
    check(sg.bin(2,5)==2*36+5,"bin layout l*ndir + m (direction fastest)");
    check(close(sg.costhf[0],std::cos(5.0*pi/180.0),1e-14) && close(sg.sinthf[35],std::sin(-5.0*pi/180.0),1e-12),"direction faces at theta_m + dtheta/2 (periodic)");
    check(sg.quad[0]==0 && sg.quad[8]==0 && sg.quad[9]==1 && sg.quad[18]==2 && sg.quad[27]==3 && sg.quad[35]==3,"quadrants [q 90, (q+1) 90) deg");

    check(!seastate_grid(1,0.04,1.0,36).valid(),"nsig < 2 rejected");
    check(!seastate_grid(32,0.5,0.1,36).valid(),"fmax <= fmin rejected");
    check(!seastate_grid(32,0.0,1.0,36).valid(),"fmin <= 0 rejected");
    check(!seastate_grid(32,0.04,1.0,2).valid(),"ndir < 4 rejected");
}

// ---------------------------------------------------------------------------------------------
// block-sparse storage
// ---------------------------------------------------------------------------------------------
static void test_store()
{
    std::cout<<"seastate_store"<<std::endl;

    // rank range with 3 ghost layers: 37 x 22 cells (not multiples of the tile size)
    const int imin=-3, jmin=-3, ni=37, nj=22, nbin=12, tile=8;
    std::vector<int> mask(ni*nj,1);

    // land in i >= 16 (global index i-imin >= 19) for all j: tiles 3 and 4 in i are partly / fully land
    for(int ii=0;ii<ni;++ii)
    for(int jj=0;jj<nj;++jj)
        if(ii+imin>=16) mask[ii*nj+jj]=-10;

    seastate_store st(imin,jmin,ni,nj,nbin,tile);
    st.build(mask.data());

    const int ntx=(ni+tile-1)/tile, nty=(nj+tile-1)/tile;
    check(st.tiles_total()==ntx*nty,"tile count covers the padded range");

    // active cells: ii = 0..18 -> tiles ti 0,1,2 (ii 0..23) contain sea; ti 3,4 (ii 24..36) all land
    check(st.tiles_allocated()==3*nty,"only tiles with sea cells are allocated");
    check(st.cells_active()==long(19)*nj,"active cell count");
    check(st.cells_allocated()==long(24)*nj,"allocated cells: 3 tile rows of 8 cells x the full width (tiles clipped at the range end)");
    check(st.bytes()<st.bytes_dense(),"sparse storage smaller than dense");

    check(st.spec(21,0)==nullptr && st.spec(33,18)==nullptr,"land tile has no storage");
    check(st.spec(17,0)!=nullptr && !st.active(17,0),"land cell in an allocated tile has storage but is inactive");
    check(st.spec(imin,jmin)!=nullptr && st.active(imin,jmin),"ghost corner cell is stored");
    check(st.spec(imin-1,0)==nullptr && st.spec(0,jmin+nj)==nullptr,"outside the range: nullptr");

    // every allocated cell gets its own non-overlapping, contiguous spectrum
    for(int i=imin;i<imin+ni;++i)
    for(int j=jmin;j<jmin+nj;++j)
    {
        float *s=st.spec(i,j);
        if(s)
            for(int b=0;b<nbin;++b) s[b]=float(1000*(i-imin)+10*(j-jmin))+0.01f*b;
    }
    bool ok=true;
    for(int i=imin;i<imin+ni;++i)
    for(int j=jmin;j<jmin+nj;++j)
    {
        const float *s=st.spec(i,j);
        if(s)
            for(int b=0;b<nbin;++b) ok = ok && s[b]==float(1000*(i-imin)+10*(j-jmin))+0.01f*b;
    }
    check(ok,"spectra do not overlap (write/read back every bin of every cell)");

    st.fill(2.5f);
    check(st.spec(0,0)[nbin-1]==2.5f && st.spec(17,0)[0]==0.0f,"fill sets active cells, keeps inactive cells at zero");

    // all land: nothing allocated
    std::vector<int> land(ni*nj,-10);
    st.build(land.data());
    check(st.tiles_allocated()==0 && st.cells_active()==0 && st.spec(0,0)==nullptr,"all land: no spectral memory");

    // tile 1: one tile per cell, allocation equals active cells
    seastate_store s1(imin,jmin,ni,nj,nbin,1);
    s1.build(mask.data());
    check(s1.cells_allocated()==s1.cells_active(),"tile size 1: allocated cells = active cells");
}

// ---------------------------------------------------------------------------------------------
// integrated parameters
// ---------------------------------------------------------------------------------------------
// JONSWAP / PM in omega (independent of the REEF3D wave library)
static double jonswap(double w, double Hs, double Tp, double gamma)
{
    const double wp=2.0*pi/Tp;
    const double s=(w<=wp) ? 0.07 : 0.09;
    double S=(5.0/16.0)*Hs*Hs*std::pow(wp,4.0)*std::pow(w,-5.0)*std::exp(-1.25*std::pow(w/wp,-4.0));
    if(gamma>1.0)
        S*=(1.0-0.287*std::log(gamma))*std::pow(gamma,std::exp(-0.5*std::pow((w-wp)/(s*wp),2.0)));
    return S;
}

// N(sig,theta) = S(sig) D(theta) / sig, D = cos^2s((theta-theta0)/2), normalised on the discrete grid
static std::vector<float> make_spectrum(const seastate_grid &sg, double Hs, double Tp, double gamma, double theta0, double s)
{
    std::vector<double> D(sg.ndir);
    double sum=0.0;
    for(int m=0;m<sg.ndir;++m)
    {
        D[m]=std::pow(std::fabs(std::cos(0.5*(sg.theta[m]-theta0))),2.0*s);
        sum+=D[m]*sg.dtheta;
    }
    std::vector<float> N(sg.nbin);
    for(int l=0;l<sg.nsig;++l)
    for(int m=0;m<sg.ndir;++m)
        N[sg.bin(l,m)]=float(jonswap(sg.sig[l],Hs,Tp,gamma)*D[m]/sum/sg.sig[l]);
    return N;
}

static void test_param()
{
    std::cout<<"seastate_param"<<std::endl;

    seastate_grid sg(40,0.03,1.0,36);
    seastate_param sp;

    // zero spectrum
    std::vector<float> zero(sg.nbin,0.0f);
    sp.compute(sg,zero.data());
    check(sp.Hs==0.0 && sp.Tp==0.0 && sp.lpeak==-1,"zero spectrum: zero parameters");
    sp.compute(sg,nullptr);
    check(sp.Hs==0.0,"nullptr spectrum: zero parameters");

    // single bin: exact
    std::vector<float> one(sg.nbin,0.0f);
    one[sg.bin(12,9)]=0.5f;
    sp.compute(sg,one.data());
    check(close(sp.Hs,4.0*std::sqrt(sg.sig[12]*0.5*sg.dsig[12]*sg.dtheta),1e-6),"single bin: Hs = 4 sqrt(sig N dsig dtheta)");
    check(close(sp.Tp,2.0*pi/sg.sig[12],1e-12) && close(sp.Tm01,sp.Tp,1e-12) && close(sp.Tm10,sp.Tp,1e-12),"single bin: Tp = Tm01 = Tm-1,0 = 2 pi/sig");
    check(close(sp.dir,90.0,1e-9) && sp.spread<1e-3,"single bin: direction 90 deg, no spread");

    // JONSWAP, Hs 2 m, Tp 10 s, gamma 3.3, main direction 30 deg (bin 3), s = 10
    const double s=10.0;
    std::vector<float> N=make_spectrum(sg,2.0,10.0,3.3,30.0*pi/180.0,s);
    sp.compute(sg,N.data());
    std::cout<<"        JONSWAP: Hs "<<sp.Hs<<" Tp "<<sp.Tp<<" Tm01 "<<sp.Tm01<<" Tm-10 "<<sp.Tm10<<" dir "<<sp.dir<<" spread "<<sp.spread<<std::endl;
    check(close(sp.Hs,2.0,0.01),"JONSWAP: Hs within 1 % of the input (discretisation and tail)");
    check(std::fabs(std::log(sp.Tp/10.0))<=0.5*std::log(sg.ratio)+1e-12,"JONSWAP: discrete Tp is the bin nearest to 10 s");
    check(close(sp.dir,30.0,1e-6),"JONSWAP: mean direction = main direction");
    const double spread_exact=std::sqrt(2.0/(s+1.0))*180.0/pi;    // cos^2s(theta/2): a1 = s/(s+1)
    check(close(sp.spread,spread_exact,0.02),"cos^2s spreading: spread = sqrt(2/(s+1)) within 2 %");

    // Pierson-Moskowitz: Tm01/Tp = 0.772, Tm-1,0/Tp = 0.857
    seastate_grid sg2(60,0.02,2.0,24);
    std::vector<float> P=make_spectrum(sg2,3.0,12.0,1.0,0.0,2.0);
    sp.compute(sg2,P.data());
    std::cout<<"        PM:      Hs "<<sp.Hs<<" Tm01/Tp "<<sp.Tm01/12.0<<" Tm-10/Tp "<<sp.Tm10/12.0<<std::endl;
    check(close(sp.Hs,3.0,0.01),"PM: Hs within 1 %");
    check(close(sp.Tm01/12.0,0.772,0.01),"PM: Tm01 = 0.772 Tp within 1 %");
    check(close(sp.Tm10/12.0,0.857,0.01),"PM: Tm-1,0 = 0.857 Tp within 1 %");
    check(sp.dir<1e-6 || sp.dir>360.0-1e-6,"PM: direction 0 deg (+x)");
}

// ---------------------------------------------------------------------------------------------
// memory budget (plan Section 6): 2000 x 2000 cells, 36 x 36 bins, 50 % land, float32
// ---------------------------------------------------------------------------------------------
static void test_memory()
{
    std::cout<<"memory budget"<<std::endl;

    // scaled-down check of the sparse saving: 200 x 200 cells, half land, tile 16
    const int n=200, nbin=36*36, tile=16;
    std::vector<int> mask(n*n,1);
    for(int i=0;i<n;++i)
    for(int j=0;j<n;++j)
        if(i>=n/2) mask[i*n+j]=-10;

    seastate_store st(0,0,n,n,nbin,tile);
    st.build(mask.data());
    const double frac=double(st.bytes())/double(st.bytes_dense());
    const double pad=double(st.cells_allocated())/double(st.cells_active());
    std::cout<<"        200x200 cells, 50 % land: "<<st.bytes()/1048576.0<<" MB of "<<st.bytes_dense()/1048576.0<<" MB dense ("<<frac*100.0<<" %), allocated/active cells "<<pad<<std::endl;
    check(st.tiles_allocated()==7*13,"half land: 7 x 13 tiles of 16 x 16 (coast tile and edge padding included)");
    check(frac<0.6 && pad<1.2,"half land: < 60 % of dense, tile overhead < 20 % of the active cells");

    const double full=2000.0*2000.0*36.0*36.0*4.0/1.0e9;
    check(close(full*0.5,10.368,1e-3),"plan example: 1 copy float32, 50 % land = 10.4 GB");
}

// ---------------------------------------------------------------------------------------------
// dispersion relation
// ---------------------------------------------------------------------------------------------
static void test_dispersion()
{
    std::cout<<"seastate_dispersion"<<std::endl;

    // deep water
    const double sd=2.0*pi/5.0;
    check(close(seastate_wavenumber(sd,1000.0),sd*sd/g,1e-12),"deep water: k = sig^2/g");
    check(close(seastate_cg(sd,sd*sd/g,1000.0),0.5*g/sd,1e-12),"deep water: cg = g/(2 sig)");

    // shallow water
    const double ss=2.0*pi/100.0, ks=seastate_wavenumber(ss,1.0);
    check(close(ks,ss/std::sqrt(g),1e-4),"shallow water: k = sig/sqrt(g d)");
    check(close(seastate_cg(ss,ks,1.0),std::sqrt(g),5e-4),"shallow water: cg = sqrt(g d) (kd = 0.02)");

    // T = 10 s, d = 10 m: L = 92.32 m
    const double s10=2.0*pi/10.0, k10=seastate_wavenumber(s10,10.0);
    check(close(2.0*pi/k10,92.374,2e-4),"T 10 s, d 10 m: L = 92.37 m (g = 9.81)");

    // residual and n over a range of depths and frequencies
    double res=0.0, nerr=0.0;
    for(double d : {0.1,1.0,5.0,20.0,100.0,4000.0})
    for(double f=0.03; f<2.0; f*=1.17)
    {
        const double sg_=2.0*pi*f, k=seastate_wavenumber(sg_,d);
        if(k*d<=30.0)
            res=std::max(res,std::fabs(g*k*std::tanh(k*d)-sg_*sg_)/(sg_*sg_));
        const double kd=std::min(k*d,30.0);
        const double n=(kd>=30.0)?0.5:0.5*(1.0+2.0*kd/std::sinh(2.0*kd));
        nerr=std::max(nerr,std::fabs(seastate_cg(sg_,k,d)-n*sg_/k)/(n*sg_/k));
    }
    check(res<1e-12,"dispersion residual < 1e-12 for d 0.1-4000 m, f 0.03-2 Hz");
    check(nerr<1e-14,"cg = n sig/k");

    check(seastate_refraction(sd,sd*sd/g,1000.0)==0.0,"refraction coefficient 0 for kd > 30");
    check(close(seastate_refraction(s10,k10,10.0),s10/std::sinh(2.0*k10*10.0),1e-14),"refraction coefficient sig/sinh(2kd)");
    check(seastate_wavenumber(1.0,0.0)==0.0,"dry: k = 0");
}

// ---------------------------------------------------------------------------------------------
// SWAN spectrum file
// ---------------------------------------------------------------------------------------------
// writes a SWAN 2D spectral file of a JONSWAP spectrum with cos^2s spreading, main direction
// theta0 (Cartesian, deg); nautical (NDIR) or Cartesian (CDIR) directions, VaDens or EnDens
static void write_spc(const std::string &file, double hs, double tp, double theta0, double s, bool naut, bool energy)
{
    const int nf=40, nd=36;
    std::vector<double> f(nf), dir(nd);
    for(int l=0;l<nf;++l) f[l]=0.03*std::pow(1.1,l);
    for(int m=0;m<nd;++m) dir[m]=5.0+10.0*m;                     // file directions (as written)
    std::vector<double> E(nf*nd);                                  // m2/Hz/degr
    double emax=0.0, sumD=0.0;
    for(int m=0;m<nd;++m)
    {
        const double thc=naut ? 270.0-dir[m] : dir[m];
        sumD+=std::pow(std::fabs(std::cos(0.5*(thc-theta0)*pi/180.0)),2*s)*10.0;
    }
    for(int l=0;l<nf;++l)
    for(int m=0;m<nd;++m)
    {
        const double thc=naut ? 270.0-dir[m] : dir[m];
        const double D=std::pow(std::fabs(std::cos(0.5*(thc-theta0)*pi/180.0)),2*s)/sumD;
        E[l*nd+m]=jonswap(2.0*pi*f[l],hs,tp,3.3)*2.0*pi*D*(energy ? 1025.0*9.81 : 1.0);
        emax=std::max(emax,E[l*nd+m]);
    }
    const double factor=1.01*emax*1.0e-4;
    std::ofstream o(file.c_str());
    o<<"SWAN   1                                Swan standard spectral file, version\n";
    o<<"$   Data produced by a test\n";
    o<<"LOCATIONS                               locations in x-y-space\n     1                                  number of locations\n      100.0000      200.0000\n";
    o<<"AFREQ                                   absolute frequencies in Hz\n"<<nf<<"                                       number of frequencies\n";
    for(int l=0;l<nf;++l) o<<f[l]<<"\n";
    o<<(naut ? "NDIR" : "CDIR")<<"                                    spectral directions in degr\n"<<nd<<"   number of directions\n";
    for(int m=0;m<nd;++m) o<<dir[m]<<"\n";
    o<<"QUANT\n     1                                  number of quantities in table\n";
    o<<(energy ? "EnDens" : "VaDens")<<"                                  densities\nm2/Hz/degr   unit\n   -0.9900E+02                          exception value\n";
    o<<"FACTOR\n"<<factor<<"\n";
    for(int l=0;l<nf;++l)
    {
        for(int m=0;m<nd;++m) o<<" "<<std::lround(E[l*nd+m]/factor);
        o<<"\n";
    }
}

static void test_swan_spc()
{
    std::cout<<"seastate_swan_spc"<<std::endl;

    seastate_grid sg(32,0.04,0.8,36);
    seastate_param sp;
    std::string err;

    struct {bool naut, energy; const char *name;} cases[3] = {{false,false,"CDIR VaDens"},{true,false,"NDIR VaDens"},{false,true,"CDIR EnDens"}};
    for(auto &c : cases)
    {
        write_spc("seastate_test.spc",1.5,9.0,40.0,8.0,c.naut,c.energy);
        seastate_swan_spc spc;
        const bool ok=spc.read("seastate_test.spc",err);
        check(ok && spc.f.size()==40 && spc.dir.size()==36,std::string(c.name)+": read 40 x 36 spectrum");
        std::vector<float> N;
        spc.to_grid(sg,N);
        sp.compute(sg,N.data());
        std::cout<<"        "<<c.name<<": Hs "<<sp.Hs<<" dir "<<sp.dir<<std::endl;
        check(close(sp.Hs,1.5,0.02),std::string(c.name)+": Hs on the spectral grid within 2 %");
        check(std::fabs(sp.dir-40.0)<0.5,std::string(c.name)+": mean direction 40 deg (Cartesian)");
        check(spc.x==100.0 && spc.y==200.0,std::string(c.name)+": location");
    }
    std::remove("seastate_test.spc");

    seastate_swan_spc bad;
    check(!bad.read("no_such_file.spc",err),"missing file reported");
}

// ---------------------------------------------------------------------------------------------
// source terms
// ---------------------------------------------------------------------------------------------
struct cell_kin
{
    std::vector<float> k, cg;
    cell_kin(const seastate_grid &sg, double d) : k(sg.nsig), cg(sg.nsig)
    {
        for(int l=0;l<sg.nsig;++l)
        {
            k[l]=float(seastate_wavenumber(sg.sig[l],d));
            cg[l]=float(seastate_cg(sg.sig[l],k[l],d));
        }
    }
};

// directionally integrated energy rate sig*S(sig) [m^2/rad], energy and action budgets
struct budget
{
    std::vector<double> s1d;
    double dE=0.0, absE=0.0, dA=0.0, absA=0.0;
    budget(const seastate_grid &sg, const std::vector<double> &S) : s1d(sg.nsig,0.0)
    {
        for(int l=0;l<sg.nsig;++l)
        for(int m=0;m<sg.ndir;++m)
        {
            const double a=S[sg.bin(l,m)]*sg.dsig[l]*sg.dtheta;
            s1d[l]+=sg.sig[l]*S[sg.bin(l,m)]*sg.dtheta;
            dE+=sg.sig[l]*a; absE+=std::fabs(sg.sig[l]*a);
            dA+=a; absA+=std::fabs(a);
        }
    }
};

static void test_source()
{
    std::cout<<"seastate_source"<<std::endl;

    seastate_grid sg(36,0.04,1.0,36);
    const int nb=sg.nbin;
    std::vector<double> P(nb), D(nb), S(nb);

    // Battjes-Janssen fraction of breaking waves: SWAN FRABRE vs. the implicit relation (1-Qb)/(-ln Qb) = (Hrms/Hm)^2
    check(seastate_source::Qb_bj(0.2,1.0)==0.0 && seastate_source::Qb_bj(1.0,1.0)==1.0 && seastate_source::Qb_bj(1.3,1.0)==1.0,"Qb = 0 for Hrms/Hm <= 0.2, 1 for Hrms/Hm >= 1");
    double qerr=0.0;
    for(double b=0.3;b<0.99;b+=0.05)
    {
        double lo=1e-300, hi=1.0-1e-15;            // bisection of the exact relation
        for(int it=0;it<200;++it)
        {
            const double q=0.5*(lo+hi);
            ((1.0-q)/(-std::log(q))<b*b ? lo : hi)=q;
        }
        qerr=std::max(qerr,std::fabs(seastate_source::Qb_bj(b,1.0)-0.5*(lo+hi)));
    }
    std::cout<<"        Qb: max. deviation from the implicit relation "<<qerr<<std::endl;
    check(qerr<0.02,"Qb within 0.02 of the implicit Battjes-Janssen relation for Hrms/Hm 0.3 - 0.95");

    // bottom friction (JONSWAP): D = Cb/g^2 (sig/sinh(kd))^2, P = 0
    {
        seastate_source_param sp; sp.friction=true;
        seastate_source src(sg,sp);
        const double d=8.0;
        cell_kin ck(sg,d);
        std::vector<float> N=make_spectrum(sg,1.0,8.0,3.3,0.0,10.0);
        src.compute(N.data(),d,ck.k.data(),ck.cg.data(),P.data(),D.data());
        double err=0.0, pmax=0.0;
        for(int l=0;l<sg.nsig;++l)
        for(int m=0;m<sg.ndir;++m)
        {
            const double ex=0.038/(g*g)*std::pow(sg.sig[l]/std::sinh(std::min(30.0,double(ck.k[l])*d)),2.0);
            err=std::max(err,std::fabs(D[sg.bin(l,m)]-ex)/ex);
            pmax=std::max(pmax,std::fabs(P[sg.bin(l,m)]));
        }
        check(err<1e-6 && pmax==0.0,"friction: D = C_b/g^2 (sig/sinh kd)^2, P = 0");
    }

    // depth-induced breaking (Battjes-Janssen)
    {
        seastate_source_param sp; sp.breaking=true;
        seastate_source src(sg,sp);
        const double d=1.5;
        cell_kin ck(sg,d);
        std::vector<float> N=make_spectrum(sg,1.5,8.0,3.3,0.0,10.0);      // Hrms = 1.06 m > Hm = 1.095 m? -> close to saturation
        src.compute(N.data(),d,ck.k.data(),ck.cg.data(),P.data(),D.data());
        const double Hrms=std::sqrt(8.0*src.Etot), Hm=0.73*d, bb=Hrms*Hrms/(Hm*Hm);
        const double Qb=seastate_source::Qb_bj(Hrms,Hm);
        const double ex=(bb<1.0 ? Qb/bb : 1.0)*src.sigm01/pi;
        std::cout<<"        breaking: Hrms/Hm "<<Hrms/Hm<<", Qb "<<src.Qb<<", D "<<D[0]<<" 1/s"<<std::endl;
        const double sbrd=ex*(1.0-Qb)/(bb-Qb);
        double nerr=0.0;
        for(int b : {0,100,nb-1}) nerr=std::max(nerr,std::fabs((P[b]-D[b]*N[b])-(-ex*N[b]))/(ex*N[b]+1e-300));
        check(close(D[0],ex+sbrd,1e-12) && close(D[nb-1],ex+sbrd,1e-12) && close(P[100],sbrd*N[100],1e-12) && src.Qb==Qb,"breaking: D = ws + sbrd, P = sbrd N (Newton linearisation, SWAN SbrD), ws = alpha/pi Qb sig_01 Hm^2/Hrms^2");
        check(nerr<1e-6,"breaking: net source P - D N = -ws N at the linearisation point");
        std::vector<float> N2=make_spectrum(sg,3.0,8.0,3.3,0.0,10.0);
        src.compute(N2.data(),d,ck.k.data(),ck.cg.data(),P.data(),D.data());
        check(src.Qb==1.0 && close(D[0],src.sigm01/pi,1e-12),"breaking: Hrms >= Hm gives Qb = 1 and D = alpha/pi sig_01");
        std::vector<float> N3=make_spectrum(sg,0.1,8.0,3.3,0.0,10.0);
        src.compute(N3.data(),20.0,ck.k.data(),ck.cg.data(),P.data(),D.data());
        check(src.Qb==0.0 && D[0]==0.0,"breaking: no dissipation for Hrms/Hm <= 0.2");
    }

    // integral parameters with the sig^-4 tail: Pierson-Moskowitz, Tm01 = 0.772 Tp
    {
        seastate_source_param sp; sp.friction=true;
        seastate_source src(sg,sp);
        cell_kin ck(sg,1000.0);
        std::vector<float> N=make_spectrum(sg,2.0,10.0,1.0,0.0,10.0);
        src.compute(N.data(),1000.0,ck.k.data(),ck.cg.data(),P.data(),D.data());
        std::cout<<"        PM Hs 2 m, Tp 10 s: Hs "<<src.Hs<<", Tm01 "<<2.0*pi/src.sigm01<<", Tm-10 "<<2.0*pi/src.sigm_10<<std::endl;
        check(close(src.Hs,2.0,0.01) && close(2.0*pi/src.sigm01,7.72,0.01) && close(2.0*pi/src.sigm_10,8.57,0.01),"moments with tail: PM Hs, Tm01 = 0.772 Tp, Tm-1,0 = 0.857 Tp within 1 %");
        const double kp=std::pow(2.0*pi/10.0,2.0)/g;
        check(src.km_wam>kp && src.km_wam<2.0*kp,"k_WAM between k_p and 2 k_p");
    }

    // van der Westhuysen et al. (2007) as SWAN GEN3 WESTH (A 732 2, Phase 9): whitecapping of SWCAP (IWCAP 7) and the
    // Yan wind input of SWIND5, recomputed here from the formulas
    {
        seastate_source_param sp; sp.komen=true; sp.westh=true; sp.wind=true; sp.U10=12.0; sp.wdir=0.3;
        seastate_source src(sg,sp);
        const double d=30.0;
        cell_kin ck(sg,d);
        std::vector<float> N=make_spectrum(sg,1.5,5.0,3.3,0.3,10.0);
        src.compute(N.data(),d,ck.k.data(),ck.cg.data(),P.data(),D.data());
        // u* of Wu (1982) as the source terms, k_WAM, sig_-10, E_tot from the source object
        const double us=src.ustar, stp=src.km_wam*std::sqrt(src.Etot)/std::sqrt(3.02e-3), ck_=3.0e-5*std::pow(stp,4.0);
        double ew=0.0, ey=0.0;
        for(int l=2;l<sg.nsig;l+=5)
        {
            double El=0.0;
            for(int m=0;m<sg.ndir;++m) El+=N[sg.bin(l,m)];
            El*=sg.sig[l]*sg.dtheta;
            const double kl=ck.k[l], B=ck.cg[l]*kl*kl*kl*El;
            const double fbr=0.5*(1.0+std::tanh(10.0*(std::sqrt(B/1.75e-3)-1.0)));
            const double pp=3.0+std::tanh(25.76*(us*kl/sg.sig[l]-0.1));
            const double fac2=std::sqrt(g*kl);
            const double wc=fbr*5.0e-5*std::pow(B/1.75e-3,0.5*pp)*std::pow(fac2/sg.sig[l],0.5*pp-1.0)*fac2+(1.0-fbr)*ck_*src.sigm_10*kl/src.km_wam;
            ew=std::max(ew,std::fabs(D[sg.bin(l,5)]-wc)/wc);
            // Yan: max(0, ((0.04 x^2 + 0.00552 x + 0.000052) cos - 0.000302) sig) N with x = u* k/sig, plus the linear growth (P without N)
            const double x=us*kl/sg.sig[l], cosd=std::cos(sg.theta[5]-0.3);
            const double yan=std::max(0.0,((0.04*x*x+0.00552*x+0.000052)*cosd-0.000302)*sg.sig[l]);
            seastate_source_param sl=sp; sl.komen=false; sl.westh=false;
            seastate_source lin(sg,sl);
            std::vector<double> P0(nb), D0(nb);
            lin.compute(N.data(),d,ck.k.data(),ck.cg.data(),P0.data(),D0.data());
            const double py=P[sg.bin(l,5)]-P0[sg.bin(l,5)];
            if(yan*N[sg.bin(l,5)]>0.0) ey=std::max(ey,std::fabs(py-yan*N[sg.bin(l,5)])/(yan*N[sg.bin(l,5)]));
        }
        std::cout<<"        Westhuysen: max. rel. difference of D "<<ew<<", of the Yan input "<<ey<<std::endl;
        check(ew<1e-9 && ey<1e-9,"van der Westhuysen whitecapping (SWAN IWCAP 7) and Yan wind input (SWIND5) as the formulas");
    }

    // whitecapping (Komen): D proportional to k^2, value
    {
        seastate_source_param sp; sp.komen=true;
        seastate_source src(sg,sp);
        cell_kin ck(sg,1000.0);
        std::vector<float> N=make_spectrum(sg,2.0,8.0,3.3,0.0,10.0);
        src.compute(N.data(),1000.0,ck.k.data(),ck.cg.data(),P.data(),D.data());
        const double stp=src.km_wam*std::sqrt(src.Etot)/std::sqrt(3.02e-3);
        const double ex=2.36e-5*std::pow(stp,4.0)*src.sigm_10*std::pow(ck.k[20]/src.km_wam,2.0);
        check(close(D[sg.bin(20,3)],ex,1e-10) && close(D[sg.bin(30,0)]/D[sg.bin(10,0)],std::pow(double(ck.k[30])/ck.k[10],2.0),1e-6),"whitecapping: D = C_ds (s/s_PM)^4 sig_-10 (k/k_WAM)^2");
    }

    // wind input
    {
        seastate_source_param sp; sp.wind=true; sp.komen=true; sp.U10=20.0; sp.wdir=0.0;
        seastate_source src(sg,sp);
        const double us=seastate_source::ustar_wu(20.0);
        check(close(us,std::sqrt((0.8+0.065*20.0)*1e-3)*20.0,1e-12) && close(seastate_source::ustar_wu(5.0),std::sqrt(1.2875e-3)*5.0,1e-12),"u* with the drag of Wu (1982)");
        cell_kin ck(sg,1000.0);
        std::vector<float> N0(nb,0.0f);
        src.compute(N0.data(),1000.0,ck.k.data(),ck.cg.data(),P.data(),D.data());
        const int l=20;
        const double sig=sg.sig[l], spm=g/(28.0*us);
        const double exA=1.5e-3/(2.0*pi*g*g*sig)*std::pow(us,4.0)*std::exp(-std::pow(std::min(2.0,spm/sig),4.0));
        check(close(P[sg.bin(l,0)],exA,1e-10) && P[sg.bin(l,18)]==0.0 && P[sg.bin(l,9)]<1e-12*P[sg.bin(l,0)],"linear growth along the wind, none across or against it");
        check(P[sg.bin(0,0)]==0.0,"no linear growth below 0.7 sig_PM");
        std::vector<float> N=make_spectrum(sg,0.5,4.0,3.3,0.0,10.0);
        std::vector<double> P1(nb), D1(nb);
        src.compute(N.data(),1000.0,ck.k.data(),ck.cg.data(),P1.data(),D1.data());
        sp.dia=false;
        const double B=std::max(0.0,0.25*1.28/1025.0*(28.0*us*ck.k[l]/sig-1.0))*sig;
        check(close(P1[sg.bin(l,0)]-P[sg.bin(l,0)],B*N[sg.bin(l,0)],1e-6),"exponential growth B N (Komen) along the wind");
        check(P1[sg.bin(l,18)]==0.0,"no exponential growth against the wind");
    }

    // quadruplets (DIA)
    {
        seastate_source_param sp; sp.komen=true; sp.dia=true;
        seastate_source src(sg,sp);
        cell_kin ck(sg,1000.0);
        std::vector<float> N=make_spectrum(sg,4.0,10.0,3.3,0.0,10.0);
        src.quadruplets(N.data(),1000.0,ck.k.data(),S.data());
        budget bu(sg,S);
        const double fp=0.1;
        auto at=[&](double f){int lb=0; for(int l=0;l<sg.nsig;++l) if(std::fabs(sg.f[l]-f)<std::fabs(sg.f[lb]-f)) lb=l; return bu.s1d[lb];};
        std::cout<<"        DIA (JONSWAP Hs 4 m, Tp 10 s, deep): energy balance "<<bu.dE/bu.absE<<", action balance "<<bu.dA/bu.absA
                 <<"; sig S(sig) at 0.85 fp "<<at(0.085)<<", 1.4 fp "<<at(0.14)<<", 2.5 fp "<<at(0.25)<<std::endl;
        check(std::fabs(bu.dE/bu.absE)<2e-3,"DIA conserves energy (relative to the transferred energy)");
        check(std::fabs(bu.dA/bu.absA)<2e-2,"DIA conserves action (within the grid interpolation)");
        check(at(0.085)>0.0 && at(0.14)<0.0 && at(0.25)>0.0,"DIA: gain below the peak, loss above it, gain in the tail");
        double asym=0.0, smax=0.0;
        for(int l=0;l<sg.nsig;++l)
        for(int m=1;m<sg.ndir;++m)
        {
            asym=std::max(asym,std::fabs(S[sg.bin(l,m)]-S[sg.bin(l,sg.ndir-m)]));
            smax=std::max(smax,std::fabs(S[sg.bin(l,m)]));
        }
        check(asym<=1e-10*smax,"DIA symmetric for a spectrum symmetric in theta");
        // diagonal derivative dS/dN against a finite difference (depth scaling and k_WAM frozen: deep water)
        std::vector<double> L(nb), Sp(nb), Sm(nb);
        src.quadruplets(N.data(),1000.0,ck.k.data(),S.data(),L.data());
        double derr=0.0;
        for(int l : {8,10,13,20})
        {
            const int b=sg.bin(l,0);
            std::vector<float> Np=N, Nm=N;
            const double h=1e-3*N[b];
            Np[b]+=float(h); Nm[b]-=float(h);
            const double hh=double(Np[b])-double(Nm[b]);
            src.quadruplets(Np.data(),1000.0,ck.k.data(),Sp.data());
            src.quadruplets(Nm.data(),1000.0,ck.k.data(),Sm.data());
            const double fd=(Sp[b]-Sm[b])/hh;
            derr=std::max(derr,std::fabs(fd-L[b])/std::fabs(fd));
        }
        std::cout<<"        DIA: diagonal derivative vs. finite difference, max. rel. deviation "<<derr<<std::endl;
        check(derr<1e-6,"DIA: diagonal derivative dS/dN (SWAN DSNL) matches a central finite difference");
        std::vector<double> S2(nb);
        src.quadruplets(N.data(),5.0,ck.k.data(),S2.data());
        check(std::fabs(S2[sg.bin(9,0)])>std::fabs(S[sg.bin(9,0)]),"DIA enhanced in shallow water (WAM depth scaling)");
        std::vector<float> Z(nb,0.0f);
        src.quadruplets(Z.data(),1000.0,ck.k.data(),S2.data());
        double z=0.0; for(double v : S2) z=std::max(z,std::fabs(v));
        check(z==0.0,"DIA: zero spectrum, zero transfer");
    }

    // triads (LTA)
    {
        seastate_source_param sp; sp.triads=true;
        seastate_source src(sg,sp);
        const double d=2.0;
        cell_kin ck(sg,d);
        std::vector<float> N=make_spectrum(sg,0.6,8.0,3.3,0.0,20.0);
        src.triads(N.data(),d,ck.k.data(),ck.cg.data(),S.data());
        budget bu(sg,S);
        int lp=0; for(int l=0;l<sg.nsig;++l) if(std::fabs(sg.f[l]-0.125)<std::fabs(sg.f[lp]-0.125)) lp=l;
        int l2=0; for(int l=0;l<sg.nsig;++l) if(std::fabs(sg.f[l]-0.25)<std::fabs(sg.f[l2]-0.25)) l2=l;
        std::cout<<"        LTA (Hs 0.6 m, Tp 8 s, d 2 m): Ursell "<<src.ursell<<", energy balance "<<bu.dE/bu.absE
                 <<"; sig S(sig) at fp "<<bu.s1d[lp]<<", 2 fp "<<bu.s1d[l2]<<std::endl;
        check(src.ursell>=0.1,"LTA case is above the Ursell limit");
        check(std::fabs(bu.dE/bu.absE)<0.03,"LTA conserves energy (within the interpolation at 2 sig)");
        check(bu.s1d[lp]<0.0 && bu.s1d[l2]>0.0,"LTA: energy from the peak to the second harmonic");
        std::vector<float> Nd=make_spectrum(sg,0.6,8.0,3.3,0.0,20.0);
        cell_kin ckd(sg,50.0);
        src.triads(Nd.data(),50.0,ckd.k.data(),ckd.cg.data(),S.data());
        double z=0.0; for(double v : S) z=std::max(z,std::fabs(v));
        check(src.ursell<0.1 && z==0.0,"LTA inactive below the Ursell limit");
    }

    // Patankar split of all terms together: P >= 0, D >= 0
    {
        seastate_source_param sp; sp.wind=true; sp.U10=15.0; sp.wdir=0.5; sp.komen=true; sp.dia=true;
        sp.breaking=true; sp.friction=true; sp.triads=true;
        seastate_source src(sg,sp);
        const double d=3.0;
        cell_kin ck(sg,d);
        std::vector<float> N=make_spectrum(sg,1.2,7.0,3.3,0.3,5.0);
        src.compute(N.data(),d,ck.k.data(),ck.cg.data(),P.data(),D.data());
        double pmin=0.0, dmin=0.0;
        for(int b=0;b<nb;++b) {pmin=std::min(pmin,P[b]); dmin=std::min(dmin,D[b]);}
        check(pmin>=0.0 && dmin>=0.0,"Patankar split: P >= 0 and D >= 0 with all source terms");
    }
}

// ---------------------------------------------------------------------------------------------
// surfbeat: single-frequency grid, Roelvink breaking, wave-group boundary generator, bound wave
// ---------------------------------------------------------------------------------------------
static void test_surfbeat()
{
    std::cout<<"seastate_surfbeat"<<std::endl;

    // single-frequency grid of the wave-group model
    {
        seastate_grid sg(0.1,24);
        check(sg.valid() && sg.nsig==1 && sg.nbin==24 && sg.dsig[0]==1.0 && close(sg.sig[0],2.0*pi*0.1,1e-15),"single-frequency grid: nsig 1, nbin = ndir, dsig 1");
        seastate_grid bad(0.0,24);
        check(!bad.valid(),"single-frequency grid: f_rep <= 0 rejected");
    }

    // Roelvink (1993) breaking: D/E = 2 alpha f_rep Qb H/h, Qb = 1 - exp(-(H/(gamma h))^n)
    {
        seastate_grid sg(0.1,24);
        seastate_source_param sp; sp.breaking=true; sp.breaking_model=2; sp.alpha=1.0; sp.gamma=0.55; sp.nroel=10.0;
        seastate_source src(sg,sp);
        const double d=2.0, H=1.0;
        cell_kin ck(sg,d);
        std::vector<float> N(sg.nbin,0.0f);
        N[0]=float(H*H/8.0/sg.dtheta/sg.sig[0]);       // E = sig N dsig dtheta = H^2/8
        std::vector<double> P(sg.nbin), D(sg.nbin);
        src.compute(N.data(),d,ck.k.data(),ck.cg.data(),P.data(),D.data());
        const double Qb=1.0-std::exp(-std::pow(H/(0.55*d),10.0));
        const double ex=2.0*0.1*Qb*H/d;
        std::cout<<"        Roelvink: Qb "<<Qb<<", D/E "<<D[5]<<" (exact "<<ex<<")"<<std::endl;
        check(close(D[0],ex,1e-5) && close(D[5],ex,1e-5) && close(src.brk_rate,ex,1e-5) && P[0]==0.0,"Roelvink breaking: D/E = 2 alpha f_rep Qb H/h in every direction, P = 0");
        std::fill(N.begin(),N.end(),0.0f);
        N[0]=float(0.01/8.0/sg.dtheta/sg.sig[0]);
        src.compute(N.data(),d,ck.k.data(),ck.cg.data(),P.data(),D.data());
        check(D[0]<1e-12,"Roelvink breaking: no dissipation for small waves");
    }

    // Herbers (1994) coefficient vs. Longuet-Higgins and Stewart (1962) for narrow-band collinear groups
    for(double kh : {0.8,1.2,2.0})
    {
        const double h=2.0, k1=kh/h;
        const double w1=std::sqrt(g*k1*std::tanh(kh));
        const double f1=w1/(2.0*pi), f2=f1*1.001;
        const double k2=seastate_wavenumber(2.0*pi*f2,h);
        double c3, th3;
        const double D=seastate_surfbeat::herbers(f1,0.0,k1,f2,0.0,k2,h,c3,th3);
        const double cg=seastate_cg(w1,k1,h), n=cg/(w1/k1);
        const double lhs=-g*(2.0*n-0.5)/(g*h-cg*cg);
        std::cout<<"        kh "<<kh<<": D "<<D<<", LHS62 "<<lhs<<", c3/cg "<<c3/cg<<std::endl;
        check(close(D,lhs,0.01) && std::fabs(th3)<1e-12,"Herbers D -> LHS62 -g(2n-1/2)/(gh-cg^2) for narrow-band collinear waves, kh "+std::to_string(kh).substr(0,3));
    }

    // bichromatic waves: envelope and bound wave time series
    {
        seastate_grid sg(0.4,24);
        std::vector<float> N0(sg.nbin,0.0f);
        seastate_surfbeat sb(sg,N0.data(),1.0,1);
        check(sb.K==0,"generator: empty spectrum gives no components");

        const double h=0.85, T=15.0, f1=6.0/15.0, f2=7.0/15.0, a1=0.09, a2=0.01;   // GLOBEX B1
        sb.trec=T; sb.df=1.0/T; sb.trep=1.0/f1; sb.frep=f1; sb.m0=0.5*(a1*a1+a2*a2);
        sb.fk={f1,f2}; sb.ak={a1,a2}; sb.thk={0.0,0.0}; sb.phk={0.3,1.1}; sb.K=2;
        sb.series({0.0},{h},true);
        const double k1=seastate_wavenumber(2.0*pi*f1,h), k2=seastate_wavenumber(2.0*pi*f2,h);
        double c3, th3;
        const double D=seastate_surfbeat::herbers(f1,0.0,k1,f2,0.0,k2,h,c3,th3);
        double eE=0.0, eZ=0.0, eQ=0.0;
        for(int n=0;n<60;++n)
        {
            const double t=n*0.25;
            double E, z, qx, qy;
            sb.at(0,t,E,z,qx,qy);
            const double ps=2.0*pi*(f2-f1)*t+1.1-0.3;
            eE=std::max(eE,std::fabs(E-0.5*(a1*a1+a2*a2+2.0*a1*a2*std::cos(ps))));
            eZ=std::max(eZ,std::fabs(z-D*a1*a2*std::cos(ps)));
            eQ=std::max(eQ,std::fabs(qx-c3*D*a1*a2*std::cos(ps))+std::fabs(qy));
        }
        std::cout<<"        bichromatic: D "<<D<<" 1/m, bound amplitude "<<D*a1*a2<<" m, max. errors E "<<eE<<", zeta "<<eZ<<", q "<<eQ<<std::endl;
        check(eE<2e-4*sb.m0 && eZ<2e-3*std::fabs(D*a1*a2) && eQ<2e-3*std::fabs(c3*D*a1*a2),"bichromatic: E = (a1^2+a2^2+2a1a2 cos)/2, zeta_b = D a1 a2 cos, q = c3 zeta_b (time interpolation)");
    }

    // random-phase generator from a directional JONSWAP spectrum
    {
        seastate_grid sg(36,0.04,1.0,36);
        std::vector<float> N=make_spectrum(sg,1.0,8.0,3.3,0.2,10.0);
        seastate_param prm; prm.compute(sg,N.data());
        seastate_surfbeat sb(sg,N.data(),1200.0,7), sb2(sg,N.data(),1200.0,7), sb3(sg,N.data(),1200.0,8);
        double var=0.0, dsum=0.0, cx=0.0, cy=0.0;
        for(int m=0;m<sb.K;++m) {var+=0.5*sb.ak[m]*sb.ak[m]; cx+=sb.ak[m]*sb.ak[m]*std::cos(sb.thk[m]); cy+=sb.ak[m]*sb.ak[m]*std::sin(sb.thk[m]);}
        for(double d : sb.Dbar) dsum+=d*sg.dtheta;
        const double mth=std::atan2(cy,cx);
        std::cout<<"        JONSWAP Hs 1 m: "<<sb.K<<" components, Hm0 "<<4.0*std::sqrt(var)<<", T_rep "<<sb.trep<<" s, mean component direction "<<mth<<" rad"<<std::endl;
        check(sb.K>100 && close(var,sb.m0,1e-12) && close(4.0*std::sqrt(sb.m0),1.0,0.01),"generator: component variance = m0 of the spectrum");
        check(close(sb.trep,prm.Tm10,0.05) && close(dsum,1.0,1e-12) && std::fabs(mth-0.2)<0.05,"generator: T_rep = Tm-1,0, mean directional distribution normalised, directions around the main direction");
        sb.series({0.0,50.0},{8.0,8.0},true);
        sb2.series({0.0,50.0},{8.0,8.0},true);
        sb3.series({0.0,50.0},{8.0,8.0},true);
        const int nt=int(std::lround(sb.trec/sb.dtbc));
        double Em=0.0, Zm=0.0, Z2=0.0, Emax=0.0, dif2=0.0, dif3=0.0, dify=0.0;
        for(int n=0;n<nt;++n)
        {
            double E,z,qx,qy, E2,z2,qx2,qy2, E3,z3,qx3,qy3, Ey,zy,qxy,qyy;
            sb.at(0,n*sb.dtbc,E,z,qx,qy);
            sb2.at(0,n*sb.dtbc,E2,z2,qx2,qy2);
            sb3.at(0,n*sb.dtbc,E3,z3,qx3,qy3);
            sb.at(1,n*sb.dtbc,Ey,zy,qxy,qyy);
            Em+=E/nt; Zm+=z/nt; Z2+=z*z/nt; Emax=std::max(Emax,E);
            dif2=std::max(dif2,std::fabs(E-E2)+std::fabs(z-z2)+std::fabs(qx-qx2));
            dif3=std::max(dif3,std::fabs(E-E3));
            dify=std::max(dify,std::fabs(E-Ey));
        }
        std::cout<<"        envelope mean "<<Em<<" (m0 "<<sb.m0<<"), max/mean "<<Emax/Em<<", bound wave mean "<<Zm<<" m, rms "<<std::sqrt(Z2)<<" m"<<std::endl;
        check(close(Em,sb.m0,1e-10) && std::fabs(Zm)<1e-12 && std::sqrt(Z2)>1e-4 && std::sqrt(Z2)<0.05,"envelope mean = m0 over T_rec, bound wave zero-mean and small");
        check(dif2==0.0 && dif3>0.1*sb.m0 && dify>0.1*sb.m0,"generator: identical for the same seed, different for another seed and along the boundary");
    }
}

// ---------------------------------------------------------------------------------------------
// external forcing: date-times, series of SWAN spectra, wind field
// ---------------------------------------------------------------------------------------------
static void test_forcing()
{
    std::cout<<"seastate_forcing"<<std::endl;

    check(seastate_datetime(19700101.000000)==0.0 && seastate_datetime(20000301.000000)==951868800.0
          && seastate_datetime(20261005.123045)==1791203445.0 && seastate_datetime(19991231.235959)==946684799.0,
          "date-time YYYYMMDD.HHMMSS -> seconds since 1970 (leap years, end of year)");

    // series of spectra: 2 locations, 3 times (1 h apart), nautical directions, energy density
    seastate_grid sg(20,0.05,0.5,24);
    const char *spc = "/tmp/seastate_test_series.spc";
    {
        std::ofstream o(spc);
        o<<"SWAN   1\n$ test\nTIME\n     1\nLOCATIONS\n     2\n   0.0  0.0\n   0.0  1000.0\nAFREQ\n     3\n 0.08\n 0.10\n 0.12\nNDIR\n     4\n 0.0\n 90.0\n 180.0\n 270.0\n";
        o<<"QUANT\n     1\nEnDens\nJ/m2/Hz/degr\n -0.9900E+02\n";
        const char *t[3] = {"20261005.000000","20261005.010000","20261005.020000"};
        for(int k=0;k<3;++k)
        {
            o<<t[k]<<"\n";
            for(int l=0;l<2;++l)
            {
                if(k==1 && l==1) {o<<"ZERO\n"; continue;}
                o<<"FACTOR\n  1.0\n";
                for(int fq=0;fq<3;++fq)
                    o<<"  "<<(k+1)*(l+1)*10.0*1025.0*9.81<<"  0  0  "<<(fq==1 ? -99 : 0)<<"\n";   // waves from north (nautical 0): propagate to -y
            }
        }
    }
    seastate_spc_series ser;
    std::string err;
    const bool ok = ser.open(spc,sg,err);
    if(!ok) std::cout<<"        "<<err<<std::endl;
    check(ok && ser.nloc==2 && ser.ys[1]==1000.0 && !ser.stationary() && ser.time(0)==seastate_datetime(20261005.0) && ser.time(1)-ser.time(0)==3600.0,
          "spectra series: header, locations, first two records");

    // reference: the same record through the stationary reader's interpolation
    seastate_swan_spc ref;
    ref.f={0.08,0.10,0.12}; ref.dir={0.0,90.0,180.0,270.0}; ref.E.assign(12,0.0);
    for(int fq=0;fq<3;++fq) ref.E[fq*4+3]=10.0;          // nautical 0 -> Cartesian 270 (to -y), sorted last
    std::vector<float> Nr; ref.to_grid(sg,Nr);
    double e0=0.0, e1=0.0;
    for(int b=0;b<sg.nbin;++b)
    {
        e0=std::max(e0,double(std::fabs(ser.N(0,0)[b]-Nr[b])+std::fabs(ser.N(0,1)[b]-2.0f*Nr[b])));
        e1=std::max(e1,double(std::fabs(ser.N(1,1)[b])));
    }
    check(e0<2e-5*Nr[sg.bin(5,18)] && e1==0.0 && Nr[sg.bin(5,18)]>0.0f,"spectra series: EnDens / rho g, nautical -> Cartesian, exception value, ZERO block");
    const double t2 = ser.time(0)+5400.0;
    const bool adv = ser.advance(t2,err);
    check(adv && ser.time(0)==seastate_datetime(20261005.01) && ser.time(1)==seastate_datetime(20261005.02) && std::fabs(ser.weight(t2)-0.5)<1e-12
          && std::fabs(ser.N(1,0)[sg.bin(5,18)]-3.0f*Nr[sg.bin(5,18)])<2e-5*Nr[sg.bin(5,18)],"spectra series: advance and linear time weight");
    ser.advance(ser.time(1)+1e5,err);
    check(ser.weight(ser.time(1)+1e5)==1.0 && ser.records==3,"spectra series: clamped after the last record");

    // wind field: linear in x and y, two times
    const char *wf = "/tmp/seastate_test_wind.dat";
    {
        std::ofstream o(wf);
        o<<"wind test\n3 2\n100.0 -50.0 200.0 300.0\n";
        for(int k=0;k<2;++k)
        {
            o<<(k==0 ? "20261005.000000" : "20261005.060000")<<"\n";
            for(int c=0;c<2;++c)
            for(int jj=0;jj<2;++jj)
            {
                for(int ii=0;ii<3;++ii)
                {
                    const double x=100.0+200.0*ii, y=-50.0+300.0*jj;
                    o<<(c==0 ? 5.0+0.01*x+0.002*y+k*4.0 : -3.0+0.004*x-0.01*y)<<" ";
                }
                o<<"$ row\n";
            }
        }
    }
    seastate_wind_series ws;
    const bool okw = ws.open(wf,err);
    if(!okw) std::cout<<"        "<<err<<std::endl;
    double werr=0.0;
    const double tw = seastate_datetime(20261005.03);
    ws.advance(tw,err);
    for(double x : {100.0,230.0,499.0}) for(double y : {-50.0,10.0,250.0})
    {
        double u,v; ws.at(tw,x,y,u,v);
        werr=std::max(werr,std::fabs(u-(5.0+0.01*x+0.002*y+2.0))+std::fabs(v-(-3.0+0.004*x-0.01*y)));
    }
    double uc,vc; ws.at(tw,-1000.0,1000.0,uc,vc);
    std::cout<<"        wind: max. interpolation error "<<werr<<std::endl;
    check(okw && ws.nx==3 && ws.ny==2 && werr<1e-12,"wind field: bilinear in space and linear in time are exact for a linear field");
    check(std::fabs(uc-(5.0+0.01*100.0+0.002*250.0+2.0))<1e-12,"wind field: clamped to the grid outside");
}

// Phase 5b: bathymetry raster (A 790 1) and the maximum energy of SWAN (A 737 1)
static void test_bathy()
{
    std::cout<<"seastate_bathy, maximum energy"<<std::endl;

    // raster of a linear bed z = -20 + 0.01 x - 0.02 y on 6 x 4 nodes, 50 m x 100 m spacing, origin (1000, 2000)
    const std::string f = "seastate_test_bathy.dat";
    {
        std::ofstream o(f);
        o<<"$ test raster\n6 4   $ nodes\n1000 2000 50 100\n";
        for(int j=0; j<4; ++j)
        {
            for(int i=0; i<6; ++i)
            o<<(-20.0 + 0.01*(1000.0+50.0*i) - 0.02*(2000.0+100.0*j))<<" ";
            o<<"\n";
        }
    }
    auto zl = [](double x, double y) {return -20.0 + 0.01*x - 0.02*y;};

    seastate_bathy b;
    std::string err;
    const bool ok = b.read(f,err);
    if(!ok) std::cout<<"        "<<err<<std::endl;
    check(ok && b.nx==6 && b.ny==4 && b.dx==50.0 && b.dy==100.0,"raster header and values read ($ comments)");

    double e1=0.0;
    for(double x : {1000.0,1033.0,1210.0,1250.0}) for(double y : {2000.0,2077.0,2300.0})
    e1 = std::max(e1,std::fabs(b.at(x,y)-zl(x,y)));
    check(e1<1e-12,"bilinear interpolation exact for a linear bed");
    check(std::fabs(b.at(0.0,0.0)-zl(1000.0,2000.0))<1e-12 && std::fabs(b.at(5000.0,9000.0)-zl(1250.0,2300.0))<1e-12,"clamped outside the raster");

    // cell [1000,1100) x [2000,2200): nodes x 1000, 1050 and y 2000, 2100 -> mean at (1025, 2050)
    check(std::fabs(b.cell(1000.0,1100.0,2000.0,2200.0)-zl(1025.0,2050.0))<1e-12,"cell value: mean of the nodes inside the cell");
    // cell finer than the raster: no node inside -> bilinear at the centre
    check(std::fabs(b.cell(1010.0,1030.0,2010.0,2030.0)-zl(1020.0,2020.0))<1e-12,"cell without nodes: bilinear at the centre");

    // majority rule with zdry: a cell with 3 dry nodes (z = +2) and 1 wet node (z = -10) is dry (bed +2),
    // with 1 dry and 3 wet nodes wet with the mean of the wet nodes (-10), without zdry the plain mean
    {
        std::ofstream o("seastate_test_bathy2.dat");
        o<<"islet\n2 2\n0 0 10 10\n2 2\n2 -10\n";
        std::ofstream o2("seastate_test_bathy3.dat");
        o2<<"islet\n2 2\n0 0 10 10\n-10 -10\n2 -10\n";
    }
    seastate_bathy b2, b3;
    b2.read("seastate_test_bathy2.dat",err);
    b3.read("seastate_test_bathy3.dat",err);
    check(std::fabs(b2.cell(-5.0,15.0,-5.0,15.0,-0.05)-2.0)<1e-12 && std::fabs(b3.cell(-5.0,15.0,-5.0,15.0,-0.05)+10.0)<1e-12,"cell: majority of wet/dry nodes, mean bed of the majority");
    check(std::fabs(b2.cell(-5.0,15.0,-5.0,15.0)+1.0)<1e-12,"cell without zdry: mean of all nodes");

    {
        std::ofstream o("seastate_test_bathy_bad.dat");
        o<<"title\n3 3\n0 0 1 1\n1 2 3\n4 5\n";
    }
    seastate_bathy bb;
    check(!bb.read("seastate_test_bathy_bad.dat",err),"incomplete raster rejected");
    check(!bb.read("does_not_exist.dat",err),"missing file rejected");

    // maximum energy: E_tot (tail included) <= (gamma d)^2/4 with Battjes-Janssen breaking
    seastate_grid sg(30,0.05,1.0,24);
    std::vector<float> N = make_spectrum(sg,3.0,8.0,3.3,0.0,10.0);
    std::vector<float> Nc = N;

    seastate_source_param sp; sp.breaking=true; sp.emax=true; sp.gamma=0.73;
    seastate_source src(sg,sp);
    const double d = 2.0;
    const bool capped = src.cap(Nc.data(),d);

    cell_kin ck(sg,d);
    seastate_param pm; pm.compute(sg,Nc.data());
    std::vector<double> P(sg.nbin), D(sg.nbin);
    src.compute(Nc.data(),d,ck.k.data(),ck.cg.data(),P.data(),D.data());
    const double emax = 0.25*(0.73*d)*(0.73*d);
    std::cout<<"        E_tot after the cap "<<src.Etot<<", (gamma d)^2/4 = "<<emax<<std::endl;
    check(capped && std::fabs(src.Etot-emax)<1e-6*emax,"cap: E_tot (with tail) = (gamma d)^2/4");
    double rat=-1.0; bool shape=true;
    for(int bn=0; bn<sg.nbin; ++bn)
    if(N[bn]>0.0f)
    {
        const double r = double(Nc[bn])/double(N[bn]);
        if(rat<0.0) rat=r; else if(std::fabs(r-rat)>1e-5*rat) shape=false;
    }
    check(shape,"cap: the spectral shape is kept (one factor)");

    std::vector<float> N2 = N;
    check(!src.cap(N2.data(),50.0) && N2==N,"no cap in deep water");
    seastate_source_param sp0; sp0.breaking=true; sp0.emax=false;
    seastate_source src0(sg,sp0);
    std::vector<float> N3 = N;
    check(!src0.cap(N3.data(),d) && N3==N,"no cap with A 737 0");
}

// ---------------------------------------------------------------------------------------------
// Phase 6b: vegetation (Dalrymple / Suzuki et al. 2011, Jacobsen et al. 2019), structures
// ---------------------------------------------------------------------------------------------
static void test_vegetation()
{
    std::cout<<"seastate_source: vegetation"<<std::endl;

    seastate_grid sg(36,0.04,1.0,36);
    const int nb=sg.nbin;
    std::vector<double> P(nb), D(nb);
    const double d=2.0, ah=1.0, bv=0.01, nv=1000.0, cd=1.0;

    // IVEG 1: the dissipation rate D E_tot equals eps/(rho g) of Mendez and Losada (2004) for Hrms = sqrt(8 E_tot),
    // eps = 1/(2 sqrt(pi)) rho Cd bv Nv (g k/(2 sig))^3 (sinh^3 k ah + 3 sinh k ah)/(3 k cosh^3 k d) Hrms^3, with the
    // k and sig of the source term (k_WAM, sig_01)
    {
        seastate_source_param sp; sp.vegetation=1; sp.vh=ah; sp.vd=bv; sp.vn=nv; sp.vcd=cd;
        seastate_source src(sg,sp);
        cell_kin ck(sg,d);
        std::vector<float> N=make_spectrum(sg,0.4,4.0,3.3,0.0,50.0);
        src.compute(N.data(),d,ck.k.data(),ck.cg.data(),P.data(),D.data());
        const double k=src.km_wam, s=src.sigm01, H=std::sqrt(8.0*src.Etot);
        const double sh=std::sinh(k*ah), ch=std::cosh(k*d);
        const double eps=1.0/(2.0*std::sqrt(pi))*cd*bv*nv*std::pow(g*k/(2.0*s),3.0)*(sh*sh*sh+3.0*sh)/(3.0*k*ch*ch*ch)*H*H*H/g;
        bool uni=true;
        for(int b=1;b<nb;++b) uni = uni && D[b]==D[0];
        std::cout<<"        D E_tot "<<D[0]*src.Etot<<", eps/(rho g) "<<eps<<std::endl;
        check(uni,"IVEG 1: the same dissipation rate in all bins");
        check(close(D[0]*src.Etot,eps,1e-10),"IVEG 1: D E_tot = eps/(rho g) of Mendez and Losada (2004)");

        // emergent vegetation: the height is limited to the depth
        seastate_source_param sp2=sp; sp2.vh=5.0;
        seastate_source src2(sg,sp2), src3(sg,[&]{seastate_source_param q=sp; q.vh=d; return q;}());
        std::vector<double> D2(nb), D3(nb);
        src2.compute(N.data(),d,ck.k.data(),ck.cg.data(),P.data(),D2.data());
        src3.compute(N.data(),d,ck.k.data(),ck.cg.data(),P.data(),D3.data());
        check(D2[0]==D3[0],"IVEG 1: emergent vegetation (height > depth) acts over the water depth");
    }

    // IVEG 2 (per frequency) for a spectrum in one frequency: int_0^ah S_u sqrt(mu) dz with S_u = (sig cosh kz/sinh kd)^2 E
    // gives sqrt(2/pi)/g Cd bv Nv sig^3 sqrt(E) (sinh^3 k ah + 3 sinh k ah)/(3 k sinh^3 kd), with the dispersion relation
    // equal to IVEG 1 at k and sig of that frequency (Simpson with 20 intervals)
    {
        cell_kin ck(sg,d);
        std::vector<float> N(nb,0.0f);
        const int l0=12;
        for(int m=0;m<sg.ndir;++m) N[sg.bin(l0,m)] = (m==0) ? 0.05f : 0.0f;
        seastate_source_param sp; sp.vegetation=2; sp.vh=ah; sp.vd=bv; sp.vn=nv; sp.vcd=cd;
        seastate_source src(sg,sp);
        src.compute(N.data(),d,ck.k.data(),ck.cg.data(),P.data(),D.data());
        const double k=ck.k[l0], s=sg.sig[l0];
        const double E=0.05*s*sg.dsig[l0]*sg.dtheta;
        const double sh=std::sinh(k*ah), sd=std::sinh(k*d);
        const double ref=std::sqrt(2.0/pi)/g*cd*bv*nv*s*s*s*std::sqrt(E)*(sh*sh*sh+3.0*sh)/(3.0*k*sd*sd*sd);
        std::cout<<"        IVEG 2 one frequency "<<D[sg.bin(l0,0)]<<", closed form "<<ref<<std::endl;
        check(close(D[sg.bin(l0,0)],ref,2e-4),"IVEG 2: one frequency, Simpson integral within 2e-4 of the closed form");
        const double ch=std::cosh(k*d);
        const double v1=std::sqrt(2.0/pi)*g*g*std::pow(k/s,3.0)*std::sqrt(E)/(3.0*k*ch*ch*ch)*cd*bv*nv*(sh*sh*sh+3.0*sh);
        check(close(ref,v1,1e-5),"IVEG 2 = IVEG 1 for one frequency (dispersion relation, k in single precision)");
    }
}

static void test_structure()
{
    std::cout<<"seastate_structure"<<std::endl;

    // d'Angremond et al. (1996), as SWAN DAM DANGREMOND
    {
        const double Tp=8.0, Hs=1.0, sl=26.565;
        const double xi=std::tan(sl*pi/180.0)/std::sqrt(Hs/(1.5613*Tp*Tp));
        const double a=-0.4*0.5 + 0.64*std::pow(4.0,-0.31)*(1.0-std::exp(-0.5*xi));
        check(close(seastate_dangremond(0.5,Hs,Tp,sl,4.0),a,1e-14),"d'Angremond: B/Hs < 8");
        const double b=-0.35*(-0.2) + 0.51*std::pow(15.0,-0.65)*(1.0-std::exp(-0.41*xi));
        check(close(seastate_dangremond(-0.2,Hs,Tp,sl,15.0),std::max(std::min(b,0.93-0.006*15.0),0.05),1e-14),"d'Angremond: B/Hs > 12");
        check(seastate_dangremond(3.0,Hs,Tp,sl,4.0)==0.075 && seastate_dangremond(-3.0,Hs,Tp,sl,4.0)==0.9,"d'Angremond: limits 0.075 and 0.9");
        const double k8=seastate_dangremond(0.3,Hs,Tp,sl,8.0), k12=seastate_dangremond(0.3,Hs,Tp,sl,12.0), k10=seastate_dangremond(0.3,Hs,Tp,sl,10.0);
        check(std::fabs(k10-0.5*(k8+k12))<0.02,"d'Angremond: 8 < B/Hs < 12 between the two branches");
    }

    // porous slab: the grid and spectrum of porous_ref.py (validation 37), Hs 1 m, depth 10 m, B 8 m, n 0.4, D50 0.5 m
    seastate_grid sg(20,0.06,0.3,36);
    std::vector<double> E(sg.nsig);
    double m0=0.0;
    for(int l=0;l<sg.nsig;++l) {E[l]=std::exp(-std::pow((sg.sig[l]-0.8)/0.12,2.0)); m0+=E[l]*sg.dsig[l];}
    for(int l=0;l<sg.nsig;++l) E[l]*=1.0/16.0/m0;
    std::vector<float> kt2(sg.nsig), kr2(sg.nsig);
    const double q=seastate_porous(sg,E.data(),10.0,8.0,0.4,0.5,kt2.data(),kr2.data());
    const int ls[4]={2,6,10,15};
    const double rt[4]={0.20019382200794292,0.16897093535502558,0.11834474190989125,0.04646329078514286};
    const double rr[4]={0.3244055517575106,0.37460072513004794,0.44970520824901283,0.44808032888830795};
    double dev=0.0;
    for(int k=0;k<4;++k) dev=std::max(dev,std::max(std::fabs(kt2[ls[k]]-rt[k])/rt[k],std::fabs(kr2[ls[k]]-rr[k])/rr[k]));
    std::cout<<"        q_rms "<<q<<" (reference 0.0783368), max. relative deviation of Kt^2, Kr^2 "<<dev<<std::endl;
    check(std::fabs(q-0.0783368)<2e-3*0.0783368,"porous: rms discharge velocity of the independent reference (porous_ref.py) within 0.2 %");
    check(dev<5e-3,"porous: Kt^2, Kr^2 at four frequencies within 0.5 % of the independent reference");

    // stronger resistance (D50 0.2 m, depth 20 m): the progressive mode by continuation from the root without resistance
    {
    std::vector<float> a2(sg.nsig), b2(sg.nsig);
    const double qq=seastate_porous(sg,E.data(),20.0,10.0,0.4,0.2,a2.data(),b2.data());
    const int lt[3]={3,8,13};
    const double at[3]={0.16103970538356113,0.06660913899566258,0.020084685873719924};
    const double bt[3]={0.35228954274549995,0.5135362697553263,0.5347080822396175};
    double dv=0.0;
    for(int k=0;k<3;++k) dv=std::max(dv,std::max(std::fabs(a2[lt[k]]-at[k])/at[k],std::fabs(b2[lt[k]]-bt[k])/bt[k]));
    std::cout<<"        D50 0.2 m: q_rms "<<qq<<" (reference 0.0367835), max. relative deviation "<<dv<<std::endl;
    check(std::fabs(qq-0.0367835)<2e-3*0.0367835 && dv<5e-3,"porous, stronger resistance: within 0.5 % of the independent reference");
    }

    // lossless limit (very coarse stones: no resistance): Kt^2 + Kr^2 = 1
    seastate_porous(sg,E.data(),10.0,8.0,0.4,1.0e6,kt2.data(),kr2.data());
    double loss=0.0;
    for(int l=0;l<sg.nsig;++l) loss=std::max(loss,std::fabs(1.0-double(kt2[l])-double(kr2[l])));
    check(loss<1e-5,"porous: Kt^2 + Kr^2 = 1 without resistance (energy conservation of the matching)");

    // long waves without resistance: Madsen (1974), gamma = n/sqrt(s), k_s = k sqrt(s)
    seastate_grid sl(3,0.005,0.01,36);
    std::vector<double> El(3,1.0e-4);
    seastate_porous(sl,El.data(),2.0,20.0,0.4,1.0e6,kt2.data(),kr2.data());
    const double s=1.0+0.34*0.6/0.4, kk=sl.sig[1]/std::sqrt(g*2.0), ga=0.4/std::sqrt(s), ph=kk*std::sqrt(s)*20.0;
    const double t2=16.0*ga*ga/(std::pow((1+ga)*(1+ga)-(1-ga)*(1-ga),2.0)*std::pow(std::cos(ph),2.0)+std::pow((1+ga)*(1+ga)+(1-ga)*(1-ga),2.0)*std::pow(std::sin(ph),2.0));
    std::cout<<"        long wave Kt^2 "<<kt2[1]<<", Madsen (1974) "<<t2<<std::endl;
    check(std::fabs(kt2[1]-t2)<2e-3,"porous: long waves within 2e-3 of Madsen (1974)");
}

// ---------------------------------------------------------------------------------------------
// Phase 8: DIA of the window of a sweep (only the interaction terms it needs, periodic directions)
// ---------------------------------------------------------------------------------------------
static void test_dia_window()
{
    std::cout<<"DIA of a sweep window (Phase 8)"<<std::endl;

    for(int ndir : {36,72})
    {
    seastate_grid sg(31,0.04,1.0,ndir);
    const int nb=sg.nbin;
    seastate_source_param sp; sp.dia=true;
    seastate_source src(sg,sp);
    const double d=15.0;
    cell_kin ck(sg,d);
    std::vector<float> N=make_spectrum(sg,1.5,6.0,3.3,0.6,4.0);
    for(int b=0;b<nb;b+=7) N[b]=0.0f;                               // some empty bins

    // the full DIA, split as in compute: P = S+ - L- N, D = -S-/N - L-
    std::vector<double> S(nb), L(nb), P(nb), D(nb);
    src.quadruplets(N.data(),d,ck.k.data(),S.data(),L.data());
    double err=0.0;
    for(int q=0;q<4;++q)
    for(int sub=0;sub<2;++sub)
    {
        const int a=q*ndir/4+sub, b=(q+1)*ndir/4-1-2*sub, la=3*sub, lb=sg.nsig-1-4*sub;
        std::fill(P.begin(),P.end(),0.0); std::fill(D.begin(),D.end(),0.0);
        src.compute(N.data(),d,ck.k.data(),ck.cg.data(),P.data(),D.data(),a,b,la,lb);
        for(int l=la;l<=lb;++l)
        for(int m=a;m<=b;++m)
        {
            const int x=sg.bin(l,m);
            const double s=S[x], v=N[x], ln=std::min(L[x],0.0);
            const double pe=std::max(s,0.0)-ln*v, de=((s<0.0 && v>0.0) ? -s/v : 0.0)-ln;
            err=std::max(err,std::fabs(P[x]-pe)/(std::fabs(pe)+std::fabs(de)*v+1e-300));
            err=std::max(err,std::fabs(D[x]-de)/(std::fabs(de)+1e-300));
        }
    }
    std::cout<<"        "<<ndir<<" directions: windows (quadrants and sub-windows) vs the full DIA, max. rel. difference "<<err<<std::endl;
    check(err<1e-10,"DIA of a window equals the full DIA in the window (to round-off), "+std::to_string(ndir)+" directions");

    }
}

int main()
{
    test_grid();
    test_store();
    test_param();
    test_dispersion();
    test_swan_spc();
    test_source();
    test_surfbeat();
    test_forcing();
    test_bathy();
    test_vegetation();
    test_structure();
    test_dia_window();
    test_memory();

    std::cout<<std::endl<<(nfail ? "FAILED: " : "all passed")<<(nfail ? std::to_string(nfail) : std::string())<<std::endl;
    return nfail ? 1 : 0;
}

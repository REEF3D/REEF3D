// Architect: Hans Bihs
// Standalone verification of the REEF3D::SEASTATE kernels: spectral grid, block-sparse action
// storage, integrated wave parameters, dispersion relation, SWAN spectrum files, source terms. No MPI, no REEF3D binary.
// Build:  g++ -O2 -std=c++20 -I../../src seastate_test.cpp ../../src/seastate_grid.cpp ../../src/seastate_store.cpp ../../src/seastate_param.cpp ../../src/seastate_dispersion.cpp ../../src/seastate_swan_spc.cpp ../../src/seastate_source.cpp -o seastate_test
// Run:    ./seastate_test
#include"seastate_grid.h"
#include"seastate_store.h"
#include"seastate_param.h"
#include"seastate_dispersion.h"
#include"seastate_swan_spc.h"
#include"seastate_source.h"
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

int main()
{
    test_grid();
    test_store();
    test_param();
    test_dispersion();
    test_swan_spc();
    test_source();
    test_memory();

    std::cout<<std::endl<<(nfail ? "FAILED: " : "all passed")<<(nfail ? std::to_string(nfail) : std::string())<<std::endl;
    return nfail ? 1 : 0;
}

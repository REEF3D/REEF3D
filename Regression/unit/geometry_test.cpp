// Standalone check of the horizontal geometry layer of the grid (grid_geometry.cpp): no MPI, no
// REEF3D solver.
// Build:  g++ -O2 -std=c++17 -I../../src geometry_test.cpp -o geometry_test
// Run:    ./geometry_test
// Checks: Cartesian grids (uniform, stretched, large UTM coordinates): 2D nodes equal XN/YN bit
// for bit, cell areas and face vectors are the exact products of the 1D spacings, the quadrilateral
// formulas agree, closure exactly zero; curvilinear grids (annulus sector, sheared and wavy grids):
// face vectors close each cell (geometric conservation), the cell areas sum to the polygon area of
// the boundary, centroids of parallelograms; bilinear subdivision (as the AMR patches) keeps the
// Cartesian nodes bit for bit and tiles each curvilinear cell exactly (children's areas and outer
// face vectors sum to the parent's).
// Architect: Hans Bihs

// grid_geometry.cpp needs only the coordinate arrays of the grid and marge of increment
#include<cmath>
#define INCREMENT_H_
class increment
{
public:
    virtual ~increment() {}
    static int marge;
};
int increment::marge = 5;
#include"../../src/grid_geometry.cpp"

#include<cstdio>
#include<iostream>
#include<string>
#include<vector>
#include<functional>

static int nfail = 0, npass = 0;
static void check(bool ok, const std::string &what)
{
    if(ok) ++npass;
    else
    {
        ++nfail;
        if(nfail<=20)
        std::cout<<"  FAIL  "<<what<<std::endl;
        if(nfail==20)
        std::cout<<"  (further failures not printed)"<<std::endl;
    }
}

// a grid with 1D nodes x(i), y(j); the 1D arrays as read_grid and gridspacing set them
struct testgrid : public grid
{
    std::vector<double> xn,yn,xp,yp,dxn,dyn,dxp,dyp;

    testgrid(int nx, int ny, std::function<double(int)> fx, std::function<double(int)> fy)
    {
        const int m = marge;
        knox = nx; knoy = ny; knoz = 1;
        imin = jmin = -margin;
        imax = knox+2*margin;
        jmax = knoy+2*margin;

        xn.resize(knox+1+4*m); xp=dxn=dxp=xn;
        yn.resize(knoy+1+4*m); yp=dyn=dyp=yn;
        for(int i=-m; i<knox+1+m; ++i) xn[i+m] = fx(i);
        for(int j=-m; j<knoy+1+m; ++j) yn[j+m] = fy(j);
        for(int i=-m; i<knox+m; ++i) { xp[i+m] = 0.5*(xn[i+m]+xn[i+m+1]); dxn[i+m] = xn[i+m+1]-xn[i+m]; }
        for(int j=-m; j<knoy+m; ++j) { yp[j+m] = 0.5*(yn[j+m]+yn[j+m+1]); dyn[j+m] = yn[j+m+1]-yn[j+m]; }

        XN = xn.data(); YN = yn.data(); XP = xp.data(); YP = yp.data();
        DXN = dxn.data(); DYN = dyn.data(); DXP = dxp.data(); DYP = dyp.data();
    }

    ~testgrid() { geometry_free(); }

    int sl(int i, int j) const { return (i-imin)*jmax + (j-jmin); }
};

static void cartesian_case(const char *name, int nx, int ny, std::function<double(int)> fx, std::function<double(int)> fy)
{
    testgrid g(nx,ny,fx,fy);
    g.geometry_alloc();
    g.geometry_cartesian_nodes();
    g.geometry_metrics();

    double err,gcl;
    const int mm = g.geometry_check(err,gcl);
    check(mm==0, std::string(name)+": 2D nodes equal XN/YN");
    check(err<1.0e-9, std::string(name)+": metrics agree with the quadrilateral formulas ("+std::to_string(err)+")");
    check(gcl==0.0, std::string(name)+": closure exactly zero");

    bool exact = true;
    for(int i=g.imin; i<g.imin+g.imax; ++i)
    for(int j=g.jmin; j<g.jmin+g.jmax; ++j)
    {
        const int q = g.sl(i,j), a = i+increment::marge, b = j+increment::marge;
        exact = exact && g.AREA2D[q]==g.DXN[a]*g.DYN[b] && g.SX1[q]==g.DYN[b] && g.SY1[q]==0.0
                      && g.SX2[q]==0.0 && g.SY2[q]==g.DXN[a] && g.XC2D[q]==g.XP[a] && g.YC2D[q]==g.YP[b];
    }
    check(exact, std::string(name)+": area, face vectors, centres are the 1D products bit for bit");

    // the general (curvilinear) branch on the same nodes: equal to round-off
    testgrid h(nx,ny,fx,fy);
    h.geo_curv = 1;
    h.geometry_alloc();
    h.geometry_cartesian_nodes();
    h.geometry_metrics();
    double d = 0.0;
    for(int q=0; q<g.geo_nslice; ++q)
    {
        const double L = std::sqrt(g.AREA2D[q]);
        d = std::max(d,std::fabs(g.AREA2D[q]-h.AREA2D[q])/g.AREA2D[q]);
        d = std::max(d,std::max(std::fabs(g.SX1[q]-h.SX1[q]),std::fabs(g.SY2[q]-h.SY2[q]))/L);
        d = std::max(d,std::max(std::fabs(g.XC2D[q]-h.XC2D[q]),std::fabs(g.YC2D[q]-h.YC2D[q]))/L);
    }
    check(d<1.0e-9, std::string(name)+": curvilinear branch agrees on a Cartesian grid ("+std::to_string(d)+")");
}

// a curvilinear grid: 2D nodes from a mapping (i,j) -> (x,y), geo_curv=1
struct curvgrid : public testgrid
{
    curvgrid(int nx, int ny, std::function<void(double,double,double&,double&)> map)
    : testgrid(nx,ny,[](int i){return double(i);},[](int j){return double(j);})
    {
        geo_curv = 1;
        geometry_alloc();
        for(int i=imin; i<=imin+imax; ++i)
        for(int j=jmin; j<=jmin+jmax; ++j)
        map(double(i),double(j),XN2D[nij(i,j)],YN2D[nij(i,j)]);
        geometry_metrics();
    }
};

static void curvilinear_case(const char *name, int nx, int ny, std::function<void(double,double,double&,double&)> map)
{
    curvgrid g(nx,ny,map);
    double err,gcl;
    g.geometry_check(err,gcl);
    check(err<1.0e-12, std::string(name)+": stored metrics = quadrilateral formulas");
    check(gcl<1.0e-12, std::string(name)+": face vectors close every cell ("+std::to_string(gcl)+")");

    bool positive = true;
    for(int q=0; q<g.geo_nslice; ++q)
    positive = positive && g.AREA2D[q]>0.0;
    check(positive, std::string(name)+": positive cell areas (i,j right-handed)");

    // sum of the interior cell areas = shoelace area of the boundary polygon
    double A = 0.0;
    for(int i=0; i<g.knox; ++i)
    for(int j=0; j<g.knoy; ++j)
    A += g.AREA2D[g.sl(i,j)];

    std::vector<double> px,py;
    for(int i=0; i<=g.knox; ++i) { px.push_back(g.XN2D[g.nij(i,0)]); py.push_back(g.YN2D[g.nij(i,0)]); }
    for(int j=1; j<=g.knoy; ++j) { px.push_back(g.XN2D[g.nij(g.knox,j)]); py.push_back(g.YN2D[g.nij(g.knox,j)]); }
    for(int i=g.knox-1; i>=0; --i) { px.push_back(g.XN2D[g.nij(i,g.knoy)]); py.push_back(g.YN2D[g.nij(i,g.knoy)]); }
    for(int j=g.knoy-1; j>=1; --j) { px.push_back(g.XN2D[g.nij(0,j)]); py.push_back(g.YN2D[g.nij(0,j)]); }
    double P = 0.0;
    for(size_t n=0; n<px.size(); ++n)
    {
        const size_t m = (n+1)%px.size();
        P += 0.5*((px[n]-px[0])*(py[m]-py[0])-(px[m]-px[0])*(py[n]-py[0]));
    }
    check(std::fabs(A-P)<=1.0e-10*std::fabs(P), std::string(name)+": cell areas sum to the boundary polygon area");

    // bilinear subdivision by 2 (as the AMR patches): children tile the parent
    double dA = 0.0, dS = 0.0;
    for(int i=0; i<g.knox; ++i)
    for(int j=0; j<g.knoy; ++j)
    {
        const double x00=g.XN2D[g.nij(i,j)], y00=g.YN2D[g.nij(i,j)];
        const double x10=g.XN2D[g.nij(i+1,j)], y10=g.YN2D[g.nij(i+1,j)];
        const double x01=g.XN2D[g.nij(i,j+1)], y01=g.YN2D[g.nij(i,j+1)];
        const double x11=g.XN2D[g.nij(i+1,j+1)], y11=g.YN2D[g.nij(i+1,j+1)];
        double X[3][3],Y[3][3];
        for(int a=0; a<3; ++a)
        for(int b=0; b<3; ++b)
        {
            const double fr=0.5*a, fs=0.5*b;
            X[a][b] = x00 + (x10-x00)*fr + (x01-x00)*fs + (x11-x10-x01+x00)*fr*fs;
            Y[a][b] = y00 + (y10-y00)*fr + (y01-y00)*fs + (y11-y10-y01+y00)*fr*fs;
        }
        double As = 0.0;
        for(int a=0; a<2; ++a)
        for(int b=0; b<2; ++b)
        As += geo2d::area(X[a][b],Y[a][b],X[a+1][b],Y[a+1][b],X[a][b+1],Y[a][b+1],X[a+1][b+1],Y[a+1][b+1]);
        const double Ap = g.AREA2D[g.sl(i,j)];
        dA = std::max(dA,std::fabs(As-Ap)/Ap);

        // east face of the parent = sum of the east faces of the two children
        const double sx = (Y[2][1]-Y[2][0]) + (Y[2][2]-Y[2][1]);
        const double sy = -(X[2][1]-X[2][0]) - (X[2][2]-X[2][1]);
        dS = std::max(dS,std::max(std::fabs(sx-g.SX1[g.sl(i,j)]),std::fabs(sy-g.SY1[g.sl(i,j)]))/std::sqrt(Ap));
    }
    check(dA<1.0e-12, std::string(name)+": 2x2 bilinear children tile the parent (area)");
    check(dS<1.0e-12, std::string(name)+": 2x2 bilinear children tile the parent (face vectors)");
}

int main()
{
    // Cartesian
    cartesian_case("uniform", 20, 12, [](int i){return 0.5*i;}, [](int j){return -3.0+0.25*j;});
    cartesian_case("stretched", 30, 17, [](int i){return 10.0*std::sinh(0.05*i);}, [](int j){return std::pow(1.07,j);});
    cartesian_case("UTM", 25, 25, [](int i){return 569000.0+1.37*i;}, [](int j){return 7034000.0+0.91*j;});
    cartesian_case("2D (one cell in y)", 40, 1, [](int i){return 0.1*i;}, [](int j){return 0.5*j;});

    // bilinear subdivision of Cartesian nodes = 1D subdivision bit for bit
    {
        testgrid g(16,9,[](int i){return 3.0+std::sinh(0.1*i);},[](int j){return 100.0+0.3*j+0.01*j*j;});
        g.geometry_alloc();
        g.geometry_cartesian_nodes();
        bool same = true;
        for(int i=0; i<g.knox; ++i)
        for(int j=0; j<g.knoy; ++j)
        for(int r=0; r<4; ++r)
        for(int s=0; s<4; ++s)
        {
            const double fr=r/4.0, fs=s/4.0;
            const double x00=g.XN2D[g.nij(i,j)], y00=g.YN2D[g.nij(i,j)];
            const double x10=g.XN2D[g.nij(i+1,j)], y10=g.YN2D[g.nij(i+1,j)];
            const double x01=g.XN2D[g.nij(i,j+1)], y01=g.YN2D[g.nij(i,j+1)];
            const double x11=g.XN2D[g.nij(i+1,j+1)], y11=g.YN2D[g.nij(i+1,j+1)];
            const double x = x00 + (x10-x00)*fr + (x01-x00)*fs + (x11-x10-x01+x00)*fr*fs;
            const double y = y00 + (y10-y00)*fr + (y01-y00)*fs + (y11-y10-y01+y00)*fr*fs;
            const double x1 = (r==0) ? g.XN[i+5] : g.XN[i+5] + (g.XN[i+6]-g.XN[i+5])*fr;
            const double y1 = (s==0) ? g.YN[j+5] : g.YN[j+5] + (g.YN[j+6]-g.YN[j+5])*fs;
            same = same && x==x1 && y==y1;
        }
        check(same, "bilinear subdivision of Cartesian nodes = 1D subdivision bit for bit");
    }

    // curvilinear
    curvilinear_case("annulus sector", 24, 10, [](double i, double j, double &x, double &y)
    {
        const double r = 80.0-2.0*j, t = 0.03*i;    // j inwards: (i,j) right-handed
        x = 1000.0 + r*std::cos(t);
        y = 2000.0 + r*std::sin(t);
    });
    curvilinear_case("sheared", 20, 15, [](double i, double j, double &x, double &y)
    {
        x = 1.0*i + 0.4*j;
        y = 0.2*i + 0.8*j;
    });
    curvilinear_case("wavy river", 60, 12, [](double i, double j, double &x, double &y)
    {
        // a meandering corridor: centreline y = 30 sin(x/40), width 2 m per cell, normal offsets
        const double s = 2.0*i, n = 2.0*(j-6.0);
        const double yc = 30.0*std::sin(s/40.0), dy = 0.75*std::cos(s/40.0);
        const double l = std::sqrt(1.0+dy*dy);
        x = s - n*dy/l;
        y = yc + n/l;
    });

    // parallelogram centroid = mean of the corners
    {
        double xc,yc;
        geo2d::centroid(0.0,0.0, 2.0,0.5, 1.0,3.0, 3.0,3.5, xc,yc);
        check(std::fabs(xc-1.5)<1.0e-14 && std::fabs(yc-1.75)<1.0e-14, "parallelogram centroid");
        geo2d::centroid(0.0,0.0, 4.0,0.0, 0.0,2.0, 1.0,2.0, xc,yc);   // trapezoid, analytic centroid
        check(std::fabs(xc-(16.0+4.0+1.0)/(3.0*5.0))<1.0e-14 && std::fabs(yc-2.0*(4.0+2.0*1.0)/(3.0*5.0))<1.0e-14, "trapezoid centroid");
    }

    std::cout<<npass<<" passed, "<<nfail<<" failed"<<std::endl;
    return nfail==0 ? 0 : 1;
}

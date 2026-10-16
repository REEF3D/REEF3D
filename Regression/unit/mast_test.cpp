// Standalone tests of REEF3D::MAST (src/mast/): no MPI, no REEF3D.
// Build:  g++ -O2 -std=c++20 -I../../ThirdParty/eigen-5.0.0 -DEIGEN_MPL2_ONLY -I../../src
//         mast_test.cpp ../../src/mast/core/mast_version.cpp -o mast_test
// Run:    ./mast_test
//
// Tests: the MAST core builds on its own (Eigen only), name and version, the core types.
// New MAST kernels add their tests here (one function per kernel, called from main).
// Architect: Hans Bihs
#include"mast/core/mast_types.h"
#include"mast/core/mast_version.h"
#include<iostream>
#include<string>
#include<cstring>

static int nfail = 0;
static void check(bool ok, const std::string &what)
{
    std::cout<<(ok ? "  PASS  " : "  FAIL  ")<<what<<std::endl;
    if(!ok) ++nfail;
}

static void test_version()
{
    std::cout<<"version"<<std::endl;
    check(std::string(mast::name())=="REEF3D::MAST", "name is REEF3D::MAST");
    check(std::strlen(mast::version())>0, "version set");
    const std::string b = mast::banner();
    check(b.find(mast::name())==0 && b.find(mast::version())!=std::string::npos
          && b.find(mast::long_name())!=std::string::npos, "banner = name version (long name): "+b);
}

static void test_types()
{
    std::cout<<"core types"<<std::endl;
    check(mast::Vec6::RowsAtCompileTime==6 && mast::Mat6::ColsAtCompileTime==6, "Vec6/Mat6 are 6x1/6x6");
    mast::Mat3 R = mast::Quat(Eigen::AngleAxisd(0.3, mast::Vec3::UnitZ())).toRotationMatrix();
    check((R.transpose()*R - mast::Mat3::Identity()).norm()<1.0e-14, "Quat -> Mat3 orthonormal");
    check(sizeof(mast::index)==sizeof(std::ptrdiff_t), "index is ptrdiff_t");
}

int main()
{
    std::cout<<mast::banner()<<std::endl;
    test_version();
    test_types();
    std::cout<<(nfail ? "FAILED: " : "all passed")<<(nfail ? std::to_string(nfail) : "")<<std::endl;
    return nfail ? 1 : 0;
}

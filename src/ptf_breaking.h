#include"ptf_fsfcd.h"
#include"sliceint4.h"

class ptf_laplace;
class field;
class fnpf_convection;

using namespace std;

#ifndef PTF_BREAKING_H_
#define PTF_BREAKING_H_

class ptf_breaking : public ptf_fsfbc, public increment 
{
public:
	ptf_breaking(lexer*, fdm*, ghostcell*);
	virtual ~ptf_breaking();
    
    virtual void breaking_algorithm(lexer*,fdm*,ghostcell*,slice&,slice&,slice&,double);

    
    
    void filter(lexer*, fdm*,ghostcell*, slice&);

    double ivel,jvel,kvel;
    
private:
    double rb3(lexer*,double);
    double rb4(lexer*,double);
    
    double dist3,dist4,db;
    
    double visc;
    
    sliceint4 bx,by;
    int count_n;
    
};

#endif

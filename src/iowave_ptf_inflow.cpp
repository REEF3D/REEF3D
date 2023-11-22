
#include"iowave.h"
#include"lexer.h"
#include"fdm_ptf.h"
#include"ghostcell.h"
#include"patchBC_interface.h"

void iowave::inflow_ptf(lexer *p, fdm_ptf* e, ghostcell* pgc, field& u, field& v, field& w)
{
    if(p->I230==0)
    {
    
	if(p->B98==3)
	dirichlet_wavegen_ptf(p,e,pgc,u,v,w);
	
	if(p->B98==4)
	active_wavegen_ptf(p,e,pgc,u,v,w);
	}
    
	if(p->B99==3||p->B99==4||p->B99==5)
	active_beach_ptf(p,e,pgc,u,v,w);
    
    pBC->patchBC_ioflow_ptf(p,e,pgc,u,v,w);
}
/* Example of a user PTO controller for REEF3D::CORAL (pto ... user <library> params <n> ...),
   used by coral_wec_test.cpp. Build:
       cc -shared -fPIC -O2 coral_pto_user_example.c -o coral_pto_user_example.so
   in:  t, q, qd, F, P, then the parameters (here K and B)
   out: K, B, F0 of the generator of the PTO
   The library counts its calls (its own state, kept over a restart with the save/load functions).
   Architect: Hans Bihs */
#include<stdlib.h>

typedef struct {double calls;} ctrl_data;

int coral_pto_control(int n_in, const double *in, int n_out, double *out, void **user)
{
    if(!*user)
        *user = calloc(1, sizeof(ctrl_data));
    if(!*user || n_in < 7 || n_out < 3)
        return 1;
    ((ctrl_data*)*user)->calls += 1.0;
    out[0] = in[5];         /* K */
    out[1] = in[6];         /* B */
    out[2] = 0.0;           /* F0 */
    return 0;
}

int coral_pto_control_save(void *user, double *buf, int n_max)
{
    if(n_max >= 1)
        buf[0] = user ? ((ctrl_data*)user)->calls : 0.0;
    return 1;
}

int coral_pto_control_load(void **user, const double *buf, int n)
{
    if(!*user)
        *user = calloc(1, sizeof(ctrl_data));
    if(!*user || n < 1)
        return 1;
    ((ctrl_data*)*user)->calls = buf[0];
    return 0;
}

void coral_pto_control_end(void *user)
{
    free(user);
}

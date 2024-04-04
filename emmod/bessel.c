#include <math.h>
#include <emmod.h>

double besj0_(double *arg)
{
    return j0(*arg);
}
double besj1_(double *arg)
{
    return j1(*arg);
}
double besj2_(double *arg)
{
    return jn(2, *arg);
}

double besjn_(int *n, double *arg)
{
    return jn(*n, *arg);
}


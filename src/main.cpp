#include "BESO.h"

int main()
{
    BESO beso(80,50, 3.0, 0.5, 0.02, 100);
    while(!beso.convergence)
    {
        beso.Optimize();
    }
    beso.FreeAll();
    return 0;
}

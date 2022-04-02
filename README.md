# BESO
This is an efficient implementation for Bi-directional Evolutionary Structural Optimization (BESO) method using C++.


## Dependencies
- [Eigen](https://gitlab.com/libeigen/eigen)

## Usage

```cpp
#include "BESO.h"

int main()
{
    BESO beso(80, 50, 3.0, 0.5, 0.02, 100);
    while(!beso.convergence)
    {
        beso.Optimize();
    }
    return 0;
}

```

## References

1. [CISM_BESO_2D](https://www.cism.org.au/tools)
2. Zuo, Z.H. and Xie, Y.M., 2015. A simple and compact Python code for complex 3D topology optimization. Advances in Engineering Software, 85, pp.1-11.
5. Huang, X. and Xie, M., 2010. Evolutionary topology optimization of continuum structures: methods and applications. John Wiley & Sons.

#include "BESO.h"

using namespace std;

int main() {
    std::cout << "Hello, World!" << std::endl;
    BESO beso(50,50, 3.0, 0.5, 0.02, 100);
    beso.Optimize();
    return 0;
}

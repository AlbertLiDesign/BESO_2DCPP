//
// Created by Albert Li on 1/4/2022.
//
#include <iostream>
#include <algorithm>
#include <iterator>
#include <vector>
#include <cmath>

#include <Eigen/Eigen>
#include <Eigen/Dense>
#include <unsupported/Eigen/KroneckerProduct>
#include <unsupported/Eigen/SparseExtra>
#include <unsupported/Eigen/src/SparseExtra/MarketIO.h>

using namespace Eigen;
using namespace std;

class BESO
{
public:
    int nelx;
    int nely;

    double rmin;
    double vf;
    double ert;
    int maxIter;

    double p;
    double Xmin;

    MatrixXd Ke;
    MatrixXd Xe;
    VectorXd dc;
    VectorXd dc_old;

    SparseMatrix<double> FltM;

    double compliance;
    double delta;
    int iter;
    double vol;
    bool convergence;
    std::vector<double> HistoryC;

    BESO(int nelx, int nely, double rmin, double vf, double ert, int maxIter);
    void Optimize();
    void FreeAll();

private:
    VectorXd F;
    VectorXd U;
    int *ik;
    int *jk;
    double *vk;
    VectorXd sh;
    VectorXi free_dofs;

    void PreFE();
    void PreFlt();
    void FE();
    void Flt();
    void ADD_DEL(double volfra);
    void GetDc();
    void GetKe();
};
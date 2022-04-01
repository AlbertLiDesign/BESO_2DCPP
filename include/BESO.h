//
// Created by Albert Li on 1/4/2022.
//
#include <iostream>
#include <math.h>
#include <Eigen/Eigen>
#include <unsupported/Eigen/KroneckerProduct>

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

    double delta;
    int iter;
    double vol;
    bool convergence;
    //std::vector<double> HistoryC;

    BESO(int nelx, int nely, double rmin, double vf, double ert, int maxIter);
    void Optimize();

    void ADD_DEL(double volfra);
    void GetDc();
    void FE();

//    void Assembly_Solve(int num_freeDofs, int num_allDofs, int num_triplets, int* free_dofs, int* ik, int* jk, double* vk, double* F, double* U)
//    {
//        std::vector<Triplet<double>> triplets;
//        triplets.reserve(num_triplets);
//
//        for (int i = 0; i < num_triplets; i++)
//        {
//            triplets.push_back(Triplet<double>(ik[i], jk[i], vk[i]));
//        }
//
//        SparseMatrix<double> K(num_allDofs, num_allDofs);
//        K.setFromTriplets(triplets.begin(), triplets.end());
//
//        std::vector<Triplet<double>> P_triplets;
//        P_triplets.reserve(num_freeDofs);
//
//        for (int i = 0; i < num_freeDofs; i++)
//        {
//            P_triplets.push_back(Triplet<double>(i, free_dofs[i], 1.0));
//        }
//
//        SparseMatrix<double> P(num_freeDofs, num_allDofs);
//        P.setFromTriplets(P_triplets.begin(), P_triplets.end());
//
//        SparseMatrix<double> K_freedof = P * K * P.transpose();
//        // saveMarket(K_freedof, "E:/test/K_freedof_cpp.mtx");
//
//        VectorXd F_freedof(num_freeDofs);
//        for (int i = 0; i < num_freeDofs; i++)
//        {
//            F_freedof(i) = F[free_dofs[i]];
//        }
//
//        VectorXd result;
//
//        //PardisoLLT<Eigen::SparseMatrix<double>> llt(K_freedof);
//        //llt.pardisoParameterArray()[59] = 0;
//        //mkl_set_num_threads(1);
//        //CholmodSimplicialLLT<Eigen::SparseMatrix<double>> llt(K_freedof);
//        //CholmodSupernodalLLT<Eigen::SparseMatrix<double>> llt(K_freedof);
//        SimplicialLLT<Eigen::SparseMatrix<double>> llt(K_freedof);
//        llt.analyzePattern(K_freedof);
//        llt.factorize(K_freedof);
//        result = llt.solve(F_freedof);
//        Eigen::VectorXd::Map(U, result.rows()) = result;
//    }

    void FreeAll();
//    void Flt(int dc_length, double* dc, double* sh);


private:
//    double *U;
    int *ik;
    int *jk;
    double *sh;
    double *vk;

    void PreFE();
    void PreFlt();
    void GetKe();
};
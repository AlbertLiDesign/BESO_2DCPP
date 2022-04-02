#include "BESO.h"

BESO::BESO(int nelx, int nely, double rmin, double vf, double ert, int maxIter)
{
    // Resolution
    this->nelx = nelx;
    this->nely = nely;
    int N = nely * nelx;

    // BESO parameters
    this->rmin = rmin;
    this->vf = vf;
    this->ert = ert;
    this->maxIter = maxIter;

    Xmin = 0.001;
    p = 3;

    // build elemental stiffness matrix
    Ke.resize(8, 8);
    GetKe();

    // prepare sensitivity numbers
    dc.resize(N);
    dc_old.resize(N);

    // prepare design variables
    Xe.resize(nely,nelx);
    for (int j = 0; j < nely; j++)
        for (int i = 0; i < nelx; i++)
            Xe(j, i) = 1.0;

    // prepare finite element method
    ik = new int[N*64];
    jk = new int[N*64];
    PreFE();
    cout<<"Prepare FEM..."<<endl;

    // prepare filtering
    PreFlt();
    cout<<"Prepare Flt..."<<endl;

    // init opt paras
    compliance = 0.0;
    delta = 1.0;
    iter = 0;
    vol = 1.0;
    convergence = false;
}

void BESO::Optimize()
{
    if (delta > 0.001 && iter < maxIter)
    {
        cout<<"====================== Iter: " + to_string(iter) + " ======================" + '\n'<<endl;
        iter += 1;
        vol = max(vf, vol * (1.0 - ert));

        // run finite element analysis
        FE();
        GetDc();
        HistoryC.emplace_back(compliance);
        Flt();
        if (iter > 1)
            for (int j = 0; j < nely; j++)
            {
                for (int i = 0; i < nelx; i++)
                {
                    dc[i * nely + j] = (dc[i * nely + j] + dc_old[i * nely + j]) * 0.5;
                }
            }
        // Record the sensitiveies in each step
        dc_old = dc;
        ADD_DEL(vol);
        // Check convergence
        if (iter > 10)
        {
            double newV = 0.0;
            double lastV = 0.0;
            for (int i = 1; i < 6; i++)
            {
                newV += HistoryC[HistoryC.size() - i];
                lastV += HistoryC[HistoryC.size() - 5 - i];
            }
            delta = abs((newV - lastV) / lastV);
        }
        cout<<"Iter: " << iter << ", Volume: " << vol << ", Compliance: " << compliance << ", Change: " << delta <<endl;
    }
    else
    {
        convergence = true;
    }
}

void BESO::PreFE()
{
    MatrixXi nodenrs(nely + 1, nelx + 1);
    int* edofVec = new int[nelx * nely];
    MatrixXi edofMat(nelx * nely, 8);
    int edofs[8] = { -1, 0, 2 * nely + 1, 2 * nely + 2, 2 * nely + 3, 2 * nely + 4, 1, 2 };

    for (size_t y = 0; y < nely + 1; y++)
    {
        for (size_t x = 0; x < nelx + 1; x++)
        {
            nodenrs(y, x) = x * (nely + 1) + y;
        }
    }

    for (size_t y = 0; y < nely; y++)
    {
        for (size_t x = 0; x < nelx; x++)
        {
            edofVec[y + x * nely] = 2 * nodenrs(y, x) + 1;
        }
    }

    for (size_t i = 0; i < nelx * nely; i++)
    {
        for (size_t j = 0; j < 8; j++)
        {
            edofMat(i, j) = edofVec[i] + edofs[j];
        }
    }

    auto a = kroneckerProduct(edofMat, MatrixXi::Ones(8, 1)).eval();
    auto za = a.transpose();
    auto b = kroneckerProduct(edofMat, MatrixXi::Ones(1, 8)).eval();
    auto zb = b.transpose();

    for (size_t i = 0; i < za.cols(); i++)
    {
        for (size_t j = 0; j < za.rows(); j++)
        {
            ik[i * za.rows() + j] = za(j, i);
        }
    }

    for (size_t i = 0; i < zb.cols(); i++)
    {
        for (size_t j = 0; j < zb.rows(); j++)
        {
            jk[i * zb.rows() + j] = zb(j, i);
        }
    }
    delete[] edofVec;
}

void BESO::PreFlt() {
    int rminf = (int)floor(rmin) - 1;

    int length = (int)(nelx * nely * pow((2 * rminf + 1), 2));
    int ih[length];
    int jh[length];
    double vh[length];
    sh.resize(nelx * nely);

    int sum = 0;
    for (int i = 0; i < nelx; i++)
    {
        for (int j = 0; j < nely; j++)
        {
            auto e1 = i * nely + j + 1;
            for (int k = std::max(i - rminf, 0); k < std::min(i + rminf + 1, nelx); k++)
            {
                for (int l = std::max(j - rminf, 0); l < std::min(j + rminf + 1, nely); l++)
                {
                    auto e2 = k * nely + l + 1;
                    ih[sum] = e1 - 1;
                    jh[sum] = e2 - 1;
                    vh[sum] = std::max(0.0, rmin - sqrt((i - k) * (i - k) + (j - l) * (j - l)));
                    sum++;
                }
            }
        }
    }

    // get the sum of each row
    vector<Triplet<double>> triplets;
    for (size_t i = 0; i < sum; i++)
        triplets.emplace_back(ih[i], jh[i], vh[i]);
    FltM.resize(nelx * nely, nelx * nely);
    FltM.setFromTriplets(triplets.begin(), triplets.end());
    VectorXd result = FltM * VectorXd::Ones(FltM.cols());
    VectorXd::Map(sh.data(), result.rows()) = result;
}


void BESO::ADD_DEL(double volfra)
{
    auto lowest = min_element(dc.data(),dc.data()+dc.size());
    auto highest = max_element(dc.data(),dc.data()+dc.size());
    double th = 0.0;
    double _vol = volfra * nelx * nely;
    while (((highest[0] - lowest[0]) / highest[0]) > 1e-5)
    {
        th = (highest[0] + lowest[0]) * 0.5;
        double sum = 0.0;
        for (int j = 0; j < nely; j++)
        {
            for (int i = 0; i < nelx; i++)
            {
                Xe(j, i) = dc(i * nely + j) > th ? 1.0 : Xmin;
                sum += Xe(j, i);
            }
        }
        if (sum - _vol > 0.0) lowest[0] = th;
        else highest[0] = th;
    }
}

void BESO::GetDc()
{
    compliance = 0.0;
    for (int ely = 0; ely < nely; ely++)
    {
        for (int elx = 0; elx < nelx; elx++)
        {
            auto n1 = (nely + 1) * elx + ely + 1;
            auto n2 = (nely + 1) * (elx + 1) + ely + 1;

            VectorXd Ue(8);
            Ue <<  U[2 * n1 - 2], U[2 * n1 - 1], U[2 * n2 - 2], U[2 * n2 - 1],
                            U[2 * n2], U[2 * n2 + 1], U[2 * n1], U[2 * n1 + 1];
            double v = (Ue.transpose() * Ke * Ue).value();
            compliance += 0.5 * pow(Xe(ely, elx), p) * v;
            dc[elx * nely + ely] = 0.5 * pow(Xe(ely, elx), p - 1) * v;
        }
    }
}

void BESO::FE()
{
    int num_allDofs = 2 * (nelx + 1) * (nely + 1);
    int num_fixedDofs = 2 * (nely + 1);
    int num_freeDofs = num_allDofs - num_fixedDofs;

    // Assemble stiffness matrix with all DOFs
    vk = new double[64 * nelx * nely];
    for (int i = 0; i < nelx; i++)
    {
        for (int j = 0; j < nely; j++)
        {
            auto ex = pow(Xe(j, i), p);
            for (int a = 0; a < 8; a++)
            {
                for (int b = 0; b < 8; b++)
                {
                    vk[i * nely * 64 + j * 64 + a * 8 + b] = ex * Ke(a * 8 + b);
                }
            }
        }
    }

    F.resize(num_allDofs);
    U.resize(num_allDofs);

    // Define force vector
    F(2 * (nelx + 1) * (nely + 1) - nely - 1) = -1.0;

    // Define fixed dofs
    VectorXi fixed_dofs(num_fixedDofs),all_dofs(num_allDofs);
    for (int i = 0; i < num_fixedDofs; i++)
        fixed_dofs[i] = i;

    for (int i = 0; i < num_allDofs; i++)
        all_dofs[i] = i;

    // Obtain free dofs (get the difference)
    free_dofs.resize(num_allDofs);
    auto it = set_difference(all_dofs.data(), all_dofs.data() + all_dofs.size(),
                   fixed_dofs.data(),fixed_dofs.data()+fixed_dofs.size(),
                   free_dofs.data());

    free_dofs.conservativeResize(std::distance(free_dofs.data(), it)); // resize the result

    // Assemble the global stiffness matrix
    VectorXd U_freedof(num_freeDofs);

    std::vector<Triplet<double>> triplets;
    int num_triplets = nelx * nely * 64;
    triplets.reserve(num_triplets);

    for (int i = 0; i < num_triplets; i++)
    {
        triplets.emplace_back(ik[i], jk[i], vk[i]);
    }

    SparseMatrix<double> K(num_allDofs, num_allDofs);
    K.setFromTriplets(triplets.begin(), triplets.end());

    std::vector<Triplet<double>> P_triplets;
    P_triplets.reserve(num_freeDofs);

    for (int i = 0; i < num_freeDofs; i++)
    {
        P_triplets.emplace_back(i, free_dofs[i], 1.0);
    }

    SparseMatrix<double> P(num_freeDofs, num_allDofs);
    P.setFromTriplets(P_triplets.begin(), P_triplets.end());

    SparseMatrix<double> K_freedof = P * K * P.transpose();
    //saveMarket(K_freedof, "K_freedof_cpp.mtx");

    VectorXd F_freedof(num_freeDofs);
    for (int i = 0; i < num_freeDofs; i++)
    {
        F_freedof(i) = F[free_dofs[i]];
    }

    // Solve FEM
    VectorXd result;
    SimplicialLLT<Eigen::SparseMatrix<double>> llt(K_freedof);
    llt.analyzePattern(K_freedof);
    llt.factorize(K_freedof);
    result = llt.solve(F_freedof);
    VectorXd::Map(U_freedof.data(), result.rows()) = result;

    for (int i = 0; i < num_freeDofs; i++)
    {
        U[free_dofs[i]] = U_freedof[i];
    }
}

void BESO::Flt()
{
    VectorXd result = (FltM.selfadjointView<Lower>() * dc).array() / sh.array();
    Eigen::VectorXd::Map(dc.data(), result.rows()) = result;
}

void BESO::GetKe()
{
    double E = 1.0;
    double nu = 0.3;
    double w = E / (1 - nu * nu);
    double k[8] =
            {
                w * (0.5 - nu / 6.0), w * (0.125 + nu / 8.0), w * (-0.25 - nu / 12.0), w * (-0.125 + 3 * nu / 8.0),
                w * (-0.25 + nu / 12.0), w * (-0.125 - nu / 8.0), w * (nu / 6.0), w * (0.125 - 3 * nu / 8.0)
            };
    Ke <<
       k[0],k[1],k[2],k[3],k[4],k[5],k[6],k[7],
       k[1],k[0],k[7],k[6],k[5],k[4],k[3],k[2],
       k[2],k[7],k[0],k[5],k[6],k[3],k[4],k[1],
       k[3],k[6],k[5],k[0],k[7],k[2],k[1],k[4],
       k[4],k[5],k[6],k[7],k[0],k[1],k[2],k[3],
       k[5],k[4],k[3],k[2],k[1],k[0],k[7],k[6],
       k[6],k[3],k[4],k[1],k[2],k[7],k[0],k[5],
       k[7],k[2],k[1],k[4],k[3],k[6],k[5],k[0];
}

void BESO::FreeAll()
{
    delete []ik;
    delete []jk;
    delete []vk;
}
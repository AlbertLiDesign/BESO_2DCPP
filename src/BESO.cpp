#include "../include/BESO.h"

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
    delta = 1.0;
    iter = 0;
    vol = 1.0;
    convergence = false;
}

void BESO::Optimize()
{
//    if (delta > 0.001 && iter < maxIter)
//    {
//        cout<<"====================== Iter: " + to_string(iter) + " ======================" + '\n'<<endl;
//        iter += 1;
//        vol = max(vf, vol * (1.0 - ert));
//
//        // run finite element analysis
//        FE();
//#region Prepare report
//        optInfo.Append("FEA:" + stopwatch.Elapsed.TotalMilliseconds + '\n');
//#endregion
//
//#region Get DC
//        stopwatch.Restart();
//        GetDc();
//        HistoryC.Add(Compliance);
//        stopwatch.Stop();
//#endregion
//#region Prepare report
//        optInfo.Append("Getting Sensitivity:" + stopwatch.Elapsed.TotalMilliseconds + '\n');
//#endregion
//
//#region Flt
//        stopwatch.Restart();
//        Flt(dc.Length, dc, sh);
//
//        if (iter > 1)
//            for (int j = 0; j < nely; j++)
//            {
//                for (int i = 0; i < nelx; i++)
//                {
//                    dc[i * nely + j] = (dc[i * nely + j] + dc_old[i * nely + j]) * 0.5;
//                }
//            }
//
//        // Record the sensitiveies in each step
//        dc_old = (double[])dc.Clone();
//        stopwatch.Stop();
//#endregion
//#region Prepare report
//        optInfo.Append("Flt:" + stopwatch.Elapsed.TotalMilliseconds + '\n');
//#endregion
//
//#region ADD & DEL
//        stopwatch.Restart();
//        ADD_DEL(vol);
//        stopwatch.Stop();
//#endregion
//#region Prepare report
//        optInfo.Append("ADD & DEL:" + stopwatch.Elapsed.TotalMilliseconds + '\n');
//#endregion
//
//
//#region Checking Convergence
//        stopwatch.Restart();
//
//        // Check convergence
//        if (iter > 10)
//        {
//            var newV = 0.0;
//            var lastV = 0.0;
//            for (int i = 1; i < 6; i++)
//            {
//                newV += HistoryC[HistoryC.Count - i];
//                lastV += HistoryC[HistoryC.Count - 5 - i];
//            }
//            delta = Math.Abs((newV - lastV) / lastV);
//        }
//#endregion
//#region Prepare report
//        optInfo.Append("Checking Convergence:" + stopwatch.Elapsed.TotalMilliseconds + '\n');
//
//        optInfo.Append("Volume: " + vol.ToString() + '\n');
//        optInfo.Append("Compliance: " + Compliance.ToString() + '\n');
//        optInfo.Append("Change: " + delta.ToString() + '\n');
//#endregion
//
//        info = "Iter: " + iter.ToString() + ", Volume: " + vol.ToString()
//               + ", Compliance: " + Compliance.ToString() + ", Change: " + delta.ToString();
//    }
//    else
//    {
//        convergence = true;
//    }
//
//    FreeAll();
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
    int *ih = new int[length];
    int *jh = new int[length];
    auto *vh = new double[length];
    sh = new double[nelx * nely];

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
        triplets.push_back(Triplet<double>(ih[i], jh[i], vh[i]));
    FltM.resize(nelx * nely, nelx * nely);
    FltM.setFromTriplets(triplets.begin(), triplets.end());
    VectorXd result = FltM * VectorXd::Ones(FltM.cols());

    VectorXd::Map(sh, result.rows()) = result;

    delete []ih;
    delete []jh;
    delete []vh;
}


void BESO::ADD_DEL(double volfra) {

}

void BESO::GetDc() {

}

//void BESO::FE()
//{
//    int num_allDofs = 2 * (nelx + 1) * (nely + 1);
//    int num_fixedDofs = 2 * (nely + 1);
//    int num_freeDofs = num_allDofs - num_fixedDofs;
//
//    // Assemble stiffness matrix with all DOFs
//    vk = new double[64 * nelx * nely];
//    for (int i = 0; i < nelx; i++)
//    {
//        for (int j = 0; j < nely; j++)
//        {
//            auto ex = pow(Xe[j, i], p);
//            for (int a = 0; a < 8; a++)
//            {
//                for (int b = 0; b < 8; b++)
//                {
//                    vk[i * nely * 64 + j * 64 + a * 8 + b] = ex * Ke[a * 8 + b];
//                }
//            }
//        }
//    }
//
//    auto F = new double[num_allDofs];
//    U = new double[num_allDofs];
//
//    // Define force vector
//    F[2 * (nelx + 1) * (nely + 1) - nely - 1] = -1.0;
//
//    if (changeSupports)
//    {
//        // Define fixed dofs
//        var fixed_dofs = new int[num_fixedDofs];
//        for (int i = 0; i < num_fixedDofs; i++)
//            fixed_dofs[i] = i;
//
//        var all_dofs = new int[num_allDofs];
//        for (int i = 0; i < num_allDofs; i++)
//            all_dofs[i] = i;
//
//        // Obtain free dofs
//        free_dofs = all_dofs.Except(fixed_dofs).ToArray();
//        changeSupports = false;
//    }
//
//    var U_freedof = new double[num_freeDofs];
//    Assembly_Solve(num_freeDofs, num_allDofs, ik.Length, free_dofs, ik, jk, vk, F, U_freedof);
//
//    for (int i = 0; i < num_freeDofs; i++)
//    {
//        U[free_dofs[i]] = U_freedof[i];
//    }
//}


//
//void BESO::Assembly_Solve(int num_freeDofs, int num_allDofs, int num_triplets, int* free_dofs, int* ik, int* jk, double* vk, double* F, double* U)
//{
//    std::vector<Triplet<double>> triplets;
//    triplets.reserve(num_triplets);
//
//    for (int i = 0; i < num_triplets; i++)
//    {
//        triplets.push_back(Triplet<double>(ik[i], jk[i], vk[i]));
//    }
//
//    SparseMatrix<double> K(num_allDofs, num_allDofs);
//    K.setFromTriplets(triplets.begin(), triplets.end());
//
//    std::vector<Triplet<double>> P_triplets;
//    P_triplets.reserve(num_freeDofs);
//
//    for (int i = 0; i < num_freeDofs; i++)
//    {
//        P_triplets.push_back(Triplet<double>(i, free_dofs[i], 1.0));
//    }
//
//    SparseMatrix<double> P(num_freeDofs, num_allDofs);
//    P.setFromTriplets(P_triplets.begin(), P_triplets.end());
//
//    SparseMatrix<double> K_freedof = P * K * P.transpose();
//    // saveMarket(K_freedof, "E:/test/K_freedof_cpp.mtx");
//
//    VectorXd F_freedof(num_freeDofs);
//    for (int i = 0; i < num_freeDofs; i++)
//    {
//        F_freedof(i) = F[free_dofs[i]];
//    }
//
//    VectorXd result;
//
//    //PardisoLLT<Eigen::SparseMatrix<double>> llt(K_freedof);
//    //llt.pardisoParameterArray()[59] = 0;
//    //mkl_set_num_threads(1);
//    //CholmodSimplicialLLT<Eigen::SparseMatrix<double>> llt(K_freedof);
//    //CholmodSupernodalLLT<Eigen::SparseMatrix<double>> llt(K_freedof);
//    SimplicialLLT<Eigen::SparseMatrix<double>> llt(K_freedof);
//    llt.analyzePattern(K_freedof);
//    llt.factorize(K_freedof);
//    result = llt.solve(F_freedof);
//    Eigen::VectorXd::Map(U, result.rows()) = result;
//}
//
//double TransposeMultiply(int rows, int cols, double* A, double* U)
//{
//    Map<Matrix<double,Dynamic,Dynamic,RowMajor>> A_(A, rows,cols);
//    Map<VectorXd> U_(U, rows);
//    auto result = U_.transpose() * A_ * U_;
//    return result.value();
//}
//

//
//void BESO::Flt(int dc_length, double* dc, double* sh)
//{
//    Map<VectorXd> dc_(dc, dc_length);
//    Map<VectorXd> sh_(sh, dc_length);
//    VectorXd result = (H.selfadjointView<Lower>() * dc_).array() / sh_.array();
//    Eigen::VectorXd::Map(dc, result.rows()) = result;
//}

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
    delete []sh;
    delete []vk;
}
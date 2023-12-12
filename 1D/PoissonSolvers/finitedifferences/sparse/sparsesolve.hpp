#include<Eigen/Sparse>
#include<Eigen/IterativeLinearSolvers>

#include <cstddef>

size_t BuildSparseLapl(Eigen::SparseMatrix<double> &A, const double dx){
    size_t status = 0;
    size_t Ax = A.rows();
    
    /*
    Check if A.rows = A.cols
    */ 

    // Poissons Equation => Tridiagonal stencil
    for (size_t iA = 0; iA < Ax; iA++){
        if (iA > 0) A.insert(iA,iA-1) = 1.0;
        A.insert(iA,iA) = -2.0;
        if (iA < Ax - 1) A.insert(iA,iA+1) = 1.0; 
    }

    // Gauge
    // for (size_t iA = 0; iA < Ax; iA++){
    //     A.coeffRef(Ax-1,iA) = 1.0;
    // }

    // Periodic Boundaries
    A.coeffRef(0,Ax-2) = 1.0;
    A.coeffRef(Ax-1,1) = 1.0;
    
    A = (-1.0 / pow(dx,2)) * A;

    A.finalize();

    return status;
}

size_t SparseFieldSolve()
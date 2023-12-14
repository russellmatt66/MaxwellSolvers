#ifndef SPARSE_HPP
#define SPARSE_HPP

#include<Eigen/Sparse>
#include<Eigen/IterativeLinearSolvers>
#include <cstddef>

#include "../../include/Grid1d.hpp"

// Finite Difference Stencil
size_t BuildSparseLaplacian(Eigen::SparseMatrix<double> &A, const double dx){
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

// Solve with Eigen::SparseLU
size_t SparseLUFieldSolve(const Eigen::SparseMatrix<double>& A, Grid1d& Grid, Eigen::VectorXd rhoEig, Eigen::VectorXd phiEig, const double dx, const size_t Nx){
    // size_t Nx = Grid.Nx();
    // double dx = Grid.dx();
    // Initialize Eigen containers
    for (size_t j = 0; j < rhoEig.rows(); j++){
        rhoEig[j] = Grid.RhoX(j); // A*phi = -rho
        phiEig[j] = 0.0; // just initialize phi to 0 for simplicity
    }

    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.compute(A);

    phiEig = solver.solve(rhoEig);
    double residual = (A * phiEig - rhoEig).norm();

    for (size_t ij = 0; ij < phiEig.rows(); ij++){
        Grid.PhiX(ij) = phiEig[ij];  
    }
    Grid.PhiX(Nx - 1) = Grid.PhiX(0); // PBC

    // Calculate grid electric field 
    for (size_t j = 0; j < Nx; j++){
        if (j == 0) // first-order central difference
            Grid.EX(j) = -(Grid.PhiX(1) - Grid.PhiX(Nx-2)) / (2.0 * dx);
        else if (j == (Nx - 1))
            Grid.EX(j) = Grid.EX(0); // E_{Nx-1} = E_{0} is second PBC
        else Grid.EX(j) = -(Grid.PhiX(j + 1) - Grid.PhiX(j - 1)) / (2.0 * dx); 
    }
    return 0;
}
#endif
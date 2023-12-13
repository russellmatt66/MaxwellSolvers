#define _USE_MATH_DEFINES

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>
#include <string>

#include "../include/Grid1d.hpp"
#include "../finitedifferences/sparse/sparsesolve.hpp"

void initDataFile(std::ofstream& datafile);
void callMatrixSolver(std::ofstream& datafile, size_t (*Solver)(const Eigen::SparseMatrix<double>&, Grid1d& , Eigen::VectorXd, Eigen::VectorXd, double, size_t), const Eigen::SparseMatrix<double>& A, Grid1d& Grid, Eigen::VectorXd rhoEig, Eigen::VectorXd phiEig);

namespace fs = std::filesystem;

// Driver code for testing finite differences
int main(int argc, char* argv[]){
    size_t Nx = std::stoi(argv[1]);

    std::string dir_path = "./data/fd/Nx" + std::to_string(Nx);
    fs::create_directory(dir_path);

    std::ofstream datafile_sparseLU;    
    datafile_sparseLU.open("./data/fd/Nx" + std::to_string(Nx) + "/sparseLU.csv");

    // Initialize datafiles
    initDataFile(datafile_sparseLU);

    // Initialize stencil, grid, and Eigen containers
    // double x_min = std::stod(argv[2]), x_max = std::stod(argv[3]);
    double x_min = -M_PI, x_max = M_PI;
    Grid1d testGrid (Nx,x_min,x_max);
    Eigen::SparseMatrix<double> A(Nx,Nx);
    A.reserve(Eigen::VectorXd::Constant(Nx,3)); // Poisson's equation so triangular
    Eigen::VectorXd rhoEig (A.rows());
    Eigen::VectorXd phiEig (A.rows());

    double dx = testGrid.dx();
    size_t RoutineFlag = BuildSparseLaplacian(A, dx);

    // Call sparse LU solve
    callMatrixSolver(datafile_sparseLU, SparseLUFieldSolve, A, testGrid, rhoEig, phiEig);

    return 0;
}

void initDataFile(std::ofstream& datafile){
    datafile << "j,rho_j,phi_j,E_j" << std::endl;
}

/* Wrapper for calling a given solver */
void callMatrixSolver(std::ofstream& datafile, size_t (*Solver)(const Eigen::SparseMatrix<double>&, Grid1d& , Eigen::VectorXd, Eigen::VectorXd, double, size_t), const Eigen::SparseMatrix<double>& A, Grid1d& Grid, Eigen::VectorXd rhoEig, Eigen::VectorXd phiEig){
    double dx = Grid.dx();
    size_t Nx = Grid.Nx();
    
    Solver(A,Grid,rhoEig,phiEig,dx,Nx);

    for (size_t j = 0; j < Nx; j++){
        datafile << j << "," << Grid.RhoX(j) << "," << Grid.PhiX(j) << "," << Grid.EX(j) << std::endl;
    }
}
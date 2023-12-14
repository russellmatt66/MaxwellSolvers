#define _USE_MATH_DEFINES

#include <cstdlib>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>
#include <string>

#include "../../include/Grid1d.hpp"
#include "../../fd/sparse/sparsesolve.hpp"

void initDataFile(std::ofstream& datafile);
void callMatrixSolver(std::ofstream& datafile, size_t (*Solver)(const Eigen::SparseMatrix<double>&, Grid1d& , Eigen::VectorXd, Eigen::VectorXd, double, size_t), const Eigen::SparseMatrix<double>& A, Grid1d& Grid, Eigen::VectorXd rhoEig, Eigen::VectorXd phiEig);
void initGridRho(Grid1d&, double (*RhoFunc)(double), size_t);

namespace fs = std::filesystem;

// Driver code for running sparse finite differences, based on proof-of-concept in test/fd-sparse.cpp
int main(int argc, char* argv[]){
    /* Parse input file */
    size_t Nx = std::stoi(argv[1]);

    /* Should there be a separate directory for the sparse data? */
    // Create data directory
    std::string dir_path = "../../data/fd/Nx" + std::to_string(Nx);
    fs::create_directory(dir_path);

    std::ofstream datafile;    
    datafile.open("../../data/fd/Nx" + std::to_string(Nx) + "/sparseLU.csv");

    // Initialize datafiles
    initDataFile(datafile);
   
    // Create stencil, grid, and Eigen containers
    double x_min = -M_PI, x_max = M_PI; /* Base values on input file */
    Grid1d testGrid (Nx,x_min,x_max);
    double dx = testGrid.dx();

    // Document the form of the stencil
    Eigen::SparseMatrix<double> A(Nx,Nx);
    A.reserve(Eigen::VectorXd::Constant(Nx,3)); // Poisson's equation so triangular
    Eigen::VectorXd rhoEig (A.rows());
    Eigen::VectorXd phiEig (A.rows());
    size_t RoutineFlag = BuildSparseLaplacian(A, dx);

    /* 
    Refactor this so it runs every k \in [ki,kf], for every Nx \in [Nxi,Nxf]
    Make sure to update where the data is going, and put the data in the correct location
    */
    // Initialize Grid.RhoX
    initGridRho(testGrid, sin, 1);

    /* Call the solver specified in fd.inp */
    callMatrixSolver(datafile, SparseLUFieldSolve, A, testGrid, rhoEig, phiEig);

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

void initGridRho(Grid1d& Grid, double (*RhoFunc)(double x), size_t n){
    const size_t Nx = Grid.Nx();
    for (size_t j = 0; j < Nx; j++){
        Grid.RhoX(j) = RhoFunc(n * Grid.Xgrid(j));
    }
}

void parseInputFileFD(std::ofstream& inputfile){
    /* Parse the input file */
}

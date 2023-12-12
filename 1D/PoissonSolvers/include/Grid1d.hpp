#ifndef GRID_1D_HPP
#define GRID_1D_HPP
#include <cstddef> 
#include <vector>

class Grid1d{
    public:
        Grid1d(size_t Nx) : 
            Nx_(Nx), 
            xgrid_(Nx,0.0), 
            rhox_(Nx,0.0), Ex_(Nx,0.0), phix_(Nx,0.0) 
        {}

    /* Write accessors */

    private:
        size_t Nx_;
        std::vector<double> xgrid_;
        std::vector<double> rhox_;
        std::vector<double> Ex_;
        std::vector<double> phix_;
};
#endif
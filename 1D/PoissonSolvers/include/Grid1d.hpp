#ifndef GRID_1D_HPP
#define GRID_1D_HPP
#include <cstddef> 
#include <vector>

class Grid1d{
    public:
        Grid1d(size_t Nx, double x_min, double x_max) : 
            Nx_(Nx), xmin_(x_min), xmax_(x_max),
            xgrid_(Nx,0.0), rhox_(Nx,0.0), Ex_(Nx,0.0), phix_(Nx,0.0) 
        {
            // Initialize xgrid_
            double L = x_max - x_min;
            double dx = L / (Nx_ - 1);
            for (size_t j = 0; j < Nx_; j++){
                xgrid_[j] = x_min + j * dx; 
            }
        }

        /* Write accessors */
        const double RhoX(size_t j) const { return rhox_[j]; }
        const double PhiX(size_t j) const { return phix_[j]; }
        const double EX(size_t j) const { return Ex_[j]; }
        const double Xgrid(size_t j) const { return xgrid_[j]; }

        double& RhoX(size_t j) { return rhox_[j]; }
        double& PhiX(size_t j) { return phix_[j]; }
        double& EX(size_t j) { return Ex_[j]; }
        double& Xgrid(size_t j) { return xgrid_[j]; }

        const size_t Nx() const { return Nx_; }
        const double L() const { return xmax_ - xmin_;}
        const double dx() const { return L() / (Nx_ - 1); }

    private:
        size_t Nx_;
        double xmin_, xmax_;
        std::vector<double> xgrid_;
        std::vector<double> rhox_;
        std::vector<double> Ex_;
        std::vector<double> phix_;
};
#endif
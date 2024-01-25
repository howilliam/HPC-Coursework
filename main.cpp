/**
 * @file main.cpp
 *
 * High-Performance Computing Coursework.
 * Solves Shallow Water Equation PDE, parallelised code with OpenMP.
 * 
 */
#include <iostream>
#include <cmath>
#include <boost/program_options.hpp>
#include <omp.h>
#include "ShallowWater.h"

using namespace std;
namespace po = boost::program_options;

/**
 * @brief Solves SWE PDE using 6th order stencil and RK4.
 */
int main(int argc, char* argv[]) {
    // Task 1: Command line input parser
    // Declare the supported options.
    po::options_description opts("Available options.");
    opts.add_options()
    ("dt", po::value<double>()->default_value(0.1), "Time-step to use")
    ("T", po::value<int>()->default_value(80), "Total integration time")
    ("Nx", po::value<int>()->default_value(100), "Number of grid points in x")
    ("Ny", po::value<int>()->default_value(100), "Number of grid points in y")
    ("ic", po::value<int>()->default_value(4), "Index of the initial condition to use (1-4)")
    ("option1", po::value<int>()->default_value(1), "For loop or BLAS option");

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, opts), vm);
    po::notify(vm);

    //Prasing arguments without need for user input on terminal
    double dt = vm["dt"].as<double>();
    double T = vm["T"].as<int>();
    int Nx = vm["Nx"].as<int>();
    int Ny = vm["Ny"].as<int>();
    int ic = vm["ic"].as<int>();
    int option1 = vm["option1"].as<int>();
    // Display current parameters
    cout << "Arguments: " << dt << ", " << T << ", " << Nx
         << ", " << Ny << ", " << ic << endl;

    if (dt <= 0.0 || T < 0.0 || Nx < 0 || Ny < 0 || ic < 1 ||
        ic > 4 || option1 > 2 || option1 < 1) {
            throw invalid_argument("Invalid parameters!");
        }
    
    int dx=1, dy=1;
    int np = 4;
    // Set number of threads for OpenMP
    omp_set_num_threads(np);

    // Store parameters in Shallow Water Object
    ShallowWater solve(dt, T, Nx, Ny, ic, dx, dy);

    // Set initial conditions
    solve.SetInitialConditions();

    // Perform integration
    solve.TimeIntegrate(option1);

    // File output, values of u,v,h at each gridpoint x,y
    solve.writeOutput();

    return 0;
}
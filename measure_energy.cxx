#include "ComplexNum.hxx"
#include "GaugeLattice.cxx"
#include "SU3.hxx"
#include "numbers.hxx"
#include "TreeWriter.hxx"

#include "argparse.hpp"

#include <TH1D.h> 
#include <TH2D.h> 
#include <TCanvas.h> 
#include <TStopwatch.h> 
#include <TString.h> 

#include <iostream>
#include <string> 
#include <stdexcept> 
#include <fstream> 

namespace {
    //dimension of the lattice 
    constexpr int D = 4; 
}

int main(int argc, char* argv[])
{   
    //parse arguments
    argparse::ArgumentParser program("measure_energy"); 

    //positional argument is the output file 
    std::string path_output;
    program.add_argument("path_output")
        .required()
        .help("Path to the output file destination. Must end with '.root'")
        .store_into(path_output);

    //lattice dimension
    program.add_argument("-s", "--");

    // lattice size 
    int lattice_size; 
    program.add_argument("-s", "--size")
        .help("Size of lattice-edge")
        .metavar("N")
        .default_value(8)
        .scan<'i', int>()
        .nargs(1)
        .store_into(lattice_size);

    // lowest Beta 
    double beta_min;
    program.add_argument("--beta-low")
        .help("lowest value of therm. Beta")
        .default_value(0.)
        .metavar("B0")
        .scan<'g', double>()
        .nargs(1)
        .store_into(beta_min); 
    
    // highest Beta
    double beta_max;
    program.add_argument("--beta-high")
        .help("highest value of therm. Beta")
        .default_value(4.5)
        .metavar("B1")
        .scan<'g', double>()
        .nargs(1)
        .store_into(beta_max); 
    
    // number of steps to measure
    int n_steps;
    program.add_argument("-n", "--n-steps")
        .help("number of measurement steps")
        .scan<'i', int>()
        .default_value((int)1e5)
        .metavar("N")
        .nargs(1)
        .store_into(n_steps); 
    
    //number of updates per temp-step
    int n_updates_per_step;
    program.add_argument("--updates-per-step")
        .help("Number of updates considered for each temperature step")
        .scan<'i', int>()
        .default_value((int)5e5) 
        .metavar("N")
        .nargs(1)
        .store_into(n_updates_per_step); 

    //choose whether heating or cooling will be executed
    program.add_argument("-m", "--mode")
        .help("Mode to operate in [heating/cooling]")
        .required()
        .metavar("MODE")
        .choices("heating", "cooling")
        .nargs(1); 
    
    try {
        program.parse_args(argc, argv);
    }
    catch (const std::exception& err) {
        std::cerr << err.what() << std::endl;
        std::cerr << program;
        return 1; 
    }

    //print help, and exit
    if (program["--help"] == true) {
        std::cout << program; 
        return 0; 
    }
    const bool cooling = (program.get<std::string>("--mode") == "cooling"); 

    //it's not strictly necessary to define this bool, but it may make the code that follows more readable
    const bool heating = !cooling; 

    std::printf(" D = %i, lattice size: %i, beta = [%.3f, %.3f] (%s). temp. steps: %i, updates per temp-step: %i\n", 
        D, 
        lattice_size, 
        beta_min, 
        beta_max, 
        program.get<std::string>("--mode").c_str(), 
        n_steps, 
        n_updates_per_step
    );
    
    TreeWriter writer("data_tree", path_output); 

    //consturct a gauge lattice, and pick the appropriate dimension
    GaugeLattice<D> lattice(lattice_size); 

    double beta = cooling ? beta_min : beta_max; 

    //change in beta between each step
    const double dBeta = 
        (cooling ? +1. : -1.)*
        (beta_max - beta_min)/((double)n_steps-1); 

    lattice.SetBeta( beta );

    if (cooling) {
        lattice.HotStart(100); //start with random lattice sites (hot) 
    } else {
        lattice.ColdStart(); //start with each lattice site as the identity (cold)
    }
    //this lattice will help us track how much the overall lattice is changing by 
    auto old_lattice = lattice; 

    const double target_prob = 0.5;
    const double theta_change = 1.15; //how much theta can be scaled by
    double theta = Nums::pi/8; 
    const double max_theta = Nums::pi; 

    int n_sites = std::pow(lattice_size, D); 

    //number of updates to consider for each site, for each update step
    const long long int n_steps_per_site = 10; 

    for (int i=0; i<n_steps; i++) {

        lattice.SetMaxTheta(theta); 
        lattice.SetBeta(beta); 

        TStopwatch timer; 
        double accept_prob = lattice.MetropolisUpdate(n_updates_per_step, n_steps_per_site);
        double real_time = timer.RealTime(); 
        double cpu_time  = timer.CpuTime(); 

        double avg_norm = lattice.GetFrobDistance(old_lattice); 

        double energy = lattice.GetEnergy(n_sites); 
    
        //write this results
        writer.WriteLine(
            beta, 
            energy, 
            accept_prob, 
            theta, 
            avg_norm,
            cpu_time/((double)n_updates_per_step*n_steps_per_site) 
        ); 

        //if the metropolis acceptance is above the target, we should scan the group more aggresively (increase theta)
        //otherwise if the accept. probability is too small, we need to scan the group less aggresively (reduce theta)
        theta = theta*( 1. + 0.5*(accept_prob-target_prob) );

        if (theta > max_theta) theta = max_theta; 

        //update beta
        beta += dBeta; 
        
        old_lattice = lattice; 

        printf("i: %-3i, beta: %.3f time %.3f (%.3f us/flip) accept. prob: % .5f, theta: % .5f/pi, avg norm: % .5f avg energy: % .5f\n", 
            i, 
            beta,
            real_time, 
            1e6*real_time/((double)n_updates_per_step*n_steps_per_site), 
            accept_prob, 
            theta/Nums::pi, 
            avg_norm, 
            energy
        );
    }
    
    writer.CloseFile(); 

    printf("done.\n"); 

    return 0; 
}
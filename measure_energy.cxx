#include "ComplexNum.hxx"
#include "GaugeLattice.cxx"
#include "SU3.hxx"
#include "numbers.hxx"

#include <TH1D.h> 
#include <TH2D.h> 
#include <TCanvas.h> 

#include <iostream>

int main(int argc, char* argv[])
{   
    using namespace std; 

    //consturct a gauge lattice 
    GaugeLattice<4> lattice(10); 

    auto hist_energy = new TH1D("h_dist", "Avg energy of lattice;#braket{E};", 200, 0., 2.); 

    //randomly roatate matrices
    lattice.SetBeta( 0. );

    auto old_lattice = lattice; 

    lattice.SetMaxTheta( Nums::pi/2. );

    for (int i=0; i<1000; i++) {

        double accept_prob = lattice.MetropolisUpdate(5e5, 10);

        double avg_norm = lattice.GetFrobDistance(old_lattice); 

        double energy = lattice.GetEnergy(1000); 

        hist_energy->Fill( energy );
        old_lattice = lattice; 

        printf("i: %-3i, accept. prob: % .5f, avg norm: % .5f avg energy: % .5f\n", i, accept_prob, avg_norm, energy);
    }

    auto c = new TCanvas; 
    hist_energy->Draw(); 

    c->SaveAs("average_energy_beta-0.pdf"); 

    printf("done.\n"); 

    return 0; 
}
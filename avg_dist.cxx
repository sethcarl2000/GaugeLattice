#include "ComplexNum.hxx"
#include "GaugeLattice.cxx"
#include "SU3.hxx"
#include "numbers.hxx"

#include <TH1D.h> 
#include <TCanvas.h> 

#include <iostream>

int main(int argc, char* argv[])
{   
    using namespace std; 

    //consturct a gauge lattice 
    GaugeLattice<2> lattice(5); 

    auto hist_dist = new TH1D("h_dist", "Avg. frobenius dist between SU(3) members; || g_{1} - g_{2} ||_{F} / 9", 200, 0.5, 1.5); 

    //randomly roatate matrices
    lattice.SetBeta( 0.0 );

    auto old_lattice = lattice; 

    lattice.SetMaxTheta( Nums::pi/2. );

    for (int i=0; i<10000; i++) {

        double accept_prob = lattice.MetropolisUpdate(1e3, 10);

        double avg_norm = lattice.GetFrobDistance(old_lattice); 

        hist_dist->Fill( sqrt(avg_norm) );
        old_lattice = lattice; 

        //printf("i: %-3i, accept. prob: % .5f, avg norm: % .5f\n", i, accept_prob, avg_norm);
    }

    auto c = new TCanvas; 
    hist_dist->Draw(); 

    c->SaveAs("average_dist.pdf"); 

    printf("done.\n"); 

    return 0; 
}
#include "ComplexNum.hxx"
#include "GaugeLattice.cxx"
#include "SU3.hxx"
#include "numbers.hxx"

#include <iostream>

int main(int argc, char* argv[])
{   
    using namespace std; 

    //consturct a gauge lattice 
    GaugeLattice<4> lattice(10); 

    lattice.SetBeta( 0.5 );

    lattice.SetMaxTheta( Nums::pi/2. );

    for (int i=0; i<1000; i++) {

        double accept_prob = lattice.MetropolisUpdate(10000, 10);
        
        printf("i: %-3i, accept. prob: % .5f\n", i, accept_prob);
    }

    printf("done.\n"); 

    return 0; 
}
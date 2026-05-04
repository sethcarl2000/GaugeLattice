#include "ComplexNum.hxx"
#include "GaugeLattice.hxx"
#include "SU3.hxx"

#include <iostream>

namespace {

    constexpr double pi = 3.1415926535; 
}


int main(int argc, char* argv[])
{   
    using namespace std; 

    std::cout << "here" << std::endl;

    SU3::Element g = SU3::Identity(); 

    std::cout << "identitiy element: \n"; 
    SU3::Print( g ); 




    std::cout << "generators: (theta = pi/4)\n";
    std::cout << "--------------------------------------------------------------\n";

    for (int i=0; i<SU3::n_generators; i++) {
        auto ti = SU3::Generator(i, pi/4.); 

        std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~";
        std::cout << "generator t" << i << "\n";
        std::cout << "trace: " << PrintComplex( SU3::Trace(ti) ) << "\n";
        SU3::Print( ti ); 
        std::cout << "adjoint: \n";
        SU3::Print( Adjoint(ti) ); 
        
        //inner product
        std::cout << " ti^dag * ti:\n";
        SU3::Print( Adjoint(ti)*ti ); 
    }

    return 0; 
}
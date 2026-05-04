#ifndef ComplexNum_HXX
#define ComplexNum_HXX

#include <TString.h> 
#include <complex> 
#include <string> 

using ComplexNum = std::complex<double>; 

using namespace std::complex_literals; 

/// @brief reutrns formatted complex number as std::string
/// @param z complex num
/// @param n_decimal_places number of decimal places to print
/// @return string of human-readable complex nubmer, in a + bi form. 
inline std::string PrintComplex(ComplexNum z, int n_decimal_places=3) {
    const char* fmt = Form("%% -4.%if %% -4.%ifi",n_decimal_places,n_decimal_places); 
    return std::string{ Form(fmt, z.real(), z.imag()) };
}

#endif
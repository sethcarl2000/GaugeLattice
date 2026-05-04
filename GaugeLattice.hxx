#ifndef GaugeLattice_HXX
#define GaugeLattice_HXX

#include "ComplexNum.hxx"
#include "SU3.hxx"

#include <vector> 

/// @brief gauge lattice
/// @tparam D dimenion of lattice 
template<int D> class GaugeLattice {
private:

    /// @brief array of all lattice-site connections 
    std::vector<SU3::Element> fConnections[D]; 

    ClassDef(GaugeLattice,1);
}; 

#endif
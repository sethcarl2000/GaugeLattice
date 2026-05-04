#include "GaugeLattice.hxx"
#include "SU3.hxx"
#include <math.h> 

namespace { 

}

//____________________________________________________________________________________________________
template<int D> GaugeLattice<D>::GaugeLattice(int side_length, double beta, double max_theta)
    : 
    fSideLength{side_length},
    fBeta{beta}, 
    fMaxTheta{max_theta}
{   
    //initialize the lattice 
    for (int i=0; i<D; i++) {
        fSideLength[i] = std::vector<SU3::Element>{std::pow(side_length, D), SU3::Identity()}; 
    }
}
//____________________________________________________________________________________________________
template<int D> double GaugeLattice<D>::MetropolisUpdate(long long int n_site_updates, int n_updates_per_site)
{
    long long int n_accepted=0; 
    
    for (long long int i_site=0; i<n_site_updates; i++) {
        
        //pick a site to update
        Index ind = RandIndex(); 

        auto& U = Site(ind);  
        //compute the product of all matrices traced with this one 
    }

    return 0.; 
}
//____________________________________________________________________________________________________
template<int D> double GaugeLattice<D>::GetEnergy() const 
{
    /*noop*/
    return 0.; 
}
//____________________________________________________________________________________________________
template<int D> void GaugeLattice<D>::Print() const 
{
    /*noop*/
}
//____________________________________________________________________________________________________
template<int D> GaugeLattice<D>::Index GaugeLattice<D>::RandIndex()
{
    Index site; 
    for (auto& i : site.pos) i = RandIdx(); 
    site.dir = fRint_direction(fTwister); 
    return site; 
}
//____________________________________________________________________________________________________
//____________________________________________________________________________________________________
//____________________________________________________________________________________________________
//____________________________________________________________________________________________________
//____________________________________________________________________________________________________
//____________________________________________________________________________________________________


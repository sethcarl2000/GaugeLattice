#ifndef GaugeLattice_CXX
#define GaugeLattice_CXX

#include "GaugeLattice.hxx"
#include "numbers.hxx"
#include "SU3.hxx"
#include <math.h> 

//____________________________________________________________________________________________________
template<int D> GaugeLattice<D>::GaugeLattice(int side_length, double beta, double max_theta)
    : 
    fSideLength{side_length},
    fBeta{beta}, 
    fMaxTheta{max_theta}
{   
    //initialize the random-index generator
    std::random_device rd; 
    fTwister = std::mt19937(rd()); 
    fRint_site = std::uniform_int_distribution<int>{0, fSideLength-1};

    //compute how many lattice sites we will have 
    int n_connections_per_direction = (int)std::pow(fSideLength, D); 
    
    //initialize the lattice 
    for (int i=0; i<D; i++) {
        fConnections[i] = std::vector<SU3::Element>(n_connections_per_direction, SU3::Identity()); 
    }
}
//____________________________________________________________________________________________________
template<int D> double GaugeLattice<D>::MetropolisUpdate(long long int n_site_updates, long long int n_updates_per_site)
{
    long long int n_accepted=0; 
    
    for (long long int i_site=0; i_site<n_site_updates; i_site++) {
        
        //pick a site to update
        Index ind_mu = RandIndex(); 
        const int mu=ind_mu.dir; 
        
        auto& U = Site(ind_mu);  
        
        //compute the product of all matrices traced with this one 
        SU3::Element U_trace{{
            0.,0.,0.,
            0.,0.,0.,
            0.,0.,0.
        }}; 

        for (int nu=0; nu<D; nu++) {
            
            //compute the loops with all (orthogonal) directions on the lattice
            if (mu==nu) continue; 

            auto ind_nu = ind_mu; 
            ind_nu.dir = nu; //get new index

            U_trace += (
                Site(next(ind_nu, nu)) *
                SU3::Adjoint(Site(next(ind_mu, mu))) *
                SU3::Adjoint(Site(ind_nu))
            ); 

            U_trace += (
                SU3::Adjoint(Site(next(prev(ind_nu, mu), nu))) *
                SU3::Adjoint(Site(prev(ind_mu, mu))) *
                Site(prev(ind_nu, mu))
            );  
        }

        double tr_old = SU3::Trace( U * U_trace ).real()/3.;

        //now, we rotate the matrix a bit

        for (long long int i=0; i<n_updates_per_site; i++) {

            //update the matrix using a random generator (group element close to the identity)
            auto U_new = ( SU3::Generator(RandGeneratorIdx(), fMaxTheta*(1. - 2.*Rand())) * U ); 

            double tr_new = SU3::Trace( U_new * U_trace ).real()/3.; 

            if ( Rand() < std::exp( fBeta*(tr_new - tr_old) ) ) {

                //printf("old trace: %.3f     new trace: %.3f     exp(-beta*dS) = %.4f\n", tr_old, tr_new, std::exp( fBeta*(tr_new - tr_old) ));

                U = U_new; 
                tr_old = tr_new; 
                ++n_accepted; 
            }
        }
    }
    
    return ((double)n_accepted)/((double)n_site_updates*n_updates_per_site); 
}
//____________________________________________________________________________________________________
template<int D> double GaugeLattice<D>::GetEnergy(int n_measurements)
{
    double ReTr=0.; 
    for (int i=0; i<n_measurements; i++) {

        //pick a random site to measure
        Index ind_mu = RandIndex(); 
        int mu = ind_mu.dir; 

        auto ind_nu = ind_mu; 
        do { ind_nu.dir = fRint_direction(fTwister); } while (ind_nu.dir == ind_mu.dir); 
        int nu = ind_nu.dir; 

        ReTr += SU3::Trace( 
            Site(ind_mu) * 
            Site(next(ind_nu, nu)) * 
            SU3::Adjoint(Site(next(ind_mu, mu))) *
            SU3::Adjoint(Site(ind_nu)) 
        ).real();
    }
    return 1. - ReTr/(3.*((double)n_measurements));
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
template<int D> double GaugeLattice<D>::GetFrobDistance(const GaugeLattice<D>& rhs) const
{
    if (rhs.GetSideLength() != GetSideLength()) {
        throw std::logic_error(Form("<GaugeLattice::GetFrobDistance>: side length of rhs (%i) does not match lhs (%i)",rhs.GetSideLength(),GetSideLength()));
        return Nums::NaN;
    }

    double avg_norm = 0.; 

    const auto& arr_L = GetArray(); 
    const auto& arr_R = rhs.GetArray(); 

    for (int i=0; i<D; i++) {

        const auto& cL = arr_L[i]; 
        const auto& cR = arr_R[i]; 

        for (int j=0; j<cL.size(); j++) avg_norm += SU3::FrobeniusNorm( cL[j] - cR[j] )/9.;
    }

    return avg_norm/std::pow(GetSideLength(), D); 
} 
//____________________________________________________________________________________________________
//____________________________________________________________________________________________________
//____________________________________________________________________________________________________
//____________________________________________________________________________________________________
//____________________________________________________________________________________________________

#endif 
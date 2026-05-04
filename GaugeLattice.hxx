#ifndef GaugeLattice_HXX
#define GaugeLattice_HXX

#include "ComplexNum.hxx"
#include "SU3.hxx"

#include <vector> 
#include <random> 

/// @brief gauge lattice
/// @tparam D dimension of lattice 
template<int D> class GaugeLattice {
private:
    //check to make sure that the lattice is at least D = 2. 
    static_assert((D==2)||(D==3)||(D==4), "Number of dimensions for lattice (template arg.) must be 2, 3, or 4."); 

    /// number of gauge connections per side 
    int fSideLength; 

    /// thermodynamic beta
    double fBeta; 

    /// max theta 
    double fMaxTheta; 

    /// @brief array of all lattice-site connections 
    std::vector<SU3::Element> fConnections[D];

    //random number generators
    std::mt19937 fTwister{std::random_device{}}; 
    std::uniform_int_distribution<int> fRint_generator{0, SU3::n_generators-1};
    std::uniform_int_distribution<int> fRint_direction{0, D-1};

    std::uniform_int_distribution<int> fRint_site{};
    inline int RandIdx() { return fRint_site(fTwister); }

    std::uniform_real_distribution<double> fRand_uniform{0., 1.};
    
    /// @return random number uniformly distributed over [0, 1)
    inline double Rand() { return fRand_uniform(fTwister); }
    
    //indexes a connection
    struct Index{ std::array<int,D> pos; int dir; };

    //check to make sure our struct is trivially constructable and copyable (thanks, claude)
    static_assert(std::is_trivially_copyable_v<Index>, "<GaugeLattice<int D>> Index is not trivially copyable");
    static_assert(std::is_trivially_move_constructible_v<Index>, "<GaugeLattice<int D>> Index is not trivially move constructable");


    /// @return random connection on the lattice 
    GaugeLattice<D>::Index RandIndex(); 
    
    /// @return connection on lattice indexed by 'ind'
    inline SU3::Element& Site(const GaugeLattice<D>::Index& ind) { 
        int idx=0;
        for (auto i : ind.pos) idx += idx*i; 
        return fConnections[ind.dir][idx]; 
    }

public: 

    GaugeLattice(int side_length, double beta=0., double max_theta=0.2*3.1415926535);
    
    /// @brief Perform metropolis update for each lattice site, randomly. 
    /// @param n_site_updates the number of individual sites to update 
    /// @param n_updates_per_site
    /// @return fraction of proposed updates accepted 
    double MetropolisUpdate(long long int n_site_updates, int n_updates_per_site); 

    /// set thermodynamic beta
    void SetBeta(double b) { fBeta=b; }

    /// set maximum generalied angle around which to rotate 
    void SetMaxTheta(double theta) { fMaxTheta=theta; }

    /// get expectation value of the energy 
    double GetEnergy() const; 

    //print current lattice configuration 
    void Print() const; 

    ClassDef(GaugeLattice,1);
}; 

#endif
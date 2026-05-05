#ifndef SU3_HXX
#define SU3_HXX

#include "ComplexNum.hxx"
#include <TString.h> 
#include <math.h> 
#include <stdexcept> 
#include <array> 

//generators specific to SU(3)
namespace SU3 {

    /// @brief element of a group (connects each lattice site) => SU(3)
    struct Element {

        //square matrix dim 
        static constexpr int dim =3; 

        //arrray of complex matrix elements
        std::array<ComplexNum, dim*dim> matrix; 
        
        //element access 
        inline ComplexNum& operator()(int i, int j) { return matrix[i*dim + j]; }

        inline const ComplexNum& get(int i, int j) const { return matrix[i*dim + j]; }

        /// matrix multiplication  
        Element& operator*=(const Element& rhs) noexcept; 

        /// matrix addition  
        Element& operator+=(const Element& rhs) noexcept;

        Element& operator-=(const Element& rhs) noexcept; 

        /// transform this matrix to its adjoint (hermitian conjugate) 
        Element& adjoint() noexcept;
    };
    
    //now, define some out-of-function body operators. this way, we can operate on 'Element' structs which 
    // are not modifiable lvalues

    // matrix addition 
    inline Element operator+(Element lhs, const Element& rhs) noexcept 
    {
        lhs += rhs; 
        return lhs; 
    }
    inline Element operator-(Element lhs, const Element& rhs) noexcept 
    {
        lhs -= rhs; 
        return lhs; 
    }

    // matrix multiplication 
    inline Element operator*(const Element& lhs, const Element& rhs) noexcept 
    {
        Element prod;  
        for (int i=0; i<Element::dim; i++) {
            for (int j=0; j<Element::dim; j++) {
                prod(i,j) 
                    = lhs.get(i,0)*rhs.get(0,j) 
                    + lhs.get(i,1)*rhs.get(1,j) 
                    + lhs.get(i,2)*rhs.get(2,j);
            }
        }
        return prod;
    }   

    /// adjoint (hermitian conjugate) 
    inline Element Adjoint(const Element& g) noexcept { 
        Element adj; 
        for (int i=0; i<Element::dim; i++) {
            for (int j=0; j<Element::dim; j++) {
                adj(i,j) = std::conj( g.get(j,i) ); 
            }
        } 
        return adj; 
    }

    inline double FrobeniusNorm(const Element& g) noexcept {
        double norm =0.; 
        for (const auto& ele : g.matrix) norm += std::norm( ele );
        return norm; 
    }

    
    /// @return Tr(g)
    inline ComplexNum Trace(const Element& g) noexcept { return g.get(0,0) + g.get(1,1) + g.get(2,2); }

    /// @return Identitiy element
    inline Element Identity() noexcept {
        return Element{
            1., 0., 0.,
            0., 1., 0., 
            0., 0., 1.
        };
    };

    // number of generators for this group 
    constexpr int n_generators = 8; 

    /// @brief build a generator, which rotates using generator 'i', with generalized angle '2*theta' 
    /// @param i generator number, indexed [0,7] 
    /// @param theta twice the generalized angle around which to rotate 
    /// @return the exponentiaed generator, with generalized angle 2*theta (2*sqrt(3)*theta for t_7) 
    //___________________________________________________________________________________________________________
    Element Generator(int i, double theta);
    //___________________________________________________________________________________________________________

    // print elements of the group member 'g' in 3x3 mat. rep
    void Print(const Element& g, int n_decimal_places=3); 
}

#endif 
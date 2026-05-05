
#include <TString.h> 

#include "SU3.hxx"
#include <stdio.h> 
#include <string> 
#include <iostream> 


namespace SU3 
{

//_____________________________________________________________________________________________________________________
Element& Element::operator*=(const Element& rhs) noexcept 
{ 
    *this = (*this) * rhs; 
    return *this; 
}
//_____________________________________________________________________________________________________________________
Element& Element::operator+=(const Element& rhs) noexcept 
{ 
    for (int i=0; i<Element::dim; i++) {
        for (int j=0; j<Element::dim; j++) {
            (*this)(i,j) += rhs.get(i,j); 
        }
    }
    return *this;
}
//_____________________________________________________________________________________________________________________
Element& Element::operator-=(const Element& rhs) noexcept 
{ 
    for (int i=0; i<Element::dim; i++) {
        for (int j=0; j<Element::dim; j++) {
            (*this)(i,j) -= rhs.get(i,j); 
        }
    }
    return *this;
}
//_____________________________________________________________________________________________________________________
Element& Element::adjoint() noexcept 
{
    ComplexNum temp;
    for (int i=0; i<Element::dim; i++) {

        //conjugate the diagonal elements
        auto& ele_ii = (*this)(i,i); 
        ele_ii = std::conj( ele_ii );
        
        //conjugate and flip the off-diagonal elements
        for (int j=0; j<i; j++) {
            auto& ele_ij = temp = (*this)(i,j); 
            auto& ele_ji = (*this)(j,i); 
            ele_ij = std::conj(ele_ji);
            ele_ji = std::conj(temp); 
        }
    }
    return *this; 
}
//_____________________________________________________________________________________________________________________
Element Generator(int i, double theta)
{
    double cos = std::cos(theta); 
    double sin = std::sin(theta); 

    /// all of these operators are computed by exponentiating the Gell-Mann representation of the 8 
    /// SU(3) group generators  
    switch (i) {

        case 0 : return Element{
            cos,     1i*sin,  0.,
            1i*sin,  cos,     0.,
            0.,      0.,      1.
        }; break; 

        case 1 : return Element{
            cos,    sin,    0.,
            -sin,   cos,    0.,
            0.,     0.,     1.
        }; break; 

        case 2 : return Element{
            cos + 1i*sin,   0.,             0.,
            0.,             cos - 1i*sin,   0.,
            0.,             0.,             1.
        }; break; 

        case 3 : return Element{
            cos,    0.,     1i*sin,
            0.,     1.,     0.,
            1i*sin, 0.,     cos
        }; break; 

        case 4 : return Element{
            cos,    0.,     sin,
            0.,     1.,     0.,
            -sin,   0.,     cos
        }; break; 

        case 5 : return Element{
            1.,     0.,     0.,            
            0.,     cos,    1i*sin,
            0.,     1i*sin, cos
        }; break; 
        
        case 6 : return Element{
            1.,     0.,     0.,            
            0.,     cos,    sin,
            0.,     -sin,   cos
        }; break; 

        case 7 : return Element{
            cos + 1i*sin,   0.,             0.,
            0.,             cos + 1i*sin,   0.,
            0.,             0.,             (cos - 1i*sin)*(cos - 1i*sin)
        }; break; 

        default :  {
            throw std::logic_error(Form("<SU3::Generator>: invalid group generator index (%i). valid range is [0,7].",i)); 
            return SU3::Identity(); 
        } 
    }
}
//_____________________________________________________________________________________________________________________
void Print(const Element& g, int n_decimal_places) 
{
    int ii=0; 
    std::cout << 
        PrintComplex(g.matrix[ii++], n_decimal_places) << "   " << PrintComplex(g.matrix[ii++], n_decimal_places) << "   " << PrintComplex(g.matrix[ii++], n_decimal_places) << "\n" <<
        PrintComplex(g.matrix[ii++], n_decimal_places) << "   " << PrintComplex(g.matrix[ii++], n_decimal_places) << "   " << PrintComplex(g.matrix[ii++], n_decimal_places) << "\n" <<
        PrintComplex(g.matrix[ii++], n_decimal_places) << "   " << PrintComplex(g.matrix[ii++], n_decimal_places) << "   " << PrintComplex(g.matrix[ii++], n_decimal_places) << "\n"; 
    
    return; 
}
//_____________________________________________________________________________________________________________________
//_____________________________________________________________________________________________________________________
//_____________________________________________________________________________________________________________________
//_____________________________________________________________________________________________________________________
//_____________________________________________________________________________________________________________________
//_____________________________________________________________________________________________________________________

//check to make sure our struct is trivially constructable and copyable (thanks, claude)
static_assert(std::is_trivially_copyable_v<Element>);
static_assert(std::is_trivially_move_constructible_v<Element>);


};
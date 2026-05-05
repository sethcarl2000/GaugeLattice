#ifndef numbers_HXX
#define numbers_HXX

#include <limits> 

namespace Nums {
    
    constexpr double pi = 3.14159265359; 

    constexpr double NaN = std::numeric_limits<double>::quiet_NaN(); 
}

#endif
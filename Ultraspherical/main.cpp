//
//  main.cpp
//  Ultraspherical
//
//  Created by Sheehan Olver on 04/06/2012.
//  Copyright (c) 2012 School of Mathematics and Statistics, The University of Sydney. All rights reserved.
//

#include <iostream>
#include <boost/numeric/ublas/banded.hpp>
#include <boost/numeric/ublas/io.hpp>

int main(int argc, const char * argv[])
{

    using namespace boost::numeric::ublas;
//    banded_matrix<double> m (6, 6, 1, 1);
//    for (signed i = 0; i < signed (m.size1 ()); ++ i)
//        for (signed j = std::max (i - 1, 0); j < std::min (i + 2, signed (m.size2 ())); ++ j)
//            m (i, j) = 3 * i + j;
//    std::cout << m << std::endl;
//    std::cout << m.lower() << std::endl;
//    std::cout << m.upper() << std::endl;
    
    banded_matrix<double> m (6, 6);
    
    m.resize(6,6);
    m(3,3) = 3;
    
    std::cout << m(3,3) << std::endl;

    return 0;
}


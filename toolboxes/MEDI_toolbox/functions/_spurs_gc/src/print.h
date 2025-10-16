#ifndef SPURS_GE_SRC_PRINT_H_
#define SPURS_GE_SRC_PRINT_H_
#include <iostream>

#define CONCAT(x,y) x ## y
#define CONCAT3(x,y,z) CONCAT(CONCAT(x,y),z)
#define PRINT_stringize(y) #y
#define PRINTVAR(a) PRINT_stringize(a) << " = " << (a) 
#define PRINT(a) std::cout << __FILE__ << "(" << __LINE__ << "): " << PRINTVAR(a) << std::endl;

#endif //SPURS_GE_SRC_PRINT_H_

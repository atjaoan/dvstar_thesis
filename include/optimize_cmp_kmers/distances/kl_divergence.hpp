#pragma once 

#include "vlmc_container.hpp"

namespace distance {
  
using vlmc_c = container::VLMC_Container;

double kl(vlmc_c left, vlmc_c right){
  return 0.0;  
}
}
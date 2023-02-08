#pragma once 

#include "vlmc_template.hpp"

namespace distance {
  
using vlmc_c = container::VLMC_template;

float kl(vlmc_c left, vlmc_c right){
  return 0.0;  
}
}
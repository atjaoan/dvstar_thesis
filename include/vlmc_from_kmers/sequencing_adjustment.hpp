#pragma once

namespace vlmc {

struct SequencingParameters {
  bool adjust_for_sequencing_errors = false;
  double depth = 1;
  double error_rate = 0;
};

double adjusted_value(SequencingParameters sequencing_parameters, uint length) {
  if (sequencing_parameters.adjust_for_sequencing_errors) {
    return sequencing_parameters.depth *
           std::pow(1 - sequencing_parameters.error_rate, double(length));
  } else {
    return 1.0;
  }
}

} // namespace vlmc
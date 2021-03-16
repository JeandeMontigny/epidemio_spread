// -----------------------------------------------------------------------------
//
// Copyright (C) Jean de Montigny.
// All Rights Reserved.
//
// -----------------------------------------------------------------------------

#ifndef SIM_PARAM_H_
#define SIM_PARAM_H_

#include "core/param/param_group.h"

namespace bdm {

/// This class defines parameters that are specific to this simulation.
struct SimParam : public ParamGroup {

  BDM_PARAM_GROUP_HEADER(SimParam, 1);

  uint64_t number_of_steps = 3600;
  double human_diameter = 25; // cm
  int map_pixel_size = 4;
};

}  // namespace bdm

#endif  // SIM_PARAM_H_

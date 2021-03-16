// -----------------------------------------------------------------------------
//
// Copyright (C) The BioDynaMo Project, Liudmila Rozanova, Alexander Temerev
// and Jean de Montigny
// All Rights Reserved.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
//
// See the LICENSE file distributed with this work for details.
// See the NOTICE file distributed with this work for additional information
// regarding copyright ownership.
//
// -----------------------------------------------------------------------------

#ifndef OPENLB_SIM_H_
#define OPENLB_SIM_H_

#include "biodynamo.h"
#include "openlb_export.h"

namespace openlb_sim {

  int main( int argc, char* argv[], std::vector<bdm::dummy_agent> bdm_agents, int init );

} // namespace openlb_sim

#endif // OPENLB_SIM_H_

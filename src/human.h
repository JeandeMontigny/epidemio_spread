// -----------------------------------------------------------------------------
//
// Copyright (C) Jean de Montigny.
// All Rights Reserved.
//
// -----------------------------------------------------------------------------

#ifndef HUMAN_H_
#define HUMAN_H_

#include "core/agent/cell.h"
//#include "core/behavior/behavior.h"

namespace bdm {

enum State { kHealthy, kIncubation, kInfected, kRecovered };

class Human : public Cell {
  BDM_AGENT_HEADER(Human, Cell, 1);

 public:
  Human() {}

  explicit Human(const Double3& position) : Base(position) {}
  virtual ~Human() {}

  // This data member stores the current state of the agent.
  int state_ = State::kHealthy;
  int incubation_counter_ = 1e9;
  int recovery_counter_ = 1e9;
  // store the destinations
  std::vector<std::pair<double, double>> destinations_list_;
  // store the path to a destination
  std::vector<std::vector<double>> path_;
  /// store the agent's orientation. by default, follow bus orientation
  std::vector<double> orientation_ = {-1,0};
  double viral_load_ = 0;
};

}  // namespace bdm

#endif // HUMAN_H_

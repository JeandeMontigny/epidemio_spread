// -----------------------------------------------------------------------------
//
// Copyright (C) Jean de Montigny.
// All Rights Reserved.
//
// -----------------------------------------------------------------------------

#ifndef BEHAVIOR_H_
#define BEHAVIOR_H_

#include "core/behavior/behavior.h"
#include "human.h"
#include "geom_construct.h"
#include "sim-param.h"
#include "a_star.h"
#include "navigation_util.h"
#include "core/scheduler.h"
//#include "core/diffusion_grid.h"
#include "core/resource_manager.h"
#include "population_creation_utils.h"

namespace bdm {

// ---------------------------------------------------------------------------
struct Navigation : public Behavior {
  BDM_BEHAVIOR_HEADER(Navigation, Behavior, 1);

  Navigation() { AlwaysCopyToNew(); }

  Navigation(std::vector<std::vector<bool>>* navigation_map) {
    navigation_map_ = navigation_map;
  }

  void Run(Agent* so) override {
    auto* sim = Simulation::GetActive();
    auto* random = sim->GetRandom();

    auto* human = bdm_static_cast<Human*>(so);
    const auto& position = human->GetPosition();
    std::vector<double> previous_position = {position[0], position[1]};

    // if human is at the bus exit
    if (IsExit(position)) {
      human->RemoveFromSimulation();
    }

    // some passengers (x%) leave if sim step == y
    // add destination to bus exit
    if ((sim->GetScheduler()->GetSimulatedSteps() % 800 == 0
          || sim->GetScheduler()->GetSimulatedSteps() % 1800 == 0
          || sim->GetScheduler()->GetSimulatedSteps() % 3000 == 0)
        && random->Uniform() < 0.2) {
      human->destinations_list_ =
        AddDestinationToList(human->destinations_list_, navigation_map_,
                             200, 200, 110, 110);
    }

    std::vector<std::vector<double>> path;

    // if agent has to calculate path to destination
    if (!path_calculated_ && !human->destinations_list_.empty()) {
      std::pair<double, double> start =
        std::make_pair(GetBDMToMapLoc(position[0]),
                       GetBDMToMapLoc(position[1]));
      std::pair<double, double> dest = human->destinations_list_[0];

      // calculate path using A*
      path = AStar((*navigation_map_), start, dest, navigation_map_->size());

      human->path_ = path;
      // remove this travel form destination_list
      human->destinations_list_.erase(human->destinations_list_.begin());
      path_calculated_ = true;
    } // end if has to calculate path

    // if agent has its path, has to move to it
    else if (path_calculated_) {
      Double3 dest = {GetMapToBDMLoc(human->path_[0][0]),
                      GetMapToBDMLoc(human->path_[0][1]), -40};
      // if human is close to destination, check if seat is still empty
      if (human->path_.size() < 25 && !IsExit(dest) &&
        SeatTaken(dest, position)) {
        dest = GetEmptySeat();
        std::vector<std::pair<double, double>> destinations_list;
        destinations_list.push_back(std::make_pair(
          GetBDMToMapLoc(dest[0]), GetBDMToMapLoc(dest[1])));
        human->destinations_list_= destinations_list;
        path_calculated_ = false;
        return;
      }

      if (!human->path_.empty()) {
        Double3 next_position = {
          GetMapToBDMLoc(human->path_[human->path_.size()-1][0]),
          GetMapToBDMLoc(human->path_[human->path_.size()-1][1]),
          position[2] };
        // navigate according to path
        human->SetPosition(next_position);

        double vec_x = human->GetPosition()[0] - previous_position[0];
        double vec_y = human->GetPosition()[1] - previous_position[1];
        double vec_norm = std::sqrt((vec_x*vec_x)+(vec_y*vec_y));
        human->orientation_ = {vec_x/vec_norm, vec_y/vec_norm};

        // if either on or after this path position
        human->path_.erase(human->path_.end());
      }
      // path is empty, so destination is reached
      else {
        // agent orientation follow bus orientation
        human->orientation_ = {-1, 0};
        // can add an other destination here
        path_calculated_ = false;
      }
    } // end has its path

  } // end Run

private:
  bool path_calculated_ = false;
  std::vector<std::vector<bool>>* navigation_map_;
}; // end Navigation

// ---------------------------------------------------------------------------
struct GetInfectedBehaviour : public Behavior {
  BDM_BEHAVIOR_HEADER(GetInfectedBehaviour, Behavior, 1);

  GetInfectedBehaviour() { AlwaysCopyToNew(); }

  void Run(Agent* so) override {
    // auto* sim = Simulation::GetActive();
    // auto* rm = sim->GetResourceManager();

    auto* human = bdm_static_cast<Human*>(so);

    // nothing to do if already infected, so remove BM
    if (human->state_ == State::kInfected) {
      human->RemoveBehavior(this);
      return;
    }

    // infection
    if (human->state_ == State::kHealthy) {
      double concentration = 0;

      //TODO: get droplet and aerosol concentrations from OpenLB

      if (concentration > 1e-6) {
        human->state_ = State::kIncubation;
      }
    } // end infection

  } // end Run
}; // end GetInfectedBehaviour

}  // namespace bdm

#endif  // BEHAVIOR_H_

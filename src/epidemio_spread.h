// ------------------------------------------------------------------------
//
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
#ifndef EPIDEMIO_SPREAD_H_
#define EPIDEMIO_SPREAD_H_

#include "biodynamo.h"
#include "core/simulation.h"
#include "human.h"
#include "behavior.h"
#include "sim-param.h"
#include "geom_construct.h"
#include "geom.h"
#include "util_methods.h"
#include "navigation_util.h"
#include "a_star.h"
#include "population_creation.h"
#include "openlb_export.h"
#include "openlb_sim.h"
#include <mpi.h>

namespace bdm {

inline int Simulate(int argc, const char** argv) {

  Param::RegisterParamGroup(new SimParam());

  Simulation simulation(argc, argv);

  auto* param = simulation.GetParam();
  auto* scheduler = simulation.GetScheduler();
  auto* sparam = param->Get<SimParam>();

  scheduler->UnscheduleOp(scheduler->GetOps("mechanical forces")[0]);

  simulation.GetRandom()->SetSeed(2975); // rand() % 10000

  // create OpenLB directories and files
  std::string openlbDir = Concat(param->output_dir, "/openlb/");  
  std::vector<std::vector<bool>> navigation_map;
  // contains bdm agents informations for openLB 
  std::vector<dummy_agent> bdm_agents;

  // Initialise MPI environment
  openlb_sim::main(0, nullptr, bdm_agents, 1);
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  // if process is the master
  if (world_rank == 0) {

    // choose if print navigation map
    bool print_navigation_map = false;
    
    //construct geom
    BuildBus();

    // construct the 2d array for navigation
    navigation_map = GetNavigationMap(-40);

    // print navigation_map
    if (print_navigation_map) {
      for (size_t map_x = 0; map_x < navigation_map.size(); map_x++) {
	for (size_t map_y = 0; map_y < navigation_map[0].size(); map_y++) {
	  std::cout << navigation_map[map_x][map_y];
	}
	std::cout << std::endl;
      }
    }

    //bus population creation
    InitialBusPopulationCreation(&navigation_map);

    // Run simulation for number_of_steps timestep
    std::cout << "simulating.." << std::endl;
    ExportOpenlbFiles(openlbDir, "stl");

  } // end if rank 0

  for (uint64_t i = 0; i <= sparam->number_of_steps; ++i) {

    // if process is the master
    if (world_rank == 0) {

      // passenger pop in at bus entrace position
      // first bus stop
      if (i == 10) {
	std::cout << "first bus strop" << std::endl;
	AddPassenger(State::kHealthy, &navigation_map);
      }
      if (i == 20) {
	AddPassenger(State::kHealthy, &navigation_map);
      }
      if (i == 30) {
	AddPassenger(State::kHealthy, &navigation_map);
      }
      if (i == 40) {
	AddPassenger(State::kHealthy, &navigation_map);
      }
      if (i == 50) {
	AddPassenger(State::kHealthy, &navigation_map);
      }
      if (i == 60) {
	AddPassenger(State::kHealthy, &navigation_map);
      }
      if (i == 70) {
	AddPassenger(State::kHealthy, &navigation_map);
      }
      // second bus stop
      if (i == 810) {
	std::cout << "second bus strop" << std::endl;
	AddPassenger(State::kInfected, &navigation_map);
      }
      if (i == 820) {
	AddPassenger(State::kHealthy, &navigation_map);
      }
      if (i == 830) {
	AddPassenger(State::kHealthy, &navigation_map);
      }
      // third bus stop
      if (i == 1810) {
	std::cout << "third bus strop" << std::endl;
	AddPassenger(State::kHealthy, &navigation_map);
      }
      if (i == 1820) {
	AddPassenger(State::kInfected, &navigation_map);
      }
      if (i == 1830) {
	AddPassenger(State::kHealthy, &navigation_map);
      }
      if (i == 1840) {
	AddPassenger(State::kHealthy, &navigation_map);
      }
      if (i == 1850) {
	AddPassenger(State::kHealthy, &navigation_map);
      }

      //run bdm simulation
      scheduler->Simulate(1);

      // if have to run an openLB simulation
      if (i % 500 == 0) {
	std::cout << "step " << i << " out of " << sparam->number_of_steps << std::endl;
	//write stl file, and export agents data to file here
	ExportOpenlbFiles(openlbDir, "stl");
      } // end if step 500 %
    } // end if rank 0

    // if step 500 %, for all processes
    if (i % 500 == 0) {
      // all processes wait here until the master is done with bdm simulation
      MPI_Comm comm;
      comm = MPI_COMM_WORLD;
      MPI_Barrier(comm);
      // get agents info for openLB sim
      bdm_agents = ReadBdmStateFile();
      //run openLB simulation
      openlb_sim::main(0, nullptr, bdm_agents, 0);

      // if master process
      if ( world_rank == 0 ) {
	printf("openLB simulation done\n");

	// read bdm_sim_state.txt file to update agents' viral load
	std::vector<double> viral_load_list;
	std::string line;
	int line_count = 0;

	std::ifstream simu_file;
	// TODO: no hard coded address
	simu_file.open("./output/openlb/bdm_sim_state.txt");

	while (getline(simu_file, line)) {
	  // if this line is not the viral load data
	  if (line_count != 8) {
	    line_count++;
	  }
	  // agent viral load
	  if (line_count == 8) {
	    line_count++;
	    viral_load_list.push_back(std::stod(line));
	  }
	  // if new agent
	  if (line == "{") {
	    line_count = 1;
	  }
	} // end while line
	simu_file.close();

	// update agents' viral load
	int agent_count = 0;
	auto update_agents = [&viral_load_list, &agent_count](bdm::Agent* so) {
	  auto* hu = bdm::bdm_static_cast<bdm::Human*>(so);
	  hu->viral_load_ = viral_load_list[agent_count];
	  agent_count++;
	};
	simulation.GetResourceManager()->ForEachAgent(update_agents);

      } // end if rank 0
    } // end if % 500, for all processes

  } // end for number_of_steps

  std::cout << "done" << std::endl;

  MPI_Finalize();

  return 0;
}

}  // namespace bdm

#endif  // EPIDEMIO_SPREAD_H_

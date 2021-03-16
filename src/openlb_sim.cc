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
#include "openlb_sim.h"
#include "olb3D.h"
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <random>
#include "biodynamo.h"
#include "openlb_export.h"
#include "sim-param.h"

namespace openlb_sim {

using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;
using namespace olb::util;

typedef double T;

// Choose your turbulent model of choice
//#define RLB
#define Smagorinsky //default
//#define ConsitentStrainSmagorinsky
//#define ShearSmagorinsky
//#define Krause

#ifdef ShearSmagorinsky
#define DESCRIPTOR ShearSmagorinskyD3Q19Descriptor
#else
#define DESCRIPTOR D3Q19<>
#endif

// ---------------------------------------------------------------------------
// Parameters for the simulation setup

// resolution of the model: number of voxel per physLength
const int N = 125;
// number of time steps per seconds 
// 300 for 1 m/s ; 4000 for 4 m/s ; 4000 for 10 m/s without seg fault, but with lost particles ; 40000++ for 10 m/s
const int M = 300;
// block profile (mode=0), power profile (mode=1)
const int inflowProfileMode = 0;
// max. simulation time
const T maxPhysT = 5;
// export visual every outputSteps steps
int outputSteps = maxPhysT*M/10; // 100

// ---------------------------------------------------------------------------
template <typename T, typename _DESCRIPTOR>
class TurbulentVelocity3D : public AnalyticalF3D<T,T> {

protected:
  // block profile (mode=0), power profile (mode=1)
  int _mode;
  T rho;
  T nu;
  T u0;
  T p0;
  T charL;
  T dx;

public:
  TurbulentVelocity3D( UnitConverter<T,_DESCRIPTOR> const& converter, int mode=0 ) : AnalyticalF3D<T,T>( 3 ) {
    _mode = mode;
    u0 = converter.getCharLatticeVelocity();
    rho = converter.getPhysDensity();
    nu = converter.getPhysViscosity();
    charL = converter.getCharPhysLength();
    p0 = converter.getCharPhysPressure();
    dx = converter.getConversionFactorLength();

    this->getName() = "turbulentVelocity3d";
  };

  bool operator()( T output[], const BaseType<T> input[] ) override
  {
    T y = input[1];
    T z = input[2];
    // block profile inititalization
    T u_calc = u0;
    // power profile inititalization
    if ( _mode==1 ) {
      T obst_y = 5.5+dx;
      T obst_z = 5.5+dx;
      T obst_r = 0.5;

      T B      = 5.5;
      T kappa  = 0.4;
      T ReTau  = 183.6;

      u_calc = u0/7.*( 2.*nu*ReTau/( charL*kappa )*log( fabs( 2.*ReTau/charL*( obst_r - sqrt( pow( y - obst_y, 2. )
                       + pow( z - obst_z, 2. ) ) )*1.5*( 1 + sqrt( pow( y - obst_y, 2. )
                           + pow( z - obst_z, 2. ) )/obst_r )/( 1 + 2.*pow( sqrt( pow( y - obst_y, 2. )
                               + pow( z - obst_z, 2. ) )/obst_r, 2. ) ) ) + B ) );
    }
    T a = -1., b = 1.;
    T nRandom = rand()/( T )RAND_MAX*( b-a ) + a;

    output[0] = u_calc+  0.15*u0*nRandom;
    output[1] = 0.15*u0*nRandom;
    output[2] = 0.15*u0*nRandom;
    return true;
  };
};

// ---------------------------------------------------------------------------
void prepareGeometry( UnitConverter<T,DESCRIPTOR> const& converter,
                      IndicatorF3D<T>& indicator, STLreader<T>& stlReader,
                      SuperGeometry3D<T>& superGeometry, std::vector<bdm::dummy_agent> bdm_agents) {

  OstreamManager clout( std::cout,"prepareGeometry" );
  clout << "Prepare Geometry ..." << std::endl;
  // Sets material number for fluid and boundary
  superGeometry.rename( 0, 2, indicator );
  superGeometry.rename( 0, 1 );
  superGeometry.rename( 2, 0, stlReader );

  superGeometry.clean();
  // -------- biodynamo -------- //
  for (size_t i = 0; i < bdm_agents.size(); i++ ) {
    bdm::dummy_agent agent = bdm_agents[i];
    // NOTE: spread only for infected agents
    if (agent.state == bdm::State::kInfected) {
      bdm::Double3 pos = agent.pos;
      auto dir = agent.dir;
      // convert radius to m
      double radius = agent.radius / 100;

      // -------- create inlets -------- //
      Vector<T,3> spread_pos_in(
        (pos[0]/100 + dir[0]*radius)*converter.getCharPhysLength(),
        (pos[1]/100 + dir[1]*radius)*converter.getCharPhysLength(),
        pos[2]/100 );
      Vector<T,3> spread_pos_out(
        (pos[0]/100 + dir[0]*(radius+0.01))*converter.getCharPhysLength(),
        (pos[1]/100 + dir[1]*(radius+0.01))*converter.getCharPhysLength(),
        pos[2]/100 );

      IndicatorCylinder3D<T> layerInflow( spread_pos_in, spread_pos_out,
        0.04*converter.getCharPhysLength());
      superGeometry.rename( 1, 3, layerInflow );
    } // if agent is infected
  } // for each agent in sim

  // Removes all not needed boundary voxels outside the surface
  superGeometry.clean();
  // Removes all not needed boundary voxels inside the surface
  superGeometry.innerClean();
  superGeometry.checkForErrors();

  clout << "Prepare Geometry ... OK" << std::endl;
}

// ---------------------------------------------------------------------------
void prepareLattice( SuperLattice3D<T,DESCRIPTOR>& sLattice,
                     UnitConverter<T,DESCRIPTOR> const& converter,
                     Dynamics<T, DESCRIPTOR>& bulkDynamics,
                     SuperGeometry3D<T>& superGeometry ) {

  OstreamManager clout( std::cout,"prepareLattice" );
  clout << "Prepare Lattice ..." << std::endl;

  const T omega = converter.getLatticeRelaxationFrequency();

  // Material=0 -->do nothing
  sLattice.defineDynamics(
    superGeometry, 0, &instances::getNoDynamics<T, DESCRIPTOR>() );

  // Material=1 -->bulk dynamics
  sLattice.defineDynamics( superGeometry, 1, &bulkDynamics );

  // Material=2 -->bounce back
  sLattice.defineDynamics(
    superGeometry, 2, &instances::getBounceBack<T, DESCRIPTOR>() );

  // Material=3 -->bulk dynamics (inflow)
  setInterpolatedVelocityBoundary<T,DESCRIPTOR>(sLattice, omega, superGeometry, 3);

  // Material=4 -->bulk dynamics (outflow)

  clout << "Prepare Lattice ... OK" << std::endl;
}

// ---------------------------------------------------------------------------
void prepareParticles(SuperLattice3D<T,DESCRIPTOR>& sLattice,
                  UnitConverter<T,DESCRIPTOR> const& converter,
                  SuperGeometry3D<T>& superGeometry,
                  SuperParticleSystem3D<T, Particle3D>& superParticleSystem,
                  SuperLatticeInterpPhysVelocity3D<T, DESCRIPTOR>& velocityProbe,
		  size_t p_nb, bool cough, std::vector<bdm::dummy_agent> bdm_agents) {

  OstreamManager clout( std::cout,"prepareParticules" );

  clout << "Prepare Particules ..." << std::endl;

  // -------- biodynamo -------- //
  for (size_t i = 0; i < bdm_agents.size(); i++ ) {
    bdm::dummy_agent agent = bdm_agents[i];
    // NOTE: spread only for infected agents
    if (agent.state == bdm::State::kInfected) {
      bdm::Double3 pos = agent.pos;
      auto dir = agent.dir;
      // convert radius to m
      double radius = agent.radius / 100;

      // create inlets
      Vector<T,3> spread_pos_in(
        (pos[0]/100 + dir[0]*(radius+0.01))*converter.getCharPhysLength(),
        (pos[1]/100 + dir[1]*(radius+0.01))*converter.getCharPhysLength(),
        pos[2]/100 );
      Vector<T,3> spread_pos_out(
        (pos[0]/100 + dir[0]*(radius+0.03))*converter.getCharPhysLength(),
        (pos[1]/100 + dir[1]*(radius+0.03))*converter.getCharPhysLength(),
        pos[2]/100 );

      IndicatorCylinder3D<T> inletCylinder(spread_pos_in, spread_pos_out,
        0.03*converter.getCharPhysLength());

      // create particles
      double p_density = 1000.; // kg/m^3
      double p_radius; // m
      double p_mass; // kg

      // gamma distribution generator. should be from 10e-6 to 1000e-6, ave ~180
      std::default_random_engine generator;
      std::gamma_distribution<double> distribution_c(1.4, 2.6);
      std::gamma_distribution<double> distribution_b(3.5, 2);

      for (size_t p = 0; p < p_nb; p++){
	if (cough) {
	  p_radius = distribution_c(generator)*1e-4;
	}
	else {
	  p_radius = distribution_b(generator)*1e-7;
	}
        p_mass = 4. / 3. * M_PI * std::pow(p_radius, 3) * p_density;
	// add one particle with p_mass and p_radius
	superParticleSystem.addParticle(inletCylinder, p_mass, p_radius, 1);
      } // end for p_nb of particles

      // add forces to particles
      auto stokesDragForce = make_shared
              < StokesDragForce3D<T, Particle3D, DESCRIPTOR>
              > ( velocityProbe, converter );
      superParticleSystem.addForce(stokesDragForce);

      superParticleSystem.setVelToFluidVel(velocityProbe);

      std::set<int> boundMaterial = {2};
      auto materialBoundary = make_shared
              < MaterialBoundary3D<T, Particle3D>
              > (superGeometry, boundMaterial);
      superParticleSystem.addBoundary(materialBoundary);

    } // if agent is infected
  } // for each agent in sim

  clout << "Prepare Particules ... OK" << std::endl;
}

// ---------------------------------------------------------------------------
void setBoundaryValues( UnitConverter<T,DESCRIPTOR> const&converter,
                        SuperLattice3D<T,DESCRIPTOR>& lattice,
                        SuperGeometry3D<T>& superGeometry, int iT ) {

  OstreamManager clout( std::cout,"setBoundaryValues" );

  if ( iT==0 ) {
    AnalyticalConst3D<T,T> rhoF( 1 );
    std::vector<T> velocity( 3,T() );
    AnalyticalConst3D<T,T> uF( velocity );

    // Seeding of random fluctuations and definition of the velocity field
    srand( time( nullptr ) );
    TurbulentVelocity3D<T,DESCRIPTOR> uSol( converter, inflowProfileMode );

    lattice.iniEquilibrium(
      // superGeometry.getMaterialIndicator({1, 2, 4}), rhoF, uF );
      superGeometry.getMaterialIndicator({1, 2}), rhoF, uF );
    lattice.iniEquilibrium( superGeometry, 3, rhoF, uSol );

    lattice.defineU( superGeometry, 3, uSol );
    // lattice.defineRho( superGeometry, 4, rhoF );

    // Make the lattice ready for simulation
    lattice.initialize();
  }
}

// ---------------------------------------------------------------------------
void getResults( SuperLattice3D<T, DESCRIPTOR>& sLattice,
                 UnitConverter<T,DESCRIPTOR> const& converter, int iT,
                 SuperGeometry3D<T>& superGeometry, Timer<T>& timer ) {

  OstreamManager clout( std::cout,"getResults" );
  SuperVTMwriter3D<T> vtmWriter( "bus_spreading" );

  if ( iT == 0 ) {
    // Writes the geometry as vti file for visualization
    SuperLatticeGeometry3D<T, DESCRIPTOR> geometry( sLattice, superGeometry );
    vtmWriter.write( geometry );
    vtmWriter.createMasterFile();
  }
}

// ---------------------------------------------------------------------------
  int main( int argc, char* argv[], std::vector<bdm::dummy_agent> bdm_agents, int init ) {

  // === 1st Step: Initialization ===
  if (init) {
    olbInit( &argc, &argv );
    return 0;
  }

  if (bdm_agents.size() == 0) {
    return 0;
  }

  std::cout << "starts olb simulation.." << std::endl;

  auto* sim = bdm::Simulation::GetActive();
  auto* random = sim->GetRandom();

  //TODO: add openlb sim number in output dir?
  singleton::directories().setOutputDir( "./output/openlb/" );
  OstreamManager clout( std::cout,"main" );
  // display messages from every single mpi process
  clout.setMultiOutput(true);

  // decide if this will be a cough or nornal breathing
  bool cough = 0;
  T max_velo = 1;

  if (random->Uniform(0, 1) < 0) { // 0.1
    cough = 1;
    max_velo = 10;
  }

  UnitConverter<T, DESCRIPTOR> const converter(
    (T) 1/N,      // physDeltaX: spacing between two lattice cells in m
    (T) 1/M,      // physDeltaT: time step in s
    (T) 1,        // charPhysLength: reference/characteristic length of simulation geometry in m
    (T) max_velo, // charPhysVelocity: maximal or highest expected velocity during simulation in m / s
    (T) 0.0002,   // physViscosity: physical kinematic viscosity in m^2 / s
    (T) 1.0       // physDensity: physical density in kg / m^3
  );

  // Prints the converter log as console output
  // Writes the converter log in a file
  converter.write("bus_spreading");

  // === 2nd Step: Prepare Geometry ===
  STLreader<T> stlReader( "./output/openlb/bus.stl",
    converter.getConversionFactorLength() );
  IndicatorLayer3D<T> extendedDomain( stlReader,
    converter.getConversionFactorLength() );
  // NOTE: this step takes about 5 minutes here.
  // Instantiation of a cuboidGeometry with weights
  CuboidGeometry3D<T> cuboidGeometry( extendedDomain,
    converter.getConversionFactorLength(), singleton::mpi().getSize() );
  // Instantiation of a loadBalancer
  HeuristicLoadBalancer<T> loadBalancer( cuboidGeometry );
  // Instantiation of a superGeometry
  SuperGeometry3D<T> superGeometry( cuboidGeometry, loadBalancer, 2 );

  prepareGeometry(converter, extendedDomain, stlReader, superGeometry, bdm_agents);

  // === 3rd Step: Prepare Lattice ===
  SuperLattice3D<T, DESCRIPTOR> sLattice( superGeometry );

  Dynamics<T, DESCRIPTOR>* bulkDynamics;
  const T omega = converter.getLatticeRelaxationFrequency();
#if defined(RLB)
  bulkDynamics = new RLBdynamics<T, DESCRIPTOR>(
    omega, instances::getBulkMomenta<T, DESCRIPTOR>() );
#elif defined(Smagorinsky)
  bulkDynamics = new SmagorinskyBGKdynamics<T, DESCRIPTOR>(
    omega, instances::getBulkMomenta<T, DESCRIPTOR>(),
      0.15);
#elif defined(ShearSmagorinsky)
  bulkDynamics = new ShearSmagorinskyBGKdynamics<T, DESCRIPTOR>(
    omega, instances::getBulkMomenta<T, DESCRIPTOR>(),
      0.15);
#elif defined(Krause)
  bulkDynamics = new KrauseBGKdynamics<T, DESCRIPTOR>(
    omega, instances::getBulkMomenta<T, DESCRIPTOR>(),
      0.15);
#else //ConsitentStrainSmagorinsky
  bulkDynamics = new ConStrainSmagorinskyBGKdynamics<T, DESCRIPTOR>(
    omega, instances::getBulkMomenta<T, DESCRIPTOR>(),
      0.05);
#endif

  prepareLattice( sLattice, converter, *bulkDynamics, superGeometry);

  // === 3b Step: Prepare particles ===
  SuperParticleSystem3D<T, Particle3D> superParticleSystem(superGeometry);
  SuperLatticeInterpPhysVelocity3D<T, DESCRIPTOR> velocityProbe(sLattice, converter);

  // initialise pvd file for particles visualisation
  SuperParticleSysVtuWriter<T, Particle3D> particleVtuWriter(superParticleSystem, "particles",
      SuperParticleSysVtuWriter<T, Particle3D>::particleProperties::velocity
        |
      SuperParticleSysVtuWriter<T, Particle3D>::particleProperties::mass
        |
      SuperParticleSysVtuWriter<T, Particle3D>::particleProperties::radius);

  size_t p_nb;
  if (cough) {
    p_nb = (int)random->Uniform(300, 2500);
  }
  else {
    p_nb = (int)random->Uniform(600, 700);
  }
  Timer<T> particleTimer(converter.getLatticeTime(maxPhysT), p_nb);

  MPI_Comm comm;
  comm = MPI_COMM_WORLD;
  MPI_Barrier(comm);

  prepareParticles(sLattice, converter, superGeometry, superParticleSystem, velocityProbe, p_nb, cough, bdm_agents);

  // === 4th Step: Main Loop with Timer ===
  Timer<T> timer( converter.getLatticeTime( maxPhysT ),
    superGeometry.getStatistics().getNvoxel() );
  timer.start();

  OstreamManager cloutsim( std::cout,"openlb simulation" );
  cloutsim << "Turbulence simulation ..." << std::endl;

  int max_step = converter.getLatticeTime( maxPhysT );
  int print_step = (int)max_step/5;

  using clock_t = std::chrono::high_resolution_clock;
  using second_t = std::chrono::duration<double, std::ratio<1> >;
  auto start_time = clock_t::now();

  for ( int iT = 0; iT <= max_step; ++iT ) {

    if (iT % print_step == 0) {
      cloutsim << "Step " << iT << " out of " << max_step << std::endl;
    }

    if (iT == 10) {
      cloutsim << "Estimated run time for the simulation: "
               << (std::chrono::duration_cast<second_t>
                    (clock_t::now() - start_time ).count() ) * max_step/10
               << " seconds" << std::endl;
    }

    // === 5ath Step: Apply filter
#ifdef ADM
    SuperLatticeADM3D<T, DESCRIPTOR> admF( sLattice, 0.01, 2 );
    admF.execute( superGeometry, 1 );
#endif

    // === 5bth Step: Definition of Initial and Boundary Conditions ===
    setBoundaryValues( converter, sLattice, superGeometry, iT );

    // === 6th Step: Collide and Stream Execution ===
    sLattice.collideAndStream();

    // === 7th Step: Computation and Output of the Results ===
    getResults( sLattice, converter, iT, superGeometry, timer );
  } // end for getLatticeTime

  cloutsim << "Particles rendering ..." << std::endl;
  // particles vtk files for visualisation
  particleTimer.start();

  for (int iT = 0; iT <= max_step; ++iT) {
      superParticleSystem.simulate(converter.getConversionFactorTime());
      if (iT % outputSteps == 0) {
          particleTimer.update(iT);
          particleVtuWriter.write(iT);
      }
  }
  particleTimer.stop();


  cloutsim << "Turbulence simulation ... OK" << std::endl;
  timer.stop();
  timer.printSummary();

  // -------- biodynamo -------- //
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  // execute only if master process
  if (world_rank == 0) {
    // update agents' viral load
    for (size_t i = 0; i < bdm_agents.size(); i++ ) {
      bdm::dummy_agent agent = bdm_agents[i];
      // NOTE: count only for healthy agents
      if (agent.state == bdm::State::kHealthy) {
	bdm::Double3 pos = agent.pos;
	auto dir = agent.dir;
	// convert radius to m
	double radius = agent.radius / 100;
	double droplet_viral_load = 1000;

	// cube min and max for particle count
	double measure_min_x = (pos[0]/100 - 0.02 + dir[0]*radius)*converter.getCharPhysLength();
	double measure_min_y = (pos[1]/100 - 0.02 + dir[1]*radius)*converter.getCharPhysLength();
	double measure_min_z = (pos[2]/100 - 0.03)*converter.getCharPhysLength();
	double measure_max_x = (pos[0]/100 + 0.04 + dir[0]*radius)*converter.getCharPhysLength();
	double measure_max_y = (pos[1]/100 + 0.04 + dir[1]*radius)*converter.getCharPhysLength();
	double measure_max_z = (pos[2]/100 + 0.03)*converter.getCharPhysLength();

	// get particles
	// for each particleSystem in superParticleSystem
	for (size_t ps_it = 0 ; ps_it < superParticleSystem.getParticleSystems().size(); ps_it ++) {
	  auto particles = superParticleSystem.getParticleSystems()[ps_it]->getParticles();
	  // for each particle in particleSystem
	  for (size_t p_it = 0; p_it < particles.size(); p_it ++) {
	    auto particle = particles[p_it];
	    auto p_pos = particle.getPos();

	    if (p_pos[0] > measure_min_x && p_pos[1] > measure_min_y && p_pos[2] > measure_min_z
		&& p_pos[0] < measure_max_x && p_pos[1] < measure_max_y && p_pos[2] < measure_max_z) {
	      agent.load += ((4 * M_PI * std::pow(particle.getRad(), 3)) / 3) * droplet_viral_load;
	    }
	  } // end for particles in particleSystem
	} // end for particleSystem in superParticleSystem
      } // end if healthy agent
    } // end for each agents

    // update bdm_sim_state.txt file with new viral load
    std::ofstream simu_file;
    simu_file.open("./output/openlb/bdm_sim_state.txt");
    
    for (size_t i = 0; i < bdm_agents.size(); i++ ) {
      bdm::dummy_agent agent = bdm_agents[i];
      simu_file << "{\n"
		<< "  " << agent.pos[0] << "\n"
		<< "  " << agent.pos[1] << "\n"
		<< "  " << agent.pos[2] << "\n"
		<< "  " << agent.dir[0] << "\n"
		<< "  " << agent.dir[1] << "\n"
		<< "  " << agent.radius * 2 << "\n"
		<< "  " << agent.state << "\n"
		<< "  " << agent.load << "\n"
		<< "}\n";
    } // end for each agents
    simu_file.close();
  } // end if master process

  delete bulkDynamics;

  return 1;
} // end main

} // namespace openlb_sim

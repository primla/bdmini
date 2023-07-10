// SPDX-FileCopyrightText: Copyright 2022-2023 Harusato Kimura
// SPDX-License-Identifier: MIT
#ifndef PSO_HPP
#define PSO_HPP
#include <cerrno>
#include <cmath>
#include <random>
#include <vector>

#define NDEBUG
#include <cassert>

struct Particle
{
  std::vector<double> x;            // position
  std::vector<double> v;            // velocity
  std::vector<double> particle_opt; // particle minimum
};

class PSO
{
public:
  double precision = 1.0e-10;
  int itr_max = 10000;
  int num_of_particle = 1000;
  int dim;
  std::vector<Particle> swarm;

  std::vector<double> global_opt;
  std::vector<double> boundary_min;
  std::vector<double> boundary_max;

  double w_max = 0.9;
  double w_min = 0.4;
  // R. C. Eberhart and Y. Shi, "Comparing inertia weights and constriction
  // factors in particle swarm optimization," Proceedings of the 2000 Congress
  // on Evolutionary Computation. CEC00 (Cat. No.00TH8512), 2000, pp. 84-88
  // vol.1, doi: 10.1109/CEC.2000.870279.

  double weight4local = 1.5;
  double weight4global = 1.5;

  double mutation_rate = 0.05;

  double seed4initialization = 1.0;

  std::mt19937_64 mt;

  double (*func)(const std::vector<double>&); // Evaluation function: R^n -> R

  //---------------------------------------------------------------------------

  PSO(std::vector<double> boundary_min,
      std::vector<double> boundary_max,
      double (*func)(const std::vector<double>&));

  void set_basic_params(int num_of_particle, int itr_max);
  void set_examination_domain(std::vector<double> boundary_min,
                              std::vector<double> boundary_max);
  void set_eval_func(double (*func)(const std::vector<double>&));
  void set_weight_params(double for_local, double for_global);
  void set_seed(double new_seed);
  void set_mrate(double new_mrate);
  void set_precision(double new_precision);

  void solve();

  std::vector<double> get_ans();
  //---------------------------------------------------------------------------
private:
  void init_particles();
};

#endif

// SPDX-FileCopyrightText: Copyright 2022-2023 Harusato Kimura
// SPDX-License-Identifier: MIT
#include "pso.hpp"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

PSO::PSO(std::vector<double> boundary_min,
         std::vector<double> boundary_max,
         double (*func)(const std::vector<double>&))
{
  if (boundary_min.size() != boundary_max.size()) {
    perror("size mismatch");
    exit(errno);
  }
  this->dim = (int)boundary_min.size();

  this->boundary_min.clear();
  this->boundary_max.clear();
  for (double v : boundary_min) {
    this->boundary_min.push_back(v);
  }
  for (double v : boundary_max) {
    this->boundary_max.push_back(v);
  }

  this->func = func;
}

void
PSO::set_basic_params(int num_of_particle, int itr_max)
{
  this->num_of_particle = num_of_particle;
  this->itr_max = itr_max;
}

void
PSO::set_weight_params(double for_local, double for_global)
{
  this->weight4local = for_local;
  this->weight4global = for_global;
}

void
PSO::set_seed(double new_seed)
{
  this->seed4initialization = new_seed;
}

void
PSO::set_mrate(double new_mrate)
{
  this->mutation_rate = new_mrate;
}

void
PSO::set_precision(double new_precision)
{
  this->precision = new_precision;
}

void
PSO::set_examination_domain(std::vector<double> boundary_min,
                            std::vector<double> boundary_max)
{
  if (boundary_min.size() != boundary_max.size()) {
    perror("size mismatch");
    exit(errno);
  }
  this->dim = (int)boundary_min.size();

  this->boundary_min.clear();
  this->boundary_max.clear();
  for (double v : boundary_min) {
    this->boundary_min.push_back(v);
  }
  for (double v : boundary_max) {
    this->boundary_max.push_back(v);
  }
}

void
PSO::set_eval_func(double (*func)(const std::vector<double>&))
{
  this->func = func;
}

void
PSO::init_particles()
{
  this->swarm.clear();
  this->mt.seed(seed4initialization);
  for (int i = 0; i < this->num_of_particle; ++i) {
    Particle new_particle;
    for (int j = 0; j < this->dim; ++j) {
      std::uniform_real_distribution<> uni(this->boundary_min[j],
                                           this->boundary_max[j]);
      double val = uni(mt);
      new_particle.x.push_back(val);
      new_particle.particle_opt.push_back(val);
      new_particle.v.push_back(0);
    }
    this->swarm.push_back(new_particle);
  }
  for (int i = 0; i < dim; ++i) {
    std::uniform_real_distribution<> uni(this->boundary_min[i],
                                         this->boundary_max[i]);
    this->global_opt.push_back(uni(mt));
  }
}

void
PSO::solve()
{
  init_particles();

  std::stringstream fname;

  int noupdate = 0;
  std::uniform_real_distribution<> uni(0.0, 1.0);
  double w = 0;
  std::vector<std::vector<double>> global_opts;
  std::uniform_real_distribution mut(0.0, 1.0);
  for (int itr = 0; itr < itr_max; ++itr) {
    w = w_max - (w_max - w_min) / (double)itr_max * (double)itr;
    for (int i = 0; i < num_of_particle; ++i) {
      // update v
      double r1 = uni(this->mt);
      double r2 = uni(this->mt);
      for (int j = 0; j < dim; ++j) {
        swarm[i].v[j] *= w;
        swarm[i].v[j] +=
          weight4local * r1 * (swarm[i].particle_opt[j] - swarm[i].x[j]);
        swarm[i].v[j] += weight4global * r2 * (global_opt[j] - swarm[i].x[j]);
      }

      // update x
      for (int j = 0; j < dim; ++j) {
        swarm[i].x[j] += swarm[i].v[j];
      }

      // check boundary
      for (int j = 0; j < dim; ++j) {
        if(swarm[i].x[j] < this->boundary_min[j]){
          swarm[i].x[j] = boundary_min[j];
        }else if(swarm[i].x[j] > this->boundary_max[j]){
          swarm[i].x[j] = boundary_max[j];
        }
        //if (swarm[i].x[j] < this->boundary_min[j] ||
        //    swarm[i].x[j] > this->boundary_max[j]) {
        //  std::uniform_real_distribution uni(boundary_min[j], boundary_max[j]);
        //  swarm[i].x[j] = uni(mt);
        //}
      }

      // mutation
      // if(mut(mt) < this->mutation_rate){
      //   for(int j = 0; j < dim; j++){
      //     std::uniform_real_distribution uni(boundary_min[j],
      //     boundary_max[j]); swarm[i].x[j] = uni(mt);
      //   }
      // }

      noupdate++;
      // update local and global opts
      if (func(swarm[i].x) < func(global_opt)) {
        noupdate = 0;
        global_opt.clear();
        swarm[i].particle_opt.clear();
        for (double v : swarm[i].x) {
          global_opt.push_back(v);
          swarm[i].particle_opt.push_back(v);
        }
      } else if (func(swarm[i].x) < func(swarm[i].particle_opt)) {
        swarm[i].particle_opt.clear();
        for (double v : swarm[i].x) {
          swarm[i].particle_opt.push_back(v);
        }
      }
      // check convergence
      if (func(global_opt) < this->precision) {
        // Judge that satisfactory results have been achieved.
        break;
      } else if (noupdate == 20) {
        // Judge that no further improvement can be expected.
        break;
      }
    }
    global_opts.push_back(global_opt);
  }
}

std::vector<double>
PSO::get_ans()
{
  return this->global_opt;
}

// SPDX-FileCopyrightText: Copyright 2022-2023 Harusato Kimura
// SPDX-License-Identifier: MIT
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "pso.hpp"

// GLOBAL VARIABLES
double offset_global = 0.0; // rad
// offset(phi) to reference direction

std::vector<std::vector<double>> ccf_re_global;
std::vector<std::vector<double>> ccf_im_global;
std::vector<double> freq_global;
std::vector<double> phv_global;
int ind_global; // index of current frequency and phv.

std::vector<double> r_ij_global;
std::vector<double> the_ij_global;

//-----------------------------------------------------------------------------
// Prototype declaration

double
test_func_re(std::vector<double> const& in);

double
test_func_im(std::vector<double> const& in);

void
read_coord_spacelag(std::string const targ_dir,
                    std::vector<std::string>& site_name,
                    std::vector<std::string> const& valid_site_name);

void
read_ccf(std::string const targ_dir, std::vector<std::string> const& site_name);

double
get_rmax();

void
solve_real_part(std::string const targ_dir,
                int const n_particle,
                int const n_itr,
                double const w4loc,
                double const w4glo);

void
solve_imag_part(std::string const targ_dir,
                int const n_particle,
                int const n_itr,
                double const w4loc,
                double const w4glo);

//-----------------------------------------------------------------------------

double
test_func_re(std::vector<double> const& in)
{
  double k = in[0];

  double ans = 0;
  for (size_t i = 0; i < ccf_re_global.size(); ++i) {
    double tmp = std::cyl_bessel_j(0, k * r_ij_global[i]);
    for (int n = 1; n < 3; ++n) {
      double core = in[2 * n - 1] * std::cos(2.0 * n * the_ij_global[i]) -
                    in[2 * n] * std::sin(2.0 * n * the_ij_global[i]);
      tmp += 2.0 * std::pow(-1.0, n) *
             std::cyl_bessel_j(2 * n, k * r_ij_global[i]) * core;
    }
    ans += std::pow(ccf_re_global[i][ind_global] - tmp, 2);
  }

  return ans;
}

//-----------------------------------------------------------------------------

double
test_func_im(std::vector<double> const& in)
{
  double k = 2.0 * M_PI * freq_global[ind_global] / phv_global[ind_global];

  double ans = 0;
  for (size_t i = 0; i < ccf_im_global.size(); ++i) {
    double tmp = 0;
    for (int n = 1; n < 3; ++n) {
      double core =
        in[2 * n - 2] * std::cos((2.0 * n - 1.0) * the_ij_global[i]) -
        in[2 * n - 1] * std::sin((2.0 * n - 1.0) * the_ij_global[i]);
      tmp += 2.0 * std::pow(-1.0, n) *
             std::cyl_bessel_j(2 * n - 1, k * r_ij_global[i]) * core;
    }
    ans += std::pow(ccf_im_global[i][ind_global] + tmp, 2);
  }

  return ans;
}

//-----------------------------------------------------------------------------

void
read_coord_spacelag(std::string const targ_dir,
                    std::vector<std::string>& site_name,
                    std::vector<std::string> const& valid_site_name)
{
  r_ij_global.clear();
  the_ij_global.clear();

  std::string line, s_x, s_y, s_name_raw, s_name;
  std::ifstream inf{ targ_dir + "array_coord.csv" };
  std::vector<double> xc, yc;
  int nsite = 0;
  while (std::getline(inf, line)) {
    std::istringstream iss(line);
    std::getline(iss, s_x, ',');
    std::getline(iss, s_y, ',');
    std::getline(iss, s_name_raw, ',');
    s_name_raw.erase(
      std::remove_if(s_name_raw.begin(), s_name_raw.end(), ::isspace),
      s_name_raw.end());
    s_name = s_name_raw.substr(0, s_name_raw.find('.'));
    for (size_t i = 0; i < valid_site_name.size(); ++i) {
      if (valid_site_name[i] == s_name) {
        xc.push_back(std::stod(s_x));
        yc.push_back(std::stod(s_y));
        site_name.push_back(s_name);
        nsite++;
      }
    }
  }
  inf.close();

  for (int i = 1; i < (int)xc.size(); ++i) {
    for (int j = 0; j < i; ++j) {
      r_ij_global.push_back(
        std::sqrt(std::pow(xc[i] - xc[j], 2) + std::pow(yc[i] - yc[j], 2)));
      the_ij_global.push_back(std::atan2(yc[i] - yc[j], xc[i] - xc[j]) -
                              offset_global);
    }
  }
}

//-----------------------------------------------------------------------------

void
read_ccf(std::string const targ_dir, std::vector<std::string> const& site_name)
{
  std::string line, val;

  freq_global.clear();
  ccf_re_global.clear();
  ccf_im_global.clear();
  for (int i = 1; i < (int)site_name.size(); ++i) {
    for (int j = 0; j < i; ++j) {
      std::string path = targ_dir + "results/statistics/CCF_UD_";
      path += site_name[j] + "-" + site_name[i] + ".csv";
      std::ifstream inf{ path };

      std::vector<double> tmp1;
      std::vector<double> tmp2;
      while (std::getline(inf, line)) {
        std::istringstream iss(line);
        std::getline(iss, val, ',');
        double freq_tmp = std::stod(val);
        if (j == 0 && i == 1) {
          if (!(freq_tmp > 30.0)) {
            freq_global.push_back(freq_tmp);
          }
        }
        // read real part of CCF
        std::getline(iss, val, ',');
        double tmp1_v;
        try {
          tmp1_v = std::stod(val);
        } catch (std::out_of_range& e) {
          std::cerr << "std::out_of_range: what(): " << e.what() << std::endl;
          tmp1_v = 0;
        }
        tmp1.push_back(tmp1_v);
        // read imag part of CCF
        std::getline(iss, val, ',');
        double tmp2_v;
        try {
          tmp2_v = std::stod(val);
        } catch (std::out_of_range& e) {
          std::cerr << "std::out_of_range: what(): " << e.what() << std::endl;
          tmp2_v = 0;
        }
        tmp2.push_back(tmp2_v);
      }
      ccf_re_global.push_back(tmp1);
      ccf_im_global.push_back(tmp2);
    }
  }
}

//-----------------------------------------------------------------------------

double
get_rmax()
{
  std::vector<double> tmp;
  for (auto v : r_ij_global) {
    tmp.push_back(v);
  }
  std::sort(tmp.begin(), tmp.end(), std::greater<double>());
  return tmp[0];
}

//-----------------------------------------------------------------------------

void
solve_real_part(std::string const targ_dir,
                int const n_particle,
                int const n_itr,
                double const w4loc,
                double const w4glo)
{

  phv_global.clear();
  std::ofstream ouf{ targ_dir + "results/dspac/result_real.csv" };
  double rmax = get_rmax();
  std::vector<double> bound_min{ 0, -1, -1, -1, -1 };
  std::vector<double> bound_max{ M_PI / rmax, 1, 1, 1, 1 };

  for (int i = 0; i < (int)freq_global.size(); ++i) {
    ind_global = i;

    PSO p1{ bound_min, bound_max, test_func_re };

    p1.set_seed(0);
    p1.set_basic_params(n_particle, n_itr);
    p1.set_weight_params(w4loc, w4glo);

    p1.solve();
    std::vector<double> ans = p1.get_ans();
    ans[0] = 2.0 * M_PI * freq_global[ind_global] /
             ans[0]; // conv to phv_global from k

    phv_global.push_back(ans[0]);

    // output section
    ouf << freq_global[ind_global];
    for (size_t k = 0; k < ans.size() - 2; ++k) {
      ouf << ", " << ans[k];
    }
    ouf << std::endl;
  }
  ouf.close();
}

//-----------------------------------------------------------------------------

void
solve_imag_part(std::string const targ_dir,
                int const n_particle,
                int const n_itr,
                double const w4loc,
                double const w4glo)
{
  std::ofstream ouf{ targ_dir + "results/dspac/result_imag.csv" };
  std::vector<double> bound_min{ -1, -1, -1, -1 };
  std::vector<double> bound_max{ 1, 1, 1, 1 };

  for (int i = 0; i < (int)freq_global.size(); ++i) {
    ind_global = i;

    PSO p1{ bound_min, bound_max, test_func_im };

    p1.set_seed(1);

    p1.set_basic_params(n_particle, n_itr);
    p1.set_weight_params(w4loc, w4glo);

    p1.solve();
    std::vector<double> ans = p1.get_ans();

    // output section
    ouf << freq_global[ind_global];
    for (size_t k = 0; k < ans.size() - 2; ++k) {
      ouf << ", " << ans[k];
    }
    ouf << std::endl;
  }
  ouf.close();
}

//-----------------------------------------------------------------------------

int
main(int argc, char** argv)
{
  // Parse argv
  std::string targ_dir = argv[1];
  int n_particle = std::stoi(argv[2]);
  int n_itr = std::stoi(argv[3]);
  double w4loc = std::stod(argv[4]);
  double w4glo = std::stod(argv[5]);
  std::vector<std::string> valid_site_name;
  for (int i = 6; i < argc; ++i) {
    valid_site_name.push_back(argv[i]);
  }

  // Read array_coord.csv, reslts/statistics/CCF_*
  std::vector<std::string> site_name;
  read_coord_spacelag(targ_dir, site_name, valid_site_name);
  read_ccf(targ_dir, site_name);

  solve_real_part(targ_dir, n_particle, n_itr, w4loc, w4glo);
  solve_imag_part(targ_dir, n_particle, n_itr, w4loc, w4glo);
}

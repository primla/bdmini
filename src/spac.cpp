// SPDX-FileCopyrightText: Copyright 2022 Harusato Kimura
// SPDX-License-Identifier: MIT

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

//-----------------------------------------------------------------------------

std::vector<double> read_vector_line(std::string const path, int const nl);

std::vector<std::vector<double>>
read_re_ccfs(std::string const targ_dir,
             std::vector<std::string> const &array_info);

std::vector<double>
calc_spac_coeff(std::vector<std::vector<double>> const &ccfs);

void write_spac_coeff(std::string const targ_dir,
                      std::vector<double> const &freq,
                      std::vector<double> const &spac_coeff,
                      std::string const array_name);

void write_phv(std::string const targ_dir, std::vector<double> const &freq,
               std::vector<double> const &phv, std::string const array_name);

double find_krmax();

double bisection_j0(double const val, double const x_min, double const x_max);

std::vector<double> calc_phv(double const radius,
                             std::vector<double> const &freq,
                             std::vector<double> const &spac_coeff);

void read_coordinate(std::string const targ_dir, std::vector<double> &x_coord,
                     std::vector<double> &y_coord,
                     std::vector<std::string> &site_name);

double calc_radius(std::string const targ_dir,
                   std::vector<std::string> const &array_info);

//-----------------------------------------------------------------------------

std::vector<double> read_vector_line(std::string const path, int const nl) {
  std::vector<double> ans;
  std::ifstream inf{path};
  std::string line, val;
  while (std::getline(inf, line)) {
    std::stringstream iss(line);
    for (int i = 0; i < nl; ++i) {
      std::getline(iss, val, ',');
    }
    std::getline(iss, val, ',');
    ans.push_back(std::stod(val));
  }
  return ans;
}

//-----------------------------------------------------------------------------

std::vector<std::vector<double>>
read_re_ccfs(std::string const targ_dir,
             std::vector<std::string> const &array_info) {
  std::vector<std::vector<double>> ans;
  for (size_t i = 0; i < array_info.size() / 2; ++i) {
    std::string path = targ_dir + "results/statistics/CCF_UD_" +
                       array_info[2 * i] + "-" + array_info[2 * i + 1] + ".csv";
    std::vector<double> re_ccf = read_vector_line(path, 1);
    ans.push_back(re_ccf);
  }
  return ans;
}

//-----------------------------------------------------------------------------

std::vector<double>
calc_spac_coeff(std::vector<std::vector<double>> const &ccfs) {
  size_t n_pairs = ccfs.size();
  size_t dim = ccfs[0].size();
  std::vector<double> ans;
  for (size_t i = 0; i < dim; ++i) {
    ans.push_back(0.0);
  }
  for (size_t i = 0; i < n_pairs; ++i) {
    for (size_t j = 0; j < dim; ++j) {
      ans[j] += ccfs[i][j];
    }
  }
  for (size_t j = 0; j < dim; ++j) {
    ans[j] /= (double)n_pairs;
  }
  return ans;
}

//-----------------------------------------------------------------------------

void write_spac_coeff(std::string const targ_dir,
                      std::vector<double> const &freq,
                      std::vector<double> const &spac_coeff,
                      std::string const array_name) {
  std::ofstream ouf{targ_dir + "results/spac/spr_" + array_name + ".csv"};
  for (size_t i = 0; i < freq.size(); ++i) {
    ouf << freq[i] << ", ";
    ouf << spac_coeff[i] << std::endl;
  }
}

//-----------------------------------------------------------------------------

void write_phv(std::string const targ_dir, std::vector<double> const &freq,
               std::vector<double> const &phv, std::string const array_name) {
  std::ofstream ouf{targ_dir + "results/spac/phv_" + array_name + ".csv"};
  for (size_t i = 0; i < freq.size(); ++i) {
    ouf << freq[i] << ", ";
    ouf << phv[i] << std::endl;
  }
}

//-----------------------------------------------------------------------------

double find_krmax() {
  double const MAX_KR_TO_SEARCH = 10;
  double const INCR = 0.001;
  int num_sample = (int)(MAX_KR_TO_SEARCH / INCR);

  double j0_bf, j0, j0_af;

  j0_bf = std::cyl_bessel_j(0, 0);
  j0 = std::cyl_bessel_j(0, INCR * 1.0);
  j0_af = std::cyl_bessel_j(0, INCR * 2.0);

  for (int i = 1; i <= num_sample; ++i) {
    if ((j0 - j0_bf) * (j0_af - j0) < 0) {
      return INCR * (double)i;
    }
    j0_bf = j0;
    j0 = j0_af;
    j0_af = std::cyl_bessel_j(0, INCR * (double)(i + 2));
  }

  perror("Error -> krmax");
  exit(errno);
}

//-----------------------------------------------------------------------------

double bisection_j0(double const val, double const x_min, double const x_max) {
  int const MAX_ITR = 1e4;
  double const ACC = 1e-5;
  double x_l, x_r;
  double y_l, y_r;
  double x_c, y_c;

  y_l = std::cyl_bessel_j(0, x_min) - val;
  y_r = std::cyl_bessel_j(0, x_max) - val;
  x_l = x_min;
  x_r = x_max;

  if (y_l * y_r > 0.0) {
    perror("Error -> Bisection");
    exit(errno);
  }

  for (int i = 0; i < MAX_ITR; ++i) {
    x_c = (x_l + x_r) / 2;
    if (x_r - x_l < ACC) {
      if (i == MAX_ITR) {
        perror("Error -> Bisection");
        exit(errno);
      }
      return x_c;
    }

    y_c = std::cyl_bessel_j(0, x_c) - val;
    y_l = std::cyl_bessel_j(0, x_l) - val;
    if (y_c * y_l <= 0.0) {
      x_r = x_c;
    } else {
      x_l = x_c;
    }
  }

  perror("Error -> Bisection");
  exit(errno);
}

//-----------------------------------------------------------------------------

std::vector<double> calc_phv(double const radius,
                             std::vector<double> const &freq,
                             std::vector<double> const &spac_coeff) {
  double krmax = find_krmax();
  std::vector<double> ans;
  for (size_t i = 0; i < freq.size(); ++i) {
    if (spac_coeff[i] >= 1) {
      ans.push_back(2.0 * M_PI * freq[i] /
                    (bisection_j0(1, 0, krmax) / radius));
    } else if (spac_coeff[i] <= std::cyl_bessel_j(0, krmax)) {
      ans.push_back(
          2.0 * M_PI * freq[i] /
          (bisection_j0(std::cyl_bessel_j(0, krmax), 0, krmax) / radius));
    } else {
      ans.push_back(2.0 * M_PI * freq[i] /
                    (bisection_j0(spac_coeff[i], 0, krmax) / radius));
    }
  }
  return ans;
}

//-----------------------------------------------------------------------------

void read_coordinate(std::string const targ_dir, std::vector<double> &x_coord,
                     std::vector<double> &y_coord,
                     std::vector<std::string> &site_name) {
  std::string line, val;
  std::ifstream inf{targ_dir + "array_coord.csv"};

  // read array_coord.csv
  while (std::getline(inf, line)) {
    std::stringstream iss{line};
    std::getline(iss, val, ',');
    x_coord.push_back(std::stod(val));
    std::getline(iss, val, ',');
    y_coord.push_back(std::stod(val));
    std::getline(iss, val, ',');
    val.erase(std::remove_if(val.begin(), val.end(), ::isspace), val.end());
    site_name.push_back(val.substr(0, val.find('.')));
  }
}

//-----------------------------------------------------------------------------

double calc_radius(std::string const targ_dir,
                   std::vector<std::string> const &array_info) {
  std::vector<double> x_coord;
  std::vector<double> y_coord;
  std::vector<std::string> site_name;
  read_coordinate(targ_dir, x_coord, y_coord, site_name);

  double radius = 0;
  for (size_t i = 0; i < array_info.size() / 2; ++i) {
    double x1 = 0;
    double y1 = 0;
    double x2 = 0;
    double y2 = 0;
    for (size_t j = 0; j < site_name.size(); ++j) {
      if (site_name[j] == array_info[2 * i]) {
        x1 = x_coord[j];
        y1 = y_coord[j];
      }
    }
    for (size_t j = 0; j < site_name.size(); ++j) {
      if (site_name[j] == array_info[2 * i + 1]) {
        x2 = x_coord[j];
        y2 = y_coord[j];
      }
    }
    if (std::sqrt(std::pow(x1 - x2, 2) + std::pow(y1 - y2, 2)) == 0) {
      perror("Error -> array_info");
      exit(errno);
    }
    radius += std::sqrt(std::pow(x1 - x2, 2) + std::pow(y1 - y2, 2));
  }
  return radius / ((double)array_info.size() / 2.0);
}

//-----------------------------------------------------------------------------
// #############################################################################

int main(int argc, char **argv) {

  //===========================================================================
  // Preparation
  std::string targ_dir = argv[1];
  std::string array_name = argv[2];
  std::vector<std::string> array_info;
  for (int i = 0; i < argc - 3; ++i) {
    array_info.push_back(argv[i + 3]);
  }
  // calc radius
  double radius = calc_radius(targ_dir, array_info);

  std::vector<double> freq =
      read_vector_line(targ_dir + "results/statistics/UD_" + array_info[0] +
                           "-" + array_info[0] + ".csv",
                       0);

  //===========================================================================

  std::vector<std::vector<double>> re_ccfs = read_re_ccfs(targ_dir, array_info);
  std::vector<double> spac_coeff = calc_spac_coeff(re_ccfs);
  write_spac_coeff(targ_dir, freq, spac_coeff, array_name);

  std::vector<double> phv = calc_phv(radius, freq, spac_coeff);
  write_phv(targ_dir, freq, phv, array_name);

  //===========================================================================
}

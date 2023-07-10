// SPDX-FileCopyrightText: Copyright 2022-2023 Harusato Kimura
// SPDX-License-Identifier: MIT

#include <algorithm>
#include <cassert>
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

struct CspecMatrix
{
  double freq;
  std::vector<std::vector<std::complex<double>>> csdmat; // matrix of CSD
  std::vector<std::vector<std::complex<double>>>
    csdinv; // inverse matrix of csdmat
};

//-----------------------------------------------------------------------------
// Prototype declaration

std::vector<std::complex<double>>
read_cspec(std::string const path);

std::vector<double>
read_vector_line(std::string const path, int const nl);

void
read_coordinate(std::string const targ_dir,
                std::vector<double>& x_coord,
                std::vector<double>& y_coord,
                std::vector<std::string>& site_name);

std::vector<std::vector<std::complex<double>>>
calc_inv(std::vector<std::vector<std::complex<double>>> in);

std::vector<CspecMatrix>
read_cspec_matrices(std::string const targ_dir,
                    std::vector<std::string> const& site_name,
                    std::vector<double> const& freq);

void
dump_fk_estimate(std::string const targ_dir,
                 std::vector<std::vector<double>>& ans,
                 std::vector<double> const& s_val,
                 std::vector<double> const& phi,
                 double const freq);

void
dump_vecs_2(std::string const targ_dir,
            std::string const fname,
            std::vector<double> const& freq,
            std::vector<double> const& vec);

double
calc_slw_peak_circ_stat(std::vector<std::vector<double>> const& ans,
                        std::vector<double> const& s_val,
                        std::vector<double> const& phi,
                        std::vector<double>& phase,
                        std::vector<double>& amp,
                        size_t const trig_mom_max);

void
do_fkanaly(std::string const targ_dir,
           std::vector<double> const& x_coord,
           std::vector<double> const& y_coord,
           std::vector<CspecMatrix> const& cspec_matrices,
           int const density,
           int const density_phi,
           double const max_phv,
           double const min_phv,
           size_t const trig_mom_max);

void
dump_fcoeffs(std::string const targ_dir,
             std::vector<double> const& freq,
             std::vector<std::vector<double>> const& phases,
             std::vector<std::vector<double>> const& amps,
             size_t const trig_mom_max);

//-----------------------------------------------------------------------------

std::vector<std::complex<double>>
read_cspec(std::string const path)
{
  std::vector<std::complex<double>> ans;
  std::ifstream inf{ path };
  std::string line, val;
  while (std::getline(inf, line)) {
    std::stringstream iss{ line };
    // ignore $1 (= freq)
    std::getline(iss, val, ',');

    // real part
    std::getline(iss, val, ',');
    double re_part = std::stod(val);
    // double re_part = 0;
    // try {
    //   re_part = std::stod(val);
    // } catch (std::out_of_range& e) {
    //   std::cerr << "std::out_of_range: what(): " << e.what() << std::endl;
    //   re_part = 0;
    // }

    // imag part
    std::getline(iss, val, ',');
    double im_part = std::stod(val);
    // double im_part = 0;
    // try {
    //   im_part = std::stod(val);
    // } catch (std::out_of_range& e) {
    //   std::cerr << "std::out_of_range: what(): " << e.what() << std::endl;
    //   im_part = 0;
    // }

    std::complex<double> tmp{ re_part, im_part };
    ans.push_back(tmp);
  }
  return ans;
}

//-----------------------------------------------------------------------------

std::vector<double>
read_vector_line(std::string const path, int const nl)
{
  std::vector<double> ans;
  std::ifstream inf{ path };
  std::string line, val;
  while (std::getline(inf, line)) {
    std::stringstream iss(line);
    for (int i = 0; i < nl; ++i) {
      std::getline(iss, val, ',');
    }
    std::getline(iss, val, ',');
    double val_d = std::stod(val);
    // double val_d = 0.0;
    // try {
    //   val_d = std::stod(val);
    // } catch (std::out_of_range& e) {
    //   std::cerr << "std::out_of_range: what(): " << e.what() << std::endl;
    //   val_d = 0.0;
    // }
    ans.push_back(val_d);
  }
  return ans;
}

//-----------------------------------------------------------------------------

void
read_coordinate(std::string const targ_dir,
                std::vector<double>& x_coord,
                std::vector<double>& y_coord,
                std::vector<std::string>& site_name)
{
  std::string line, val;
  std::ifstream inf{ targ_dir + "array_coord.csv" };

  // read array_coord.csv
  while (std::getline(inf, line)) {
    std::stringstream iss{ line };

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

std::vector<std::vector<std::complex<double>>>
calc_inv(std::vector<std::vector<std::complex<double>>> in)
{
  // Preparation
  size_t dim = in.size();
  std::complex<double> zero{ 0.0, 0.0 };
  std::vector<std::complex<double>> zeros(dim, zero);
  std::vector<std::vector<std::complex<double>>> ans(dim, zeros);

  std::complex<double> one{ 1.0, 0.0 };
  for (size_t i = 0; i < dim; ++i) {
    ans[i][i] = one;
  }

  // Gauss-Jordan elimination
  for (size_t i = 0; i < dim; ++i) {
    size_t pivot_row = i;
    for (size_t j = i + 1; j < dim; ++j) {
      if (std::norm(in[j][i]) > std::norm(in[i][i])) {
        pivot_row = j;
      }
    }
    if (pivot_row != i) {
      std::vector<std::complex<double>> tmp_1 = in[i];
      std::vector<std::complex<double>> tmp_2 = ans[i];
      in[i] = in[pivot_row];
      ans[i] = ans[pivot_row];
      in[pivot_row] = tmp_1;
      ans[pivot_row] = tmp_2;
    }

    // do sweep
    std::complex<double> key = in[i][i];
    in[i][i] = one;
    for (size_t j = i + 1; j < dim; ++j) {
      in[i][j] /= key;
    }
    for (size_t j = 0; j < dim; ++j) {
      ans[i][j] /= key;
    }
    // for l-th row
    for (size_t l = 0; l < dim; ++l) {
      if (l != i) {
        std::complex<double> s = in[l][i];
        if (s != zero) {
          for (size_t k = i; k < dim; ++k) {
            in[l][k] -= s * in[i][k];
          }
          for (size_t k = 0; k < dim; ++k) {
            ans[l][k] -= s * ans[i][k];
          }
        }
      }
    }
  }

  return ans;
}

//-----------------------------------------------------------------------------

std::vector<CspecMatrix>
read_cspec_matrices(std::string const targ_dir,
                    std::vector<std::string> const& site_name,
                    std::vector<double> const& freq)
{
  size_t nsites = site_name.size();

  std::vector<std::vector<std::vector<std::complex<double>>>> cspecmat_raw(
    nsites, std::vector<std::vector<std::complex<double>>>(nsites));

  for (size_t i = 0; i < nsites; ++i) {
    for (size_t j = 0; j < nsites; ++j) {
      // std::string path = targ_dir + "results/statistics/UD_";
      // path += site_name[i] + "-" + site_name[j] + ".csv";
      // cspecmat_raw[i][j] = read_cspec(path);

      // ----------------------------------------------------------------------
      // use ccf instead of cspec
      std::vector<std::complex<double>> ones(freq.size(),
                                             std::complex<double>{ 1.0, 0.0 });
      if (i != j) {
        std::string path = targ_dir + "results/statistics/CCF_UD_";
        path += site_name[i] + "-" + site_name[j] + ".csv";
        cspecmat_raw[i][j] = read_cspec(path);
      } else {
        cspecmat_raw[i][j] = ones;
      }
      // ----------------------------------------------------------------------
    }
  }

  // Conv to CspecMatrix
  std::vector<CspecMatrix> ans(freq.size());
  for (size_t k = 0; k < freq.size(); ++k) {
    CspecMatrix tmp;
    tmp.csdmat.resize(nsites);
    for (size_t i = 0; i < nsites; ++i) {
      tmp.csdmat[i].resize(nsites);
    }
    for (size_t i = 0; i < nsites; ++i) {
      for (size_t j = 0; j < nsites; ++j) {
        tmp.csdmat[i][j] = cspecmat_raw[i][j][k];
      }
    }
    // std::vector<std::vector<std::complex<double>>> inv =
    // calc_inv(tmp.csdmat); tmp.csdinv = inv;
    tmp.csdinv = calc_inv(tmp.csdmat);
    tmp.freq = freq[k];
    ans[k] = tmp;
  }
  return ans;
}

//-----------------------------------------------------------------------------

void
dump_fk_estimate(std::string const targ_dir,
                 std::vector<std::vector<double>>& ans,
                 std::vector<double> const& s_val,
                 std::vector<double> const& phi,
                 double const freq)
{
  std::stringstream fname;
  fname << "results/fk/FK_";
  fname << std::setw(8) << std::setfill('0') << std::fixed
        << std::setprecision(5);
  fname << freq;
  std::string path = fname.str();
  size_t pos;
  while ((pos = path.find(".")) != std::string::npos) {
    path.replace(pos, 1, "p");
  }
  path += "_Hz.csv";

  // normalize fk spectrum
  std::vector<double> rowmax(s_val.size());
  std::vector<double> rowmin(s_val.size());
  for (size_t i = 0; i < s_val.size(); ++i) {
    std::vector<double> tmp = ans[i];
    std::sort(tmp.begin(), tmp.end(), std::greater<double>{});
    rowmax[i] = tmp[0];
    rowmin[i] = tmp[tmp.size() - 1];
  }
  std::sort(rowmax.begin(), rowmax.end(), std::greater<double>{});
  std::sort(rowmin.begin(), rowmin.end(), std::greater<double>{});
  double maxval = rowmax[0];
  double minval = rowmin[rowmin.size() - 1];

  for (size_t i = 0; i < s_val.size(); ++i) {
    for (size_t j = 0; j < phi.size(); ++j) {
      ans[i][j] = (ans[i][j] - minval) / (maxval - minval);
    }
  }

  // write fk spectra
  std::ofstream ouf{ targ_dir + path };
  for (size_t i = 0; i < s_val.size(); ++i) {
    for (size_t j = 0; j < phi.size(); ++j) {
      ouf << s_val[i] * std::cos(phi[j]) << ", ";
      ouf << s_val[i] * std::sin(phi[j]) << ", ";
      ouf << ans[i][j] << std::endl;
    }
    ouf << s_val[i] * std::cos(phi[0]) << ", ";
    ouf << s_val[i] * std::sin(phi[0]) << ", ";
    ouf << 0.5 * (ans[i][phi.size() - 1] + ans[i][0]) << std::endl;
    ouf << std::endl;
  }
}

//-----------------------------------------------------------------------------

void
dump_vecs_2(std::string const targ_dir,
            std::string const fname,
            std::vector<double> const& freq,
            std::vector<double> const& vec)
{
  std::ofstream ouf{ targ_dir + "results/fk/" + fname };
  for (size_t i = 0; i < freq.size(); ++i) {
    ouf << freq[i] << ", " << vec[i] << std::endl;
  }
}

//-----------------------------------------------------------------------------

void
dump_fcoeffs(std::string const targ_dir,
             std::vector<double> const& freq,
             std::vector<std::vector<double>> const& phases,
             std::vector<std::vector<double>> const& amps,
             size_t const trig_mom_max)
{
  std::ofstream ouf{ targ_dir + "results/fk/phases.csv" };
  for (size_t i = 0; i < freq.size(); ++i) {
    ouf << freq[i] << ", ";
    for (size_t j = 0; j < trig_mom_max - 1; ++j) {
      ouf << phases[i][j] << ", ";
    }
    ouf << phases[i][trig_mom_max - 1] << std::endl;
  }
  ouf.close();

  // ---

  ouf.open(targ_dir + "results/fk/amps.csv");
  for (size_t i = 0; i < freq.size(); ++i) {
    ouf << freq[i] << ", ";
    for (size_t j = 0; j < trig_mom_max - 1; ++j) {
      ouf << amps[i][j] << ", ";
    }
    ouf << amps[i][trig_mom_max - 1] << std::endl;
  }
  ouf.close();

  // ---

  ouf.open(targ_dir + "results/fk/re_and_im_coeff.csv");
  for (size_t i = 0; i < freq.size(); ++i) {
    ouf << freq[i] << ", ";
    for (size_t j = 0; j < trig_mom_max - 1; ++j) {
      ouf << amps[i][j] * cos(phases[i][j]) << ", ";
      ouf << amps[i][j] * sin(phases[i][j]) << ", ";
    }
    ouf << amps[i][trig_mom_max - 1] * cos(phases[i][trig_mom_max]) << ", ";
    ouf << amps[i][trig_mom_max - 1] * sin(phases[i][trig_mom_max])
        << std::endl;
  }
  ouf.close();
}

//-----------------------------------------------------------------------------

double
calc_slw_peak_circ_stat(std::vector<std::vector<double>> const& ans,
                        std::vector<double> const& s_val,
                        std::vector<double> const& phi,
                        std::vector<double>& phase,
                        std::vector<double>& amp,
                        size_t const trig_mom_max)
{
  size_t peak_ind = 0;
  bool flag = false;
  for (size_t i = 0; i < s_val.size(); ++i) {
    for (size_t j = 0; j < phi.size(); ++j) {
      if (ans[i][j] == 1.0) {
        peak_ind = i;
        flag = true;
        break;
      }
    }
    if (flag) {
      break;
    }
  }

  phase.resize(trig_mom_max);
  amp.resize(trig_mom_max);
  for (size_t m = 0; m < trig_mom_max; ++m) {
    double cos_m = 0.0;
    double sin_m = 0.0;
    double w_sum = 0.0;
    for (size_t j = 0; j < phi.size(); ++j) {
      cos_m += ans[peak_ind][j] * std::cos(m * phi[j]);
      sin_m -= ans[peak_ind][j] * std::sin(m * phi[j]);
      w_sum += ans[peak_ind][j];
    }
    cos_m /= w_sum;
    sin_m /= w_sum;
    phase[m] = std::atan2(sin_m, cos_m);
    amp[m] = std::sqrt(sin_m * sin_m + cos_m * cos_m);
  }

  return s_val[peak_ind];
}

//-----------------------------------------------------------------------------

void
do_fkanaly(std::string const targ_dir,
           std::vector<double> const& x_coord,
           std::vector<double> const& y_coord,
           std::vector<CspecMatrix> const& cspec_matrices,
           int const density,
           int const density_phi,
           double const max_phv,
           double const min_phv,
           size_t const trig_mom_max)
{
  // prepare slowness vector
  std::vector<double> phi(density_phi);
  std::vector<double> c_val(density);
  for (int i = 0; i < density_phi; ++i) {
    phi[i] = (double)i / (double)density_phi * 2.0 * M_PI - M_PI;
  }
  for (int i = 0; i < density; ++i) {
    c_val[i] = min_phv + (max_phv - min_phv) * (double)i / (double)density;
  }
  std::vector<double> s_val(density);
  for (int i = 0; i < density; ++i) {
    s_val[i] = 1.0 / c_val[i];
  }

  std::sort(s_val.begin(), s_val.end());

  // estimate FWD spectral density
  std::vector<double> freq;
  std::vector<double> phv;
  std::vector<std::vector<double>>
    phases; // phase of first- to (trig_mom_max)th-order coeff.
  std::vector<std::vector<double>>
    amps; // amp of first- to (trig_mom_max)th-order coeff.
  std::complex<double> imunit{ 0.0, 1.0 };
  int nsites = (int)x_coord.size();
  for (int f = 0; f < (int)cspec_matrices.size(); ++f) {
    if (cspec_matrices[f].freq > 30) {
      break;
    }
    if (f % 10 != 0) {
      continue;
    }

    freq.push_back(cspec_matrices[f].freq);

    std::vector<std::vector<double>> ans(density,
                                         std::vector<double>(density_phi));
    // for each bin
    for (int i = 0; i < density; ++i) {
      for (int j = 0; j < density_phi; ++j) {
        // for each csd
        std::complex<double> val{ 0.0, 0.0 };
        for (int k = 0; k < nsites; ++k) {
          for (int l = 0; l < nsites; ++l) {
            double lag_x = x_coord[k] - x_coord[l];
            double lag_y = y_coord[k] - y_coord[l];
            double circ_freq = 2.0 * M_PI * cspec_matrices[f].freq;
            std::complex<double> tmp =
              exp(imunit * circ_freq *
                  (lag_x * s_val[i] * std::cos(phi[j]) +
                   lag_y * s_val[i] * std::sin(phi[j])));
            tmp *= cspec_matrices[f].csdinv[k][l];
            val += tmp;
          }
        }
        ans[i][j] = 1.0 / val.real();
      }
    }
    // !CAUTION! ans is passed by reference
    dump_fk_estimate(targ_dir, ans, s_val, phi, cspec_matrices[f].freq);

    std::vector<double> phase;
    std::vector<double> amp;
    double slw_peak =
      calc_slw_peak_circ_stat(ans, s_val, phi, phase, amp, trig_mom_max);
    phv.push_back(1.0 / slw_peak); // m/s
    phases.push_back(phase);
    amps.push_back(amp);
  }

  dump_vecs_2(targ_dir, "phv_fk.csv", freq, phv);
  dump_fcoeffs(targ_dir, freq, phases, amps, trig_mom_max);
}

//-----------------------------------------------------------------------------

//#############################################################################

int
main(int argc, char** argv)
{
  // Preparation
  if (argc != 6) {
    perror("Error: Invalid argument");
    exit(errno);
  }
  std::string targ_dir = argv[1];
  int density = std::stoi(argv[2]);
  int density_phi = std::stoi(argv[3]);
  double min_phv = std::stoi(argv[4]);
  double max_phv = std::stoi(argv[5]);

  size_t trig_mom_max = 20;

  // Read coordinate
  std::vector<double> x_coord;
  std::vector<double> y_coord;
  std::vector<std::string> site_name;
  read_coordinate(targ_dir, x_coord, y_coord, site_name);

  // Read cspec matrices and calc inverse of cspec matrices.
  std::string path = targ_dir + "results/statistics/UD_" + site_name[0] + "-" +
                     site_name[0] + ".csv";
  std::vector<double> freq = read_vector_line(path, 0);
  std::vector<CspecMatrix> cspec_matrices =
    read_cspec_matrices(targ_dir, site_name, freq);

  //===========================================================================

  do_fkanaly(targ_dir,
             x_coord,
             y_coord,
             cspec_matrices,
             density,
             density_phi,
             max_phv,
             min_phv,
             trig_mom_max);

  //===========================================================================
}

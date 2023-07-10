// SPDX-FileCopyrightText: Copyright 2023 Harusato Kimura
// SPDX-License-Identifier: MIT

#include "obs_site_group.hpp"
#include <filesystem>

baidu::core::ObsSiteGroup::ObsSiteGroup(std::string new_label,
                                        std::string new_targ_dir,
                                        int new_seg_len,
                                        char delimiter,
                                        std::string suffix)
{
  this->label = new_label;
  this->targ_dir = new_targ_dir;
  this->seg_len = new_seg_len;

  baidu::io::read_coordinate(
    this->targ_dir, this->x_coord, this->y_coord, this->site_name);
  this->comp = baidu::io::get_comp(targ_dir, this->site_name[0]);

  this->time = baidu::io::read_vector_line(
    this->targ_dir + this->site_name[0] + suffix, 0, delimiter);
  this->time_seg.resize(this->seg_len);
  for (int i = 0; i < this->seg_len; ++i) {
    this->time_seg[i] = this->time[i];
  }

  this->nsites = site_name.size();
  for (int n = 0; n < this->nsites; ++n) {
    std::vector<std::vector<double>> input = baidu::io::read_vector_multilines(
      targ_dir + site_name[n] + suffix, 1, this->comp, delimiter);
    baidu::core::ObsSite tmp_site{ site_name[n], seg_len, time, input };
    this->sites.push_back(tmp_site);
  }

  for (int n = 0; n < this->nsites; ++n) {
    this->sites[n].creat_segments();
    this->sites[n].calc_fft();
  }

  std::vector<std::string> comps{ "UD", "NS", "EW" };
  for (int i = 0; i < this->nsites; ++i) {
    for (int j = 0; j < this->comp; ++j) {
      std::string fname = targ_dir + "results/inputs/" + this->sites[i].name +
                          "_" + comps[j] + ".csv";
      baidu::io::write_vector(fname, ", ", this->time, sites[i].timeseries[j]);
    }
  }
}

// ----------------------------------------------------------------------------

void
baidu::core::ObsSiteGroup::write_segments()
{
  std::vector<std::string> comps{ "UD", "NS", "EW" };
  for (int i = 0; i < this->nsites; ++i) {
    for (int j = 0; j < this->comp; ++j) {
      std::string segdir =
        targ_dir + "results/segments/" + this->sites[i].name + "_" + comps[j];
      std::filesystem::create_directory(segdir);
      for (int k = 0; k < this->sites[i].num_of_segment; ++k) {
        std::stringstream fname;
        fname << segdir + "/";
        fname << std::setw(6) << std::setfill('0') << k;
        fname << ".csv";
        baidu::io::write_vector(
          fname.str(), ", ", this->time_seg, sites[i].segments[j][k]);
      }
    }
  }
}

// ----------------------------------------------------------------------------

void
baidu::core::ObsSiteGroup::write_fspectra()
{
  size_t fspec_size = this->sites[0].fspectra[0][0].size();
  double seg_total_time = this->seg_len * (this->time[1] - this->time[0]);
  std::vector<std::string> comps{ "UD", "NS", "EW" };
  for (int i = 0; i < this->nsites; ++i) {
    for (int j = 0; j < this->comp; ++j) {
      std::string segdir =
        targ_dir + "results/spectra/" + this->sites[i].name + "_" + comps[j];
      std::filesystem::create_directory(segdir);
      for (int k = 0; k < this->sites[i].num_of_segment; ++k) {
        std::stringstream fname;
        fname << segdir + "/";
        fname << std::setw(6) << std::setfill('0') << k;
        fname << ".csv";

        std::vector<double> amp(fspec_size);
        std::vector<double> phase(fspec_size);
        for (size_t l = 0; l < fspec_size; ++l) {
          amp[l] =
            std::abs(this->sites[i].fspectra[j][k][l]) * 4.0 * seg_total_time;
          // 4.0 = 2.0 * 2.0 (for minus part and hann window)
          phase[l] = std::arg(this->sites[i].fspectra[j][k][l]);
        }
        baidu::io::write_vector(
          fname.str(), ", ", this->sites[i].freq, amp, phase);
      }
    }
  }
}

// ----------------------------------------------------------------------------

void
baidu::core::ObsSiteGroup::scrutinize_segments(int allowance)
{
  if (std::filesystem::exists(this->targ_dir + "valid_segments.csv")) {
    std::vector<std::string> valseg = baidu::io::read_valid_seg(targ_dir);
    this->valid_seg_ind.resize(valseg.size());
    for (size_t i = 0; i < valseg.size(); ++i) {
      size_t pos = valseg[i].find_first_of(".");
      valid_seg_ind[i] = std::stoi(valseg[i].substr(0, pos));
    }
    return;
  } else {
    // calc normalized RMS
    for (size_t i = 0; i < this->sites.size(); ++i) {
      this->sites[i].calc_rms();
    }

    int n_class = 20;
    std::vector<std::vector<std::vector<std::pair<int, int>>>> histos;
    // site, comp, class (key, frequency i.e. dosuu)
    histos.resize(sites.size());
    for (size_t i = 0; i < sites.size(); ++i) {
      histos[i].resize(sites[i].comp);
      for (int j = 0; j < sites[i].comp; ++j) {
        histos[i][j].resize(n_class);
        for (int k = 0; k < n_class; ++k) {
          histos[i][j][k].first = k;
          histos[i][j][k].second = 0;
        }
      }
    }

    int nseg = this->sites[0].num_of_segment;
    // create histogram
    for (size_t i = 0; i < this->sites.size(); ++i) {
      for (int j = 0; j < this->sites[i].comp; ++j) {
        for (int k = 0; k < nseg; ++k) {
          // scan all the classes
          for (int l = 0; l < n_class; ++l) {
            if (0.05 + (double)l * 0.1 <= sites[i].n_rms_seg[j][k] &&
                0.05 + (double)(l + 1) * 0.1 > sites[i].n_rms_seg[j][k]) {
              histos[i][j][l].second += 1;
            }
          }
        }
      }
    }

    // sort histogram
    for (size_t i = 0; i < this->sites.size(); ++i) {
      for (int j = 0; j < this->sites[i].comp; ++j) {
        std::sort(histos[i][j].begin(),
                  histos[i][j].end(),
                  [](std::pair<int, int> o1, std::pair<int, int> o2) {
                    return o1.second > o2.second;
                  });
      }
    }

    // get index of valid segment
    this->valid_seg_ind.clear();
    for (int i = 0; i < nseg; ++i) {
      bool valid = true;
      for (size_t j = 0; j < this->sites.size(); ++j) {
        for (int k = 0; k < this->sites[j].comp; ++k) {
          if (std::abs(sites[j].n_rms_seg[k][i] -
                       ((double)histos[j][k][0].first * 0.1 + 0.1)) >
              allowance * 0.1 + 0.05) {
            valid = false;
          }
        }
      }
      if (valid == true) {
        this->valid_seg_ind.push_back(i);
      }
    }

    std::vector<std::string> valseg(this->valid_seg_ind.size());
    for (size_t i = 0; i < this->valid_seg_ind.size(); ++i) {
      std::stringstream fname;
      fname << std::setw(6) << std::setfill('0') << this->valid_seg_ind[i];
      fname << ".csv";
      valseg[i] = fname.str();
    }
    baidu::io::write_vector(targ_dir + "valid_segments.csv", valseg);
  }
}

// ----------------------------------------------------------------------------

void
baidu::core::ObsSiteGroup::calc_statistics(int nsmooth)
{
  if (this->valid_seg_ind.empty()) {
    perror("Error: No valid segment");
    exit(errno);
  }
  std::vector<std::string> comps{ "UD", "NS", "EW" };
  double seg_total_time = this->seg_len * (this->time[1] - this->time[0]);
  std::vector<double> zeros(this->sites[0].fspectra[0][0].size(), 0.0);

  this->pspec.resize(this->comp);
  for (int i = 0; i < this->comp; ++i) {
    for (size_t j = 0; j < this->sites.size(); ++j) {
      std::vector<double> pspec_j(this->sites[j].fspectra[i][0].size());
      for (int k : this->valid_seg_ind) {
        for (size_t l = 0; l < pspec_j.size(); ++l) {
          pspec_j[l] += (std::conj(this->sites[j].fspectra[i][k][l]) *
                         this->sites[j].fspectra[i][k][l])
                          .real();
        }
      }
      for (size_t k = 0; k < pspec_j.size(); ++k) {
        pspec_j[k] /= (double)this->sites[j].num_of_segment;
        pspec_j[k] *= seg_total_time / 2.0 * 3.0; // /8.0*3.0 : for hann window
      }
      baidu::util::do_smoothing(pspec_j, nsmooth);
      std::string label =
        this->targ_dir + "results/statistics/" + comps[i] + "_";
      label += this->site_name[j] + "-" + this->site_name[j] + ".csv";
      baidu::io::write_vector(label, ", ", this->sites[j].freq, pspec_j, zeros);
      this->pspec[i].push_back(pspec_j);
    }
  }

  this->cspec.resize(this->comp);
  for (int i = 0; i < this->comp; ++i) {
    for (size_t j = 1; j < this->sites.size(); ++j) {
      for (size_t k = 0; k < j; ++k) {
        // CSPEC
        std::vector<std::complex<double>> cspec_jk(
          this->sites[j].fspectra[i][0].size());
        for (int l : this->valid_seg_ind) {
          for (size_t m = 0; m < cspec_jk.size(); ++m) {
            cspec_jk[m] += std::conj(this->sites[j].fspectra[i][l][m]) *
                           this->sites[k].fspectra[i][l][m];
          }
        }
        for (size_t m = 0; m < cspec_jk.size(); ++m) {
          cspec_jk[m] /= (double)this->sites[j].num_of_segment;
          cspec_jk[m] *=
            seg_total_time / 2.0 * 3.0; // /8.0*3.0 : for hann window
        }
        baidu::util::do_smoothing(cspec_jk, nsmooth);
        std::string label_1 =
          this->targ_dir + "results/statistics/" + comps[i] + "_";
        std::string label_2 =
          this->targ_dir + "results/statistics/" + comps[i] + "_";
        label_1 += this->site_name[j] + "-" + this->site_name[k] + ".csv";
        label_2 += this->site_name[k] + "-" + this->site_name[j] + ".csv";
        baidu::io::write_vector(label_1, ", ", this->sites[0].freq, cspec_jk);
        baidu::io::write_vector_conj(
          label_2, ", ", this->sites[0].freq, cspec_jk);
        this->cspec[i].push_back(cspec_jk);
        // CCF
        for (size_t m = 0; m < cspec_jk.size(); ++m) {
          cspec_jk[m] /= std::sqrt(this->pspec[i][j][m] * this->pspec[i][k][m]);
        }
        label_1 = targ_dir + "results/statistics/CCF_" + comps[i] + "_";
        label_2 = targ_dir + "results/statistics/CCF_" + comps[i] + "_";
        label_1 += this->site_name[j] + "-" + this->site_name[k] + ".csv";
        label_2 += this->site_name[k] + "-" + this->site_name[j] + ".csv";
        baidu::io::write_vector(label_1, ", ", this->sites[0].freq, cspec_jk);
        baidu::io::write_vector_conj(
          label_2, ", ", this->sites[0].freq, cspec_jk);
      }
    }
  }
}

// ----------------------------------------------------------------------------

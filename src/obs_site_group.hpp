// SPDX-FileCopyrightText: Copyright 2023 Harusato Kimura
// SPDX-License-Identifier: MIT

#ifndef OBS_SITE_GROUP_HPP
#define OBS_SITE_GROUP_HPP

#include "io.hpp"
#include "obs_site.hpp"

namespace baidu::core {

class ObsSiteGroup
{
public:
  std::string label;
  std::string targ_dir;

  std::vector<ObsSite> sites;
  std::vector<std::string> site_name;
  std::vector<double> x_coord;
  std::vector<double> y_coord;
  int comp;
  int nsites;

  int seg_len;

  std::vector<int> valid_seg_ind;

  std::vector<double> time;
  std::vector<double> time_seg;
  std::vector<double> freq;

  std::vector<std::vector<std::vector<std::complex<double>>>> cspec;
  // comp, combination ij, freq
  std::vector<std::vector<std::vector<double>>> pspec;
  // comp, nsite, freq

  // ---
  ObsSiteGroup(std::string new_label,
               std::string new_targ_dir,
               int new_seg_len,
               char delimiter,
               std::string suffix);

  void write_segments();
  void write_fspectra();

  void scrutinize_segments(int allowance);

  void calc_statistics(int nsmooth);
};

}

#endif

// SPDX-FileCopyrightText: Copyright 2023 Harusato Kimura
// SPDX-License-Identifier: MIT

#ifndef OBS_SITE_HPP
#define OBS_SITE_HPP

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

#include "util.hpp"

namespace baidu::core {

class ObsSite
{
public:
  std::string name;

  int comp; // 1 (UD) or 3 (UD & NS & EW)

  int seg_len = 2048; // sample length for segmenting (default value is 2048)

  int total_sample;

  int num_of_segment; // = floor(total_sample / (0.5*seg_len) )

  // Time series without DC component
  std::vector<std::vector<double>> timeseries; // cmp, time
  std::vector<double> time;

  // RMS of timeseries
  std::vector<double> rms_all; // cmp

  // segments
  std::vector<std::vector<std::vector<double>>> segments; // cmp, seg_n, time
  std::vector<double> time_seg;

  // RMS of segment normalized by rms_all
  std::vector<std::vector<double>> n_rms_seg; // cmp, seg_n

  // fourier spectra corresponding to each segment
  std::vector<std::vector<std::vector<std::complex<double>>>>
    fspectra; // cmp, seg_n, freq
  std::vector<double> freq;

  // coordinate
  double x_coord; // m
  double y_coord; // m

  // --------------------------------------------------------------------------

  ObsSite(std::string new_name,
          int new_seg_len,
          std::vector<double> new_time,
          std::vector<std::vector<double>> new_timeseries);

  // --------------------------------------------------------------------------

  void remove_offset();

  void creat_segments();

  void calc_rms();

  double calc_rms_kernel(std::vector<double> const& vec);

  void calc_fft();
};

}

#endif

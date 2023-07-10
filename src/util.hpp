// SPDX-FileCopyrightText: Copyright 2022-2023 Harusato Kimura
// SPDX-License-Identifier: MIT

#ifndef UTIL_HPP
#define UTIL_HPP

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

namespace baidu::util {

void
apply_hann_window(std::vector<std::complex<double>>& timeseries);

void
FFTCore(std::vector<std::complex<double>>& vec, int length);

std::vector<std::complex<double>>
SortBitReverse(std::vector<std::complex<double>> in);

std::vector<std::complex<double>>
FFT(int const& total_sample,
    std::vector<std::complex<double>> const& time_series);

std::vector<std::complex<double>>
FFT(int const& total_sample, std::vector<double> const& time_series);

void
zero_padding(std::vector<std::complex<double>>& time_series);

void
creat_freq_vec(std::vector<double>& freq,
               int const& total_sample,
               double const& df);

void
do_smoothing(std::vector<std::complex<double>>& inspec, int const n_smooth);

void
do_smoothing(std::vector<double>& inspec, int const n_smooth);

}

#endif

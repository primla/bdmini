// SPDX-FileCopyrightText: Copyright 2022-2023 Harusato Kimura
// SPDX-License-Identifier: MIT

#include "util.hpp"

//-----------------------------------------------------------------------------

void
baidu::util::do_smoothing(std::vector<double>& inspec, int const n_smooth)
{
  size_t dim = inspec.size();
  for (int n = 0; n < n_smooth; ++n) {
    std::vector<double> tmp = inspec;
    for (size_t i = 1; i < dim - 1; ++i) {
      inspec[i] = 0.25 * (tmp[i - 1] + tmp[i + 1]) + 0.5 * tmp[i];
    }
  }
}

//-----------------------------------------------------------------------------

void
baidu::util::do_smoothing(std::vector<std::complex<double>>& inspec,
                          int const n_smooth)
{
  size_t dim = inspec.size();
  for (int n = 0; n < n_smooth; ++n) {
    std::vector<std::complex<double>> tmp = inspec;
    for (size_t i = 1; i < dim - 1; ++i) {
      inspec[i] = 0.25 * (tmp[i - 1] + tmp[i + 1]) + 0.5 * tmp[i];
    }
  }
}

//-----------------------------------------------------------------------------

void
baidu::util::creat_freq_vec(std::vector<double>& freq,
                            int const& total_sample,
                            double const& df)
{
  double fk = 0.0;
  freq.resize(total_sample / 2 + 1);
  for (int i = 0; i < total_sample / 2 + 1; ++i) {
    freq[i] = fk;
    fk += df;
  }
}

//-----------------------------------------------------------------------------

void
baidu::util::apply_hann_window(std::vector<std::complex<double>>& timeseries)
{
  std::vector<double> hann;
  int dim = (int)timeseries.size();
  hann.resize(dim);
  for (int i = 0; i < dim; ++i) {
    hann[i] = 0.5 * (1.0 - std::cos(2.0 * M_PI * (double)i / (dim - 1)));
  }

  for (int i = 0; i < dim; ++i) {
    timeseries[i] *= hann[i];
  }
}

//-----------------------------------------------------------------------------

std::vector<std::complex<double>>
baidu::util::FFT(int const& total_sample,
                 std::vector<std::complex<double>> const& time_series)
{
  std::vector<std::complex<double>> fourier_coeffs;
  fourier_coeffs.resize(total_sample);
  for (int i = 0; i < total_sample; ++i) {
    fourier_coeffs[i] = time_series[i];
  }

  // use hann window
  baidu::util::apply_hann_window(fourier_coeffs);

  FFTCore(fourier_coeffs, total_sample);
  for (int i = 0; i < total_sample; ++i) {
    fourier_coeffs[i] /= (double)total_sample;
  }
  return SortBitReverse(fourier_coeffs);
}

//-----------------------------------------------------------------------------

std::vector<std::complex<double>>
baidu::util::FFT(int const& total_sample,
                 std::vector<double> const& time_series)
{
  std::vector<std::complex<double>> fourier_coeffs;
  fourier_coeffs.resize(total_sample);
  for (int i = 0; i < total_sample; ++i) {
    std::complex<double> tmp{ time_series[i], 0.0 };
    fourier_coeffs[i] = tmp;
  }

  // use hann window
  baidu::util::apply_hann_window(fourier_coeffs);

  FFTCore(fourier_coeffs, total_sample);
  for (int i = 0; i < total_sample; ++i) {
    fourier_coeffs[i] /= (double)total_sample;
  }
  return SortBitReverse(fourier_coeffs);
}

//-----------------------------------------------------------------------------

void
baidu::util::FFTCore(std::vector<std::complex<double>>& vec, int length)
{
  std::vector<std::complex<double>> first;
  std::vector<std::complex<double>> last;
  std::vector<std::complex<double>> w;
  std::complex<double> imunit(0.0, 1.0);
  if (length > 1) {
    for (int i = 0; i < length / 2; ++i) {
      //w.push_back(std::exp(imunit * -2.0 * M_PI * (double)i / (double)length));
      w.push_back(std::exp(imunit * 2.0 * M_PI * (double)i / (double)length));
    }
    for (int i = 0; i < length / 2; ++i) {
      first.push_back((vec[i] + vec[i + length / 2]));
      last.push_back((vec[i] - vec[i + length / 2]) * w[i]);
    }
    vec.clear();
    FFTCore(first, length / 2);
    FFTCore(last, length / 2);
    for (int i = 0; i < length / 2; ++i) {
      vec.push_back(first[i]);
    }
    for (int i = 0; i < length / 2; ++i) {
      vec.push_back(last[i]);
    }
  } else {
    return;
  }
}

//-----------------------------------------------------------------------------

void
baidu::util::zero_padding(std::vector<std::complex<double>>& time_series)
{
  size_t dim = time_series.size();
  int log2_paddedsize;
  for (int i = 1; i < 30; ++i) {
    if (std::pow(2, i) >= dim) {
      log2_paddedsize = i;
      break;
    }
  }
  int padsize = std::pow(2, log2_paddedsize) - dim;
  std::complex<double> zeros{ 0.0, 0.0 };
  for (int i = 0; i < padsize; ++i) {
    time_series.push_back(zeros);
  }
}

//-----------------------------------------------------------------------------

std::vector<std::complex<double>>
baidu::util::SortBitReverse(std::vector<std::complex<double>> in)
{
  std::vector<int> p;
  int log2_N;
  log2_N = log2(in.size());
  int log2_N_is_odd;
  int m_max;
  if (log2_N % 2 == 1) {
    m_max = (log2_N - 1) / 2;
    log2_N_is_odd = 0;
  } else {
    m_max = log2_N / 2;
    log2_N_is_odd = -1;
  }
  int L;
  L = std::pow(2, m_max);
  p.push_back(0);
  for (int m = 0; m < m_max; ++m) {
    for (int j = 0; j < std::pow(2, m); ++j) {
      p.push_back(std::pow(2, log2_N - m - 1) + p[j]);
    }
  }

  std::complex<double> tmp;
  // log2_N_is_odd
  if (log2_N_is_odd == 0) {
    for (int j = 0; j < L; ++j) {
      for (int k = 0; k < j; ++k) {
        tmp = in[j + p[k]];
        in[j + p[k]] = in[k + p[j]];
        in[k + p[j]] = tmp;
        tmp = in[j + L + p[k]];
        in[j + L + p[k]] = in[k + L + p[j]];
        in[k + L + p[j]] = tmp;
      }
    }
    // log2_N_is_even
  } else {
    for (int j = 0; j < L; ++j) {
      for (int k = 0; k < j; ++k) {
        tmp = in[j + p[k]];
        in[j + p[k]] = in[k + p[j]];
        in[k + p[j]] = tmp;
      }
    }
  }
  return in;
}

//-----------------------------------------------------------------------------

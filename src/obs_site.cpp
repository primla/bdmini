// SPDX-FileCopyrightText: Copyright 2023 Harusato Kimura
// SPDX-License-Identifier: MIT

#include "obs_site.hpp"

//-----------------------------------------------------------------------------

baidu::core::ObsSite::ObsSite(std::string new_name,
                              int new_seg_len,
                              std::vector<double> new_time,
                              std::vector<std::vector<double>> new_timeseries)
{
  this->name = new_name;
  this->comp = new_timeseries.size();

  this->seg_len = new_seg_len;

  // Adjust total sample to smallest one
  std::vector<int> sample_len(this->comp);
  for (int i = 0; i < this->comp; ++i) {
    sample_len[i] = new_timeseries[i].size();
  }
  std::sort(sample_len.begin(), sample_len.end());
  this->total_sample = sample_len[0];

  this->num_of_segment =
    std::floor((double)total_sample / std::floor(0.5 * (double)seg_len)) - 1;

  // set time vec
  this->time.resize(this->total_sample);
  for (int i = 0; i < this->total_sample; ++i) {
    this->time[i] = new_time[i];
  }
  this->time_seg.resize(this->seg_len);
  for (int i = 0; i < this->seg_len; ++i) {
    this->time_seg[i] = this->time[i];
  }

  // set timeseries
  this->timeseries.resize(comp);
  for (int i = 0; i < this->comp; ++i) {
    this->timeseries[i].resize(this->total_sample);
    for (int j = 0; j < this->total_sample; ++j) {
      this->timeseries[i][j] = new_timeseries[i][j];
    }
  }

  remove_offset();
}

//-----------------------------------------------------------------------------

void
baidu::core::ObsSite::remove_offset()
{
  // Get DC component
  std::vector<double> ave{ 0, 0, 0 };

  for (int i = 0; i < this->comp; ++i) {
    for (int j = 0; j < total_sample; ++j) {
      ave[i] += this->timeseries[i][j];
    }
    ave[i] /= (double)this->total_sample;
  }

  // Romove offset
  for (int i = 0; i < this->comp; ++i) {
    for (int j = 0; j < this->total_sample; ++j) {
      this->timeseries[i][j] -= ave[i];
    }
  }
}

//-----------------------------------------------------------------------------

void
baidu::core::ObsSite::creat_segments()
{
  this->segments.resize(this->comp);
  for (int i = 0; i < this->comp; ++i) {
    this->segments[i].resize(this->num_of_segment);
    for (int j = 0; j < this->num_of_segment; ++j) {
      this->segments[i][j].resize(this->seg_len);
      for (int k = 0; k < this->seg_len; ++k) {
        this->segments[i][j][k] =
          this->timeseries[i][j * std::floor(seg_len * 0.5) + k];
      }
    }
  }
}

//-----------------------------------------------------------------------------

void
baidu::core::ObsSite::calc_rms()
{
  // rms all
  this->rms_all.resize(this->comp);
  for (int i = 0; i < this->comp; ++i) {
    this->rms_all[i] = calc_rms_kernel(this->timeseries[i]);
  }

  // normalized rms of seg
  this->n_rms_seg.resize(this->comp);
  for (int i = 0; i < this->comp; ++i) {
    this->n_rms_seg[i].resize(this->num_of_segment);
    for (int j = 0; j < this->num_of_segment; ++j) {
      n_rms_seg[i][j] = calc_rms_kernel(this->segments[i][j]);
      n_rms_seg[i][j] /= this->rms_all[i];
    }
  }
}

//-----------------------------------------------------------------------------

double
baidu::core::ObsSite::calc_rms_kernel(std::vector<double> const& vec)
{
  double ave = 0;
  for (size_t i = 0; i < vec.size(); ++i) {
    ave += vec[i];
  }
  ave /= (double)vec.size();

  double ans = 0;
  for (size_t i = 0; i < vec.size(); ++i) {
    ans += (vec[i] - ave) * (vec[i] - ave);
  }
  ans /= (double)vec.size();
  return std::sqrt(ans);
}

//-----------------------------------------------------------------------------

void
baidu::core::ObsSite::calc_fft()
{
  double dt = this->time_seg[1] - this->time_seg[0];
  double df = 1.0 / (double)this->seg_len / dt;
  baidu::util::creat_freq_vec(this->freq, this->seg_len, df);

  // apply fft for each segment
  this->fspectra.resize(this->comp);
  for (int i = 0; i < this->comp; ++i) {
    this->fspectra[i].resize(this->num_of_segment);
    for (int j = 0; j < this->num_of_segment; ++j) {
      this->fspectra[i][j] =
        baidu::util::FFT(this->seg_len, this->segments[i][j]);
    }
  }
}

//-----------------------------------------------------------------------------

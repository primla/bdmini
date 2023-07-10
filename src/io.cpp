// SPDX-FileCopyrightText: Copyright 2022-2023 Harusato Kimura
// SPDX-License-Identifier: MIT

#include "io.hpp"

//-----------------------------------------------------------------------------

void
baidu::io::read_coordinate(std::string const targ_dir,
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

std::vector<std::string>
baidu::io::read_valid_seg(std::string const targ_dir)
{
  std::ifstream inf{ targ_dir + "valid_segments.csv" };
  std::vector<std::string> valid_seg;
  std::string line;
  while (std::getline(inf, line)) {
    valid_seg.push_back(line);
  }
  return valid_seg;
}

//-----------------------------------------------------------------------------

int
baidu::io::get_comp(std::string const targ_dir, std::string site0)
{
  std::ifstream inf{ targ_dir + site0 + ".csv" };
  int comp = 0;
  std::string line, val;
  std::getline(inf, line);
  std::stringstream iss{ line };
  while (std::getline(iss, val, ',')) {
    comp++;
  }
  return comp - 1;
}

//-----------------------------------------------------------------------------

std::vector<double>
baidu::io::read_vector_line(std::string const path, int nl, char delimiter)
{
  std::vector<double> ans;
  std::ifstream inf{ path };
  std::string line, val;
  while (std::getline(inf, line)) {
    std::stringstream iss(line);
    for (int i = 0; i < nl; ++i) {
      std::getline(iss, val, delimiter);
    }
    std::getline(iss, val, delimiter);
    double val_d = 0.0;
    try {
      val_d = std::stod(val);
    } catch (std::out_of_range& e) {
      // do nothing i.e. set val_d = 0
    }
    ans.push_back(val_d);
  }
  return ans;
}

//-----------------------------------------------------------------------------

std::vector<std::vector<double>>
baidu::io::read_vector_multilines(std::string const path,
                                  int il,
                                  int nl,
                                  char delimiter)
{
  std::vector<std::vector<double>> ans;
  ans.resize(nl);
  std::ifstream inf{ path };
  std::string line, val;
  while (std::getline(inf, line)) {
    std::stringstream iss(line);
    for (int i = 0; i < il; ++i) {
      std::getline(iss, val, delimiter);
    }
    double val_d;
    for (int i = 0; i < nl; ++i) {
      std::getline(iss, val, delimiter);
      val_d = 0;
      try {
        val_d = std::stod(val);
      } catch (std::out_of_range& e) {
        // do nothing i.e. set val_d = 0
      }
      ans[i].push_back(val_d);
    }
  }
  return ans;
}

//-----------------------------------------------------------------------------

void
baidu::io::write_vector(std::string const path,
                        std::vector<std::string> const& v1)
{
  std::ofstream ouf{ path };
  for (size_t i = 0; i < v1.size(); ++i) {
    ouf << v1[i] << std::endl;
  }
}

//-----------------------------------------------------------------------------

void
baidu::io::write_vector(std::string const path,
                        std::string const separator,
                        std::vector<double> const& v1,
                        std::vector<double> const& v2)
{
  std::ofstream ouf{ path };
  std::vector<size_t> len{ v1.size(), v2.size() };
  std::sort(len.begin(), len.end());
  for (size_t i = 0; i < len[0]; ++i) {
    ouf << v1[i] << separator;
    ouf << v2[i] << std::endl;
  }
}

void
baidu::io::write_vector(std::string const path,
                        std::string const separator,
                        std::vector<double> const& v1,
                        std::vector<double> const& v2,
                        std::vector<double> const& v3)
{
  std::ofstream ouf{ path };
  std::vector<size_t> len{ v1.size(), v2.size(), v3.size() };
  std::sort(len.begin(), len.end());
  for (size_t i = 0; i < len[0]; ++i) {
    ouf << v1[i] << separator;
    ouf << v2[i] << separator;
    ouf << v3[i] << std::endl;
  }
}

void
baidu::io::write_vector(std::string const path,
                        std::string const separator,
                        std::vector<double> const& v1,
                        std::vector<std::complex<double>> const& v2)
{
  std::ofstream ouf{ path };
  std::vector<size_t> len{ v1.size(), v2.size() };
  std::sort(len.begin(), len.end());
  for (size_t i = 0; i < len[0]; ++i) {
    ouf << v1[i] << separator;
    ouf << v2[i].real() << separator;
    ouf << v2[i].imag() << std::endl;
  }
}

void
baidu::io::write_vector_conj(std::string const path,
                             std::string const separator,
                             std::vector<double> const& v1,
                             std::vector<std::complex<double>> const& v2)
{

  std::ofstream ouf{ path };
  std::vector<size_t> len{ v1.size(), v2.size() };
  std::sort(len.begin(), len.end());
  for (size_t i = 0; i < len[0]; ++i) {
    ouf << v1[i] << separator;
    ouf << v2[i].real() << separator;
    ouf << -v2[i].imag() << std::endl;
  }
}
//-----------------------------------------------------------------------------

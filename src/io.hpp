// SPDX-FileCopyrightText: Copyright 2022-2023 Harusato Kimura
// SPDX-License-Identifier: MIT

#ifndef IO_HPP
#define IO_HPP

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

namespace baidu::io {

//-----------------------------------------------------------------------------

void
read_coordinate(std::string const targ_dir,
                std::vector<double>& x_coord,
                std::vector<double>& y_coord,
                std::vector<std::string>& site_name);

//-----------------------------------------------------------------------------

std::vector<std::string>
read_valid_seg(std::string const targ_dir);

//-----------------------------------------------------------------------------

int
get_comp(std::string const targ_dir, std::string site0);

//-----------------------------------------------------------------------------

std::vector<double>
read_vector_line(std::string const path, int nl, char delimiter);

//-----------------------------------------------------------------------------

std::vector<std::vector<double>>
read_vector_multilines(std::string const path, int il, int nl, char delimiter);

//-----------------------------------------------------------------------------

void
write_vector(std::string const path, std::vector<std::string> const& v1);

//-----------------------------------------------------------------------------

void
write_vector(std::string const path,
             std::string const separator,
             std::vector<double> const& v1,
             std::vector<double> const& v2);

//-----------------------------------------------------------------------------

void
write_vector(std::string const path,
             std::string const separator,
             std::vector<double> const& v1,
             std::vector<double> const& v2,
             std::vector<double> const& v3);

//-----------------------------------------------------------------------------

void
write_vector(std::string const path,
             std::string const separator,
             std::vector<double> const& v1,
             std::vector<std::complex<double>> const& v2);

//-----------------------------------------------------------------------------

void
write_vector_conj(std::string const path,
                  std::string const separator,
                  std::vector<double> const& v1,
                  std::vector<std::complex<double>> const& v2);

//-----------------------------------------------------------------------------
}
#endif

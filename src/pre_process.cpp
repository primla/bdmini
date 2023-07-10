// SPDX-FileCopyrightText: Copyright 2022-2023 Harusato Kimura
// SPDX-License-Identifier: MIT

#include "obs_site_group.hpp"

int
main(int argc, char** argv)
{
  if (argc != 5) {
    perror("Error: Invalid argument");
    exit(errno);
  }

  std::string targ_dir = argv[1];
  int seg_len = std::stoi(argv[2]);
  std::string label = argv[3];
  int nsmooth = std::stoi(argv[4]);
  char delimiter = ',';

  baidu::core::ObsSiteGroup osg{ label, targ_dir, seg_len, delimiter, ".csv" };

  osg.scrutinize_segments(2);
  osg.calc_statistics(nsmooth);
}

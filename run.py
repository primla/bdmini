#!/usr/bin/python3

# SPDX-FileCopyrightText: Copyright 2022-2023 Harusato Kimura
# SPDX-License-Identifier: MIT

# ==============================================================================

import sys
import os
import shutil
import json
import math
import subprocess
import numpy as np

# ==============================================================================

path2paramfile = sys.argv[1]
targ_dir = os.path.dirname(path2paramfile) + "/"
case_name = os.path.splitext(os.path.basename(path2paramfile))[0]

res_dir = targ_dir + "results"
if os.path.exists(res_dir):
    shutil.rmtree(res_dir)

# creating a directory to restore results
os.makedirs(res_dir)

# creating subdirectories
subdirs = [
    "inputs",
    "statistics",
    "spac",
    "fk",
    "dspac",
]
for i in range(0, len(subdirs), 1):
    os.makedirs(res_dir + "/" + subdirs[i])

# ==============================================================================

# load params.json
json_open = open(path2paramfile, "r")
params = json.load(json_open)
seg_len = int(params["seg_len"])
if "n_smoothing" in params:
    nsmooth = str(params["n_smoothing"])
else:
    nsmooth = "0"

# ==============================================================================

print("Pre process:")

# pre_process: segmenting, scrutinizing, fft, calc PSD, CSD, and CCF
argv = " " + targ_dir
argv += " " + str(params["seg_len"])
argv += " " + case_name
argv += " " + nsmooth
subprocess.run(["./bin/pre_process" + argv], shell=True)

print("Pre process done.")

print("")
print("------------------------------")
print("")

# ==============================================================================

print("SPAC analysis:")

# SPAC analysis
if "SPAC" in params:
    if "arrays" in params["SPAC"]:
        for j in range(0, len(params["SPAC"]["arrays"])):
            array_name = params["SPAC"]["arrays"][j]
            array_info = " "
            for j in range(0, len(params["SPAC"][array_name])):
                array_info += params["SPAC"][array_name][j]
                array_info += " "
            argv = targ_dir + " " + array_name + array_info
            subprocess.run(["./bin/spac " + argv], shell=True)

print("SPAC analysis done.")

print("")
print("------------------------------")
print("")

# ==============================================================================

print("FK analysis:")

# FK analysis (capon)
if "FK" in params:
    if "bounds" in params["FK"] and "density" in params["FK"]:
        argv = targ_dir + " "
        argv += str(params["FK"]["density"][0]) + " "
        argv += str(params["FK"]["density"][1]) + " "
        argv += str(params["FK"]["bounds"][0]) + " "
        argv += str(params["FK"]["bounds"][1])
        subprocess.run(["./bin/fkanaly " + argv], shell=True)

print("FK analy done.")

print("")
print("------------------------------")
print("")

# ==============================================================================

print("DSPAC analysis:")
print("This may take hours.")

if "DSPAC" in params:
    calc_param_pso = " " + str(params["DSPAC"]["n_particle"])
    calc_param_pso += " " + str(params["DSPAC"]["n_itr"])
    calc_param_pso += " " + str(params["DSPAC"]["w4loc"])
    calc_param_pso += " " + str(params["DSPAC"]["w4glo"])
    array_conf = ""
    for sname in params["DSPAC"]["array"]:
        array_conf += " " + sname
    subprocess.run(
        ["./bin/dspac " + targ_dir + calc_param_pso + array_conf], shell=True
    )


print("DSPAC analy done.")

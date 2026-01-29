#!/bin/bash

targ=document

lualatex $targ.tex
biber $targ
lualatex $targ.tex
lualatex $targ.tex
evince $targ.pdf &

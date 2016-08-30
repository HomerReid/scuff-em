#!/bin/bash

ARGS="-setnumber L 1.0 -setnumber N 10"
gmsh ${ARGS} -2 Square_N.geo -o Square.msh; RenameMesh Square.msh

ARGS="-setnumber L 1.0 -setnumber N 20"
gmsh ${ARGS} -2 Square_N.geo -o Square.msh; RenameMesh Square.msh

#!/bin/bash

DIR=.
if [ $# == 1 ]
then
  DIR=$1
fi

if [ ! -r ${DIR}/Fresnel/TestFresnel.sh ]
then
  echo "Usage: $0 [scuff_tests_path]"
  echo " where scuff_tests_path defaults /usr/local/share/scuff-em/tests"
  exit
fi

cd ${DIR}/Fresnel
TestFresnel.sh

#!/bin/bash

echo "the start time is"
date

###############################################################################

inputDir="../input/homogeneousStrain"
caseName="homogeneousStrain"          # was "channel"

###############################################################################

runCase () {

    rm -rf "../data/$caseName" > /dev/null 2>&1
    mkdir  "../data/$caseName"
    mkdir  "../data/$caseName/data"
    mkdir  "../data/$caseName/input"
    mkdir  "../data/$caseName/runtime"
    cp     "$inputDir/"*        "../data/$caseName/input/" > /dev/null 2>&1
    cp -r  "$inputDir/restart"* "../data/$caseName/input/" > /dev/null 2>&1

    echo "*** RUNNING ***"
    echo "Output is being written to ../data/$caseName/runtime and ../data/$caseName/data"
    ./odt.x $caseName 0
}

###############################################################################

rebuild () {
  echo '*** REBUILDING ***'
  cd ../build
  make -j8
  if [ $? -ne 0 ] ; then
    echo ; echo 'FATAL: error in the build' ; echo
    exit 0
  fi
  echo '*** DONE REBUILDING ***'
  cd ../run
}

###############################################################################

if [ "$1" == "-r" ]; then rebuild; fi

runCase "$caseName"

echo
echo "the end simulation time is"
date

exit 0
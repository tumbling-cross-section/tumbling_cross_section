#!/bin/bash


../src/crossArea -n 2000 -g 2.5 -i miniTest.pqr -r ../lib/sc_radii.lib -v > testC.log 2>testC.elog 


../src/crossArea_old -a 80 -g 2.5 -i miniTest.pqr -r ../lib/sc_radii.lib -v  > testC_old.log



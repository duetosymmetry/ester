#!/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
BINDIR=${DIR}/../../bin
ESTER=${BINDIR}/ester

MASSES=(02.0 02.5 03.0 03.5 04.0 04.5 05.0 05.5 06.0 06.5 07.0 07.5 08.0 08.5 09.0 09.5 10.0)
OMEGAS=(0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9)

for mass in ${MASSES[*]}; do
  mkdir -p M-${mass}
  ${ESTER} 1d -M ${mass} -o M-${mass}/1d.dat > M-${mass}/1d.log 2>&1
  for omega in ${OMEGAS[*]}; do
    ${ESTER} 2d -i M-${mass}/1d.dat -Omega_bk ${omega} -o M-${mass}/Om-${omega}.dat > M-${mass}/Om-${omega}.log 2>&1
 done
done

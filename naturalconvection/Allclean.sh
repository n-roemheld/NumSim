#! /usr/bin/env bash

fluid_solver="numsim_fluid"
solid_solver="numsim_solid"

fluidsolver_name="FluidSolver"
solidsolver_name="SolidSolver"

rm -rf "precice-run/"
rm -rf "out_${fluidsolver_name}" "out_${solidsolver_name}"
rm -f "${fluid_solver}.log" "${solid_solver}.log" 
rm -f precice-*.log
rm -f precice-*.json


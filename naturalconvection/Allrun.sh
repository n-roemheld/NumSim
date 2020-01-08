#! /usr/bin/env bash

fluid_input="naturalconvection_fluid.input"
solid_input="naturalconvection_solid.input"

fluid_solver="numsim_fluid"
solid_solver="numsim_solid"

rm -rf "precice-run/"

./${fluid_solver} ${fluid_input} > ${fluid_solver}.log 2>&1 &
PIDFluid=$!
./${solid_solver} ${solid_input} > ${solid_solver}.log 2>&1 &
PIDSolid=$!

echo "Waiting for the participants to exit..."
echo "(you may run 'tail -f ${fluid_solver}.log' or 'tail -f ${solid_solver}.log' in another terminal to check the progress)"

wait ${PIDFluid}
wait ${PIDSolid}

if [ $? -ne 0 ] || [ "$(grep -c -E "error:" ${fluid_solver}.log)" -ne 0 ] || [ "$(grep -c -E "error:" ${solid_solver}.log)" -ne 0 ]; then
    echo ""
    echo "Something went wrong... See the log files for more."
else
    echo ""
    echo "The simulation completed!"
fi


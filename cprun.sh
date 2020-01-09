
testcase="naturalconvection"
#testcase="forcedconvection"
#testcase="dummy"


mkdir -p run

if [ "$testcase" == "forcedconvection" ]
then 
    cp -a ./forcedconvection/. ./run
    cd run 
    ./Allclean.sh
    cp ../build/src/numsim_fluid .
    cp ../build/srcSolid/numsim_solid .
    ./Allrun.sh
fi
if [ "$testcase" == "naturalconvection" ]
then 
    cp -a ./naturalconvection/. ./run
    cd run 
    ./Allclean.sh
    cp ../build/src/numsim_fluid .
    cp ../build/srcSolid/numsim_solid .
    ./Allrun.sh


fi
if [ "$testcase" == "dummy" ]
then
    cp ./ldchv/SolidSolverDummySrc/SolidSolverDummy ./run
    cp -a ./ldchv/. ./run
    cd run
    ./Allclean.sh
    cp ../build/src/numsim_fluid .

    fluid_input="ldchp.input"
    solid_input="precice-config.xml"

    fluid_solver="numsim_fluid"
    solid_solver="SolidSolverDummy"

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
fi

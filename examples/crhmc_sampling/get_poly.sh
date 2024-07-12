#!/bin/bash

# dim facets 0=skinny, 1=order
./simple_crhmc $1 $2 $3
python3 to_matlab.py $1 $3
mv matlab_poly.mat ../../../../idk/PolytopeSamplerMatlab/matlab_poly.mat
cd ../../../../idk/PolytopeSamplerMatlab/
let "a = $1 * 60 / 4 + 1"
matlab -batch "run(run_thing($a))"
mv matlab_results.txt ../../volesti_dev/volesti_fork/examples/crhmc_sampling/matlab_results.txt
cd ../../volesti_dev/volesti_fork/examples/crhmc_sampling


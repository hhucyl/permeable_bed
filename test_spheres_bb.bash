#! /bin/bash

echo "bbtype, tau, Nproc"
for tau  in '2.0' '1.8' '1.6' '1.4' '1.2' '1.0' '0.8' '0.6'
do
    ./test_spheres_bb 3 200 $tau 12 
done
for tau  in '2.0' '1.8' '1.6' '1.4' '1.2' '1.0' '0.8' '0.6'
do
    # for h in '10' '30' '50' '70' '90' '110'
    # do
        
        ./test_spheres_bb 0 50 $tau 12
        ./test_spheres_bb 1 50 $tau 12
        ./test_spheres_bb 2 50 $tau 12
        ./test_spheres_bb 3 50 $tau 12
        ./test_spheres_bb 4 50 $tau 12
        # ./test_spheres_ga -1 $h $tau 12
        # ./test_spheres_ibm -2 $h $tau 12
        # ./test_spheres_ga -1 30 $tau 12
        # ./test_spheres_ibm -2 30 $tau 12
        # ./test_spheres_ibm -3 30 $tau 12
        # ./test_spheres_ga -4 30 $tau 12
    # done
done
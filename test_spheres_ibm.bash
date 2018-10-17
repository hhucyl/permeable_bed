#! /bin/bash

echo "bbtype, tau, Nproc"

for tau  in '2.0' '1.8' '1.6' '1.4' '1.2' '1.0' '0.8' '0.6'
do
    # for h in '10' '30' '50' '70' '90' '110'
    # do
        
        ./test_spheres_ibm 0 90 $tau 12
        ./test_spheres_ibm 1 90 $tau 12
        ./test_spheres_ibm 2 90 $tau 12
        ./test_spheres_ibm 3 90 $tau 12
        ./test_spheres_ibm 4 90 $tau 12
        ./test_spheres_ibm 0 90 $tau 12
        ./test_spheres_ibm 1 90 $tau 12
        ./test_spheres_ibm 2 90 $tau 12
        ./test_spheres_ibm 3 90 $tau 12
        ./test_spheres_ibm 4 90 $tau 12
        # ./test_spheres_ga -1 $h $tau 12
        # ./test_spheres_ibm -2 $h $tau 12
        # ./test_spheres_ga -1 30 $tau 12
        # ./test_spheres_ibm -2 30 $tau 12
        # ./test_spheres_ibm -3 30 $tau 12
        # ./test_spheres_ga -4 30 $tau 12
    # done
done
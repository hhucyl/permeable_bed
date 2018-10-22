#! /bin/bash

echo "bbtype, tau, Nproc"

for tau  in '2.0' '1.8' '1.6' '1.4' '1.2' '1.0' '0.8' '0.6'
do
    # for h in '10' '30' '50' '70' '90' '110'
    # do
        
        ./test_spheres_ga -1 90 $tau 12 1
        ./test_spheres_ga1 -4 90 $tau 12 1
        ./test_spheres_ga -1 90 $tau 12 0
        ./test_spheres_ga1 -4 90 $tau 12 0
    # done
done

 for h in '30' '50' '70' '90' '110' '130'
     do
        
        ./test_spheres_ga -1 $h 1.6 12 1
        ./test_spheres_ga1 -4 $h 1.6 12 1
        ./test_spheres_ga -1 $h 1.6 12 0
        ./test_spheres_ga1 -4 $h 1.6 12 0
     done

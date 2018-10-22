#! /bin/bash

echo "bbtype, tau, Nproc"

for tau  in '2.0' '1.8' '1.6' '1.4' '1.2' '1.0' '0.8' '0.6'
do
    # for h in '10' '30' '50' '70' '90' '110'
    # do
        
        ./test_spheres_ibm -2 90 $tau 12 1
        ./test_spheres_ibm_s -2 90 $tau 12 0
	./test_spheres_ibm1 -3 90 $tau 12 1
        ./test_spheres_ibm_s1 -3 90 $tau 12 0 
    # done
done

for h in '30' '50' '90' '110' '130' '200'
do
	./test_spheres_ibm -2 $h 1.6 12 1
	./test_spheres_ibm_s -2 $h 1.6 12 0
	./test_spheres_ibm1 -3 $h 1.6 12 1
	./test_spheres_ibm_s1 -3 $h 1.6 12 0
done

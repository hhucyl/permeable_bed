#! /bin/bash

echo "bbtype, tau, Nproc"
#for tau  in '2.0' '1.8' '1.6' '1.4' '1.2' '1.0' '0.8' '0.6'
#do
#    ./test_spheres_bb 3 200 $tau 12 
#done
for tau  in '2.0' '1.8' '1.6' '1.4' '1.2' '1.0' '0.8'
do
    # for h in '10' '30' '50' '70' '90' '110'
    # do
        
       ./test_spheres_bb 0 90 $tau 12 >> log.txt 2>&1
       ./test_spheres_bb 1 90 $tau 12 >> log.txt 2>&1
       ./test_spheres_bb 2 90 $tau 12 >> log.txt 2>&1
       ./test_spheres_bb 3 90 $tau 12 >> log.txt 2>&1
       ./test_spheres_bb 4 90 $tau 12 >> log.txt 2>&1
        
    # done
done
#  for h in '30' '50' '70' '90' '110' '130' '200'
#      do
        
#     #    ./test_spheres_bb 0 $h 1.6 12
#         ./test_spheres_bb 1 $h 1.6 12
#   #      ./test_spheres_bb 2 $h 1.6 12
#  #       ./test_spheres_bb 3 $h 1.6 12
# #        ./test_spheres_bb 4 $h 1.6 12
        
#      done

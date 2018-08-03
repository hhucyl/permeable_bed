#! /bin/bash
# make test_sc_bb
# make test_sc_ga
echo "bbtype, h, tau, Nproc"
for tau  in '2.0' '1.8' '1.6' '1.4' '1.2' '1.0' '0.8' '0.6' '0.53' '0.51' '0.503'
do
    ./test_sc_bb 0 30 $tau 12 >> log.txt 2>&1
    ./test_sc_bb 1 30 $tau 12 >> log.txt 2>&1
    ./test_sc_bb 2 30 $tau 12 >> log.txt 2>&1
    ./test_sc_bb 3 30 $tau 12 >> log.txt 2>&1
    ./test_sc_bb 4 30 $tau 12 >> log.txt 2>&1
    ./test_sc_ga -1 30 $tau 12 >> log.txt 2>&1
    ./test_sc_ibm -2 30 $tau 12 >> log.txt 2>&1
done

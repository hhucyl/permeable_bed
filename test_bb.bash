#! /bin/bash
make test_bb
echo "bbtype, tau, Nproc"
./test_bb 1 0.8 12 >> log.txt 2>&1
./test_bb 2 0.8 12 >> log.txt 2>&1
./test_bb 3 0.8 12 >> log.txt 2>&1
./test_bb 4 0.8 12 >> log.txt 2>&1
# ./test_bb 1 1.0 12
# ./test_bb 2 1.0 12
# ./test_bb 3 1.0 12
# ./test_bb 4 1.0 12
# ./test_bb 1 1.2 12
# ./test_bb 2 1.2 12
# ./test_bb 3 1.2 12
# ./test_bb 4 1.2 12
# ./test_bb 1 0.6 12
# ./test_bb 2 0.6 12
# ./test_bb 3 0.6 12
# ./test_bb 4 0.6 12
# ./test_bb 1 0.53 12
# ./test_bb 2 0.53 12
# ./test_bb 3 0.53 12
# ./test_bb 4 0.53 12
# ./test_bb 1 0.51 12
# ./test_bb 2 0.51 12
# ./test_bb 3 0.51 12
# ./test_bb 4 0.51 12
# ./test_bb 1 1.4 12
# ./test_bb 2 1.4 12
# ./test_bb 3 1.4 12
# ./test_bb 4 1.4 12
# ./test_bb 1 1.6 12
# ./test_bb 2 1.6 12
# ./test_bb 3 1.6 12
# ./test_bb 4 1.6 12
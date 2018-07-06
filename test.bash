#! /bin/bash
make test_cd_bb
echo "collide type, bbtype, h, tau, Nproc"
./test_cd_bb 1 1 100 0.8 12
./test_cd_bb 1 2 100 0.8 12
./test_cd_bb 1 3 100 0.8 12
./test_cd_bb 0 1 100 0.8 12
./test_cd_bb 0 2 100 0.8 12
./test_cd_bb 0 3 100 0.8 12
./test_cd_bb 0 0 100 0.8 12
./test_cd_bb 1 0 100 0.8 12
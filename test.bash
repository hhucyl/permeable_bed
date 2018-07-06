#! /bin/bash
make test_cd_bb
echo "collide type, bbtype, h, tau, Nproc"
./test_cd_bb 1 1 100 0.8 12

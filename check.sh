g++-14 -O2 -Wall himenoBMTxpa.cpp    -o cpp.out
gcc-14 -O2 -Wall himenoBMTxpa.c      -o c.out
gcc-14 -O2 -Wall himenoBMTxpa-copy.c -o copy.out
echo "\nC Language:"
./c.out    L
echo "\nC++ Language:"
./cpp.out  L
echo "\nC with copy:"
./copy.out L


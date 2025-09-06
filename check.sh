#CC=gcc-14
#CXX=g++-14
CC=clang-21
CXX=clang-21
FLAGS="-Wall -O2"
$CXX $FLAGS himenoBMTxpa.cpp    -o cpp.out
$CC  $FLAGS himenoBMTxpa.c      -o c.out
$CC  $FLAGS himenoBMTxpa-copy.c -o copy.out
echo "========================================\nC Language:"
./c.out    M
echo "========================================\nC++ Language:"
./cpp.out  M
echo "========================================\nC with copy:"
./copy.out M


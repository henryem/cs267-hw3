cat pgen-p*t${1}.out | sort | diff -q serial-sorted.out /dev/stdin

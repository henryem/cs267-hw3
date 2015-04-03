cat pgen-p*t${1}.out | sort | diff serial-sorted.out /dev/stdin

#!/bin/bash

# ./testit.sh PROGRAM OPTS PROCESSES TESTNUM

mkdir -p tmptest/

make $1 > tmptest/maketmp 2>&1;

CMD="mpiexec -n $3 ./$1 $2"

if [[ ! -f output/$1.test$4 ]]; then
    echo "FAIL: Test #$4 of $1"
    echo "       command = '$CMD'"
    echo "       OUTPUT MISSING"

else

    $CMD > tmptest/tmp

    diff output/$1.test$4 tmptest/tmp > tmptest/difftmp

    if [[ -s tmptest/difftmp ]] ; then
       echo "FAIL: Test #$4 of $1  ($5)"
       echo "       command = '$CMD'"
       echo "       diffs follow:"
       cat tmptest/difftmp
    else
       echo "PASS: Test #$4 of $1  ($5)"
    fi

fi


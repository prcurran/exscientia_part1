#!/bin/sh
#   run_tests.sh    -   run tests from within a distribution

for t in test_*.py ; do
    python $t
done
# Actual development corrections

The first local SymPy run stopped at its removable-line check: the test had
interchanged the signs of V on the two lines. Direct endpoint substitution
gives V=L at t=alpha and V=-L at s=alpha. The two assertions were corrected;
the candidate D formula and rectangle source were not changed. No result
JSON had been written by that failed run. Subsequent successful records
refer to the corrected script hash.

The shell's default sandbox helper intermittently failed at process setup.
The same local reads and bounded CAS execution were run with the reviewed
execution mode. No repository or toolchain configuration was changed.

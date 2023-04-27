To compile the code use the following commands:

```
  gcc trap.c -o traps.exe
  mpicc mpi_trap0.c -o trap0.exe
  mpicc mpi_trap4_do.c -o trap.exe
```

And to run the programs use the following commands:

```
  ./traps.exe
  mpirun --use-hwthread-cpus --oversubscribe -np 4 ./trap0.exe
  mpirun --use-hwthread-cpus --oversubscribe -np 4 ./trap.exe
```

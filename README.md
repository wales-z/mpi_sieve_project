# MPI sieve experiment

## Brief Description

This is the code for distributed and parallel computing class's MPI experiment in UESTC (电子科技大学，分布式并行计算MPI实验)。

In this project, I implemented sieve of Eratosthenes based on MPI, with 3 step improvement. Also, experiment's data are included.

## 3-step improvement

step1, remove all even numbers before sieve.

step2, remove MPI communication code.

step3, redesign the loop in order to get higher cache hit-rate.

## Code Files

base.cpp: The original version, implementing sieve of Eratosthenes in a serial way

optimizer1.cpp: corresponding to improvement step1

optimizer2.cpp: corresponding to improvement step2

optimizer3.cpp: corresponding to improvement step3
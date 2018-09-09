# Adaptive Model-based Transfer for MOO

This repository contains the Matlab implementation of our work [Curbing Negative Influences Online for Seamless Transfer Evolutionary Optimization](https://www.researchgate.net/publication/326846571_Curbing_Negative_Influences_Online_for_Seamless_Transfer_Evolutionary_Optimization).

The `runZDT.m` file contains examples on ZDT test sets. When solving any one of the problems in the test sets using the AMTEA, all the other problems act as source tasks providing source probabilistic models for transfer.

To run the code:
```matlab
addpath('../')
runZDT
```
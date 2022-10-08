# Accumulation Network for Circular Economy
The code in this repository allows to analyse Circular Economy with a Linear Network Flow Dynamics approach. CE is strongly related to reuse, share, repaire and recycle materials. This circularity may lead to the accumulation of some unwanted compounds at some stages of the cycle that may harm human health or the environment on the long run. A specific code is developped for checking whether this accumulation occurs and veryfing if it is robust with a Monte Carlo approach.

The main folder contains all functions needed for measuring accumulation and robustness. 

The Aanalysis(A, ver; graph=false) function takes as inputs:
  - A that is the normalised weight matrix of the graph
  - ver that is a list containing the names of all vertices
  - graph that is a Boolean variable. By default is set to false and no graph representation is provided in the output. If the graph representation is needed, set it to true.

This function print as outputs the list of all accumulation vertices and their respective net value.

The Mcanalysis(A, ver; iter = 1000, var = 0.1, plot = true, deb = false, mdvar = true) function takes as inputs:
  - A that is the normalised weight matrix of the graph
  - ver that is a list containing the names of all vertices
  - iter that is the number of changes done by the algorithm in a single simulation
  - var that is the default variance of the truncated normal distribution
  - plot that is a Boolean variable. By default is set to true and this allows to plot the accumulation of a single simulation
  - deb that allows the debug when set true
  - mdvar that allows to decide whether the variance is related (true) or not (false) to the mean value (big values can change more than smaller one) of the truncated normal       distribution
  
The MC(A, ver; itmc = 10000, itr = 1000, var = 0.1, deb = false, mdvar = true, savef = true) function takes as inputs:
  - A that is the normalised weight matrix of the graph
  - ver that is a list containing the names of all vertices
  - itmc that is the number of single MC simulation
  - itr that is the number of changes done by the algorithm in a single MC simulation
  - var that is the default variance of the truncated normal distribution
  - deb that allows the debug when set true
  - mdvar that allows to decide whether the variance is related (true) or not (false) to the mean value (big values can change more than smaller one) of the truncated normal       distribution
  - saverf that is a Boolean variable. By default is set to true and this allows to plot the mean values and the variance of the average net accumulation of all the itmc           simulations

This function print as outputs the control matrices (a matrix that states how much is changed on average the matrix. All entries are basically equal due to Law of large number) and all the plots for veryfing accumulation.


# Erasmus Genetic Programming - Problem

A problem in Erasmus is defined by input data, expected output data for the input data and a fitness function. The fitness function calculates a score between 0.0 and 1.0 for an individual GC that was tasked with producing the expected output from the input data. Because Erasmus generates executable code that could, in theory, be manipulated to have adverse effects the input data, expected output and fitness function are defined by a problem definition hash. A change to any component (even one that does not change the fitness score) results in a different problem hash.

Verified problems 
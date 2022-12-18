# Dimension_Function
Computation of the dimension of the Markov Spectrum

The file treewalk contains fortran code that computes a set of forbidden words up to length n for a given value s and writes them in a file setup????.day.  To run treewalk.out 3.005 10 1000. It gives a file fruits1000.dat which contains the set of forbidden words, alphabet, number of forbidden words, and their length. This has to be preprocessed for the next step 
(head -n 2 fruits1000.dat ; tail -n 3 fruits1000.dat | head -n 2; tail -n+3 fruits1000.dat | head -n -3; tail -n 1 fruits1000.dat) > setup1000.day 

The file dimsetup uses the data file setup1000.day to compute a Markov dynamical system whose limit set is the set of continued fractions corresponding the to set of forbidden words from the file setup1000.day. To run: dimsetup.out 1000. Finally, the function computedim.c computes the Hausdorff dimension of the limit set. 

The shell script has been used to compute the dimension function on the second interval (3.0117, 3.0171). It took a few hours on a 24 threads machine. 

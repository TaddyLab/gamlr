A python wrapper for Matt Taddy's C code in the gamlr package.

To do:

1. Decide on what type of python-c linkage	
DONE: ctypes

2. Familiarize self with ctypes:	
DONE: github.com/nelsonauner/ctypes

3. Familiarize self with gamlr.R wrapper	
DONE: I included error checking, but some options (precalc of XX) is left out

4. Simple python wrapper for gamlr.c
TODO:
- Error: cannot find Rmath.h (even though I copied it from Rcpp)..this shouldn't be necessary?
- Goal: Maintain original C code
- Numerous examples for passing double * and int * to c function


5. Error-checking in python wrapper	
- Run same example in R & Python and check result

6. Additional features missing in original wrapper	
- Precalculation of XX

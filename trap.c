/* File:     mpi_trap1.c
 * Purpose:  Use MPI to implement a parallel version of the trapezoidal 
 *           rule.  In this version the endpoints of the interval and
 *           the number of trapezoids are hardwired.
 *
 * Input:    None (hardcoded).
 * Output:   Estimate of the integral from a to b of f(x)
 *           using the trapezoidal rule and n trapezoids.
 *
 * Compile:  mpicc mpi_trap1.c -o mpi_trap1
 * Run:      mpirun -n <number of processes> ./mpi_trap1
 *
 * Algorithm:
 *    1.  Each process calculates "its" interval of
 *        integration.
 *    2.  Each process estimates the integral of f(x)
 *        over its interval using the trapezoidal rule.
 *    3a. Each process != 0 sends its integral to 0.
 *    3b. Process 0 sums the calculations received from
 *        the individual processes and prints the result.
 *
 */
#include <stdio.h>
#include <time.h>
/* We'll be using MPI routines, definitions, etc. */

/* Calculate local integral  */
double Trap(double left_endpt, double right_endpt, long int trap_count, 
   double base_len);    

/* Function we're integrating */
double f(double x); 

int main(void) {
   long int n = 1316134912;   
   double a = -100000000000, b = 100000000000, h;
   double total_int;

   clock_t begin = clock();
   h = (b-a)/n;          /* h is the same for all processes */
   total_int = Trap(a, b, n, h);
   clock_t end = clock();

   double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
   printf("Time to make operations: %lf\n", time_spent);
   /* Print the result */
   printf("With n = %ld trapezoids, our estimate\n", n);
   printf("of the integral from %f to %f = %.15f\n", a, b, total_int);

   return 0;
} /*  main  */


/*------------------------------------------------------------------
 * Function:     Trap
 * Purpose:      Serial function for estimating a definite integral 
 *               using the trapezoidal rule
 * Input args:   left_endpt
 *               right_endpt
 *               trap_count 
 *               base_len
 * Return val:   Trapezoidal rule estimate of integral from
 *               left_endpt to right_endpt using trap_count
 *               trapezoids
 */
double Trap(
      double left_endpt  /* in */, 
      double right_endpt /* in */, 
      long int    trap_count  /* in */, 
      double base_len    /* in */) {
   double estimate, x; 
   int i;

   estimate = (f(left_endpt) + f(right_endpt))/2.0;
   for (i = 1; i <= trap_count-1; i++) {
      x = left_endpt + i*base_len;
      estimate += f(x);
   }
   estimate = estimate*base_len;

   return estimate;
} /*  Trap  */


/*------------------------------------------------------------------
 * Function:    f
 * Purpose:     Compute value of function to be integrated
 * Input args:  x
 */
double f(double x) {
   return x*x;
} /* f */

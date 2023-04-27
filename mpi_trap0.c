#include <stdio.h>
#include <mpi.h>

#define A 1.0
#define B 100.0
#define N 2048

double Trap(double left_endpt, double right_endpt, int trap_count, double base_len);    

double f(double x); 

void Get_input(int my_rank, int comm_sz, double* a_p, double* b_p, int* n_p);

int main(void) {
   int my_rank, comm_sz, source, local_n, n; 
   double h, local_a, local_b, local_int, total_int, a, b;

   double tstart, tend;
   double p_a, p_b;
   int p_n; 

   

   /* Let the system do what it needs to start up MPI */
   MPI_Init(NULL, NULL);
   MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
   MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);

   if (my_rank == 0) {
      printf("Enter a, b and n:\n");
      scanf("%lf %lf %d", &p_a, &p_b, &p_n);
      tstart = MPI_Wtime();

      printf("[Master]: Sending messages...\n");
      for (int dest = 1; dest < comm_sz; dest++) {
         MPI_Send(&p_a, 1, MPI_DOUBLE, dest, 0, MPI_COMM_WORLD);
         MPI_Send(&p_b, 1, MPI_DOUBLE, dest, 0, MPI_COMM_WORLD);
         MPI_Send(&p_n, 1, MPI_INT, dest, 0, MPI_COMM_WORLD);
      }
      printf("[Master]: Sent all messages\n");
   } else {
      MPI_Recv(&p_a, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(&p_b, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(&p_n, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

      printf("[Worker %i]: Recived all messages as: a=%lf b=%lf n=%d \n", my_rank, p_a, p_b, p_n);
   }

   a = p_a;
   b = p_b;
   n = p_n;

   h = (b-a)/n;
   local_n = n/comm_sz;

   local_a = a + my_rank*local_n*h;
   local_b = local_a + local_n*h;
   local_int = Trap(local_a, local_b, local_n, h);

   /* Add up the integrals calculated by each process */
   if (my_rank != 0) { 
      MPI_Send(&local_int, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD); 
   } else {
      total_int = local_int;
      for (source = 1; source < comm_sz; source++) {
         MPI_Recv(&local_int, 1, MPI_DOUBLE, source, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
         total_int += local_int;
      }
   } 
   tend = MPI_Wtime();

   /* Print the result */
   if (my_rank == 0) {
      printf("\nTook %f s to run\n", (tend-tstart));
      printf("With n = %d trapezoids, our estimate\n", n);
      printf("of the integral from %f to %f = %.15f\n", a, b, total_int);
   }

   /* Shut down MPI */
   MPI_Finalize();

   return 0;
} /*  main  */

void Get_input (
   int my_rank,
   int comm_sz,
   double* p_a,
   double* p_b,
   int* p_n
) {
   int dest;

   if (my_rank == 0) {
      printf("Enter a:\n");
      scanf("%lf", p_a);
      printf("Enter b:\n");
      scanf("%lf", p_b);
      printf("Enter n:\n");
      scanf("%d", p_n);

      printf("[Master]: Sending messages...\n");
      // for (int dest = 1; dest < comm_sz; dest++) {
      //    MPI_Send(p_a, 1, MPI_DOUBLE, dest, 0, MPI_COMM_WORLD);
      //    MPI_Send(p_b, 1, MPI_DOUBLE, dest, 0, MPI_COMM_WORLD);
      //    MPI_Send(p_n, 1, MPI_INT, dest, 0, MPI_COMM_WORLD);
      // }
      printf("[Master]: Sent all messages\n");
   } else {
      // MPI_Recv(p_a, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      // MPI_Recv(p_b, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      // MPI_Recv(p_n, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

      printf("[Worker %i]: Recived all messages\n", my_rank);
      // a = *p_a;
      // b = *p_b;
      // n = *p_n;
   }

}

double Trap(
   double left_endpt  /* in */, 
   double right_endpt /* in */, 
   int    trap_count  /* in */, 
   double base_len    /* in */
) {
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

double f(double x) { return x*x; }

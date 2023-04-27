#include <stdio.h>
#include <mpi.h>

/* Construir tipo derivado */
void Build_mpi_type(double* a_p, double* b_p, int* n_p,
   MPI_Datatype* input_mpi_t_p);

/* Valores de entrada */
void Get_input(int my_rank, int comm_sz, double* tstart, double* a_p, double* b_p,
   int* n_p);

/* Funcion de Trapezoides  */
double Trap(double left_endpt, double right_endpt, int trap_count, 
   double base_len);    

/* Funcion f(x) a integrar */
double f(double x); 

int main(void) {
   int my_rank, comm_sz, n, local_n;   
   double a, b, h, local_a, local_b;
   double local_int, total_int;
   double tstart, tend;

   MPI_Init(NULL, NULL);

   MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

   MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);

   Get_input(my_rank, comm_sz, &tstart, &a, &b, &n);

   h = (b-a)/n;          
   local_n = n/comm_sz; 

   /* calculo de valores locales*/
   local_a = a + my_rank*local_n*h;
   local_b = local_a + local_n*h;
   local_int = Trap(local_a, local_b, local_n, h);

   /* reducir resultados */
   MPI_Reduce(&local_int, &total_int, 1, MPI_DOUBLE, MPI_SUM, 0,
         MPI_COMM_WORLD);
   tend = MPI_Wtime();

   /* resultado */
   if (my_rank == 0) {
      printf("\nTook %f s to run\n", (tend-tstart));
      printf("With n = %d trapezoids, our estimate\n", n);
      printf("of the integral from %f to %f = %.15e\n", a, b, total_int);
   }

   /* Shut down MPI */
   MPI_Finalize();

   return 0;
} /*  main  */

//-----------------------------------------------------------

void Build_mpi_type(
      double*        a_p            /* in  */, 
      double*        b_p            /* in  */, 
      int*           n_p            /* in  */,
      MPI_Datatype*  input_mpi_t_p  /* out */) {

   /*
    * Los enunciados están en desorden.  Colóquelos en el orden
    * correcto para que Build_mpi_type cree el tipo derivado
    */
   MPI_Aint array_of_displacements[3] = {0};
   MPI_Aint a_addr, b_addr, n_addr;

   MPI_Get_address(a_p, &a_addr);
   MPI_Get_address(b_p, &b_addr);
   MPI_Get_address(n_p, &n_addr);

   array_of_displacements[1] = n_addr-b_addr;
   array_of_displacements[2] = a_addr-b_addr; 

   int array_of_blocklengths[3] = {1, 1, 1};
   MPI_Datatype array_of_types[3] = {MPI_DOUBLE, MPI_DOUBLE, MPI_INT};
   
   MPI_Type_create_struct(3, array_of_blocklengths, 
         array_of_displacements, array_of_types,
         input_mpi_t_p);
   MPI_Type_commit(input_mpi_t_p);
   
}  /* Build_mpi_type */

//------------------------------------------------ -100000000000 100000000000 10000000000000

void Get_input(
      int      my_rank  /* in  */, 
      int      comm_sz  /* in  */,
      double*   tstart  /* out */,
      double*  a_p      /* out */, 
      double*  b_p      /* out */,
      int*     n_p      /* out */) {
   MPI_Datatype input_mpi_t;

   Build_mpi_type(a_p, b_p, n_p, &input_mpi_t);

   if (my_rank == 0) {
      printf("Enter a, b, and n\n");
      scanf("%lf %lf %d", a_p, b_p, n_p);
      *tstart = MPI_Wtime();
   }

   MPI_Bcast(b_p, 1, input_mpi_t, 0, MPI_COMM_WORLD);
   MPI_Type_free(&input_mpi_t);
}  /* Get_input */

//---------------------------------------------

double Trap(
      double left_endpt  /* in */, 
      double right_endpt /* in */, 
      int    trap_count  /* in */, 
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

double f(double x /* in */) {
   return x*x;
} /* f */

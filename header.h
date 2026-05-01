#ifndef HEADER
#define HEADER

#include <memory>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <fenv.h>
#include <fstream>
#include <mutex>
#include <pthread.h>
#include <sched.h>
#include <string>
#include <sys/resource.h>
#include <sys/sysinfo.h>
#include <sys/time.h>
#include <time.h>
#include <unistd.h>

int feenableexcept (int excepts);
int fedisableexcept (int excepts);
int fegetexcept (void);

//main
void reduce_sum (int p, int *a = nullptr, int n = 0);
void reduce_sum(int p, double *a, int n);
void *thread_fun (void *pa);
double get_full_time ();
double get_cpu_time ();

//fill_msr_matrix
void get_my_rows (int nx, int ny, int p, int k, int &i1, int &i2);
void ij_to_l (int nx, int /*ny*/, int i, int j, int &l);
void l_to_ij (int nx, int ny, int &i, int &j, int l);

void fill_I_diag (int nx, int ny, double hx, double hy, int *I);
int fill_IA (int nx, int ny, double hx, double hy,
             int *I, double *A, int p, int k);
int IA_ij (int nx, int ny, double hx, double hy, int i, int j,
           int is, int js, int s, int */*I*/, double *A);
int get_diag (int nx, int ny, double hx, double hy,
              int i, int j, int */*I*/, double *A);
int get_off_diag (int nx, int ny, double hx, double hy,
                  int i, int j, int *I, double *A);
int get_len_msr (int nx, int ny);
int get_len_msr_off_diag (int nx, int ny, double hx, double hy, int p, int k);

//fill_right_side
void fill_right_side (int nx, int ny, double a, double b, double c, double d,
                      double *B, int p, int k, int func_id,
                      double (*f) (int, double, double));
double F_ij (int nx, int ny, double hx, double hy,
             double x0, double y0, int i, int j,
             int func_id, double (*f) (int, double, double));

//residual
double f (int k, double x, double y);
double approximation_fun (int nx, int ny, double hx, double hy,
                          double *x, double *y,
                          double x0, double y0, double a, double c,
                          double *solution);
void calculate_residual (int nx, int ny, double hx, double hy,
                         double *x, double *y, double a, double c,
                         double &r1, double &r2, double &r3, double &r4,
                         int func_id, double (*f) (int, double, double),
                         double *solution, int p, int k);

//check_alg
int check_row_sum (int nx, int ny, double hx, double hy,
                   int *I, double *A, int p, int k);
int check_symm (int nx, int ny, double hx, double hy,
                int *I, double *A, int p, int k);
int get_triangles (int nx, int ny, int i, int j);

//solve_msr_system
void mult_msr_matrix_vector (const double *A, const int *I, int n,
                             const double *x, double *b, int p, int k);
double scalar_product (int n, const double *x, const double *y,
                       int p, int k, double *sp);
void mult_sub_vector (int n, double *x, const double *y,
                      double t, int p, int k);
double norm_vector (int n, const double *b, int p, int k, double *sp);
void thread_rows (int n, int p, int k, int &i1, int &i2);

int minimal_errors_msr_martix_full (int n, const double *A, const int *I,
                                    const double *b, double *x, double *r,
                                    double *u, double *v, double eps,
                                    int maxit, int max_steps, int p, int k, double *sp);
int minimal_errors_msr_martix (int n, const double *A, const int *I,
                               const double *b, double *x, double *r,
                               double *u, double *v, double eps,
                               int maxit, int p, int k, double *sp);

void apply_preconditioner_msr_matrix (int n, const double *A, const int *I,
                                      double *v, double *u, double *r, int p, int k);

class my_args_t
{
public:
  int nx      = 0;
  int ny      = 0;
  int func_id = 0;
  int mi      = 0;
  int p       = 0;
  int k       = 0;
  int it      = 0;

  double a   = 0;
  double b   = 0;
  double c   = 0;
  double d   = 0;
  double eps = 0;
  double t1  = 0;
  double t2  = 0;

  double r1 = -1;
  double r2 = -1;
  double r3 = -1;
  double r4 = -1;

  double *x = nullptr;
  double *y = nullptr;
  double *A = nullptr;
  int    *I = nullptr;
  double *solution = nullptr;
  double *right_side = nullptr;
  double *u = nullptr;
  double *v = nullptr;
  double *r = nullptr;
  double *sp = nullptr;
};

#endif
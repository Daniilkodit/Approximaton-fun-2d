#include "../header.h"

void *
thread_fun (void *pa)
{
  my_args_t *my_args = (my_args_t *)pa;

  int nx = my_args->nx;
  int ny = my_args->ny;
  int func_id = my_args->func_id;
  int mi = my_args->mi;
  int p = my_args->p;
  int k = my_args->k;
  int it = my_args->it;

  double a = my_args->a;
  double b = my_args->b;
  double c = my_args->c;
  double d = my_args->d;
  double eps = my_args->eps;

  double r1 = my_args->r1;
  double r2 = my_args->r2;
  double r3 = my_args->r3;
  double r4 = my_args->r4;

  double *x = my_args->x;
  double *y = my_args->y;
  double *A = my_args->A;
  double *right_side = my_args->right_side;
  double *solution = my_args->solution;
  int *I = my_args->I;
  double *r = my_args->r;
  double *u = my_args->u;
  double *v = my_args->v;
  double *sp = my_args->sp;

  cpu_set_t cpu;
  CPU_ZERO (&cpu);

  int n_cpus = get_nprocs ();
  int cpu_id = n_cpus - 1 - (k % n_cpus);

  CPU_SET (cpu_id, &cpu);
  pthread_t tid = pthread_self ();
  pthread_setaffinity_np (tid, sizeof (cpu_set_t), &cpu);

  int N = (nx + 1) * (ny + 1);
  double hx = (b - a) / nx;
  double hy = (d - c) / ny;
  //fill msr matrix
  if (k == 0)
    fill_I_diag (nx, ny, hx, hy, I);

  reduce_sum (p);
  int err = fill_IA (nx, ny, hx, hy, I, A, p, k);
  if (err < 0)
    {
      if (k == 0)
        printf ("It's bug in fill_IA err = %d\n", err);
      return nullptr;
    }

  reduce_sum (p);
  fill_right_side (nx, ny, a, b, c, d,
                   right_side, p, k, func_id, f);

  err = check_row_sum (nx, ny, hx, hy, I, A, p, k);
  if (err < 0)
    {
      if (k == 0)
        printf ("It's bug in check_row_sum err = %d\n", err);
      return nullptr;
    }

  err = check_symm (nx, ny, hx, hy, I, A, p, k);
  if (err < 0)
    {
      if (k == 0)
        printf ("It's bug in check_symm err = %d\n", err);
      return nullptr;
    }

  // solve msr system
  int max_steps = 100000;
  double t1 = get_full_time ();
  it = minimal_errors_msr_martix_full (N, A, I, right_side, solution,
                                          r, u, v, eps, mi, max_steps, p, k, sp);
  t1 = get_full_time () - t1;

  if (it < 0)
    {
      if (k == 0)
        printf ("Number of steps exceeded max_steps = %d\n", max_steps);
      return nullptr;
    }

  // calculation residual
  reduce_sum (p);
  double t2 = get_full_time ();
  calculate_residual (nx, ny, hx, hy, x, y, a, c,r1, r2, r3, r4, func_id, f, solution, p, k);
  t2 = get_full_time () - t2;
  reduce_sum (p);

  my_args->r1 = r1;
  my_args->r2 = r2;
  my_args->r3 = r3;
  my_args->r4 = r4;
  my_args->t1 = t1;
  my_args->t2 = t2;
  my_args->it = it;

  return nullptr;
}

#include "../header.h"

int
main (int argc, char **argv)
{
  int task = 6, nx, ny, func_id, mi, p, it = 0;
  double a, b, c, d, eps;
  double r1 = -1, r2 = -1, r3 = -1, r4 = -1;
  //pthread_barrier_t barrier;

  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW | FE_UNDERFLOW);

  if (!(argc == 11
        && sscanf (argv[1], "%lf",  &a)   == 1
        && sscanf (argv[2], "%lf",  &b)   == 1
        && sscanf (argv[3], "%lf",  &c)   == 1
        && sscanf (argv[4], "%lf",  &d)   == 1
        && sscanf (argv[5],  "%d",  &nx)  == 1
        && sscanf (argv[6],  "%d",  &ny)  == 1
        && sscanf (argv[7],  "%d",  &func_id)   == 1
        && sscanf (argv[8],  "%lf", &eps) == 1
        && sscanf (argv[9],  "%d",  &mi)  == 1
        && sscanf (argv[10], "%d",  &p)   == 1))
    {
      printf ("Usage: a b c d nx ny func_id eps mi p\n");
      return 0;
    }

  if (!(a < b && c < d && nx > 0 && ny > 0
        && func_id >= 0 && eps > 0 && mi > 0 && p > 0))
    {
      printf ("Usage corrected data\n");
      return 0;
    }
  func_id %= 8;
  // in thread
  int N = (nx + 1) * (ny + 1); // всего функций куранта
  int len_matrix = N + 1 + get_len_msr (nx, ny); // кол-во не 0 скалярных произведений + 1;

  std::unique_ptr<double[]> x(new double[nx + 1]());
  std::unique_ptr<double[]> y(new double[ny + 1]());
  std::unique_ptr<double[]> A(new double[len_matrix]());
  std::unique_ptr<int[]>    I(new int[len_matrix]());
  std::unique_ptr<double[]> right_side(new double[N]());
  std::unique_ptr<double[]> solution(new double[N]());
  std::unique_ptr<double[]> u(new double[N]());
  std::unique_ptr<double[]> v(new double[N]());
  std::unique_ptr<double[]> r(new double[N]());
  std::unique_ptr<double[]> sp(new double[p]());
  double hx = (b - a) / nx;
  double hy = (d - c) / ny;
  for (int i = 0; i < nx + 1; i++)
    x[i] = a + i * hx;
  for (int j = 0; j < ny + 1; j++)
    y[j] = c + j * hy;

  my_args_t *ap = new my_args_t[p];
  for (int i = 0; i < p; i++)
    {
      ap[i].nx = nx;
      ap[i].ny = ny;
      ap[i].func_id = func_id;
      ap[i].mi = mi;
      ap[i].p = p;
      ap[i].k = i;
      ap[i].it = it;

      ap[i].a = a;
      ap[i].b = b;
      ap[i].c = c;
      ap[i].d = d;
      ap[i].eps = eps;

      ap[i].r1 = r1;
      ap[i].r2 = r2;
      ap[i].r3 = r3;
      ap[i].r4 = r4;

      ap[i].x = x.get ();
      ap[i].y = y.get ();
      ap[i].A = A.get ();
      ap[i].I = I.get ();
      ap[i].u = u.get ();
      ap[i].v = v.get ();
      ap[i].r = r.get ();
      ap[i].sp = sp.get ();
      ap[i].right_side = right_side.get ();
      ap[i].solution = solution.get ();
    }

  pthread_t *tid = new pthread_t[p];
  for (int i = 1; i < p; i++)
    {
      if(pthread_create (tid + i, 0, thread_fun ,ap + i))
        {
          printf("Thread %d not create\n",i);
          return -1;
        }
    }
  thread_fun (ap+0);

  double t1 = ap[0].t1, t2 = ap[0].t2;
  r1 = ap[0].r1;
  r2 = ap[0].r2;
  r3 = ap[0].r3;
  r4 = ap[0].r4;
  it = ap[0].it;

  printf (
  "%s : Task = %d R1 = %e R2 = %e R3 = %e R4 = %e T1 = %.2f T2 = %.2f\
 It = %d E = %e K = %d Nx = %d Ny = %d P = %d\n",
   argv[0], task, r1, r2, r3, r4, t1, t2, it, eps, func_id, nx, ny, p);
  return 0;
}

#include "../header.h"

#define F(IS, JS, S) (IA_ij (nx, ny, hx, hy, i, j, (IS), (JS), (S), I, A))

int
fill_IA (int nx, int ny, double hx, double hy,
         int *I, double *A, int p, int k)
{
  int l1, l2;
  int N = (nx + 1) * (ny + 1);
  int err = 0, len = 0;
  get_my_rows (nx, ny, p, k, l1, l2);

  for (int l = l1; l < l2; l++)
    {
      int i, j;
      l_to_ij (nx, ny, i, j, l);
      if (get_diag (nx, ny, hx, hy, i, j, I, A + l) < 0)
        {
          err = -1;
          break;
        }
      int r = I[l];
      int s = I[l + 1] - I[l];
      int t = get_off_diag (nx, ny, hx, hy, i, j, I + r, A + r);
      if (s != t)
        {
          err = -2;
          break;
        }
      len += s;
    }

  reduce_sum (p, &err, 1);
  if (err < 0)
    return -1;

  reduce_sum (p, &len, 1);
  if (I[N] != (N + 1) + len)
    return -2;

  return 0;
}

int
get_diag (int nx, int ny, double hx, double hy,
          int i, int j, int */*I*/, double *A)
{
  return IA_ij (nx, ny, hx, hy, i, j, i, j, 0, nullptr, A);
}

int
get_off_diag (int nx, int ny, double hx, double hy,
              int i, int j, int *I, double *A)
{
  if (i > 0 && i < nx && j > 0 && j < ny)
    {
      if (I && A)
        {
          F (i + 1, j    , 0);
          F (i    , j - 1, 1);
          F (i - 1, j - 1, 2);
          F (i - 1, j    , 3);
          F (i    , j + 1, 4);
          F (i + 1, j + 1, 5);
        }
      return 6;
    }

  if (j == 0 && i > 0 && i < nx)
    {
      if (I && A)
        {
          F (i + 1, j    , 0);
          F (i - 1, j    , 1);
          F (i    , j + 1, 2);
          F (i + 1, j + 1, 3);
        }
      return 4;
    }

  if (j == ny && i > 0 && i < nx)
    {
      if (I && A)
        {
          F (i + 1, j    , 0);
          F (i    , j - 1, 1);
          F (i - 1, j - 1, 2);
          F (i - 1, j    , 3);
        }
      return 4;
    }

  if (i == 0 && j > 0 && j < ny)
    {
      if (I && A)
        {
          F (i + 1, j    , 0);
          F (i    , j - 1, 1);
          F (i    , j + 1, 2);
          F (i + 1, j + 1, 3);
        }
      return 4;
    }

  if (i == nx && j > 0 && j < ny)
    {
      if (I && A)
        {
          F (i    , j - 1, 0);
          F (i - 1, j - 1, 1);
          F (i - 1, j    , 2);
          F (i    , j + 1, 3);
        }
      return 4;
    }

  if (i == 0 && j == 0)
    {
      if (I && A)
        {
          F (i + 1, j    , 0);
          F (i    , j + 1, 1);
          F (i + 1, j + 1, 2);
        }
      return 3;
    }

  if (i == nx && j == ny)
    {
      if (I && A)
        {
          F (i    , j - 1, 0);
          F (i - 1, j - 1, 1);
          F (i - 1, j    , 2);
        }
      return 3;
    }

  if (i == 0 && j == ny)
    {
      if (I && A)
        {
          F (i + 1, j    , 0);
          F (i    , j - 1, 1);
        }
      return 2;
    }

  if (i == nx && j == 0)
    {
      if (I && A)
        {
          F (i - 1, j    , 0);
          F (i    , j + 1, 1);
        }
      return 2;
    }

  return -99999999;
}

// does not parallelize
void
fill_I_diag (int nx, int ny, double hx, double hy, int *I)
{
  int N = (nx + 1) * (ny + 1);
  int r = N + 1;
  for (int l = 0; l < N; l++)
    {
      int i, j;
      l_to_ij (nx, ny, i, j, l);
      int s = get_off_diag (nx, ny, hx, hy, i, j, nullptr, nullptr);
      I[l] = r;
      r += s;
    }
  I[N] = r;
}

int
get_len_msr (int nx, int ny)
{
  return (nx - 1) * (ny - 1) * 6
          + 2 * 4 * (nx - 1)
          + 2 * 4 * (ny - 1)
          + 2 * 3 + 2 * 2;
}

// it's check correct get_of_diag
int
get_len_msr_off_diag (int nx, int ny, double hx, double hy, int p, int k)
{
  int i1, i2;
  i1 = (nx + 1) * k;
  i2 = (nx + 1) * (k + 1);
  i1 /= p;
  i2 /= p;
  int res = 0;
  for (int i = i1 ; i < i2; i++)
    {
      for (int j = 0; j < ny; j++)
        {
          res += get_off_diag (nx, ny, hx, hy, i, j, nullptr, nullptr);
        }
    }
  reduce_sum (p, &res, 1);
  return res;
}

void
get_my_rows (int nx, int ny, int p, int k, int &i1, int &i2)
{
  i1 = k * (nx + 1) * (ny + 1);
  i2 = (k + 1) * (nx + 1) * (ny + 1);
  i1 /= p;
  i2 /= p;
}

void
ij_to_l (int nx, int /*ny*/, int i, int j, int &l)
{
  l = i + j * (nx + 1);
}

void
l_to_ij (int nx, int /*ny*/, int &i, int &j, int l)
{
  j = l / (nx + 1);
  i = l - j * (nx + 1);
}
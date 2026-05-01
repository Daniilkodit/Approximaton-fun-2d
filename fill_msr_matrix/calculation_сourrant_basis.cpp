#include "../header.h"

int
IA_ij (int nx, int ny, double hx, double hy, int i, int j,
       int is, int js, int s, int *I, double *A)
{
  int l, ls;
  ij_to_l (nx, ny, i, j, l);;
  ij_to_l (nx, ny, is, js, ls);
  double S = hx * hy / 2.;

  if (I && l != ls)
    I[s] = ls; //(f_l, f_ls) l - строка, ls - столбец

  if (i > 0 && i < nx && j > 0 && j < ny)
    {
      if (l == ls)
        {
          A[s] = 6.0 * S / 6.0;
          return 0;
        }
      else if (s >= 0 && s <= 5)
        {
          A[s] = 2.0 * S / 12.0;
          return 0;
        }
      else
        return -1;
    }

  if (j == 0 && i > 0 && i < nx)
    {
      if (l == ls)
        {
          A[s] = 3.0 * S / 6.0;
          return 0;
        }
      else if (s == 0 || s == 1)
        {
          A[s] = 1.0 * S / 12.0;
          return 0;
        }
      else if (s == 2 || s == 3)
        {
          A[s] = 2.0 * S / 12.0;
          return 0;
        }
      else
        return -1;
    }

  if (j == ny && i > 0 && i < nx)
    {
      if (l == ls)
        {
          A[s] = 3.0 * S / 6.0;
          return 0;
        }
      else if (s == 0 || s == 3)
        {
          A[s] = 1.0 * S / 12.0;
          return 0;
        }
      else if (s == 1 || s == 2)
        {
          A[s] = 2.0 * S / 12.0;
          return 0;
        }
      else
        return -1;
    }

  if (i == 0 && j > 0 && j < ny)
    {
      if (l == ls)
        {
          A[s] = 3.0 * S / 6.0;
          return 0;
        }
      else if (s == 0 || s == 3)
        {
          A[s] = 2.0 * S / 12.0;
          return 0;
        }
      else if (s == 1 || s == 2)
        {
          A[s] = 1.0 * S / 12.0;
          return 0;
        }
      else
        return -1;
    }

  if (i == nx && j > 0 && j < ny)
    {
      if (l == ls)
        {
          A[s] = 3.0 * S / 6.0;
          return 0;
        }
      else if (s == 0 || s == 3)
        {
          A[s] = 1.0 * S / 12.0;
          return 0;
        }
      else if (s == 1 || s == 2)
        {
          A[s] = 2.0 * S / 12.0;
          return 0;
        }
      else
        return -1;
    }

  if (i == 0 && j == 0)
    {
      if (l == ls)
        {
          A[s] = 2.0 * S / 6.0;
          return 0;
        }
      else if (s == 0 || s == 1)
        {
          A[s] = 1.0 * S / 12.0;
          return 0;
        }
      else if (s == 2)
        {
          A[s] = 2.0 * S / 12.0;
          return 0;
        }
      else
        return -1;
    }

  if (i == nx && j == ny)
    {
      if (l == ls)
        {
          A[s] = 2.0 * S / 6.0;
          return 0;
        }
      else if (s == 0 || s == 2)
        {
          A[s] = 1.0 * S / 12.0;
          return 0;
        }
      else if (s == 1)
        {
          A[s] = 2.0 * S / 12.0;
          return 0;
        }
      else
        return -1;
    }

  if (i == 0 && j == ny)
    {
      if (l == ls)
        {
          A[s] = 1.0 * S / 6.0;
          return 0;
        }
      else if (s == 0 || s == 1)
        {
          A[s] = 1.0 * S / 12.0;
          return 0;
        }
      else
        return -1;
    }

  if (i == nx && j == 0)
    {
      if (l == ls)
        {
          A[s] = 1.0 * S / 6.0;
          return 0;
        }
      else if (s == 0 || s == 1)
        {
          A[s] = 1.0 * S / 12.0;
          return 0;
        }
      else
        return -1;
    }

  return -99999999;
}
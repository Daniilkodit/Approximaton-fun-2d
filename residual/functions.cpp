#include "../header.h"

double
f (int k, double x, double y)
{
  switch (k)
  {
  case 0:
    return 1.0;
  case 1:
    return x;
  case 2:
    return y;
  case 3:
    return x + y;
  case 4:
    return std::sqrt (x * x + y * y);
  case 5:
    return x * x + y * y;
  case 6:
    return std::exp (x * x - y * y);
  case 7:
    return 1.0 / (25.0 * (x * x + y * y) + 1.0);
  default:
    return 0.0;
  }
}

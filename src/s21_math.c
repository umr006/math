#include "s21_math.h"

long double s21_cos(double x) {
  x = s21_fmod(x, 2 * s21_M_PI);
  long double result = 1;
  int sign = -1;
  for (int n = 2; n < 200; n += 2) {
    result += sign * tempow(x, n) / factorial(n);
    sign *= -1;
  }
  return result;
}

long double s21_sin(double x) {
  x = s21_fmod(x, 2 * s21_M_PI);
  long double result = x;
  int sign = -1;
  for (int n = 3; n < 200; n += 2) {
    result += sign * tempow(x, n) / factorial(n);
    sign *= -1;
  }
  return result;
}

long double s21_tan(double x) { return s21_sin(x) / s21_cos(x); }

long double atan_helper(double x) {
  long double result = x;
  long double temp = x;
  double mem_pow = x;
  int i = 1;
  int sign = -1;
  while (s21_fabs(temp) > 1e-10) {
    i += 2;
    mem_pow = mem_pow * x * x;
    temp = sign * mem_pow / i;
    result += temp;
    sign *= -1;
  }
  return result;
}

long double s21_atan(double x) {
  long double result = 0;
  if (x < 1.0 && x > -1.0) {
    result = atan_helper(x);
  } else if (x == 1.0) {
    result = s21_M_PI / 4;
  } else if (x == -1.0) {
    result = -s21_M_PI / 4;
  } else if (x == 0) {
    result = 0;
  } else if (x > 1.0) {
    result = s21_M_PI / 2 - atan_helper(1 / x);
  } else if (x < -1.0) {
    result = -s21_M_PI / 2 - atan_helper(1 / x);
  }
  return result;
}

long double s21_acos(double x) {
  long double result = 0;
  if (x >= 0 && x <= 1)
    result = s21_atan(s21_sqrt(1 - x * x) / x);
  else if (x < 0 && x >= -1)
    result = s21_M_PI + s21_atan(s21_sqrt(1 - x * x) / x);
  else if (x > 1 || x < -1)
    result = s21_NAN;
  return result;
}

long double s21_asin(double x) {
  long double result = 0;
  if (x > 1 || x < -1)
    result = s21_NAN;
  else
    result = s21_atan(x / s21_sqrt(1 - x * x));
  return result;
}

long double s21_fmod(double x, double y) {
  long double result = 0;
  if (y == 0)
    result = s21_NAN;
  else {
    result = x / y;
    result = result - (int)result;
    result *= y;
  }
  return result;
}

long double s21_ceil(double x) {
  long double r = 0;
  long double l = 0;
  x = (long double)x;
  r = x - (long long)x;
  l = (long long)x;

  if (x < 0 && l == 0) {
    printf("%c", '-');
  }

  if (r > 0)
    x = (x + 1) - r;
  else if (r < 0) {
    x = x - r;
  }
  return x;
}

long double s21_floor(double x) {
  long double r = s21_fabs(x - (long long)x);
  if (x < 0 && r > 0) x = x - 1;
  return (long long)x;
}

long double s21_exp(double x) {
  long double result = 1.0;
  long double term = 1.0;
  int i;
  if (x < -14.5)
    result = 0;
  else {
    for (i = 1; term > 1e-17 || term < -1e-17; i++) {
      term *= x / i;
      result += term;
      if (result > DBL_MAX) {
        result = s21_INF;
        break;
      }
    }
  }
  return result;
}

long double s21_sqrt(double x) {
  long double result = 0;
  if (x == s21_NAN) {
    result = s21_NAN;
  } else if (x < 0) {
    result = s21_NAN;
  } else {
    long double l = -1 * x;
    long double r = x + 2;
    long double m = 0;
    long double t = x;
    while (1) {
      m = (r + l) / 2;
      if (m == t || m * m == x) break;
      if (m * m > x)
        r = m;
      else
        l = m;
      t = m;
    }
    result = m;
  }
  return result;
}

long double s21_log(double x) {
  long double res = 0;
  if (x < 0)
    res = s21_NAN;
  else if (x == 0)
    res = -s21_INF;
  else {
    int cnt = 0;
    while ((int)x >= 2) {
      x = x / 2;
      cnt++;
    }
    x = x - 1;
    int sign = 1;
    for (int i = 1; i < 999; i++) {
      res += sign * tempow(x, i) / i;
      sign *= -1;
    }
    res = res + cnt * s21_ln2;
  }
  return res;
}

long double s21_fabs(double x) {
  if (x < 0) x *= -1;
  return x;
}

int s21_abs(int x) {
  if (x < 0) {
    x *= -1;
  }
  return x;
}

long double factorial(int n) {
  if (n == 0) return 1;
  return n * factorial(n - 1);
}

int s21_isinf(double x) {
  int res = 0;
  if ((x == s21_INF) || (x == -s21_INF)) {
    res = 1;
  }
  return res;
}

int isnn(double x) {
  int res = 0;
  if (((x >= 0) && (x < 0))) {
    res = 1;
  }
  return res;
}

long double s21_pow(double base, double exp) {
  long double res = 0;
  long double r_exp = 0.0;
  if (!s21_isinf(exp)) r_exp = exp - (long long)exp;
  if (r_exp <= 0 && !isnn(base) && !isnn(exp) && !s21_isinf(base) &&
      !s21_isinf(exp)) {
    res = tempow(base, exp);
  } else {
    res = fractional_exp(base, exp);
  }
  return res;
}

long double fractional_exp(double base, double exp) {
  long double res;
  if (base == 0 && !isnn(exp))
    res = exp < 0 ? s21_INF : 0;
  else if (base == -1 && s21_isinf(exp))
    res = 1;
  else if (base == 1)
    res = 1;
  else if (s21_isinf(base) && !s21_isinf(exp))
    res = base;
  else if (exp == 0)
    res = 1;
  else if (isnn(base))
    res = s21_NAN;
  else if (exp == s21_INF)
    res = s21_fabs(base) < 1 ? 0 : s21_INF;
  else if (exp == -s21_INF)
    res = s21_fabs(base) < 1 ? s21_INF : 0;
  else if (base == -s21_INF && !isnn(exp))
    res = exp < 0 ? 0 : s21_INF;
  else
    res = s21_exp(exp * s21_log(base));
  return res;
}

long double tempow(double base, double exp) {
  long double result = 1.0;
  if (exp > 0) {
    for (int i = 0; i < exp; i++) {
      result *= base;
    }
  } else {
    exp *= -1;
    for (int i = 0; i < exp; i++) {
      result *= base;
    }
    result = 1 / result;
  }
  return result;
}

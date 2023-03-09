#ifndef SRC_S21_MATH_H_
#define SRC_S21_MATH_H_

#include <float.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>

#define s21_M_PI 3.1415926535897932
#define s21_EPS 1e-17
#define s21_INF 1.0 / 0.0
#define s21_NAN 0.0 / 0.0
#define s21_ln10 2.30258509299404590109
#define s21_MAX_double 1.7976931348623158e308
#define s21_ln2 0.6931471805599453

long double s21_cos(double x);
long double s21_sin(double x);
long double s21_tan(double x);

long double s21_acos(double x);
long double s21_asin(double x);
long double s21_atan(double x);
long double atan_helper(double x);

long double s21_fmod(double x, double y);
long double s21_ceil(double x);
long double s21_floor(double x);

long double s21_pow(double base, double exp);
long double s21_exp(double x);
long double s21_sqrt(double x);
long double s21_log(double x);
long double factorial(int n);
long double s21_fabs(double x);
int s21_abs(int x);
int isnn(double x);
int s21_isinf(double x);
long double fractional_exp(double base, double exp);
long double tempow(double base, double exp);

#endif  // SRC_S21_MATH_H_

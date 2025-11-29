/**
 * @file common_fun.h
 * @brief declaration of general functions
 * @author Liweining (liwn@cau.edu.cn)
 * @version 1.1.1
 * @date 2023-04-03
 *
 * @copyright Copyright (c) {2023}  Liu-Team
 *
 * @par change log:
 * <table>
 * <tr><th>Date       <th>Version <th>Author  <th>Description
 * <tr><td>2022-09-04 <td>1.0.0     <td>liwn     <td>content
 * <tr><td>2023-03-09 <td>1.0.2     <td>liwn     <td>content
 * <tr><td>2023-04-03 <td>1.1.0     <td>liwn     <td>content
 * </table>
 */

#include <errno.h>
#include <fcntl.h>  // count the columns
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/mman.h>  // count the columns
#include <sys/stat.h>  // count the columns
#include <sys/time.h>
#include <unistd.h>

#ifdef linux
#include <sys/types.h>
#include <sys/wait.h>
#endif

#ifndef _COMMMON_FUN_H
#define _COMMMON_FUN_H

#define FILE_BLOCK_SIZE 100 * 1024 * 1024
#define PI 3.14159265358979323846
#define MAX_LINE_LEN 2048 /* Maximum length of rows in all input file */
#define NCHAR_MAX 20
#define FILE_NOT_FOUND 1
#define NUM_DIFF 2
#define MEM_ALLOC 3
#define ERROR_MATH 4
#define UNDER_DEVELOP 5
#define ERROR_PARA 6
#define ERROR_INPUT 7
#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))

    void rand_rank(int *array, int n, int iseed);
int progress(int finish, int all);
void sort(int *array, int n);
int get_column_count(char *filename);
char *read_whole(char *file);
double rnorm01();
void rnorm(int num, double mu, double sigma, double *out);
void rnormMV(int num, int d, double *mu, double *Sigma, double *sample);
void output_ddataframe(double *df, int row, int col, char *file, double scale);
void output_d2i_dataframe(double *df, int row, int col, char *file);
void output_idataframe(int *df, int row, int col, char *file);
void read_dataframe(char *file, double *df, double scale);
long getTimeUsec();
void printf_array(char *head, double *array, int n, int len);
void var(double *mat, int row, int col, double miss, int layout, double *var);
void subset(int row, int col, double* mat, int* row_index, int trans_row, int* col_index, int trans_col, double* out);
int in_char(char *a, char **b, int n);
int in(int *a, int *b, int row, int col);
int read_dataframe_char(char *file, int rows, int cols, char **df, int cols_read, int *cols_sel, int header);
void read_dataframe_char_old(char *file, int rows, int cols, char **df, int cols_read, int *cols_sel, int header);
int read_dataframe_double(char *file, int rows, int cols, double *df, int cols_read, int *cols_sel, int header);
char **unique_char(char **vec, int n, int inc, int *num_uniq, int **num_each, int *level);
char **unique_char_old(char **vec, int n, int inc, int *n_uniq, int *level);
int *string_to_array(char *str, int *num, char *sep);
void string_to_darray(char *str, int n, double *array, int inc, char *sep);
int countLines(char *file);
double sumd(double* array, int n, int inc);
int sumi(int *array, int n, int inc);
void fill_d(double value, double *array, int n, int inc);
void fill_i(int value, int *array, int n, int inc);
void copy_d(double *a, int inca, double *b, int incb, int n, int times);
void copy_i(int *a, int inca, int *b, int incb, int n, int times);
void solve_3d(double *A, double *B);
void solve_2d(double *A, double *B);
void solve_direct(double *mat, int n, double *inv);
void vec_square_sum(double *mat, int row, int col, double *result, int margin);
void sqrt_mat(double *input, int n, double *output);
void shrink(double *vec, int n, int inci, double scale, double *out, int inco);
void add(int n, double *A, double *B, double alpha, double *C);
void add_diag(int order, double* mat, double* out, double add);
void copy_mat(int row, int col, double* a, int lda, double* b, int ldb, double alpha);
void mat_vec_product(int row, int col, double* A, int ldA, double* a, int inc, double* b);
double ddot2(int num, double *x, int incx, double *y, int incy);
void kronecker(int m, int n, double *A, int a, int b, double *B, int incy, double *K);
void cholesky(double* A, int n, double* L);
void mat_to_vec(int row, int col, double* mat);
int countd(int n, double *array, int inc, double value, int type);
int char_match(int nsub, char** subset, int nsuper, char** superset, int* index);
void char_subset(int *index, int nsuper, char **superset);
void output_ddataframe_names(double *df, int row, int col, int ld, char *file, char *header, char **rowname);
int which_int(int value, int *vec, int num);
void self_dgemm(int Layout, int TransA, int TransB, int M, int N, int K, double alpha, double *A, int lda, double *B, int ldb, double beta, double *C, int ldc);
void append_dfile(FILE *fp, int num, double *vec, int inc, char *sep);
char *check_file(char *filename, char *suffix);
int check_header(char *file_chr, int col, char **out_chr);
void error(char *message, int code);
void negate(int *array, int inc, int num);
int min_i(int a, int b);
double min_d(double a, double b);
int max_i(int a, int b);
void strcpy_l(char **to, char *from);

#endif

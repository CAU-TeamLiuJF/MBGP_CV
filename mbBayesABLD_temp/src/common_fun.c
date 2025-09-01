/**
 * @file common_fun.c
 * @brief Definition of Model Dependent Functions
 * @details If not specified, the matrices in the function are all stored in the form of one-dimensional vectors,
 * arow-major order
 * @author Liweining (liwn@cau.edu.cn)
 * @version 1.1.1
 * @date 2023-04-03
 *
 * @copyright Copyright (c) {2023}  Liu-Team
 *
 * @par change log:
 * <table>
 * <tr><th>Date       <th>Version <th>Author  <th>Description
 * <tr><td>2022-09-04 <td>1.0.0     <td>liwn     <td>first version
 * <tr><td>2023-03-09 <td>1.0.2     <td>liwn     <td>
 * <tr><td>2023-04-03 <td>1.1.0     <td>liwn     <td>
 * </table>
 */

#include "common_fun.h"

/**
 * @brief Random ordering of input array
 * @param array        Input array
 * @param n            Number of elements in the array
 * @param iseed        Random number seed
 */
void rand_rank(int* array, int n, int iseed) {
  int temp = 0;
  int j = 0;

  /* Random Number Seed */
  srand(iseed);
  for (int i = 0; i < n; i++) {
    j = (rand() + n) % n;
    temp = array[i];
    array[i] = array[j];
    array[j] = temp;
  }
}

/**
 * @brief Print progress bar
 * @param  finish           Number of times completed
 * @param  all              Total required times
 * @warning                 It will have a great impact on the efficiency of the program
 */
int progress(int finish, int all) {
  char buffer[10] = "";

  /* Report progress */
  float finish_per = (float)finish / all * 100;
  float finish_per_last = (float)(finish - 1) / all * 100;

  sprintf(buffer, "%.2f", finish_per);
  finish_per = atof(buffer);
  sprintf(buffer, "%.2f", finish_per_last);
  finish_per_last = atof(buffer);

  float dif = finish_per - finish_per_last + 0.0001;

  if (dif < 0.01) {
    return 0;
  } else if (finish_per > 100.0) {
    finish_per = 100.0;
  }

  if (isatty(fileno(stdout))) {
    int bar_out = (float)finish / all * 20; /* There are a total of 20 bins */
    if (bar_out > 20) {
      bar_out = 20;
    }
    char bar[21] = {'\0'};                          /* There are a total of 21 digits from 0-20, including '\ 0' */
    char* lable = (char*)"|/-\\";                   /* Only the escape character '\' can represent a '\'*/
    for (int i = 0; i < bar_out; i++) bar[i] = '>'; /* Output one more at a time */
    printf("[%-20s] [%.2f%%][%c]\r", bar, finish_per, lable[finish % 4]);
    fflush(stdout); /* Refresh output buffer */
  } else if ((int)finish_per != (int)finish_per_last) {
    printf("%d%%\n", (int)finish_per);
    fflush(stdout);
  }

  return 0;
}

/**
 * @brief Sort the input integer array
 * @param array       Input array
 * @param n           Number of elements in the array
 */
void sort(int* array, int n) {
  int temp = 0;
  /* Sorting subject */
  for (int i = 0; i < n - 1; i++) {
    for (int j = i + 1; j < n; j++) {
      if (array[i] > array[j]) /* If the front is larger than the back, exchange */
      {
        temp = array[i];
        array[i] = array[j];
        array[j] = temp;
      }
    }
  }
}

/**
 * @brief sums the elements in the array
 * @param  array            Input array
 * @param  n                Number of elements in the array
 * @param  inc              Increment of the array
 * @return Sum of the elements in the array
 */
double sumd(double* array, int n, int inc) {
  {
    double sum = 0;
    for (int i = 0; i < n; i++) {
      sum += array[i * inc];
    }
    return sum;
  }
}

/**
 * @brief sums the elements in the array
 * @param  array            Input array
 * @param  n                Number of elements in the array
 * @param  inc              Increment of the array
 * @return Sum of the elements in the array
 */
int sumi(int* array, int n, int inc) {
  {
    int sum = 0;
    for (int i = 0; i < n; i++) {
      sum += array[i * inc];
    }
    return sum;
  }
}

/**
 * @brief   Get the number of columns in the space separated text file (only judged by the first line, that is, the
 * number of columns in all lines is the same by default)
 * @param filename  File name
 * @return  Number of columns
 */
int get_column_count(char* filename)
{
  int fd = open(filename, O_RDONLY);
  if (fd == -1) {
    printf("Failed to open file: %s\n", filename);
    exit(FILE_NOT_FOUND);
  }

  struct stat sb;
  if (fstat(fd, &sb) == -1) {
    printf("Failed to get file size: %s\n", filename);
    exit(MEM_ALLOC);
  }

  char* file_data = mmap(NULL, sb.st_size, PROT_READ, MAP_PRIVATE, fd, 0);
  if (file_data == MAP_FAILED) {
    printf("Failed to map file into memory: %s\n", filename);
    exit(MEM_ALLOC);
  }

  int column_count = 0;
  int is_space = 1;
  int i = 0;
  while (i < sb.st_size && file_data[i] != '\n') {
    if (file_data[i] == ' ' || file_data[i] == '\t' || file_data[i] == '\r') {
      is_space = 1;
    }
    else if (is_space) {
      is_space = 0;
      column_count++;
    }
    i++;
  }

  munmap(file_data, sb.st_size);
  close(fd);

  return column_count;
}

/**
 * @brief Read and store the entire file in a string vector
 * @param file  File name
 * @return  String vector
 */
char* read_whole(char* file) {
  /* Read in binary */
  FILE* pFile = fopen(file, "rb");
  if (pFile == NULL) {
    printf("Can not load file: %s\n", file);
    exit(FILE_NOT_FOUND);
  }

  /* get file size */
  fseek(pFile, 0, SEEK_END);
  long lSize = ftell(pFile);
  rewind(pFile);

  /* allocate memory */
  char* buffer = (char*)calloc(lSize, sizeof(char));
  if (buffer == NULL) {
    fputs("Memory allocate error.", stderr);
    exit(MEM_ALLOC);
  }

  /* read whole file to buffer */
  long result = fread(buffer, 1, lSize, pFile);
  if (result != lSize) {
    fputs("fread error", stderr);
    exit(FILE_NOT_FOUND);
  }

  fclose(pFile);
  return buffer;
}

/**
 * @brief Random sampling from standard normal distribution
 * @return a double value
 */
double rnorm01() {
  static double U = 0.0, V;
  int runif = 0;
  double Z1;
  // double Z2;

  runif = rand();
  while (runif == 0) {
    runif = rand();
  }

  U = (double) runif / (RAND_MAX + 1.0);
  V = (double) rand() / (RAND_MAX + 1.0);
  Z1 = sqrt(-2.0 * log(U)) * sin(2.0 * PI * V);
  // Z2 = sqrt(-2.0 * log(U)) * cos(2.0 * PI * V);

  return Z1;
}

/**
 * @brief Random sampling from one-dimensional normal distribution
 * @param num     Number of samples
 * @param mu      Mean of the normal distribution
 * @param sigma   Standard deviation of the normal distribution
 * @param out     Output array
 */
void rnorm(int num, double mu, double sigma, double* out) {
  for (int i = 0; i < num; i++) {
    out[i] = rnorm01() * sqrt(sigma) + mu;
  }
}

/**
 * @brief Multivariate normal distribution sampling
 * @author  Jicai Jiang
 * @param  n                The order of the matrix Sigma; n ≥ 0.
 * @param  mu               This is mean vector mu with length k.
 * @param  Sigma            This is the k x k covariance matrix Sigma.
 * @param  sample           random deviates
 * @param  stream           MKL random number generation
 * @include                 mkl.h
 */
void rnormMV(int num, int d, double* mu, double* Sigma, double* sample) {
  double snRand[d];
  double VchoL[d * d];
  cholesky(Sigma, d, VchoL);

  for (size_t i = 0; i < num; i++) {
    rnorm(d, 0, 1, snRand);

    /* Convert to a normal distribution with specified mean and variance (Y = (L')^{-1}X + μ ) */
    for (size_t j = 0; j < d; j++) {
      sample[i * d + j] = mu[j];
      for (size_t k = 0; k < (j + 1); k++) {
        sample[i * d + j] += VchoL[j * d + k] * snRand[k];
      }
    }
  }
}

/**
 * @brief print double matrix to target file
 * @param  df               matrix (row order)
 * @param  row              number of lines
 * @param  col              number of columns
 * @param  file             output file name
 * @param  scale            shrink factor for elements in matrix, e.g. df[i]/scale
 */
void output_ddataframe(double* df, int row, int col, char* file, double scale) {
  FILE* fp = fopen(file, "wt");
  for (size_t i = 0; i < row; i++) {
    for (size_t j = 0; j < col; j++) {
      if (j == col - 1) {
        fprintf(fp, "%E", df[col * i + j] / scale);
      } else {
        fprintf(fp, "%E\t", df[col * i + j] / scale);
      }
    }
    putc('\n', fp);
  }
  fclose(fp);
}

/**
 * @brief print int matrix to target file
 * @param  df               matrix (row order)
 * @param  row              number of lines
 * @param  col              number of columns
 * @param  file             output file name
 * @param  scale            shrink factor for elements in matrix, e.g. df[i]/scale
 */
void output_idataframe(int* df, int row, int col, char* file)
{
  FILE* fp = fopen(file, "wt");
  for (size_t i = 0; i < row; i++) {
    for (size_t j = 0; j < col; j++) {
      if (j == col - 1) {
        fprintf(fp, "%d", df[col * i + j]);
      }
      else {
        fprintf(fp, "%d\t", df[col * i + j]);
      }
    }
    putc('\n', fp);
  }
  fclose(fp);
}

/**
 * @brief print double matrix (Keep integer part) to target file
 * @param  df               matrix (row order)
 * @param  row              number of lines
 * @param  col              number of columns
 * @param  file             output file name
 */
void output_d2i_dataframe(double* df, int row, int col, char* file) {
  FILE* fp = fopen(file, "wt");
  for (size_t i = 0; i < row; i++) {
    for (size_t j = 0; j < col; j++) {
      if (j == col - 1) {
        fprintf(fp, "%0.0f", df[col * i + j]);
      } else {
        fprintf(fp, "%0.0f\t", df[col * i + j]);
      }
    }
    putc('\n', fp);
  }
  fclose(fp);
}

/**
 * @brief  Negate a vector containing logical variables (0/1)
 * @param array Vector to be negated
 * @param inc   Increment of values to be subsetted
 * @param num   Number of values to be subsetted
 */
void negate(int* array, int inc, int num) {
  for (size_t i = 0; i < num; i++) {
    array[i] = 1 - array[i];
  }
}

/**
 * @brief Read the file and store it in a float array
 * @param file  File name
 * @param df    Float array to store the data
 * @param scale Scale factor
 */
void read_dataframe(char* file, double* df, double scale) {
  /* Check if the file is readable */
  FILE* fp;
  if ((fp = fopen(file, "rt")) == NULL) {
    printf("cannot read the file: %s\n", file);
    exit(FILE_NOT_FOUND);
  }

  /* read file */
  int index = 0;
  double tmp;
  while (fscanf(fp, "%lf", &tmp) == 1) {
    df[index] = tmp / scale;
    index++;
  }
  fclose(fp);
}

/**
 * @brief timer
 * @return
 */
long getTimeUsec() {
  struct timeval t;
  gettimeofday(&t, 0);
  return (long)((long)t.tv_sec * 1000 * 1000 + t.tv_usec);
}

/**
 * @brief Print an float array
 * @param head    What appears at the beginning of the line
 * @param array   The array to be printed
 * @param n       The length of the array
 * @param len     Print the significant digits of a floating point number
 */
void printf_array(char* head, double* array, int n, int len) {
  printf("%s", head);
  for (size_t i = 0; i < n; i++) {
    printf("%.*lf\t", len, array[i]);
  }
  printf("\n");
}

/**
 * @brief Calculate the column vector variance of the matrix
 * @param  mat              input matrix
 * @param  row              number of rows in input matrix
 * @param  col              number of columns in input matrix
 * @param  miss             missing value
 * @param  layout           1 for row-major and 0 for column-major
 * @param  var              output variance
 */
void var(double* mat, int row, int col, double miss, int layout, double* var)
{
  // two-dimensional array storage format
  int index;

  // mean
  int i, j, count;
  double sum;
  double* mean = (double*)calloc(col, sizeof(double));
  for (i = 0; i < col; i++) {
    sum = 0.0;
    count = 0;
    for (j = 0; j < row; j++) {
      if (layout == 1) {
        index = j * col + i;
      } else {
        index = i * row + j;
      }

      // skip missing value
      if (mat[index] == miss) {
        continue;
      }

      sum += mat[index];
      count++;
    }
    mean[i] = sum / count;
  }

  memset(var, 0, col * sizeof(double));
  for (i = 0; i < col; i++) {
    count = 0;
    for (j = 0; j < row; j++) {
      if (layout == 1) {
        index = j * col + i;
      } else {
        index = i * row + j;
      }
      // skip missing value
      if (mat[index] == miss) {
        continue;
      }
      var[i] += pow(mat[index] - mean[i], 2);
      count++;
    }
    var[i] /= (count - 1);
  }

  free(mean);
}

/**
 * @brief Get a subset of a matrix
 * @param row            number of rows in input matrix
 * @param col            number of columns in input matrix
 * @param mat            input matrix
 * @param row_index      row index of subset
 * @param trans_row      whether to transpose the row index, 1 for yes and 0 for no
 * @param col_index      column index of subset
 * @param trans_col      whether to transpose the column index, 1 for yes and 0 for no
 * @param out            output matrix
 */
void subset(int row, int col, double* mat, int* row_index, int trans_row, int* col_index, int trans_col, double* out) {
  int index = 0;
  int rowi, colj;

  for (size_t i = 0; i < row; i++) {
    rowi = trans_row ? 1 - row_index[i] : row_index[i];
    if (!rowi)
      continue;
    for (size_t j = 0; j < col; j++) {
      colj = trans_col ? 1 - col_index[j] : col_index[j];
      if (!colj)
        continue;
      out[index++] = mat[i * col + j];
    }
  }
}

/**
 * @brief check whether the string appears in a string array
 * @param  a                input string
 * @param  b                string array
 * @param  n                Number of string elements in b
 * @return position of string a in b, -1 if not found
 */
int in_char(char* a, char** b, int n) {
  if (n < 1) {
    return -1;
  }
  for (size_t i = 0; i < n; i++) {
    if (strcmp(a, b[i]) == 0) {
      return i;
    }
  }
  return -1;
}

/**
 * @brief check whether vector a is in b matrix
 * @param  a              input vector
 * @param  b              matrix
 * @param  row            number of rows in b
 * @param  col            number of columns in b
 * @return  1 if a is in b, 0 if not
 */
int in(int* a, int* b, int row, int col) {
  int eqnum, logi;
  for (size_t i = 0; i < row; i++) {
    eqnum = 0;
    for (size_t j = 0; j < col; j++) {
      logi = a[j] == b[i * col + j];
      eqnum += logi;
    }

    if (eqnum == col) {
      return 1;
    }
  }
  return 0;
}

/**
 * @brief Determine whether the file contains a header line based on whether the specified position in the first line
 * can be resolved to a float value
 * @param  file_chr   A string that contains at least the first line of the file, ending with a line break
 * @param  col        The column number of the file to be checked, 0 for checking all columns
 * @param  out_chr    output string
 * @return            1 for that file has a header line and 0 for no header line
 */
int check_header(char *file_chr, int col, char **out_chr) {
  int i = 0;
  float r1 = 0;
  int index = 0;
  int to_float = 0;
  char* value = "";
  char* first_rest;
  char* first_line = strtok_r(file_chr, "\n", out_chr);

  /* a copy of the first row */
  int first_line_len = strlen(first_line);
  char* first_line_copy = (char*)malloc((first_line_len + 1) * sizeof(char));
  strcpy(first_line_copy, first_line);

  while (1) {
    value = strtok_r(first_line_copy, " \t", &first_rest);
    if (value == NULL) break;
    to_float = sscanf(value, "%f", &r1);
    if (col != 0 && i == col) {
      index = to_float;
    } else {
      if (to_float > 0) {
        index++;
      };
    }
    i++;
    first_line_copy = NULL;
  }
  if (index > 0) {
    /* no header line */
    value = strchr(file_chr, '\0');
    *value = '\n';
    i = 0;
    *out_chr = file_chr;
  }
  else {
    /* with header line */
    i = 1;
  }

  free(first_line_copy);
  return i;
}

/**
 * @brief Read the specified columns in the text file and store it as a strings array
 * @param  file             file name
 * @param  rows             number of rows
 * @param  cols             number of columns
 * @param  df               strings array to store the data
 * @param  cols_read        number of columns to read
 * @param  cols_sel         position of columns to read
 * @param  header           whether the file has a header, 1 for yes and 0 for no
 */
void read_dataframe_char_old(char* file, int rows, int cols, char** df, int cols_read, int* cols_sel, int header)
{
  /* Read the entire file (note the file size) */
  char* file_char = read_whole(file);

  int stored_coli = 0;
  char* save_ptr = NULL;
  char* buff = NULL;

  /* If the number of rows and columns is not provided, get them from the file */
  if (rows == 0) {
    rows = countLines(file);
    /* alloc memory */
    df = (char**)calloc(rows * cols_read, sizeof(char*));
    for (size_t i = 0; i < rows * cols_read; i++) {
      df[i] = (char*)calloc(NCHAR_MAX, sizeof(char));
    }
  }
  if (cols == 0) {
    cols = get_column_count(file);
  }

  if (header) {
    buff = strtok_r(file_char, "\n", &save_ptr);
  }
  else {
    save_ptr = file_char;
  }

  /* Read by row */
  for (size_t i = 0; i < rows; i++) {
    stored_coli = 0;
    for (int j = 1; j <= cols; j++) {
      if (j < cols) {
        buff = strtok_r(NULL, " ", &save_ptr);
      }
      else {
        buff = strtok_r(NULL, "\n", &save_ptr);
      }

      if (in(&j, cols_sel, cols_read, 1)) {
        if (buff == NULL) {
          printf("Error: error in reading file: %s.\n", file);
          exit(2);
        }

        strcpy(df[i * cols_read + stored_coli], buff);
        stored_coli++;
      }
    }
  }
}

/**
 * @brief Read the specified columns in the text file and store it as a strings array
 * @param  file             file name
 * @param  rows             number of rows
 * @param  cols             number of columns
 * @param  df               strings array to store the data
 * @param  cols_read        number of columns to read
 * @param  cols_sel         position of columns to read
 * @param  header           whether the file has a header, 1 for yes and 0 for no
 * @return                  Number of file lines except header line
 */
int read_dataframe_char(char* file, int rows, int cols, char** df, int cols_read, int* cols_sel, int header) {
  /* If the number of rows and columns is not provided, get them from the file */
  if (rows == 0) {
    rows = countLines(file);
    
    /* alloc memory */
    df = (char**)calloc(rows * cols_read, sizeof(char*));
    if (df == NULL) {
      printf("Failed to allocate memory when reading character type columns\n");
      exit(MEM_ALLOC);
    }
    for (size_t i = 0; i < rows * cols_read; i++) {
      df[i] = (char*)calloc(NCHAR_MAX, sizeof(char));
    }
  }
  if (cols == 0) {
    cols = get_column_count(file);
  }

  /* Read the entire file (note the file size) */
  char* file_char = read_whole(file);

  /* Header */
  char *buff = NULL, *save_ptr = NULL;
  if (header == 1) {
    /* with header line */
    buff = strtok_r(file_char, "\n", &save_ptr);
  }
  else if (header < 0) {
    /* not sure */
    header = check_header(file_char, cols_sel[0], &save_ptr);
  }
  else {
    /* no header line */
    save_ptr = file_char;
  }
  rows -= header;

  /* Read by row */
  int stored_coli = 0;
  for (size_t i = 0; i < rows; i++) {
    stored_coli = 0;
    for (int j = 1; j <= cols; j++) {
      buff = strtok_r(NULL, " \t", &save_ptr);
      if (in(&j, cols_sel, cols_read, 1)) {
        if (buff == NULL) {
          printf("Error: error in reading file: %s.\n", file);
          exit(FILE_NOT_FOUND);
        }

        strcpy(df[i * cols_read + stored_coli], buff);
        stored_coli++;
        if (stored_coli == cols_read) {
          buff = strtok_r(NULL, "\n", &save_ptr);
          break;
        }
      }
    }
  }

  return rows;
}

/**
 * @brief Read the specified columns in the text file and store it as a double array
 * @param  file             file name
 * @param  rows             number of rows
 * @param  cols             number of columns
 * @param  df               double array to store the data
 * @param  cols_read        number of columns to read
 * @param  cols_sel         position of columns to read
 * @param  header           whether the file has a header line, 1 for yes, 0 for no and -1 for not sure
 * @return                  Number of file lines except header line
 */
int read_dataframe_double(char* file, int rows, int cols, double* df, int cols_read, int* cols_sel, int header) {
  /* Load the entire file (pay attention to the file size) */
  char* file_char = read_whole(file);

  /* If the number of lines and columns of the file is not provided, it is obtained from the file */
  if (rows == 0) {
    rows = countLines(file);
  }
  if (cols == 0) {
    cols = get_column_count(file);
  }

  /* If no columns to read is specified, the entire file is read */
  if (cols_read <= 0 || cols_read > cols) {
    cols_read = cols;
    cols_sel = (int*)calloc(cols_read, sizeof(int));
    for (size_t i = 0; i < cols_read; i++) {
      cols_sel[i] = i + 1;
    }
  }

  int stored_coli = 0;
  char* save_ptr = NULL;
  char* buff = NULL;

  /* Header */
  if (header == 1) {
    /* with header line */
    buff = strtok_r(file_char, "\n", &save_ptr);
  }
  else if (header < 0) {
    /* not sure */
    header = check_header(file_char, cols_sel[0], &save_ptr);
  } else {
    save_ptr = file_char;
  }
  rows-=header;

  /* Read line by line */
  for (size_t i = 0; i < rows; i++) {
    stored_coli = 0;
    for (int j = 1; j <= cols; j++) {
      if (j < cols) {
        buff = strtok_r(NULL, " \t", &save_ptr);
      } else {
        buff = strtok_r(NULL, "\n", &save_ptr);
      }

      if (in(&j, cols_sel, cols_read, 1)) {
        if (buff == NULL) {
          printf("Error: error in reading file: %s\n", file);
          exit(FILE_NOT_FOUND);
        }

        df[i * cols_read + stored_coli] = atof(buff);
        stored_coli++;
      }
    }
  }

  return rows;
}

/**
 * @brief Find the unique string and the number of unique strings in the string vector
 * @param vec       string vector
 * @param n         number of strings
 * @param inc       Increment when reading strings in vec
 * @param num_uniq  number of unique strings
 * @param num_each  Number of occurrences of each unique string
 * @param level     Index of elements in `vec` in return unique strings
 * @return          unique strings
 */
char** unique_char(char** vec, int n, int inc, int *num_uniq, int **num_each, int* level)
{
  /* Temporary string container */
  char** tmp_uniq_chr = (char**)calloc(n, sizeof(char*));
  for (size_t i = 0; i < n; i++) {
    tmp_uniq_chr[i] = (char*)calloc(NCHAR_MAX, sizeof(char));
  }
  int* tmp_each_num = (int*)calloc(n, sizeof(int));

  /* Retrieve new characters */
  int index;
  int uniq_num = 0;
  for (size_t i = 0; i < n; i++) {
    /* TODO:Starting from i may cause bugs, please pay attention to comparing the results */
    index = in_char(vec[i * inc], tmp_uniq_chr, uniq_num);
    if (index < 0) {
      /* new character */
      strcpy(tmp_uniq_chr[uniq_num], vec[i * inc]);
      index = uniq_num;
      uniq_num++;
    }
    level[i] = index;
    tmp_each_num[index]++;
  }

  /* Create a new small container to put unique strings */
  char** uniq = (char**)calloc(uniq_num, sizeof(char*));
  int *num_each_out = (int*)calloc(uniq_num, sizeof(int));
  for (size_t i = 0; i < uniq_num; i++) {
    num_each_out[i] = tmp_each_num[i];
    uniq[i] = (char*)calloc(NCHAR_MAX, sizeof(char));
    strcpy(uniq[i], tmp_uniq_chr[i]);
  }
  *num_each = num_each_out;

  /* Free the requested memory */
  for (size_t i = 0; i < n; i++) {
    free(tmp_uniq_chr[i]);
  }
  free(tmp_uniq_chr);
  free(tmp_each_num);

  *num_uniq = uniq_num;
  return uniq;
}

/**
 * @brief Find the unique string and the number of unique strings in the string vector
 * @param vec      string vector
 * @param n        number of strings
 * @param inc      Increment when reading strings in vec
 * @param n_uniq   Number of occurrences of each unique string
 * @param level    unique strings
 * @return         number of unique strings
 */
char** unique_char_old(char** vec, int n, int inc, int* n_uniq, int* level)
{
  /* Temporary string container */
  char** tmp = (char**)calloc(n, sizeof(char*));
  for (size_t i = 0; i < n; i++) {
    tmp[i] = (char*)calloc(NCHAR_MAX, sizeof(char));
  }

  int index;
  int num = 0;
  /* Add the first character */
  strcpy(tmp[0], vec[0]);

  /* Retrieve new characters */
  for (size_t i = 1; i < n; i++) {
    index = in_char(vec[i * inc], tmp, num);
    if (index < 0) {
      /* For new characters */
      strcpy(tmp[num], vec[i * inc]);
      index = num;
      num++;
    }
    level[i] = index;
  }

  /* Create a new small container to put unique strings */
  char** uniq = (char**)calloc(num, sizeof(char*));
  for (size_t i = 0; i < num; i++) {
    uniq[i] = (char*)calloc(NCHAR_MAX, sizeof(char));
    strcpy(uniq[i], tmp[i]);
  }

  /* Free the requested memory */
  for (size_t i = 0; i < n; i++) {
    free(tmp[i]);
  }
  free(tmp);

  *n_uniq = num;
  return uniq;
}

/**
 * @brief Parses the string separated by the specified delimiter into an integer array
 * @param  str              An array of strings
 * @param  num              Number of elements in the result array
 * @param  sep              delimiter of 'str'
 * @return Parsed integer array
 */
int* string_to_array(char* str, int* num, char* sep)
{
  int num_tmp = 0;
  char* tmp_str;
  int* tmp_values = (int*)calloc(MAX_LINE_LEN, sizeof(int));

  // first value
  tmp_str = strtok(str, sep);
  tmp_values[0] = atoi(tmp_str);
  (num_tmp)++;

  // next values
  tmp_str = strtok(NULL, sep);
  while (tmp_str) {
    tmp_values[num_tmp] = atoi(tmp_str);
    (num_tmp)++;
    tmp_str = strtok(NULL, sep);
  }

  // storage
  int num_cols = num_tmp;
  int* values = (int*)calloc(num_cols, sizeof(int));
  for (size_t i = 0; i < num_cols; i++) {
    values[i] = tmp_values[i];
  }

  free(tmp_values);

  *num = num_tmp;
  return values;
}

/**
 * @brief Parses the string separated by the specified delimiter into an doule array
 * @param  str              An array of strings
 * @param  n                Number of elements in the input array
 * @param  array            Parsed output array
 * @param  inc              Interval of elements when saving results to array
 * @param  sep              delimiter of 'str'
 */
void string_to_darray(char* str, int n, double* array, int inc, char* sep) {
  char* buf;
  int j = 0;

  for (size_t i = 0; i < n; i++) {
    buf = strtok(str, sep);
    array[j] = atof(buf);
    str = NULL;
    j += inc;
  }
}

/**
 * @brief Fill a setting value in the specified position of the double array
 * @param value   Value to be filled
 * @param array   input array
 * @param n       Number of fills
 * @param inc     Interval of elements when saving results to array
 */
void fill_d(double value, double* array, int n, int inc) {
  for (size_t i = 0; i < n; i++) {
    array[i * inc] = value;
  }
}

/**
 * @brief Fill a setting value in the specified position of the int array
 * @param value   Value to be filled
 * @param array   input array
 * @param n       Number of fills
 * @param inc     Interval of elements when saving results to array
 */
void fill_i(int value, int* array, int n, int inc) {
  for (size_t i = 0; i < n; i++) {
    array[i * inc] = value;
  }
}

/**
 * @brief copy values, e.g. sub(b) = sub(a)
 * @param  a            Input array
 * @param  inca         Interval of elements when reading a
 * @param  b            target array
 * @param  incb         Interval of elements when saving results to b
 * @param  n            Number of elements in the input array
 * @param  times        Number of times to copy
 */
void copy_d(double* a, int inca, double* b, int incb, int n, int times) {
  for (size_t i = 0; i < times; i++) {
    for (size_t j = 0; j < n; j++) {
      b[i * n + j * incb] = a[j * inca];
    }
  }
}

/**
 * @brief copy values, e.g. sub(b) = sub(a)
 * @param  a            Input array
 * @param  inca         Interval of elements when reading a
 * @param  b            target array
 * @param  incb         Interval of elements when saving results to b
 * @param  n            Number of elements in the input array
 * @param  times        Number of times to copy
 */
void copy_i(int* a, int inca, int* b, int incb, int n, int times) {
  for (size_t i = 0; i < times; i++) {
    for (size_t j = 0; j < n; j++) {
      b[i * n + j * incb] = a[j * inca];
    }
  }
}

/**
 * @brief Direct inversion of a 3x3 matrix
 * @param mat          Input matrix (3x3)
 * @param inv          Inverted matrix (3x3)
 */
void solve_3d(double* mat, double* inv) {
  inv[0] = mat[4] * mat[8] - mat[5] * mat[7];
  inv[1] = mat[2] * mat[7] - mat[1] * mat[8];
  inv[2] = mat[1] * mat[5] - mat[2] * mat[4];
  inv[3] = mat[5] * mat[6] - mat[3] * mat[8];
  inv[4] = mat[0] * mat[8] - mat[2] * mat[6];
  inv[5] = mat[2] * mat[3] - mat[0] * mat[5];
  inv[6] = mat[3] * mat[7] - mat[4] * mat[6];
  inv[7] = mat[6] * mat[1] - mat[0] * mat[7];
  inv[8] = mat[0] * mat[4] - mat[1] * mat[3];
  double det = mat[0] * inv[0] + mat[1] * inv[3] + mat[2] * inv[6];
  if (det == 0) {
    printf("Noninvertible Matrix.\n");
    exit(EXIT_FAILURE);
  }
  for (int i = 0; i < 9; i++) {
    inv[i] /= det;
  }
}

/**
 * @brief Direct inversion of a 2x2 matrix
 * @param mat          Input matrix (2x2)
 * @param inv          Inverted matrix (2x2)
 */
void solve_2d(double* mat, double* inv) {
  inv[0] = mat[3];
  inv[1] = -1.0 * mat[1];
  inv[2] = -1.0 * mat[2];
  inv[3] = mat[0];
  double det = mat[0] * inv[3] - mat[1] * mat[2];
  if (det == 0) {
    printf("Noninvertible Matrix.\n");
    exit(EXIT_FAILURE);
  }
  for (int i = 0; i < 4; i++) {
    inv[i] /= det;
  }
}

/**
 * @brief Applying the definition of inverse matrix to solve the inverse directly
 * @param  mat               input matrix
 * @param  n                 dimensions of the matrix
 * @param  inv               inverse of the matrix
 */
void solve_direct(double* mat, int n, double* inv) {
  int is[n], js[n], i, j, k, l, u, v;
  double d, p;
  int nn = n * n;
  for (k = 0; k < nn; k++) inv[k] = mat[k];
  for (k = 0; k <= n - 1; k++) {
    d = 0.0;
    for (i = k; i <= n - 1; i++) {
      for (j = k; j <= n - 1; j++) {
        l = i * n + j;
        p = fabs(inv[l]);
        if (p > d) {
          d = p;
          is[k] = i;
          js[k] = j;
        }
      }
    }
    if (is[k] != k) {
      for (j = 0; j <= n - 1; j++) {
        u = k * n + j;
        v = is[k] * n + j;
        p = inv[u];
        inv[u] = inv[v];
        inv[v] = p;
      }
    }
    if (js[k] != k) {
      for (i = 0; i <= n - 1; i++) {
        u = i * n + k;
        v = i * n + js[k];
        p = inv[u];
        inv[u] = inv[v];
        inv[v] = p;
      }
    }
    l = k * n + k;
    inv[l] = 1.0 / inv[l];
    for (j = 0; j <= n - 1; j++) {
      if (j != k) {
        u = k * n + j;
        inv[u] = inv[u] * inv[l];
      }
    }
    for (i = 0; i <= n - 1; i++) {
      if (i != k) {
        for (j = 0; j <= n - 1; j++) {
          if (j != k) {
            u = i * n + j;
            inv[u] = inv[u] - inv[i * n + k] * inv[k * n + j];
          }
        }
      }
    }
    for (i = 0; i <= n - 1; i++) {
      if (i != k) {
        u = i * n + k;
        inv[u] = -inv[u] * inv[l];
      }
    }
  }
  for (k = n - 1; k >= 0; k--) {
    if (js[k] != k) {
      for (j = 0; j <= n - 1; j++) {
        u = k * n + j;
        v = js[k] * n + j;
        p = inv[u];
        inv[u] = inv[v];
        inv[v] = p;
      }
    }
    if (is[k] != k) {
      for (i = 0; i <= n - 1; i++) {
        u = i * n + k;
        v = i * n + is[k];
        p = inv[u];
        inv[u] = inv[v];
        inv[v] = p;
      }
    }
  }
}

/**
 * @brief The horizontal/vertical (0/1) inner product of the matrix, that is, the sum of the squares of the elements
 * of each column/row
 * @param  mat              input matrix
 * @param  row              number of rows of the matrix
 * @param  col              number of columns of the matrix
 * @param  result           output array
 * @param  margin           0 for horizontal and 1 for vertical
 */
void vec_square_sum(double* mat, int row, int col, double* result, int margin) {
  int vector, value;
  if (margin) {
    vector = col;
    value = row;
  } else {
    vector = row;
    value = col;
  }

  for (size_t i = 0; i < vector; i++) {
    result[i] = 0;  // 20221025
    for (size_t j = 0; j < value; j++) {
      if (margin) {
        result[i] += mat[j * col + i] * mat[j * col + i];
      } else {
        result[i] += mat[i * col + j] * mat[i * col + j];
      }
    }
  }
}

/**
 * @brief Get the square of the vector
 * @param input input vector
 * @param n length of the vector
 * @param output output square of the vector
 */
void sqrt_mat(double* input, int n, double* output) {
  for (int i = 0; i < n; i++) {
    if (input[i] < 0) {
      printf("Error: Cannot get square root of negative number!\n");
      exit(ERROR_MATH);
    }

    output[i] = input[i];
  }
}

/**
 * @brief scales each element in the vector
 * @param  vec              input vector
 * @param  n                number of elements to scale
 * @param  inci             increment between elements of the input vector
 * @param  scale            scale factor, e.g. vec[i * inc] * scale
 * @param  out              output vector
 * @param  inco             increment between elements of the output vector
 */
void shrink(double* vec, int n, int inci, double scale, double* out, int inco)
{
  for (int i = 0; i < n; i++) {
    out[i * inco] = vec[i * inci] * scale;
  }
}

/**
 * @brief C = A + B * alpha
 * @param  n                number of elements in arrays
 * @param  A                array A
 * @param  B                array B
 * @param  alpha            scale factor
 * @param  C                array C
 */
void add(int n, double* A, double* B, double alpha, double* C) {
  for (int i = 0; i < n; i++) {
    C[i] = A[i] + B[i] * alpha;
  }
}

/**
 * @brief out = mat; diag(out) = diag(mat) + add
 * @param  order            order of the matrix
 * @param  mat              an square matrix
 * @param  out              output matrix
 * @param  add              value to add to the diagonal
 */
void add_diag(int order, double* mat, double* out, double add) {
  for (size_t i = 0; i < order; i++) {
    for (size_t j = 0; j < order; j++) {
      if (i == j) {
        out[i * order + j] = mat[i * order + j] + add;
      } else {
        out[i * order + j] = mat[i * order + j];
      }
    }
  }
}

/**
 * @brief sub(b) = sub(a) * alpha
 * @param  row             number of rows of matrix a
 * @param  col             number of columns of matrix a
 * @param  a                matrix a
 * @param  lda              leading dimension of matrix a
 * @param  b                matrix b
 * @param  ldb              leading dimension of matrix b
 * @param  alpha            scalar alpha
 */
void copy_mat(int row, int col, double* a, int lda, double* b, int ldb, double alpha) {
  for (int i = 0; i < row; i++) {
    for (int j = 0; j < col; j++) {
      b[i * ldb + j] = a[i * lda + j] * alpha;
    }
  }
}

/**
 * @brief b = A %*% a
 * @param  row              number of rows of matrix A
 * @param  col              number of columns of matrix A
 * @param  A                A row x col matrix
 * @param  ldA              leading dimension of matrix A
 * @param  a                a vector of length col
 * @param  inc              increment between elements of a
 * @param  b                a vector of length row
 */
void mat_vec_product(int row, int col, double* A, int ldA, double* a, int inc, double* b) {
  memset(b, 0, row * sizeof(double));
  for (int i = 0; i < row; i++) {
    for (int j = 0; j < col; j++) {
      b[i] += A[i * ldA + j] * a[j * inc];
    }
  }
}

/**
 * @brief x*y
 * @param num
 * @param x
 * @param incx
 * @param y
 * @param incy
 */
double ddot2(int num, double* x, int incx, double* y, int incy) {
  double tmp = 0;
  for (size_t i = 0; i < num; i++) {
    tmp += x[i * incx] * y[i * incy];
  }
  return tmp;
}

/**
 * @brief M = kronecker(A, B)
 * @param rowA          row of A
 * @param colA          col of A
 * @param A             matrix A
 * @param rowB          row of B
 * @param colB          col of B
 * @param B             matrix B
 * @param incB          leading dimension of matrix B
 * @param out           output matrix
 */
void kronecker(int rowA, int colA, double* A, int rowB, int colB, double* B, int incB, double* out) {
  for (size_t i = 0; i < rowA; i++) {
    for (size_t j = 0; j < colA; j++) {
      copy_mat(rowB, colB, B, incB, out + i * rowB * colA * colB + j * colB, colA * colB, A[i * colA + j]);
    }
  }
}

/**
 * @brief Cholesky decomposition is performed on the matrix A of n x n to obtain the lower triangular matrix
 * @param  A                Matrix stored in vector form
 * @param  n                The order of the matrix A; n ≥ 0.
 * @param  L                lower triangular matrix
 */
void cholesky(double* A, int n, double* L) {
  if (L == NULL) {
    exit(EXIT_FAILURE);
  }

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      L[i * n + j] = 0;
    }
  }

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < (i + 1); j++) {
      double s = 0;
      for (int k = 0; k < j; k++) s += L[i * n + k] * L[j * n + k];
      L[i * n + j] = (i == j) ? sqrt(A[i * n + i] - s) : (1.0 / L[j * n + j] * (A[i * n + j] - s));
    }
  }
}

/**
 * @brief
 * @param row          number of rows of matrix mat
 * @param col          number of columns of matrix mat
 * @param mat          input matrix mat
 */
void mat_to_vec(int row, int col, double* mat) {
  double* vec = (double*)malloc(row * col * sizeof(double));
  for (int i = 0; i < col; i++) {
    for (int j = 0; j < row; j++) {
      vec[i * row + j] = mat[j * col + i];
    }
  }
  copy_d(vec, 1, mat, 1, row * col, 1);
  free(vec);
}

/**
 * @brief count the number of elements equal to specific value in vector
 * @param n         length of vector
 * @param array     input vector
 * @param inc       increment between elements of array
 * @param value     specific value
 * @param type      1 for count the number of elements equal to value and 0 for doing the opposite
 * @return  the number of elements (not) equal to value
 */
int countd(int n, double* array, int inc, double value, int type) {
  int count = 0;
  for (int i = 0; i < n; i++) {
    if (array[i * inc] == value) {
      count++;
    }
  }

  if (type == 1) {
    return count;
  } else {
    return n - count;
  }
}

/**
 * @brief get the index of
 * @param nsub      length of subset
 * @param subset    subset
 * @param nsuper    length of superset
 * @param superset  superset
 * @param index     index of subset in superset
 * @return  the number of elements in superset blong to subset
 */
int char_match(int nsub, char** subset, int nsuper, char** superset, int* index)
{
  int nmatch = 0;
  for (int i = 0; i < nsuper; i++) {
    if (in_char(superset[i], subset, nsub) >= 0) {
      index[i] = 1;
      nmatch++;
    }
  }

  return nmatch;
}

/**
 * @brief Extract a subset of the string array based on the index, overwriting the original string
 * @param index     index of subset in superset
 * @param nsuper    number of elements in superset
 * @param superset  string array to be subsetted
 */
void char_subset(int* index, int nsuper, char** superset)
{
  int count = 0;
  for (int i = 0; i < nsuper; i++) {
    if (index[i]) {
      if (i != count) {
        // superset[count] = superset[i]; // May cause bugs
        strcpy(superset[count], superset[i]);
      }
      count++;
    }
  }
}

/**
 * @brief  output a data frame with header and row names
 * @param df       data frame
 * @param row      number of rows
 * @param col      number of columns
 * @param ld       leading dimension of matrix df
 * @param file     output file
 * @param header   header
 * @param rowname  row namess
 */
void output_ddataframe_names(double* df, int row, int col, int ld, char* file, char *header, char** rowname)
{
  FILE* fp = fopen(file, "wt");
  // header
  fprintf(fp, "%s\n", header);
  for (size_t i = 0; i < row; i++) {
    // rowname
    fprintf(fp, "%s\t", rowname[i]);
    for (size_t j = 0; j < col; j++) {
      if (j == col - 1) {
        fprintf(fp, "%E", df[ld * i + j]);
      }
      else {
        fprintf(fp, "%E\t", df[ld * i + j]);
      }
    }
    putc('\n', fp);
  }
  fclose(fp);
}

/**
 * @brief Find whether the value is in the vec, and return the position if it is
 * @param value  value
 * @param vec    vector
 * @param num    length of vector
 * @return  the position of value in vec, -1 if not found
 */
int which_int(int value, int *vec, int num)
{
  for (int i = 0; i < num; i++) {
    if (vec[i] == value) {
      return i;
    }
  }
  return -1;
}

/**
 * @brief  Append a double array to an existing file
 * @param fp   File to write
 * @param num  Length of array
 * @param vec  Array to write
 * @param inc  Increment between elements of vec
 * @param sep  Separator
 */
void append_dfile(FILE *fp, int num, double *vec, int inc, char *sep)
{
  for (int i = 0; i < num; i++) {
    if (i == num - 1) {
      fprintf(fp, "%E", vec[i * inc]);
    }
    else {
      fprintf(fp, "%E%s", vec[i * inc], sep);
    }
  }
  putc('\n', fp);
}

/**
 * @brief  Custom self_dgemm function, in order to compare with the mkl result
 * @param Layout Specifies whether two-dimensional array storage is row-major (CblasRowMajor) or column-major
 * (CblasColMajor)
 * @param TransA Specifies the form of op(A) used in the matrix multiplication
 * @param TransB Specifies the form of op(B) used in the matrix multiplication
 * @param M Specifies the number of rows of the matrix op(A) and of the matrix C. The value of m must be at least zero.
 * @param N Specifies the number of columns of the matrix op(B) and the number of columns of the matrix C. The value of
 * n must be at least zero.
 * @param K Specifies the number of columns of the matrix op(A) and the number of rows of the matrix op(B). The value of
 * k must be at least zero
 * @param alpha Specifies the scalar alpha.
 * @param A Array
 * @param lda Specifies the leading dimension of a as declared in the calling (sub)program
 * @param B Array
 * @param ldb Specifies the leading dimension of b as declared in the calling (sub)program
 * @param beta Specifies the scalar beta. When beta is equal to zero, then c need not be set on input.
 * @param C Array
 * @param ldc Specifies the leading dimension of c as declared in the calling (sub)program
 */
void self_dgemm(int Layout, int TransA, int TransB, int M, int N, int K, double alpha, double* A, int lda, double* B,
                int ldb, double beta, double* C, int ldc)
{
  if (Layout == 101) {
    if (TransA == 111 && TransB == 111) {
      for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {
          C[i * ldc + j] *= beta;
          for (int k = 0; k < K; k++) {
            C[i * ldc + j] += alpha * A[i * lda + k] * B[k * ldb + j];
          }
        }
      }
    }
  } else {
    printf("Error: Not implemented yet.\n");
    exit(UNDER_DEVELOP);
  }
}

/**
 * @brief  Check whether the file is readable
 * @param filename  file name
 * @param suffix    suffix of file name
 * @return  full file name
 */
char *check_file(char *filename, char *suffix) {
  // file name
  int prefix_len = strlen(filename);
  int suffix_len = strlen(suffix);
  char *full_filename = (char *)malloc((prefix_len + suffix_len + 1) * sizeof(char));
  sprintf(full_filename, "%s%s", filename, suffix);

  // read file
  FILE *fp = fopen(full_filename, "r");
  if (fp == NULL) {
    printf("Error: error in reading file: %s\n", full_filename);
    exit(FILE_NOT_FOUND);
  }
  fclose(fp);
  return full_filename;
}

/**
 * @brief  Output error information
 * @param message  error message
 * @param code     error code
 */
void error(char* message, int code)
{
  printf("%s\n", message);
  exit(code);
}

int min_i(int a, int b) {
  
  int min = (a > b) ? b : a;
  return min;
}

int max_i(int a, int b)
{
  int max = (a > b) ? a : b;
  return max;
}

double min_d(double a, double b)
{
  double min = (a < b) ? a : b;
  return min;
}

void strcpy_l(char** to, char* from)
{
  int len = strlen(from);
  *to = (char*)malloc((len + 1) * sizeof(char));
  if (to == NULL) {
    printf("error in copy strings!\n");
    exit(MEM_ALLOC);
  }
  strcpy(*to, from);
}

#ifdef _WIN32
    /**
     * @brief Count the number of lines in the file (Windows only)
     * @param  file             filename
     * @return int number of lines
     */
    int countLines(char* file) {
  FILE* fp;
  int rows = 0;
  char* BufContent = (char*)calloc(FILE_BLOCK_SIZE, sizeof(char));
  size_t BufContentSz;

  if ((fp = fopen(file, "rb")) == NULL) {
    printf("Can not load file: %s\n", file);
    exit(EXIT_FAILURE);
  }
  if (BufContent == NULL) {
    printf("The memory of size %d cannot be allocated in the count function.\n", FILE_BLOCK_SIZE);
    fclose(fp);
    return -2;
  } else {
    while ((BufContentSz = fread(BufContent, sizeof(unsigned char), FILE_BLOCK_SIZE, fp)) > 0) {
      for (int i = 0; i < BufContentSz; i++) {
        if (BufContent[i] == '\n') {
          rows++;
        }
      }
    }
    free(BufContent);
  }

  free(BufContent);
  return rows;
}
#endif

#ifdef linux
/**
 * @brief Count the number of lines in the file (Linux only)
 * @param  file             filename
 * @return int number of lines
 */
int countLines(char* file) {
  /* Check file readability */
  FILE* fp;
  if ((fp = fopen(file, "rt")) == NULL) {
    printf("can not open file %s\n", file);
    exit(FILE_NOT_FOUND);
  }
  fclose(fp);

  int rows = 0;
  char cmdstring[210];
  memset(cmdstring, 0, sizeof(cmdstring));
  char buff[1024];
  memset(buff, 0, sizeof(buff));
  FILE* fstream = NULL;

  /* Prepare bash command to output file line count */
  strcpy(buff, file);
  strcat(buff, " | awk '{print $1}'");
  strcat(cmdstring, "wc -l ");
  strcat(cmdstring, buff);

  if ((fstream = popen(cmdstring, "r")) == NULL) {
    fprintf(stderr, "execute command failed: %s", strerror(errno));
    return -1;
  }

  if (fgets(buff, sizeof(buff), fstream) != NULL) {
    rows = atoi(buff);
  } else {
    rows = -1;
  }

  pclose(fstream);
  return rows;
}
#endif

/**
 * @file multi_pop_fun.c
 * @brief declaration of functions for model evaluation
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
 * <tr><td>2023-04-03 <td>1.1.0     <td>liwn     <td>Add genotype QC code
 * <tr><td>2023-05-12 <td>1.1.1     <td>liwn     <td>frequency calculation within populations
 * </table>
 */

#include "multi_pop_fun.h"

#include "common_fun.h"

// -------------- definations of functions -------------- //

/**
 * @brief Generate incidence matrix according to the level of effect
 * @param  phe              Pheno structure
 * @param  model            Model structure
 * @param  output           Output structure
 */
void make_incidence(struct Pheno* phe, struct Model* model, struct Output* output)
{
  model->fixLs = (int*)calloc(model->nfix, sizeof(int));
  int** levels = (int**)calloc(model->nfix, sizeof(int*));      /* fixed effect level of individual */
  char*** label = (char***)calloc(model->nfix, sizeof(char**)); /* labels for each effect */
  int** num_lev = (int**)calloc(model->nfix, sizeof(int*));     /* Number of levels in each fixed effect */

  char** labels = (char**)calloc(phe->rows, sizeof(char*)); /* effect labels of an effect, nind x 1 */
  for (size_t j = 0; j < phe->rows; j++) {
    labels[j] = (char*)calloc(NCHAR_MAX, sizeof(char));
  }

  /* Number of levels of acquisition effect */
  model->fixL = 0;
  for (size_t i = 0; i < model->nfix; i++) {
    /* get the label of each individual effect i */
    read_dataframe_char(phe->filename, phe->rows, phe->cols, labels, 1, phe->fixCols + i, -1);

    /* remove fixed effect label in individuals with missing phenotypic values */
    char_subset(phe->pheIDInx, phe->rows, labels);

    /* get the level of each individual effect i */
    levels[i] = (int*)calloc(phe->nindPhe, sizeof(int));
    num_lev[i] = (int*)calloc(phe->nindPhe, sizeof(int));
    label[i] = unique_char(labels, phe->nindPhe, 1, model->fixLs + i, num_lev + i, levels[i]);

    /* Total number of effect levels */
    model->fixL += model->fixLs[i];
  }

  /* Save fixed effect labels */
  int count = 0;
  output->fix_labels = (char**)calloc(model->fixL, sizeof(char*));
  model->fix_nobs = (int*)calloc(model->fixL, sizeof(int));
  for (size_t i = 0; i < model->nfix; i++) {
    for (size_t j = 0; j < model->fixLs[i]; j++) {
      /* number of observation in i fix effect and j levels */
      model->fix_nobs[count] = num_lev[i][j];

      output->fix_labels[count] = (char*)calloc(NCHAR_MAX, sizeof(char));
      strcpy(output->fix_labels[count], label[i][j]);
      count++;

      /* free memory */
      free(label[i][j]);
    }
    free(label[i]);
  }
  free(label);

  /* Generate incidence matrix */
  int fixi;
  int col = model->fixL * model->npop; /* number of columns in X */
  model->X = (double*)calloc(phe->nindPhe * model->npop * col, sizeof(double));
  int rowi, coli;
  for (size_t i = 0; i < model->npop; i++) {
    fixi = 0;
    for (size_t j = 0; j < model->nfix; j++) {
      for (size_t k = 0; k < phe->nindPhe; k++) {
        rowi = k + i * phe->nindPhe;
        coli = i * model->fixL + fixi + levels[j][k];
        // It is assumed that all breed's fixed effects are the same currently
        model->X[rowi * col + coli] = 1.0;
      }
      fixi += model->fixLs[j];
    }
  }

  /* total levels of fixed effects (assume all breeds has the same fixed effects) */
  model->fixLa = model->npop * model->fixL;

  /* Release the requested memory */
  for (size_t i = 0; i < model->nfix; i++) free(levels[i]);
  free(levels);
  for (size_t i = 0; i < model->nfix; i++) free(num_lev[i]);
  free(num_lev);
  for (size_t i = 0; i < phe->rows; i++) free(labels[i]);
  free(labels);
}

/**
 * @brief Generate incidence matrix according to the level of effect
 * @param  phe              Pheno structure
 * @param  model            Model structure
 * @param  output           Output structure
 */
void make_incidence_old(struct Pheno* phe, struct Model* model, struct Output* output)
{
  model->fixLs = (int*)calloc(model->nfix, sizeof(int));
  int** levels = (int**)calloc(model->nfix, sizeof(int*));      /* fixed effect level of individual */
  char*** label = (char***)calloc(model->nfix, sizeof(char**)); /* labels for each effect */

  char** labels = (char**)calloc(phe->rows, sizeof(char*)); /* effect levels of an effect, nind x 1 */
  for (size_t j = 0; j < phe->rows; j++) {
    labels[j] = (char*)calloc(NCHAR_MAX, sizeof(char));
  }

  /* Number of levels of acquisition effect */
  model->fixL = 0;
  for (size_t i = 0; i < model->nfix; i++) {
    /* get the label of each individual effect i */
    read_dataframe_char_old(phe->filename, phe->rows, phe->cols, labels, 1, phe->fixCols + i, phe->header);

    /* remove fixed effect label in individuals with missing phenotypic values */
    char_subset(phe->pheIDInx, phe->rows, labels);

    /* get the level of each individual effect i */
    levels[i] = (int*)calloc(phe->nindPhe, sizeof(int));
    label[i] = unique_char_old(labels, phe->nindPhe, 1, model->fixLs + i, levels[i]);

    /* Total number of effect levels */
    model->fixL += model->fixLs[i];
  }

  /* Save fixed effect labels */
  int count = 0;
  output->fix_labels = (char**)calloc(model->fixL, sizeof(char*));
  for (size_t i = 0; i < model->nfix; i++) {
    for (size_t j = 0; j < model->fixLs[i]; j++) {
      output->fix_labels[count] = (char*)calloc(NCHAR_MAX, sizeof(char));
      strcpy(output->fix_labels[count], label[i][j]);
      count++;

      /* free memory */
      free(label[i][j]);
    }
    free(label[i]);
  }
  free(label);

  /* Generate incidence matrix */
  int fixi;
  int col = model->fixL * model->npop; /* number of columns in X */
  model->X = (double*)calloc(phe->nindPhe * model->npop * col, sizeof(double));
  int rowi, coli;

  for (size_t i = 0; i < model->npop; i++) {
    fixi = 0;
    for (size_t j = 0; j < model->nfix; j++) {
      for (size_t k = 0; k < phe->nindPhe; k++) {
        rowi = k + i * phe->nindPhe;
        coli = i * model->fixL + fixi + levels[j][k];
        // It is assumed that all breed's fixed effects are the same currently
        model->X[rowi * col + coli] = 1.0;
      }
      fixi += model->fixLs[j];
    }
  }

  /* total levels of fixed effects */
  model->fixLa = model->npop * model->fixL;

  /* Release the requested memory */
  for (size_t i = 0; i < model->nfix; i++) free(levels[i]);
  free(levels);
  for (size_t j = 0; j < phe->rows; j++) free(labels[j]);
  free(labels);
}

/**
 * @brief sampling residual effects for individuals with missing phenotypes
 * GIANOLA D, FERNANDO R L. A Multiple-Trait Bayesian Lasso for Genome-Enabled Analysis and Prediction of Complex Traits
 * [J]. Genetics, 2020,214(2): 305-331.
 * @param  model            Model structure
 * @param  stream           MKL random number generator
 * @include                 mkl.h
 */
void sampleMissRes(struct Model* model, struct RseSample* resmp, VSLStreamStatePtr stream)
{
  int k = 0;
  int ld = 0, inc = 0;
  int first_address = 0;
  int obs_loc;
  int ind_start = 0;
  int npop = model->npop;
  int nmiss = npop - 1;

  double R0_oo;
  // double* sample = (double*)calloc(model->nind * nmiss, sizeof(double));
  // double* tmp = (double*)calloc(nmiss, sizeof(double));
  // double* E_ym = (double*)calloc(nmiss, sizeof(double));
  // double* V_ym = (double*)calloc(nmiss * nmiss, sizeof(double));
  // double* R0_mo = (double*)calloc(nmiss, sizeof(double));
  // double* R0_mm = (double*) calloc(nmiss * nmiss, sizeof(double));

  // double* V_ymL = (double*)calloc(nmiss * nmiss, sizeof(double));  // debug 2

  // char buffer[100];  // debug 2

  for (size_t popi = 0; popi < npop; popi++) {
    obs_loc = which_int(1, model->obs + popi * npop, npop);

    // Re sub-matrix
    // memset(resmp->R0_mo, 0, nmiss * sizeof(double));
    // memset(resmp->R0_mm, 0, nmiss * nmiss * sizeof(double));
    subset(npop, npop, model->Re, model->obs + popi * npop, 1, model->obs + popi * npop, 1, resmp->R0_mm);
    subset(npop, npop, model->Re, model->obs + popi * npop, 0, model->obs + popi * npop, 0, &R0_oo);
    subset(npop, npop, model->Re, model->obs + popi * npop, 1, model->obs + popi * npop, 0, resmp->R0_mo);

    for (size_t missi = 0; missi < nmiss; missi++) {
      // R_{0}^{[miss, obs]}(R_{0}^{[obs, obs]})^{-1}
      resmp->tmp[missi] = resmp->R0_mo[missi] * (1 / R0_oo);
      for (size_t k = 0; k <= missi; k++) {
        // R_{0}^{[miss, obs]}(R_{0}^{[obs, obs]})^{-1}R_{0}^{[obs, miss]}
        resmp->V_ym[missi * nmiss + k] = resmp->tmp[missi] * resmp->R0_mo[k];
        // R_{0}^{[miss, miss]}-R_{0}^{[miss, obs]}(R_{0}^{[obs, obs]})^{-1}R_{0}^{[obs, miss]}
        resmp->V_ym[missi * nmiss + k] = resmp->R0_mm[missi * nmiss + k] - resmp->V_ym[missi * nmiss + k];
        resmp->V_ym[k * nmiss + missi] = resmp->V_ym[missi * nmiss + k];
      }
    }

    // sample
    // vdRngNor(model->npops[popi], nmiss, E_ym, V_ym, sample, stream); // last
    vdRngNor(model->npops[popi], nmiss, resmp->E_ym, resmp->V_ym, resmp->sample, stream);
    // memset(resmp->sample, 0, model->npops[popi] * nmiss * sizeof(double));  // debug 2
    // copy_d(resmp->V_ym, 1, V_ymL, 1, nmiss * nmiss, 1);                       // debug 2
    // k = LAPACKE_dpotrf(LAPACK_ROW_MAJOR, 'U', nmiss, V_ymL, nmiss);         // debug 2
    // for (size_t i = 0; i < nmiss; i++) {                                    // debug 2
    //   for (size_t j = 0; j < i; j++) {                                      // debug 2
    //     V_ymL[i * nmiss + j] = 0;                                           // debug 2
    //   }                                                                     // debug 2
    // }                                                                       // debug 2
    // self_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, model->npops[popi], nmiss, nmiss, 1.0,  // debug 2
    //            model->randn + model->randidx, nmiss, V_ymL, nmiss, 1.0, resmp->sample, nmiss);    // debug 2
    // model->randidx += model->npops[popi] * nmiss;                                        // debug 2
    // sprintf(buffer, "ReS%zuJ.txt", i + 1);                                    // debug
    // read_dataframe_double(buffer, model->npops[i], nmiss, resmp->sample, 0, &i, 0);  // debug
    // sprintf(buffer, "./data/ReS%d.txt", i + 1);                               // debug
    // output_ddataframe(resmp->sample, model->npops[i], nmiss, buffer, 1.0);           // debug

    k = 0;
    for (size_t popj = 0; popj < npop; popj++) {
      if (model->obs[popi * npop + popj]) {
        continue;
      }

      if (0) {
        ld = npop;
        inc = npop;
        first_address = popj * inc + ind_start * ld;
      }
      else {
        ld = 1;
        inc = model->nind;
        first_address = popj * inc + ind_start * ld;
      }

      // add the mean
      cblas_daxpy(model->npops[popi], resmp->tmp[k], model->yc + obs_loc * model->nind + ind_start, 1,
                  resmp->sample + k, nmiss);

      // assign
      // cblas_dcopy(model->npops[popi], resmp->sample + k, nmiss, model->yc + first_address, ld);
      first_address = popj * model->inc + ind_start * model->ld;
      cblas_dcopy(model->npops[popi], resmp->sample + k, nmiss, model->yc + first_address, model->ld);

      k++;
    }

    ind_start += model->npops[popi];
  }

  // free(sample);
  // free(tmp);
  // free(E_ym);
  // free(V_ym);
  // free(R0_mo);
  // free(R0_mm);
  // free(V_ymL); // debug 2
}

/**
 * @brief Make the mean value of each column of the gene content matrix equal to 0
 * @param  M                gene content matrix
 * @param  nmrk             number of marker in matrix
 * @param  nind             number of individuals in matrix
 * @param  trans            Matrix format, 1 for SNP x nind and 0 for nind x SNP
 */
void gtcenter(double* M, int nmrk, int nind, int trans)
/*
 * Centralized gene content matrix
 */
{
  int incx, start;
  if (trans == 1) {
    /* SNP x nind */
    incx = 1;
    start = nind;
  }
  else {
    /* nind x SNP */
    incx = nmrk;
    start = 1;
  }
  double mean;
  for (size_t i = 0; i < nmrk; i++) {
    mean = cblas_dasum(nind, M + i * start, incx) / nind;
    cblas_daxpy(nind, -1.0, &mean, 0, M + i * start, incx);
  }
}

/**
 * @brief Generating variance component matrix according to phenotype variance
 * @param  n                Dimension of square matrix
 * @param  pvar             Phenotypic variance vector, length n
 * @param  per              Proportion of generated variance components in phenotypic variance
 * @param  cor              Correlation coefficient used to generate covariance
 * @param  out              Output variance covariance component, quantity n ^ 2
 * @param  num              Not enabled yet
 * @note It is temporarily defaulted that the variance components of each bin are the same, that is, they do not
 * change with the number of snps
 */
void priorMat(int npop, double* pvar, double per, double cor, double* out, int nbin, double *sum2pq)
{
  int n_pp = npop * npop;
  double sum2pqi = 1.0;

  for (size_t i = 0; i < nbin; i++) {
    // variance
    for (size_t j = 0; j < npop; j++) {
      if (nbin > 1) sum2pqi = sum2pq[j * nbin + i];
      out[i * n_pp + j * npop + j] = pvar[j] * per / nbin / sum2pqi; /* same in all bins */
    }

    // covariance
    for (size_t j = 0; j < npop; j++) {
      for (size_t k = j + 1; k < npop; k++) {
        out[i * n_pp + j * npop + k] = cor * sqrt(out[i * n_pp + j * npop + j] * out[i * n_pp + k * npop + k]);
        out[i * n_pp + k * npop + j] = out[i * n_pp + j * npop + k];
      }
    }
  }
}

/**
 * @brief Direct inversion of matrix
 * @param  mat              Input matrix
 * @param  n                Dimension of matrix
 * @param  inv              Inverse of input matrix
 */
void inverse(double* mat, int n, double* inv)
{
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
 * @brief random number generation for the scaled Wishart distribution.
 * @param  S                This is the symmetric, positive-semidefinite p x p scale matrix S.
 * @param  df               This is the scalar degrees of freedom, nu.
 * @param  p                Dimension of matrix
 * @param  C                sample
 * @param  stream           random number generation
 * @include                 mkl.h
 */
void wishart(double* S, double df, int p, double* C, VSLStreamStatePtr stream)
{
  double L[p * p];
  cholesky(S, p, L);
  double At[p * p];
  double ranChi = 0.0;
  int status = 0;

  /* At is the upper triangular array */
  for (int i = 0; i < p; i++) {
    /* Diagonal elements with decreasing degrees of freedom (df-i) */
    status = vdRngChiSquare(METHOD_C, stream, 1, &ranChi, df - i);

    if (status != 0) {
      printf("Error in vdRngChiSquare, status = %d", status);
    }

    At[i * p + i] = sqrt(ranChi);
    /* Non diagonal elements, standard normal distribution */
    for (int j = i + 1; j < p; j++) {
      vdRngGaussian(METHOD_N, stream, 1, At + i * p + j, 0, 1);
    }
  }

  /* Initialize the triangular array under the B-matrix initialization */
  double B[p * p];
  for (int i = 0; i < p; i++)
    for (int j = 0; j < (i + 1); j++) B[i * p + j] = 0;

  /* B = chol(S) x At */
  for (int i = 0; i < p; i++)
    for (int j = 0; j < (i + 1); j++)
      for (int k = j; k < (i + 1); k++) B[i * p + j] += L[i * p + k] * At[j * p + k];

  /* C matrix initialization */
  for (int i = 0; i < p; i++)
    for (int j = 0; j < (i + 1); j++) C[i * p + j] = 0;

  /* C = B x B */
  for (int i = 0; i < p; i++)
    for (int j = 0; j < (i + 1); j++)
      for (int k = 0; k < (j + 1); k++) C[i * p + j] += B[i * p + k] * B[j * p + k];

  /* Fill C into a full storage matrix */
  for (int i = 0; i < p; i++)
    for (int j = i + 1; j < p; j++) C[i * p + j] = C[j * p + i];
}

/**
 * @brief get the scale inverse Wishart distribution samples based on the scale Wishart distribution samples
 * @param  Sigma            the symmetric, positive-semidefinite p x p scale matrix S.
 * @param  df               the scalar degrees of freedom, nu.
 * @param  p                Dimension of matrix
 * @param  C                sample
 * @param  stream           random number generation
 * @include                 mkl.h
 */
void invWishart(double* Sigma, double df, int p, double* C, VSLStreamStatePtr stream)
{
  double invV[p * p];
  inverse(Sigma, p, invV);
  double invC[p * p];
  /* Obtain samples of scale Wishart distribution */
  wishart(invV, df, p, invC, stream);
  inverse(invC, p, C);
}

/**
 * @brief  inversion of the symmetric (Hermitian) positive-definite matrix using the Cholesky factorization.
 * @param  n                The order of the matrix mat; n ≥ 0.
 * @param  mat              Input matrix
 * @param  inv              Overwritten by the upper triangle of the inverse of mat.
 * @include                 mkl.h
 */
void solve(int n, double* mat, double* inv)
{
  int status = 0;
  cblas_dcopy(n * n, mat, 1, inv, 1);
  status = LAPACKE_dpotrf(LAPACK_ROW_MAJOR, 'U', n, inv, n);
  status = LAPACKE_dpotri(LAPACK_ROW_MAJOR, 'U', n, inv, n);

  if (status) {
    printf(
        "An error occurred in the inverse operation of the matrix using the "
        "Cholesky factorization.\n");
    exit(status);
  }
  else {
    /* Fill in the lower triangle element */
    for (size_t row = 0; row < n; row++) {
      for (size_t col = 0; col < row; col++) {
        inv[row * n + col] = inv[col * n + row];
      }
    }
  }
}

/**
 * @brief random number generation for the inverse scaled Wishart distribution.
 * @param  S                This is the symmetric, positive-semidefinite p x p scale matrix S (row major).
 * @param  df               This is the scalar degrees of freedom, nu.
 * @param  n                Dimension of matrix
 * @param  sample           sample, n x n row major matrix
 * @param  stream           random number generation
 * @include                 mkl.h
 * @ref                     Gelman, A., Carlin, J., Stern, H., and Rubin, D. (2004). "Bayesian Data Analysis, Texts in
 * Statistical Science, 2nd ed.". Chapman and Hall, London; Wishart, J. (1928). "The Generalised Product Moment
 * Distribution in Samples from a Normal Multivariate Population". Biometrika, 20A(1-2), p. 32–52.
 */
void vdRngInvWishart(double* S, double df, int n, double* sample, VSLStreamStatePtr stream)
{
  int status;
  double* Z = (double*)calloc(n * n, sizeof(double));
  double* cholU = (double*)calloc(n * n, sizeof(double));

  /* Cholesky decomposition */
  status = LAPACKE_dpotrf(LAPACK_ROW_MAJOR, 'U', n, S, n);

  if (status) {
    printf(
        "An error occurred in the inverse operation of the matrix using the "
        "Cholesky factorization.\n");

    /* Release the requested memory */
    free(Z);
    free(cholU);

    exit(status);
  }
  else {
    /* Inverse Cholesky Matrix */
    status = LAPACKE_dpotri(LAPACK_ROW_MAJOR, 'U', n, S, n);

    /* Perform Cholesky decomposition on the inverse matrix again */
    status = LAPACKE_dpotrf(LAPACK_ROW_MAJOR, 'U', n, S, n);

    /* Determine whether A is positive definite (square matrix, symmetric, eigenvalues are all greater than 0) */

    /* Diagonal element: chi-square distribution sampling */
    for (size_t i = 0; i < n; i++) {
      // status = vdRngChiSquare(METHOD_C, stream, 1, Z + i * n + i, df - i +
      // 1);
      status = vdRngChiSquare(METHOD_C, stream, 1, Z + i * n + i, df - i);
      Z[i * n + i] = sqrt(Z[i * n + i]);
    }

    /* Non diagonal elements: normal distribution sampling (top triangle) */
    for (size_t col = 0; col < n; col++) {
      for (size_t row = 0; row < col; row++) {
        vdRngGaussian(METHOD_N, stream, 1, Z + row * n + col, 0, 1);
      }
    }

    /* Set the lower triangle element of A to 0 */
    for (size_t row = 0; row < n; row++) {
      for (size_t col = 0; col < row; col++) {
        S[row * n + col] = 0.0;
      }
    }

    /* Multiplying by Cholesky matrix (upper triangle) */
    memset(sample, 0, n * n * sizeof(double));
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1.0, Z, n, S, n, 1.0, sample, n);

    /* Inverse Wishart samples */
    status = LAPACKE_dpotri(LAPACK_ROW_MAJOR, 'U', n, sample, n);

    /* Supplement the triangular elements under the matrix (fully stored) */
    for (size_t row = 0; row < n; row++) {
      for (size_t col = 0; col < row; col++) {
        sample[row * n + col] = sample[col * n + row];
      }
    }

    /* Release the requested memory */
    free(Z);
    free(cholU);
  }
}

/**
 * @brief Obtain phenotype missing type, such as 1.23, - 99, - 99, then the type is 1, 0, 0, which is used to
 * impute in the residual effect of the deletion
 * @note At present, individuals in the default phenotype file are sorted by varieties, and each variety has
 * only one deletion mode, that is, if there are three varieties, there are only three types (0 0 1; 0 1 0; 1 0 0)
 * @param  npop             number of breeds
 * @param  ninds            number of individuals in each breed
 * @param  y                phenotype matrix/vector
 * @param  miss             Missing phenotypic value
 * @param  type             format of phenotype stored, 1 for matrix (nind x nbreeds) and 0 for vectors
 * @param  obs_index        index of observed individuals, e.g. 1, 0, 0, 0, 1, 0, 0, 0, 1
 */
void obs_type(int npop, int* ninds, double* y, double miss, int type, int* obs_index)
{
  int first_ind = 0;
  int index;
  int nind = sumi(ninds, npop, 1);
  for (size_t i = 0; i < npop; i++) {
    if (i > 0) {
      first_ind += ninds[i - 1];
    }
    for (size_t j = 0; j < npop; j++) {
      if (type) {
        index = first_ind * npop + j;
      }
      else {
        index = first_ind + j * nind;
      }

      obs_index[i * npop + j] = y[index] != miss;
    }
  }
}

/**
 * @brief Grouping the snps on each chromosome based on the fixed number of snps in each bins
 * @param  geno      genotype parameters structure
 * @param  bin     return: Number of bins defined
 */
void make_bins_fix(Gmat* M, Bin* bin)
{
  /* Number of bins */
  double num = bin->nsnp;
  bin->nbins = (int*)calloc(M->nchr, sizeof(int));
  bin->nbin = 0;
  for (size_t i = 0; i < M->nchr; i++) {
    bin->nbins[i] = round(M->nsnps[i] / num);
    bin->nbin += bin->nbins[i];
  }

  /* Store results */
  int bini = 0;
  bin->nsnps = (double*)calloc(bin->nbin, sizeof(double));
  for (size_t i = 0; i < M->nchr; i++) {
    for (size_t j = 0; j < bin->nbins[i]; j++) {
      if (j == (bin->nbins[i] - 1)) {
        bin->nsnps[bini] = (double)M->nsnps[i] - num * (bin->nbins[i] - 1);
      }
      else {
        bin->nsnps[bini] = num;
      }
      bini++;
    }
  }
}

/**
 * @brief  Calculate the r2 value between SNPs
 * @param M       gene content matrix
 * @param chr     Chromosome requiring SNP grouping
 * @param start   Index of the first SNP of chromosome chr in the M
 * @param nmrk    number of SNPs in chromosome chr
 * @param win     R2 will be computed between SNPs that are at most this far apart [100]
 * @param mean    R2 mean value within a wedge triangle
 * @return
 */
void compute_r2(Gmat* M, int chr, int start, int nmrk, int win, double* mean)
{
  /* format of matrix */
  gmat_format(M, 1);

  /* Whether the M has been centralized */
  double* frq_dummy;
  if (M->center) {
    frq_dummy = (double*)calloc(nmrk, sizeof(double));
  }
  else {
    frq_dummy = M->frq + start;
  }

  int snp_end;
  double r, sum;
  int row_start, col_end, num_all;
  // int k;
  // int num = (nmrk - win - 1) * win + win * (win + 1) / 2;
  double* r2mat = (double*)calloc(nmrk * win, sizeof(double));
  for (int i = 0; i < nmrk - 1; i++) { /* The last snp has only diagonal elements */
    snp_end = (i + win < nmrk - 1) ? i + win : nmrk - 1;
    for (int j = i + 1, k = 0; j <= snp_end; j++, k++) {
      r = gene_content_correlation(M->M + (i + start) * M->ld, frq_dummy[i], M->inc, M->M + (j + start) * M->ld,
                                   frq_dummy[j], M->inc, M->nind);
      r2mat[i * win + k] = pow(r, 2);
    }

    /* Calculate the mean r2 in oblique triangle */
    num_all = 0;
    sum = 0.0;
    row_start = max_i(i - win + 1, 0);
    for (int j = i, k = 0; j >= row_start; j--, k++) {
      col_end = min_i(nmrk - j - 1, win);
      for (int l = k; l < col_end; l++) {
        /* code */
        sum += r2mat[j * win + l];
        num_all++;
        // col_start--;
      }
    }
    mean[i] = sum / num_all;
  }

  if (M->center) free(frq_dummy);
  free(r2mat);
}

/**
 * @brief   Calculate the correlation coefficient between two centered vectors
 * @param va    vector a
 * @param inca  increment for the elements of a
 * @param vb    vector a
 * @param incb  increment for the elements of b
 * @param n     number of elements
 * @return
 */
double gene_content_correlation(double* va, double frqa, int inca, double* vb, double frqb, int incb, int n)
{
  double cov = 0.0, vara = 0.0, varb = 0.0, deva = 0.0, devb = 0.0;
  for (int i = 0; i < n; i++) {
    deva = va[i * inca] - 2 * frqa;
    devb = vb[i * incb] - 2 * frqb;
    vara += pow(deva, 2);
    varb += pow(devb, 2);
    cov += deva * devb;
  }

  if (vara * varb == 0) {
    return 0;
  }
  else {
    return cov / sqrt(vara * varb);
  }
}

/**
 * @brief Grouping the snps on each chromosome based on LD information
 * @param  M         genotype matrix
 * @param  win       R2 will be computed between SNPs that are at most this far apart [100]
 * @param  min       minimum number of SNPs in each bin
 * @return Number of SNPs in each bin
 */
void make_bins_LD(Gmat* M, Bin* bin)
{
  /* container */
  int num_max = (M->nmrk + bin->min) / bin->min;
  double* nsnps_tmp = (double*)calloc(num_max, sizeof(double));

  /* bins define */
  int bini = 0;
  int snp_index = 0;
  double* r2mean = (double*)calloc(M->nmrk, sizeof(double));
  for (size_t i = 0; i < M->nchr; i++) {
    if (M->nsnps[i] <= 2 * bin->min) {
      nsnps_tmp[bini] = M->nsnps[i];
      bini++;
    }
    else {
      memset(r2mean, 0, M->nmrk);
      compute_r2(M, i, snp_index, M->nsnps[i], bin->win, r2mean);
      // output_ddataframe(r2mean, M->nsnps[i], 1, "/public/home/liujf/liwn/code/GitHub/data/ld_mean_tmp.txt", 1.0);
    }
    snp_index += M->nsnps[i];
  }

  printf("Under development\n");
  exit(UNDER_DEVELOP);

  bin->nsnps = (double*)calloc(bini, sizeof(double));
  for (size_t i = 0; i < bini; i++) {
    bin->nsnps[i] = nsnps_tmp[i];
  }

  free(nsnps_tmp);
  free(r2mean);
}

/**
 * @brief  Used to group snps on the genome into different snp groups
 * @param  geno    genotype parameters structure
 * @param  M       gene content matrix structure
 * @param  bin     bin information structure
 */
void make_bins(struct Genotype* geno, Gmat* M, Bin* bin)
{
  /* need snp information */
  if (M->nchr == 0) {
    if (strlen(geno->mapf) > 0) {
      read_bim(geno->prefix, geno->mapf, M);
    }
    else {
      printf("Error: Please provide bins definition file or plink *.map (*.bim) file.\n");
      exit(ERROR_PARA);
    }
  }

  if (strcmp(bin->type, "LD") == 0) {
    make_bins_LD(M, bin);
  }
  else if (strcmp(bin->type, "fix") == 0) {
    make_bins_fix(M, bin);
  }
  else {
    printf("Error: Please provide bins definition file or plink *.map (*.bim) file.\n");
    exit(ERROR_PARA);
  }
}

/**
 * @brief Change the single column phenotypic values to be stored in multiple traits form
 * @param  vec              input phenotypic values (nind x 1)
 * @param  npop             number of breeds
 * @param  npops            number of individuals in each breeds
 * @param  miss             missing value
 * @param  multi            output phenotypes, mat
 * @param  type             stored type of output phenotypes, mat/vec
 */
void one_col_to_multi(double* vec, int npop, int* npops, double miss, double* multi, int type)
{
  int nind = sumi(npops, npop, 1);
  int indi = 0;
  int index;
  // save phenotypes

  // number of individula with phenotype
  // int nphe = count(nind, pheno, miss, 0);

  // impute missing values
  fill_d(miss, multi, nind * npop, 1);

  // fill
  for (size_t i = 0; i < npop; i++) {
    indi = sumi(npops, i, 1);
    for (size_t j = 0; j < npops[i]; j++) {
      if (type == 1) {
        index = (j + indi) * npop + i;
      }
      else if (type == 0) {
        index = j + indi + i * nind;
      }
      else {
        printf("Error: type should be mat or vec");
        exit(ERROR_PARA);
      }

      multi[index] = vec[j + indi];
    }
  }
}

/**
 * @brief Change the single column phenotypic values to be stored in multiple traits form
 * @param  y                input phenotypic values (nind x 1)
 * @param  npop             number of breeds
 * @param  npops            number of individuals in each breeds
 * @param  miss             missing value
 * @param  type             stored type of output phenotypes, mat/vec
 */
void one_phe_to_multi_new(double* y, int npop, int* npops, double miss, int type)
{
  int nind = sumi(npops, npop, 1);
  int indi = 0;
  int index;

  /*  */
  double *y1c = (double*)calloc(nind, sizeof(double));
  copy_d(y1c, 1, y, 1, nind, 1);

  /*  */
  y = realloc(y, nind * npop * sizeof(double));

  // impute missing values
  fill_d(miss, y, nind * npop, 1);

  // fill
  for (size_t i = 0; i < npop; i++) {
    indi = sumi(npops, i, 1);
    for (size_t j = 0; j < npops[i]; j++) {
      if (type == 1) {
        index = (j + indi) * npop + i;
      }
      else if (type == 0) {
        index = j + indi + i * nind;
      }
      else {
        printf("Error: type should be mat or vec");
        exit(ERROR_PARA);
      }

      y[index] = y1c[j + indi];
    }
  }
}

/**
 * @brief update left and right hand side in MME
 * @param  iter                iteration number
 * @param  model               model parameters structure
 */
void left_right_hand_update(struct Model* model, int iter)
{
  int np = model->npop;
  int nind = model->nind;
  const int npind = np * nind;
  int fixL = model->fixL;

  // initialize
  double* XpRi = (double*)calloc(np * fixL * np * nind, sizeof(double));
  double* Ri = (double*)calloc(npind * npind, sizeof(double));
  memset(model->Rhs, 0, model->fixLa * sizeof(double));
  memset(model->Lhs, 0, model->fixLa * model->fixLa * sizeof(double));

  if (iter < 2) {
    int ind_s = 0;
    int row = 0;
    for (size_t i = 0; i < model->npop; i++) {
      if (i > 0) {
        ind_s += model->npops[i - 1];
      }

      for (size_t j = 0; j < model->npops[i]; j++) {
        row = i * nind + ind_s + j;
        Ri[row * npind + row] = model->Re_inv[i * model->npop + i];
      }
    }

    cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, np * fixL, npind, npind, 1.0, model->X, np * fixL, Ri, npind,
                1.0, XpRi, npind);

    // left hand side
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, np * fixL, np * fixL, npind, 1.0, XpRi, npind, model->X,
                np * fixL, 1.0, model->Lhs, np * fixL);

    // right hand side
    cblas_dgemv(CblasRowMajor, CblasNoTrans, np * fixL, npind, 1.0, XpRi, npind, model->yc, 1, 1.0, model->Rhs, 1);
  }
  else {
    // left hand side
    kronecker(np, np, model->Re_inv, model->fixL, model->fixL, model->X1pX1, model->fixL, model->Lhs);
    // right hand side (X'*Ri*yc = R0 ⊗ X1p * yc)
    kronecker(np, np, model->Re_inv, nind, fixL, model->X, np * fixL, XpRi);
    cblas_dgemv(CblasRowMajor, CblasTrans, np * nind, np * fixL, 1.0, XpRi, np * fixL, model->yc, 1, 1.0, model->Rhs,
                1);
  }

  free(XpRi);
  free(Ri);
}

/**
 * @brief sample fixed effects
 * @param  Lhs              left hand side in MME
 * @param  Rhs              right hand side in MME
 * @param  theta            estimated parameters
 * @param  n                number of parameters
 * @param  stream           stream for random number generation
 */
void fix_sample(struct Model* model, VSLStreamStatePtr stream)
{
  double mu, inv, dot, Lhs_ii;
  for (int i = 0; i < model->fixLa; i++) {
    Lhs_ii = model->Lhs[i * model->fixLa + i];
    if (Lhs_ii != 0.0) {
      inv = 1.0 / Lhs_ii;
      // dot = cblas_ddot(n, Lhs + i, n, theta, 1);
      dot = ddot2(model->fixLa, model->Lhs + i, model->fixLa, model->theta, 1);
      mu = inv * (model->Rhs[i] - dot) + model->theta[i];
      // rnorm(1, &mu, &inv, theta + i, stream);
      // vdRngNor(1, 1, &mu, &inv, model->theta + i, stream);
      model->theta[i] = rnorm01() * sqrt(inv) + mu;
      // model->theta[i] = model->randn[model->randidx] * sqrt(inv) + mu;  // debug 2
      // model->randidx++;                                                 // debug 2
    }
  }
}

/**
 * @brief update corrected phenotypic values
 * @param  yc               corrected phenotypic values
 * @param  mat              incidence matrix of fixed effects
 * @param  vec              estimated parameters of fixed effects
 * @param  m                number of rows in incidence matrix
 * @param  n                number of columns in incidence matrix
 * @param  sign             -1 for minus and 1 for plus
 */
void yc_local_update(double* yc, double* X, double* sol, const MKL_INT npind, const MKL_INT fixLa, const double sign)
{
  cblas_dgemv(CblasRowMajor, CblasNoTrans, npind, fixLa, sign, X, fixLa, sol, 1, 1.0, yc, 1);
}

/**
 * @brief update corrected phenotypic values
 * @param  model            model parameters structure
 * @param  snp_pos          snp index
 * @param  sign             -1 for minus and 1 for plus
 */
void yc_alpha_update(struct Model* model, double* snp_last, int snp_pos)
{
  double dvi;
  for (size_t i = 0; i < model->npop; i++) {
    dvi = snp_last[i] - model->theta[model->fixLa + snp_pos * model->npop + i];
    cblas_daxpy(model->nind, dvi, model->Mr + snp_pos * model->nind, 1, model->yc + i * model->nind, 1);
  }
}

/**
 * @brief
 * @param  model            model parameters structure
 * @param  stream           My Param doc
 */
void variance_scale(struct Model* model, VSLStreamStatePtr stream)
{
  int index;
  int n_pp = model->npop * model->npop;
  double* wishartScale = (double*)calloc(n_pp, sizeof(double));
  memset(model->tmp_pp, 0, n_pp * sizeof(double));

  for (size_t i = 0; i < model->npop; i++) {
    /* priori is a unit matrix */
    model->tmp_pp[i * model->npop + i] = 1.0;
    for (size_t j = 0; j < (i + 1); j++) {
      index = i * model->npop + j;
      /* Sum of SNP effect variance */
      model->tmp_pp[index] += cblas_ddot(model->nbin, model->nsnps, 1, model->VaInv + index, n_pp);
      model->tmp_pp[j * model->npop + i] = model->tmp_pp[index];
    }
  }
  solve(model->npop, model->tmp_pp, wishartScale);
  wishart(wishartScale, model->df_va * model->nmrk + model->npop, model->npop, model->priVaS, stream);

  free(wishartScale);
}

/**
 * @brief calculate the gene frequency and sum of the 2pq
 * @param  M       gene content matrix
 * @param  format  Storage format of M, 1 for mrk x ind and 0 for ind x mrk
 * @param  nind    number of individuals
 * @param  nmrk    number of markers
 * @param  frq     minor allele frequency
 * @param  center  Whether to centralize each column in the gene content matrix, 1 for yes and 0 for no
 * @include        mkl.h
 */
void gene_frq(Gmat* M)
{
  // Matrix storage format
  gmat_format(M, 1);

  /* calculate gene frequency */
  double frqi, mean, count;
  int indj;
  for (size_t i = 0; i < M->nmrk; i++) {
    for (size_t j = 0, indi = 0; j < M->npop; j++) {
      mean = 0.0;
      for (size_t k = 0; k < M->ninds[j]; k++, indi++) {
        /* allele content */
        count = M->M[i * M->ld + indi * M->inc];
        if (count != MISS_GENE_CODE) {
          mean += count;
        }
      }
      mean /= (M->ninds[j] - M->mrk_miss[i + j * M->nmrk]);
      frqi = mean / 2.0;

      /* Centralization of gene content */
      indj = sumi(M->ninds, j, 1);
      if (M->center) {
        cblas_daxpy(M->ninds[j], -1.0, &mean, 0, M->M + i * M->ld + indj * M->inc, M->inc);
      }

      M->frq[i + j * M->nmrk] = frqi;
    }
  }
}

/**
 * @brief Point multiplication of matrix column vectors
 * @details mat = [e_1, e_2, e_3, ... , e_row], ei = [e_i1, e_i2, ... , e_col], S =
 * \sum_{i=1}^{row}\boldsymbol e_i\boldsymbol e_i\prime
 * @param  mat              input matrix
 * @param  row              number of rows in input matrix
 * @param  col              number of columns in input matrix
 * @param  S                inner product of matrix, col x col matrix
 * @include                 mkl.h
 */
void mat_cblas_ddot(double* mat, int row, int col, double* S)
{
  for (size_t i = 0; i < col; i++) {
    for (size_t j = 0; j <= i; j++) {
      /* cblas_ddot There is a difference in the order of 1e-12 between manual loop addition and manual loop addition */
      S[i * col + j] = cblas_ddot(row, mat + i, col, mat + j, col);
      S[j * col + i] = S[i * col + j];
    }
  }
}

/**
 * @brief Sampling of marking effects
 * @param model            model parameters structure
 * @param snp_pos          snp position
 * @param bin_i            bin index
 * @param C                covariance matrix
 * @param wRinv            inverse of weighted residual
 */
void alpha_sample(struct Model* model, int snp_pos, int bin_i, double* C, double* wRinv, VSLStreamStatePtr stream)
{
  double wRinvRinvj;
  double mu;
  double sigma;
  double dvi;
  int npop = model->npop;

  double mpm = model->effsqs[model->fixLa + snp_pos];
  double* alpha = model->theta + model->fixLa + snp_pos * npop;
  double* Ginv = model->VaInv + bin_i * npop * npop;

  // double snRand[npop];
  memset(model->tmp_p, 0, npop * sizeof(double));

  /* Convert to a normal distribution with specified mean and variance (Y = (L')^{-1}X + μ ) */
  for (size_t i = 0; i < npop; i++) {
    wRinvRinvj = 0.0;
    for (size_t j = 0; j < npop; j++) {
      wRinvRinvj += wRinv[j] * model->Re_inv[i * npop + j];
      if (j == i) {
        continue;
      }
      else {
        model->tmp_p[i] += (Ginv[i * npop + j] + mpm * model->Re_inv[i * npop + j]) * alpha[j];
      }
    }

    // Distributed parameter
    sigma = 1.0 / C[i * npop + i];
    mu = (wRinvRinvj - model->tmp_p[i]) * sigma;

    // Last Effect Value
    dvi = alpha[i];

    // sample from normal distribution
    // vdRngNor(1, 1, &mu, &sigma, alpha + i, stream);
    alpha[i] = mu + rnorm01() * sqrt(sigma);
    // alpha[i] = mu + model->randn[model->randidx] * sqrt(sigma);  // debug 2
    // model->randidx++;                                            // debug 2

    // update adjusted phenotypes
    dvi -= alpha[i];
    cblas_daxpy(model->nind, dvi, model->Mr + snp_pos * model->nind, 1, model->yc + i * model->nind, 1);
    // for (size_t k = 0; k < (i + 1); k++) alpha[k] += VchoL[k * npop + k] * snRand[k];
  }
}

/**
 * @brief Obtain random samples from normal distribution
 * @param num             number of samples
 * @param d               dimension of distribution
 * @param mu              mean of distribution
 * @param Sigma           variance of distribution
 * @param sample          random samples
 * @param stream          MKL random number generator
 */
void vdRngNor(int num, int d, double* mu, double* Sigma, double* sample, VSLStreamStatePtr stream)
{
  int status = 0;
  if (d > 1) {
    double* VchoL = (double*)calloc(d * d, sizeof(double));
    copy_d(Sigma, 1, VchoL, 1, d * d, 1);
    status = LAPACKE_dpotrf(LAPACK_ROW_MAJOR, 'L', d, VchoL, d);
    if (status != 0) {
      printf("Error in cholesky decomposition!\n");
      exit(status);
    }

    status = vdRngGaussianMV(METHOD_N, stream, num, sample, d, STOR_N, mu, VchoL);
    if (status != 0) {
      printf("Error in normal distribution sampling!\n");
      exit(status);
    }
    free(VchoL);
  }
  else {
    status = vdRngGaussian(METHOD_N, stream, num, sample, *mu, sqrt(*Sigma));
  }
}

/**
 * @brief  Read phenotypes from the provided file
 * @param pheno    phenotype parameters structure
 * @param model    model parameters structure
 * @param phe_pos  Location of phenotypic information in the file, e.g. "2 3"
 * @param miss     Less than or equal to this value will be considered a missing phenotypic value
 */
void load_phenotype(struct Pheno* pheno, struct Model* model, char* phe_pos, double miss)
{
  pheno->type = 0; /* 0 for vector and 1 for matrix */

  /* Location of phenotypic information in the file */
  if (strlen(phe_pos) == 0) {
    // By default, the last n columns of the phenotype file are phenotype
    // pheno->triCols = (int*)calloc(model->npop, sizeof(int));
    // for (size_t i = 0; i < model->npop; i++) {
    //   pheno->triCols[i] = pheno->cols - model->npop + i + 1;
    // }
    // pheno->ntriCol = model->npop;

    // By default, the last 1 columns of the phenotype file are phenotype
    pheno->triCols = (int*)calloc(1, sizeof(int));
    pheno->triCols[0] = pheno->cols;
    pheno->ntriCol = 1;
  }
  else {
    pheno->triCols = string_to_array(phe_pos, &pheno->ntriCol, " ");
  }

  /* prepare container */
  pheno->y = (double*)calloc(pheno->nind * model->npop, sizeof(double));
  pheno->pheID = (char**)calloc(pheno->nind, sizeof(char*));  // surplus
  read_dataframe_double(pheno->filename, pheno->rows, pheno->cols, pheno->y, pheno->ntriCol, pheno->triCols, -1);
  if (model->npop > 1) {
    if (pheno->ntriCol == 1) {
      single_phe_to_vec(pheno, model, miss);
    }
    else if (pheno->ntriCol == model->npop) {
      multi_phe_to_vec(pheno, model, miss);
    }
    else {
      printf("Error: The number of phenotypes is not equal to 1 or the number of populations!\n");
      exit(ERROR_PARA);
    }
  }
  else {
    single_phe(pheno, model, miss);
  }

  /* storage format of phenotype */
  if (pheno->type) {
    /* nind x npop -- popj * inc + ind_start * ld */
    model->ld = model->npop;
    model->inc = 1;
  }
  else {
    /* (nind x npop) x 1 */
    model->ld = 1;
    model->inc = pheno->nindPhe;
  }
}

/**
 * @brief
 * @param phe
 * @param model
 * @param geno
 * @param miss
 */
void multi_phe_to_vec(struct Pheno* phe, struct Model* model, double miss)
{
  double value;
  int nind_phe = 0;

  double* y = (double*)calloc(phe->rows, sizeof(double));
  model->npops = (int*)calloc(model->npop, sizeof(int));

  for (size_t i = 0; i < model->npop; i++) {
    for (size_t j = 0; j < phe->rows; j++) {
      value = phe->y[j * model->npop + i];
      if (value > miss) {
        phe->pheID[nind_phe] = phe->ID[j];
        y[nind_phe] = value;
        nind_phe++;
        model->npops[i]++;
      }
    }
  }
  phe->nindPhe = nind_phe;

  // index of individuals with phenotype
  phe->pheIDInx = (int*)calloc(phe->rows, sizeof(int));
  char_match(phe->nindPhe, phe->pheID, phe->rows, phe->ID, phe->pheIDInx);

  // phenotype vector
  model->yc = (double*)calloc(model->npop * nind_phe, sizeof(double));
  one_col_to_multi(y, model->npop, model->npops, miss, model->yc, phe->type);

  free(y);
}

/**
 * @brief
 * @param phe
 * @param model
 * @param geno
 * @param miss
 */
void single_phe_to_vec(struct Pheno *phe, struct Model *model, double miss)
{
  Gmat gmat = model->geno->Ma;
  double value;
  int nind_phe = 0;

  double* y = (double*)calloc(phe->rows, sizeof(double));
  model->npops = (int*)calloc(model->npop, sizeof(int));

  /* number of individuals */
  if (gmat.nind != phe->rows) {
    printf("Error: The number of individuals in genotype files and phenotype files is not equal.\n");
    exit(ERROR_INPUT);
  }

  int ind_idx = 0;
  for (size_t i = 0; i < model->npop; i++) {
    for (size_t j = 0; j < gmat.ninds[i]; j++) {
      value = phe->y[ind_idx];
      if (value > miss) {
        phe->pheID[nind_phe] = phe->ID[j];
        y[nind_phe] = value;
        nind_phe++;
        model->npops[i]++;
      }
      ind_idx++;
    }
  }
  phe->nindPhe = nind_phe;

  // index of individuals with phenotype
  phe->pheIDInx = (int*)calloc(phe->rows, sizeof(int));
  char_match(phe->nindPhe, phe->pheID, phe->rows, phe->ID, phe->pheIDInx);

  // phenotype vector
  model->yc = (double*)calloc(model->npop * nind_phe, sizeof(double));
  one_col_to_multi(y, model->npop, model->npops, miss, model->yc, phe->type);

  free(y);
}

/**
 * @brief
 * @param phe
 * @param model
 * @param geno
 * @param miss
 */
void single_phe(struct Pheno* phe, struct Model* model, double miss)
{
  phe->nindPhe = 0;
  double* y = (double*)calloc(phe->nind, sizeof(double));
  model->npops = (int*)calloc(1, sizeof(int));
  phe->pheIDInx = (int*)calloc(phe->nind, sizeof(int));

  // Number of phenotypic individuals
  for (size_t j = 0; j < phe->nind; j++) {
    if (phe->y[j] > miss) {
      phe->pheID[phe->nindPhe] = phe->ID[j];
      y[phe->nindPhe] = phe->y[j];
      phe->nindPhe++;
      model->npops[0]++;
      phe->pheIDInx[j] = 1;
    }
  }

  // phenotype vector
  model->yc = (double*)calloc(phe->nindPhe, sizeof(double));
  copy_d(y, 1, model->yc, 1, phe->nindPhe, 1);

  free(y);
}

/**
 * @brief  Extract the gene content matrix of the reference populations from the complete gene content matrix M
 * @param  nref      number of individuals in reference population
 * @param  refID     ID of individuals with phenotype
 * @param  ID        ID of individuals with genotype
 * @param  nmrk      Number of SNPs
 * @param  nind      Number of individuals with genotype
 * @param  M         genotype content matrix
 * @param  type_M    Storage format of M, 1 for mrk x ind and 0 for ind x mrk
 * @param  Mref      genotype content matrix of individuals with phenotype
 * @param  type_Mr   Storage format of Mref, 1 for mrk x ind and 0 for ind x mrk
 */
void get_ref_M(int nref, char** refID, char** ID, int nmrk, int nind, double* M, int type_M, double* Mref, int type_Mr)
{
  /* storage format */
  int ld_M, inc_M, ld_Mref, inc_Mref;
  if (type_M) {
    ld_M = 1;
    inc_M = nind;
  }
  else {
    ld_M = nmrk;
    inc_M = 1;
  }

  if (type_Mr) {
    ld_Mref = 1;
    inc_Mref = nref;
  }
  else {
    ld_Mref = nmrk;
    inc_Mref = 1;
  }

  int ind_pos;
  for (size_t i = 0; i < nref; i++) {
    ind_pos = in_char(refID[i], ID, nind);
    if (ind_pos >= 0) {
      copy_d(M + ind_pos * ld_M, inc_M, Mref + i * ld_Mref, inc_Mref, nmrk, 1);
      // if (trans) {
      //   // storage as nmrk x nind
      //   copy_d(M + ind_pos * nmrk, 1, Mref + i * 1, nref, nmrk, 1);
      // }
      // else {
      //   // storage as nind x nmrk
      //   copy_d(M + ind_pos * nmrk, 1, Mref + i * nmrk, 1, nmrk, 1);
      // }
    }
  }
}

/**
 * @brief Centralization of gene content matrix (012)
 * @param M        genotype content matrix
 * @param nind     number of individuals
 * @param nmrk     number of markers
 * @param type     type of genotype storage, 1 for nind x nmrk and 0 for nmrk x nind
 */
void center(double* M, int nind, int nmrk, int type)
{
  // Matrix storage format
  int inrow, incol;
  if (type) {
    inrow = 1;
    incol = nmrk;
  }
  else {
    inrow = nmrk;
    incol = 1;
  }

  double mean;
  for (size_t snpi = 0; snpi < nmrk; snpi++) {
    mean = cblas_dasum(nind, M + snpi * inrow, incol) / (double)nind;
    /* Centralization of gene content */
    cblas_daxpy(nind, -1.0, &mean, 0, M + snpi * inrow, incol);
  }
}

/**
 * @brief parse plink fam file
 * @param  gene   Structure storing genotype information
 */
void read_fam(char* prefix, Gmat* M, char*** iid)
{
  // Check whether the file is readable
  char* fam_file_name = check_file(prefix, ".fam");

  // read file
  FILE* fp = fopen(fam_file_name, "r");

  // number of individuals
  M->nind = countLines(fam_file_name);

  // Prepare container
  char** fids = (char**)calloc(MAX_POPS, sizeof(char*));
  int* ninds = (int*)calloc(MAX_POPS, sizeof(int));
  M->fam = (Fam*)calloc(M->nind, sizeof(Fam));
  *iid = (char**)calloc(M->nind, sizeof(char*));

  // parse file
  int i = 0;
  char line[MAX_LINE_LEN];
  char* buff;
  while (fgets(line, MAX_LINE_LEN, fp) != NULL) {
    // family ID
    strcpy_l(&(M->fam[i].fid), strtok_r(line, " \t", &buff));
    strcpy_l(&(M->fam[i].iid), strtok_r(NULL, " \t", &buff));
    strcpy_l(&(M->fam[i].father), strtok_r(NULL, " \t", &buff));
    strcpy_l(&(M->fam[i].mother), strtok_r(NULL, " \t", &buff));
    M->fam[i].sex = atoi(strtok_r(NULL, " \t", &buff));
    M->fam[i].phe = atof(strtok_r(NULL, " \t", &buff));
    (*iid)[i] = M->fam[i].iid;

    /* save unique family IDs */
    if (in_char(M->fam[i].fid, fids, M->npop) < 0) {
      fids[M->npop] = M->fam[i].fid;
      M->npop++;
      if (M->npop > MAX_POPS) {
        printf("Error: Too many populations!\n");
        exit(ERROR_PARA);
      }
    }

    /* Record the number of individuals of each breed */
    ninds[M->npop - 1]++;

    i++;
  }

  /* storage information */
  M->fids = (char**)calloc(M->npop, sizeof(char*));
  M->ninds = (int*)calloc(M->npop, sizeof(int));
  for (size_t i = 0; i < M->npop; i++) {
    M->fids[i] = fids[i];
    M->ninds[i] = ninds[i];
  }

  free(fids);
  free(ninds);
  fclose(fp);
}

/**
 * @brief parse plink bim file
 * @param  gene   Structure storing genotype information
 */
void read_bim(char* prefix, char* mapf, Gmat* M)
{
  // Check whether the file is readable
  char* bim_file_name;
  if (strlen(prefix) != 0) {
    bim_file_name = check_file(prefix, ".bim");
  }
  else {
    bim_file_name = check_file(prefix, "");
  }

  // read file
  FILE* fp = fopen(bim_file_name, "r");

  // number of individuals
  M->nmrk = countLines(bim_file_name);

  // Prepare container
  M->bim = (Bim*)calloc(M->nmrk, sizeof(Bim));

  // parse file
  int i = 0;
  char line[MAX_LINE_LEN];
  char* buff;
  M->chr = (char**)calloc(M->nmrk, sizeof(char*));
  while (fgets(line, MAX_LINE_LEN, fp) != NULL) {
    // parse lines
    strcpy_l(&(M->bim[i].chr), strtok_r(line, " \t", &buff));
    strcpy_l(&(M->bim[i].id), strtok_r(NULL, " \t", &buff));
    M->bim[i].cm = atof(strtok_r(NULL, " \t", &buff));
    M->bim[i].pos = atoi(strtok_r(NULL, " \t", &buff));

    /* allele */
    read_allele(M->bim[i].A1, strtok_r(NULL, " \t", &buff));
    read_allele(M->bim[i].A2, strtok_r(NULL, " \t", &buff));
    // strcpy(M->bim[i].A2, strtok_r(NULL, " \t", &buff));  // only snp

    M->chr[i] = M->bim[i].chr;
    i++;
  }

  /* Number of SNPs in each chromosome */
  M->chr_pos = (int*)calloc(M->nmrk, sizeof(int));
  M->chr = unique_char(M->chr, M->nmrk, 1, &M->nchr, &M->nsnps, M->chr_pos);

  fclose(fp);
}

/**
 * @brief parse plink bed file, 0/1/2/3 for Homozygous of A2/Heterozygous/Homozygous of A1/missing
 * @param  gene   Structure storing genotype information
 */
void read_bed(char* prefix, Gmat* M)
{
  // Check whether the file is readable
  char* bed_file_name = check_file(prefix, ".bed");

  // read file
  FILE* fp = fopen(bed_file_name, "rb");
  unsigned char magic[3];
  fread(magic, sizeof(unsigned char), 3, fp);

  // Check the validity of the bed file.
  if (magic[0] != 108 || magic[1] != 27) {
    printf("Invalid magic number in BED file\n");
    exit(ERROR_PARA);
  }
  if (magic[2] != 1) {
    if (magic[2] == 0) {
      printf("bed file is in individual-major format\n");
      exit(ERROR_PARA);
    }
    else {
      printf("bed file-format specifier is not valid\n");
      exit(ERROR_PARA);
    }
  }

  int snp_len = (M->nind + 3) / 4;  // number of bytes used to store each SNP
  int val;
  unsigned char buffer[snp_len];
  unsigned char snp_index[] = {2, MISS_GENE_CODE, 1, 0};  // hom1, miss, het, hom2
  M->M = (double*)calloc(M->nind * M->nmrk, sizeof(double));
  M->mrk_miss = (int*)calloc(M->npop * M->nmrk, sizeof(int));

  /* store format */
  gmat_format(M, 1);

  /* get allele content */
  for (size_t i = 0; i < M->nmrk; i++) {
    fread(buffer, sizeof(unsigned char), snp_len, fp);
    for (size_t p = 0, indij = 0, k = 0; p < M->npop; p++) {
      for (size_t j = 0; j < M->ninds[p]; j++, indij++) {
        val = (buffer[k] >> (2 * (indij % 4))) & 0x03;

        /* count the number of missing SNPs */
        if (val == 1) {
          M->mrk_miss[i + p * M->nmrk]++;
          // M->bim[i].nmiss++;
          M->fam[indij].nmiss++;
        }

        /* Number of reference bases */
        M->M[i * M->ld + indij * M->inc] = snp_index[val];

        if (indij % 4 == 3) k++;
      }
    }
  }

  fclose(fp);
}

void plink_parse(char* rawf, char* prefix, struct Genotype* gene)
{
  if (strlen(rawf) != 0) {
    plink_raw_parse(rawf, &gene->Ma, &gene->iid);
  }
  else if (strlen(prefix) != 0) {
    plink_binary_parse(prefix, gene);
  }
  else {
    printf("Error: No genotype file provided!\n");
    exit(ERROR_PARA);
  }
}

/**
 * @brief Parse the Plink raw format file
 * @param  rawf             path of the raw file
 * @param  gene             Structure for storing genotype parameters
 */
void plink_raw_parse(char* rawf, Gmat* M, char*** iid)
{
  /* number of individuals and markers */
  M->nind = countLines(rawf) - 1;
  M->nmrk = get_column_count(rawf) - 6;
  M->M = (double*)calloc(M->nmrk * M->nind, sizeof(double));
  char* raw = read_whole(rawf);

  char *buff, *whole_line, *tmp, *header;  // *pat, *mat, *sex, *phe

  /* get snp ID */
  M->bim = (Bim*)calloc(M->nmrk, sizeof(Bim));
  buff = strtok_r(raw, "\n", &whole_line);
  for (size_t i = 0; i < M->nmrk + 6; i++) {
    tmp = strtok_r(buff, " ", &header);
    if (i > 5) {
      M->bim[i].id = (char*)calloc(strlen(tmp), sizeof(char));
      strcpy(M->bim[i].id, tmp);
    }
    if (i == 0) buff = NULL;
  }

  /* Prepare container */
  char** fids = (char**)calloc(MAX_POPS, sizeof(char*));
  int* ninds = (int*)calloc(MAX_POPS, sizeof(int));
  M->fam = (Fam*)calloc(M->nind, sizeof(Fam));
  *iid = (char**)calloc(M->nind, sizeof(char*));

  /* storage format */
  gmat_format(M, 0);

  /* Parse file strings line by line */
  M->npop = 0;
  double gene_count = 0;
  for (size_t i = 0; i < M->nind; i++) {
    /* first 6 columns */
    strcpy_l(&(M->fam[i].fid), strtok_r(NULL, " \t", &whole_line));
    strcpy_l(&(M->fam[i].iid), strtok_r(NULL, " \t", &whole_line));
    strcpy_l(&(M->fam[i].father), strtok_r(NULL, " \t", &whole_line));
    strcpy_l(&(M->fam[i].mother), strtok_r(NULL, " \t", &whole_line));
    M->fam[i].sex = atoi(strtok_r(NULL, " \t", &whole_line));
    M->fam[i].phe = atof(strtok_r(NULL, " \t", &whole_line));
    (*iid)[i] = M->fam[i].iid;

    /* save unique family IDs */
    if (in_char(M->fam[i].fid, fids, M->npop) < 0) {
      // fids[M->npop] = (char*)calloc(strlen(M->fam[i].fid) + 1, sizeof(char));
      // strcpy(fids[M->npop], M->fam[i].fid);
      fids[M->npop] = M->fam[i].fid;
      M->npop++;
      if (M->npop > MAX_POPS) {
        printf("Error: Too many populations!\n");
        exit(ERROR_PARA);
      }
    }

    /* Record the number of individuals of each breed */
    ninds[M->npop - 1]++;

    /* Read the storage SNP content matrix */
    buff = strtok_r(NULL, "\n", &whole_line); /* allele count of individual i */
    for (size_t j = 0; j < M->nmrk; j++) {
      // tmp = strtok(buff, " ");

      if (*buff == 78) {
        /* missing, NA, ASCII codes for "N" is 78 */
        gene_count = MISS_GENE_CODE;
        /* count the number of missing */
        M->mrk_miss[j + (M->npop - 1) * M->nmrk]++;
        // M->bim[j].nmiss++;
        M->fam[i].nmiss++;
        buff += 3;
      }
      else if ((*buff > 47) & (*buff < 51)) {
        /* The ASCII codes for "0", "1", and "2" are 48, 49, and 50 */
        gene_count = *buff - 48;
        buff += 2;
      }
      else {
        printf("Error: Unknown genotype code in plink raw file: %c!\n", *buff);
        exit(ERROR_PARA);
      }
      M->M[i * M->ld + j * M->inc] = gene_count;
    }
  }

  /* store information */
  M->fids = (char**)calloc(M->npop, sizeof(char*));
  M->ninds = (int*)calloc(M->npop, sizeof(int));
  for (size_t i = 0; i < M->npop; i++) {
    M->fids[i] = fids[i];
    M->ninds[i] = ninds[i];
  }

  free(fids);
  free(ninds);
  free(raw);
}

/**
 * @brief Parse the Plink raw format file
 * @param  rawf             path of the binary file prefix
 * @param  gene             Structure for storing genotype parameters
 * @param  trans            set to 1 for transposing the matrix, that is, it is stored in the form of nmrk x nind;
 * set to 0 for storing in the form of nind x nmrk.
 */
void plink_binary_parse(char* prefix, struct Genotype* gene)
{
  gene->prefix = prefix;
  read_fam(gene->prefix, &gene->Ma, &gene->iid);
  read_bim(gene->prefix, gene->mapf, &gene->Ma);
  read_bed(gene->prefix, &gene->Ma);
}

// /**
//  * @brief
//  * @param  gene             Structure for storing genotype parameters
//  */
// void genotype_qc_old(Gmat* org, Gmat* new, double maf, double mind, double geno)
// {
//   /* init */
//   copy_gmat(new, org);
//   init_gmat(new);

//   /* filter SNP */
//   int maf_del = 0, geno_del = 0, loc = 0, keep = 0;
//   double miss_rate, mafi;
//   for (size_t i = 0; i < org->nchr; i++) {
//     for (size_t j = 0; j < org->nsnps[i]; j++) {
//       mafi = min_d(org->frq[loc], 1 - org->frq[loc]);
//       if (org->bim[loc].nmiss > 0) {
//         /* missing exist in SNP i */
//         miss_rate = (double)org->bim[loc].nmiss / (double)org->nind;
//         if (miss_rate > geno) {
//           remove_SNP(new, i, loc);
//           geno_del++;
//         }
//         else {
//           new->frq[keep] = org->frq[loc];
//           keep++;
//         }
//       }
//       else if (mafi <= maf) {
//         remove_SNP(new, i, loc);
//         maf_del++;
//       }
//       else {
//         new->frq[keep] = org->frq[loc];
//         keep++;
//       }
//       if (new->nsnps[i] == 0) new->nchr--; /* all SNPs in chr i has beed removed */
//       loc++;
//     }
//   }

//   /* filter sample */
//   loc = 0;
//   int mind_del = 0;
//   for (size_t i = 0; i < org->npop; i++) {
//     for (size_t j = 0; j < org->ninds[i]; j++) {
//       if (org->fam[loc].nmiss > 0) {
//         /* missing exist in individual i */
//         miss_rate = (double)org->fam[loc].nmiss / (double)org->nmrk;
//         if (miss_rate > mind) {
//           remove_individual(new, i, loc);
//           mind_del++;
//         }
//       }
//       else {
//         new->ninds[i]++;
//       }
//       if (new->ninds[i] == 0) new->npop--; /* all samples in pop i has beed removed */

//       loc++;
//     }
//   }

//   /* matrix storage format */
//   gmat_format(org, 1);
//   gmat_format(new, 1);

//   /* imputation */
//   double fill_value;
//   int num_miss;
//   for (size_t i = 0; i < org->nmrk; i++) {
//     num_miss = 0;
//     fill_value = (org->center) ? 0.0 : 2 * org->frq[i];                 /* mean */
//     if ((org->bim[i].nmiss == 0) | (new->snp_index[i] == 0)) continue;  // Not missing or removed
//     for (size_t j = 0; j < org->nind; j++) {
//       if ((org->fam[j].nmiss == 0) | (new->ind_index[j] == 0)) continue;
//       org->M[i * org->ld + j * org->inc] = fill_value;
//       num_miss++;
//       if (num_miss > org->bim[i].nmiss) break;  // All missing SNPs have been imputed
//     }
//   }

//   /* new gene content matrix */
//   new->M = (double*)calloc(new->nind* new->nmrk, sizeof(double));
//   subset(org->row, org->col, org->M, new->row_index, 0, new->col_index, 0, new->M);

//   /* report */
//   printf("%d variants and %d people pass filters and QC.\n", new->nmrk, new->nind);

//   /* Delete the original matrix to save memory */
//   free(org->M);
// }

/**
 * @brief
 * @param  gene             Structure for storing genotype parameters
 */
void genotype_qc(Gmat* org, Gmat* new, double maf, double mind, double geno)
{
  /* init */
  copy_gmat(new, org);
  init_gmat(new);

  /* filter SNP */
  int maf_del = 0, geno_del = 0;
  double miss_rate, mafi;
  for (size_t i = 0; i < org->npop; i++) {
    for (size_t j = 0, keep = 0, loc = 0; j < org->nchr; j++) {
      for (size_t k = 0; k < org->nsnps[j]; k++, loc++) {
        mafi = min_d(org->frq[loc + i * org->nmrk], 1 - org->frq[loc + i * org->nmrk]);
        if (org->mrk_miss[loc + i * org->nmrk] > 0) {
          /* missing exist in SNP i */
          miss_rate = (double)org->mrk_miss[loc + i * org->nmrk] / (double)org->ninds[i];
          if (miss_rate > geno) {
            remove_SNP(new, j, loc);
            geno_del++;
          }
          else {
            new->frq[keep + i * org->nmrk] = org->frq[loc + i * org->nmrk];
            keep++;
          }
        }
        else if (mafi <= maf) {
          remove_SNP(new, j, loc);
          maf_del++;
        }
        else {
          new->frq[keep + i * org->nmrk] = org->frq[loc + i * org->nmrk];
          keep++;
        }
        if (new->nsnps[j] == 0) new->nchr--; /* all SNPs in chr i has beed removed */
      }
    }
  }

  /* filter sample */
  int mind_del = 0;
  for (size_t i = 0, loc = 0; i < org->npop; i++) {
    for (size_t j = 0; j < org->ninds[i]; j++, loc++) {
      if (org->fam[loc].nmiss > 0) {
        /* missing exist in individual i */
        miss_rate = (double)org->fam[loc].nmiss / (double)org->nmrk;
        if (miss_rate > mind) {
          remove_individual(new, i, loc);
          mind_del++;
        }
      }
      else {
        new->ninds[i]++;
      }
      if (new->ninds[i] == 0) new->npop--; /* all samples in pop i has beed removed */
    }
  }

  /* matrix storage format */
  gmat_format(org, 1);
  gmat_format(new, 1);

  /* imputation */
  double fill_value;
  int num_miss;
  for (size_t i = 0; i < org->nmrk; i++) {
    for (size_t j = 0, indij = 0; j < org->npop; j++) {
      num_miss = 0;
      if ((org->mrk_miss[i + j * org->nmrk] == 0) | (new->snp_index[i] == 0)) continue;  // Not missing or removed
      fill_value = (org->center) ? 0.0 : 2 * org->frq[i + j * org->nmrk];                /* mean */
      for (size_t k = 0; k < org->ninds[j]; k++, indij++) {
        if ((org->fam[indij].nmiss == 0) | (new->ind_index[indij] == 0)) continue;
        org->M[i * org->ld + indij * org->inc] = fill_value;
        num_miss++;
        if (num_miss > org->mrk_miss[i + j * org->nmrk]) break;  // All missing SNPs have been imputed
      }
    }
  }

  /* new gene content matrix */
  new->M = (double*)calloc(new->nind* new->nmrk, sizeof(double));
  subset(org->row, org->col, org->M, new->row_index, 0, new->col_index, 0, new->M);

  /* report */
  printf("%d variants removed due to minor allele threshold (--maf).\n", maf_del);
  printf("%d variants removed due to missing genotype data (--geno).\n", geno_del);
  printf("%d samples removed due to missing genotype data (--mind).\n", mind_del);
  printf("%d variants and %d people pass filters and QC.\n", new->nmrk, new->nind);

  /* Delete the original matrix to save memory */
  free(org->M);
}

/**
 * @brief Sum gene frequency according to grouping
 * @param  nbin    Number of bins
 * @param  nsnps   Number of SNPs in each bin
 * @param  frq     Frequency of each snp
 * @param  sum2pq  Sum of 2pq
 */
void sum2pq_group(int npop, int nmrk, int nbin, double* nsnps, double* frq, double* sum2pq)
{
  for (size_t i = 0; i < npop; i++) {
    for (size_t j = 0, snpjk = 0; j < nbin; j++) {
      sum2pq[j + i * nbin] = 0.0;
      for (size_t k = 0; k < nsnps[j]; k++, snpjk++) {
        sum2pq[j + i * nbin] += 2.0 * frq[snpjk + i * nmrk] * (1.0 - frq[snpjk + i * nmrk]);
      }
    }
  }
}

/**
 * @brief update leading dimension and increment
 * @param mat   Gmat structure
 * @param aim   0 for subset individual's SNPs and 1 for subset SNPs of population
 */
void gmat_format(Gmat* mat, int aim)
{
  if (mat->format) {
    /* mrk x nind */
    if (aim) {
      mat->ld = mat->nind;
      mat->inc = 1;
    }
    else {
      mat->ld = 1;
      mat->inc = mat->nind;
    }
    mat->row = mat->nmrk;
    mat->col = mat->nind;
    mat->row_index = mat->snp_index;
    mat->col_index = mat->ind_index;
  }
  else {
    /* nind x nmrk */
    if (aim) {
      mat->ld = 1;
      mat->inc = mat->nmrk;
    }
    else {
      mat->ld = mat->nmrk;
      mat->inc = 1;
    }
    mat->row = mat->nind;
    mat->col = mat->nmrk;
    mat->row_index = mat->ind_index;
    mat->col_index = mat->snp_index;
  }
}

/**
 * @brief Synchronize some elements in two structures
 * @param to    Target structure
 * @param from  Source Structure
 */
void copy_gmat(Gmat* to, Gmat* from)
{
  to->npop = from->npop;
  to->nchr = from->nchr;
  to->nmrk = from->nmrk;
  to->nind = from->nind;
  to->format = from->format;
  to->center = from->center;

  to->nsnps = (int*)calloc(from->nchr, sizeof(int));
  copy_i(from->nsnps, 1, to->nsnps, 1, from->nchr, 1);

  to->ninds = (int*)calloc(from->npop, sizeof(int));
  copy_i(from->ninds, 1, to->ninds, 1, from->npop, 1);
}

/**
 * @brief Allocate memory for some variables in the matrix structure
 * @param M    Target structure
 */
void init_gmat(Gmat* M)
{
  M->snp_index = (int*)calloc(M->nmrk, sizeof(int));
  fill_i(1, M->snp_index, M->nmrk, 1);

  M->ind_index = (int*)calloc(M->nind, sizeof(int));
  fill_i(1, M->ind_index, M->nind, 1);

  M->frq = (double*)calloc(M->nmrk * M->npop, sizeof(double));
}

/**
 * @brief  remove the specified individual
 * @param M         gene content matrix
 * @param popi      The population index of deleted individual
 * @param location  Index of deleted individuals
 */
void remove_individual(Gmat* M, int popi, int location)
{
  M->nind--;
  M->ninds[popi]--;
  M->ind_index[location] = 0;
}

/**
 * @brief  remove the specified SNP
 * @param M         gene content matrix
 * @param popi      The chromosome index of deleted SNP
 * @param location  Index of deleted SNP
 */
void remove_SNP(Gmat* M, int chri, int location)
{
  M->nsnps[chri]--;
  M->nmrk--;
  M->snp_index[location] = 0;
}

void read_allele(char* to, char* from)
{
  if ((from[1] != 10) & (from[1] != 0)) {
    printf("Genotype information can only contain SNP!\n");
    exit(ERROR_INPUT);
  }

  to[0] = from[0];
}

void fid_order(Gmat *M)
{
  int fid_rank;
  int last_fid_rank = 0;
  for (size_t i = 0; i < M->nind; i++)
  {
    fid_rank = in_char(M->fam[i].fid, M->fids, M->npop);

    if (fid_rank < 0) {
      printf("Error: unknown error occurred while checking the sorting of variety identifiers.\n");
      exit(ERROR_PARA);
    } else if (fid_rank < last_fid_rank) {
      printf("Error: Individuals in genotype files need to be sorted according to breed identifiers.\n");
      exit(ERROR_PARA);
    }
    last_fid_rank = fid_rank;
  }
}

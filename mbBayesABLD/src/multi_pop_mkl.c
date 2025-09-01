/**
 * @file multi_pop_mkl.c
 * @brief mbBayesAS -- Multi-breed Bayesian Models
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

#include <argp.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <time.h>
#include <unistd.h>

#include "common_fun.h"
#include "multi_pop_fun.h"

/* Program introduction */
static char doc[] =
    "mbBayesAS -- BayesA model modeling heterogeneous (co)variances from "
    "SNP group specified by the user";
const char* argp_program_version = "mbBayesAS 1.1.1";
const char* argp_program_bug_address = "<liwn@cau.edu.cn>";

/* parse function */
static error_t parse_opt(int key, char* arg, struct argp_state* state)
{
  struct arguments* arguments = state->input;
  struct Model* model = arguments->model;
  struct Genotype* gene = arguments->gene;
  struct Pheno* phe = arguments->phe;
  Bin* bin = arguments->bin;

  switch (key) {
    case 'g':
      arguments->raw = arg;
      arguments->require_para--;
      break;
    case 'l':
      arguments->prefix = arg;
      arguments->require_para--;
      break;
    case 'p':
      phe->filename = arg;
      arguments->require_para--;
      break;
    case 'y':
      phe->phePos = arg;
      arguments->require_para--;
      break;
    case 'm':
      phe->miss = atof(arg);
      break;
    case 'N':
      bin->nsnp = atoi(arg);
      break;
    case 'b':
      arguments->binf = arg;
      break;
    case 'e':
      bin->type = arg;
      break;
    case 'f':
      gene->mapf = arg;
      break;
    case 'j':
      arguments->binPos = atoi(arg);
      break;
    case 'a':
      phe->idPos = atoi(arg);
      break;
    case 'd':
      model->df_va = atof(arg);
      break;
    case 'D':
      model->df_ve = atof(arg);
      break;
    case 'G':
      arguments->A = atoi(arg);
      break;
    case 'E':
      arguments->effOutf = arg;
      break;
    case 'A':
      arguments->varOutf = arg;
      break;
    case 'S':
      arguments->iterAll = atoi(arg);
      break;
    case 'B':
      arguments->burnin = atoi(arg);
      break;
    case 't':
      arguments->thin = atoi(arg);
      break;
    case 's':
      arguments->seed = atol(arg);
      break;
    case 'r':
      arguments->report_sep = atoi(arg);
      break;
    case 'v':
      arguments->varf = arg;
      break;
    case 'V':
      arguments->gebvOutf = arg;
      break;
    case 'c':
      arguments->chain_var = atoi(arg);
      break;
    case 'M':
      arguments->mcmcOut = arg;
      break;
    case 'L':
      arguments->logf = arg;
      break;
    case 'i':
      arguments->prioSa = arg;
      break;
    case 'I':
      arguments->prioSe = arg;
      break;
    case 'F':
      phe->fixPos = arg;
      break;
    case 'R':
      phe->ranPos = arg;
      break;
    case 'J':
      arguments->maf = atof(arg);
      break;
    case 'K':
      arguments->geno = atof(arg);
      break;
    case 'P':
      arguments->mind = atof(arg);
      break;
    case 'o':
      arguments->outdir = arg;
      break;
    case 'C':
      gene->Ma.center = 0;
      break;
    case 'T':
      arguments->constrain = 1;
      break;
    case 'n':
      arguments->thread = atoi(arg);
      break;
    case 'H':
      phe->header = 1;
      break;
    case ARGP_KEY_END:
    {
      // After processing all parameters, the processing is carried out
      printf("\n");
      if (arguments->require_para > 0) {
        printf("Missing required parameters.\n\n");
        argp_state_help(state, state->err_stream, ARGP_HELP_STD_HELP);
      }
    } break;

    default:
      return ARGP_ERR_UNKNOWN;
  }

  return 0;
}

/* Parameter description */
static struct argp_option options[] = {
    {0, 0, 0, 0, "Required parameters:", 1},
    {"bfile", 'l', "FILE", 0, "plink binary file prefix, e.g. /path/plink"},
    {"phef", 'p', "FILE", 0, "phenotype file"},
    {0, 0, 0, 0, "Optional parameters:", 2},
    {"raw", 'g', "FILE", 0, "plink raw file, e.g. /path/plink.raw"},
    {"thread", 'T', "INT", 0, "number of threads for the OpenMP run-time library [1]"},
    {"phe", 'y', "CHR", 0, "phenotypes position, e.g. '4 5' [last npop column(s)]"},
    {"fix", 'F', "CHR", 0, "fixed effect position, e.g. '2 3' [NULL]"},
    {"ran", 'R', "CHR", 0, "environmental random effect position [NULL]"},
    {"id", 'a', "INT", 0, "id position at phenotype file [1]"},
    {"varf", 'v', "FILE", 0, "SNP variance prior [varp * 0.5]"},
    {"region", 'e', "CHR", 0, "Method for dividing snp into different groups, LD/fix [LD]"},
    {"binf", 'b', "FILE", 0, "file contains number of SNPs per interval [NULL]"},
    {"bincol", 'j', "INT", 0, "position of nSNPs in each bins at binf [1]"},
    {"mapf", 'f', "FILE", 0, "plink .map (.bim) use to define bins [NULL]"},
    {"nsnp", 'N', "FILE", 0, "Number of SNPs per interval [100]"},
    {"effOut", 'E', "FILE", 0, "posteriori mean of location effect [theta.txt]"},
    {"varOut", 'A', "FILE", 0, "posteriori mean of dispersion matrix [var.txt]"},
    {"gebvOut", 'V', "FILE", 0, "breeding value output file [gebv.txt]"},
    {"varC", 'c', "INT", 0, "output first c predictions in MCMC chain [NULL]"},
    {"mcmcOut", 'M', "FILE", 0, "output predictions in MCMC chain [MCMC_var.txt]"},
    {"miss", 'm', "DOUBLE", 0, "missing value in phenotype [-99]"},
    {"dfvara", 'd', "DOUBLE", 0, "degrees of freedom Va [4]"},
    {"dfvare", 'D', "DOUBLE", 0, "degrees of freedom Ve [4]"},
    {"diagG", 'G', "INT", 0, "a integer large enough [100000]"},
    {"iter", 'S', "INT", 0, "length of MCMC chain [100]"},
    {"maf", 'J', "INT", 0, "SNPs with maf less than this value will be excluded [-0.01]"},
    {"geno", 'K', "INT", 0, "Allow the missing proportion of the SNP [0.99]"},
    {"mind", 'P', "INT", 0, "Allow the missing genotype of the individual [0.99]"},
    {"burnin", 'B', "INT", 0, "burnin steps [10]"},
    {"thin", 't', "INT", 0, "thinning rate [10]"},
    {"seed", 's', "LONG", 0, "random number seed [NULL]"},
    {"outdir", 'o', "CHR", 0, "output directory [.]"},
    {"report", 'r', "INT", 0, "MCMC sampling report interval [100]"},
    {"logf", 'L', "FILE", 0, "redirect standard output to this file"},
    {"prioSa", 'i', "CHR", 0, "distribution of scale of additive effect prior"},  // [NULL/wishart/gamma]
    {"prioSe", 'I', "CHR", 0, "distribution of scale of residual effect prior"},  // [NULL/wishart]
    {0, 0, 0, 0, "Flag parameters:", 3},
    {"nocov", 'n', 0, 0, "constrain the covariance between the residuals to zero"},
    {"dontcenter", 'C', 0, 0, "centering gene content matrix"},
    {"header", 'H', 0, 0, "phenotype and bins files has header, [automatic detection]"},
    {0}};

/* argp parser, parameters of argp_parse */
static struct argp argps = {options, parse_opt, 0, doc};

int main(int argc, char* argv[])
{
  /* time record */
  time_t t;
  char time_buf[1024];
  time(&t);
  ctime_r(&t, time_buf);

  /* ----------- Structures that store variables of different types ----------- */
  Bin bin = {.min = 100,
            .nsnp = 100,
            .win = 50,
            .type = "fix",
            .thresh = 0.2};

  struct Output out;
  struct Pheno phe = {
      .header = -1,
      .phePos = "",  // Column in phenotype file where phenotype data is located
      .miss = -99,   // Missing phenotype identifier
      .idPos = 1,    // Column in phenotype file where ID is located
      .ranPos = "",  // Column in phenotype file for individual random effects; fitting non-genetic random effects is currently not supported
      .type = 0,     // Phenotype file storage mode: 1 for matrix storage, 0 for vector storage
      .fixPos = ""   // Columns in phenotype file for individual fixed effects, e.g., 2,3
  };
  struct Genotype geno = {.Ma.format = 0,
                          .Ma.center = 1,  // Center the gene content matrix
                          .mapf = "",      // Map file of marker loci; required when partitioning binf file is not provided
                          .Mq.format = 0};

  struct Model model = {.df_va = 4.0,
                        .df_ve = 4.0,
                        .h2 = 0.5,
                        .rg = 0.0,
                        .nfix = 0,
                        .nRan = 0,
                        .geno = &geno};

  struct arguments para = {.raw = "",     // Genotype file in PLINK raw format
                           .prefix = "",  // Prefix of PLINK binary genotype files (fam/bim/bed)
                           .thread = 1,   // Number of threads for OpenMP library functions
                           .A = 100000,   // A sufficiently large value; used when prioSa is gamma
                           .nsnp_bin = 100,  // When partitioning information is not provided, divide into intervals with approximately nsnp_bin markers each
                           .iterAll = 100,     // Total number of iterations
                           .burnin = 10,       // Number of burn-in cycles
                           .thin = 10,         // Sampling interval after burn-in
                           .seed = 0,          // Random seed
                           .report_sep = 100,  // Interval for reporting time and certain variable values
                           .chain_var = 4,     // Output the first chain_var variables to observe changes during iterations
                           .effOutf = "theta.txt",     // Output file for posterior means of effects (all effects)
                           .varOutf = "var.txt",       // Output file for posterior means of variance components
                           .gebvOutf = "EBV",          // Output file for genomic estimated breeding values
                           .mcmcOut = "MCMC_var.txt",  // Variable values of the first chain_var variables
                           .logf = "",                 // Log file
                           .binf = "",                 // File indicating the number of SNPs per interval
                           .binPos = 1,                // Column in binf file indicating the number of SNPs within intervals
                           .outdir = ".",              // Output file path
                           .varf = "",    // Prior variance components; each line represents a component; the last line is residual variance
                           .prioSa = "",  // Distribution of prior variance for additive effects; can be wishart/gamma
                           .prioSe = "",  // Distribution of prior variance for residual effects; can be wishart
                           .require_para = 2,  // Number of required parameters
                           .constrain = 0,     // Constraint on residual covariance
                           .maf = -0.01,       // Allowed minimum allele frequency
                           .geno = 0.99,       // Allowed missing rate for a SNP
                           .mind = 0.99,       // Allowed missing rate of markers for an individual
                           .gene = &geno,      // Command-line arguments for genotype-related parameters
                           .bin = &bin,        // Command-line arguments for genotype-related parameters
                           .phe = &phe,        // Command-line arguments for phenotype-related parameters
                           .model = &model};   // Command-line arguments for model-related parameters

  /* ----------- Parse command line arguments ----------- */
  argp_parse(&argps, argc, argv, ARGP_IN_ORDER, 0, &para);

  /* number of threads for the OpenMP run-time library */
  if (para.thread > 0) {
    mkl_set_num_threads(para.thread);
  }

  if (para.seed == 0) {
    para.seed = time(NULL);
  }
  /* Initializing the random number generator */
  srand(para.seed);

  // Check whether the output folder exists
  struct stat st = {0};
  if (stat(para.outdir, &st) == -1) {
    mkdir(para.outdir, 0777);
  }

  /* redirect standard output to log file */
  char filename[MAX_LINE_LEN];
  if (strlen(para.logf) != 0) {
    sprintf(filename, "%s/%s", para.outdir, para.logf);
    freopen(filename, "w", stdout);
  }

  /* time report */
  printf("Program mbBayesABLD\n");
  printf("Run started at %s", time_buf);
  fflush(stdout);

  /* Random number generator */
  VSLStreamStatePtr stream;
  /***** Initialize *****/
  vslNewStream(&stream, BRNG, para.seed);

  /* ----------- Get the information needed from the file ----------- */
  /* ### Individual ID with phenotype ### */
  phe.rows = countLines(phe.filename);
  phe.cols = get_column_count(phe.filename);
  phe.ID = (char**)calloc(phe.rows, sizeof(char*));
  for (int i = 0; i < phe.rows; i++) phe.ID[i] = (char*)calloc(MAX_ID_LEN, sizeof(char));
  phe.nind = read_dataframe_char(phe.filename, phe.rows, phe.cols, phe.ID, 1, &phe.idPos, phe.header);
  phe.header = phe.rows - phe.nind;

  /* Parse plink file into matrix(array) format */
  printf("Loading genotypes file...\n");
  fflush(stdout);
  plink_parse(para.raw, para.prefix, &geno);
  // output_d2i_dataframe(geno.Ma.M, geno.Ma.nind, geno.Ma.nmrk, "genoMa.txt");

  /* Check if individuals in genotype files are sorted according to breed identifiers */
  fid_order(&geno.Ma);

  /* read phenotypes */
  model.npop = geno.Ma.npop;
  // model.npops = (int*)calloc(model.npop, sizeof(int));
  // copy_d(geno.Ma.ninds, 1, model.npops, 1, model.npop, 1);
  load_phenotype(&phe, &model, phe.phePos, phe.miss);

  /* gene frequency */
  geno.Ma.frq = (double*)calloc(geno.Ma.npop * geno.Ma.nmrk, sizeof(double));
  gene_frq(&geno.Ma);

  /* filtering */
  printf("Filtering SNPs or samples...\n\n");
  fflush(stdout);
  genotype_qc(&geno.Ma, &geno.Mq, para.maf, para.mind, para.geno);
  Gmat M = geno.Mq;

  /* variance prior distribution degrees of freedom */
  model.df_va += M.npop;
  model.df_ve += M.npop;

  /* square of the number of breeds */
  int n_pp = pow(M.npop, 2);

  /* phenotype missing pattern */
  model.obs = (int*)calloc(n_pp, sizeof(int));
  obs_type(M.npop, model.npops, model.yc, phe.miss, 0, model.obs);

  /* bins definition */
  int row_binf = 0;
  if (strlen(para.binf) != 0) {
    /* load from file */
    row_binf = countLines(para.binf);
    model.nsnps = (double*)calloc(row_binf, sizeof(double));
    model.nbin = read_dataframe_double(para.binf, row_binf, 0, model.nsnps, 1, &para.binPos, -1);
  }
  else {
    make_bins(&geno, &geno.Mq, &bin);
    model.nsnps = bin.nsnps;
    model.nbin = bin.nbin;
  }

  /* Check that the sum of all bins markers is equal to the number of markers in the genotype file */
  int nsnp_binf = sumd(model.nsnps, model.nbin, 1);
  if (M.nmrk != nsnp_binf) {
    printf("Error: The number of snps in binf is not equal to the number of snps in the genotype file.\n");
    printf("%d ≠ %d\n", M.nmrk, nsnp_binf);
    exit(ERROR_PARA);
  }

  /* Gene content matrix of validation population and reference population */
  model.nind = phe.nindPhe;
  int npind = model.nind * M.npop;
  model.Mr = (double*)calloc(M.nmrk * model.nind, sizeof(double));
  get_ref_M(model.nind, phe.pheID, geno.iid, M.nmrk, M.nind, M.M, 0, model.Mr, 1);

  /* Build incidence matrix X */
  if (strlen(phe.fixPos) != 0) {
    phe.fixCols = string_to_array(phe.fixPos, &model.nfix, " ");
    make_incidence(&phe, &model, &out);
  }
  else {
    /* If no fixed effects are specified, only the intercept term fitted */
    model.nfix = model.fixL = 1;
    model.fixLa = M.npop;
    model.X = (double*)calloc(npind * model.fixLa, sizeof(double));
    for (int i = 0; i < M.npop; i++) {
      fill_d(1.0, model.X + (i * model.nind) * model.fixLa + i, model.nind, model.fixLa);
    }
  }
  // output_d2i_dataframe(model.X, npind, model.fixLa, "Xc.txt"); // debug

  /* dot multiplication of column vectors in fixed incidence matrix */
  model.X1pX1 = (double*)calloc(model.fixL * model.fixL, sizeof(double));
  cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, model.fixL, model.fixL, model.nind, 1.0, model.X, model.fixLa,
              model.X, model.fixLa, 1.0, model.X1pX1, model.fixL);

  /* non-genetic random effects */
  model.ran_level = 0;
  if (strlen(phe.ranPos) != 0) {
    printf("Error: Fitting non-genetic random effects is not supported now.\n");
    exit(UNDER_DEVELOP);
  }

  /* Variance component prior */
  model.Va = (double*)calloc(n_pp * model.nbin + 1, sizeof(double)); /* The last vector stores R0 */
  model.Re = (double*)calloc(n_pp, sizeof(double));
  if (strlen(para.varf) != 0) {
    /* Check the number of variance components */
    int nprior = countLines(para.varf);
    if (nprior != (model.nRan + model.nbin + 1)) {
      printf("Error: Lines in the varf file provided is not equal to the number of random effects.\n");
      printf("Please note that there should be no header lines in the varf.\n");
      exit(ERROR_PARA);
    }
    read_dataframe_double(para.varf, 0, 0, model.Va, 0, 0, 0);  // need modification
    /* The last matrix stores R0 */
    copy_d(model.Va + model.nbin * n_pp, 1, model.Re, 1, n_pp, 1);
  }
  else {
    /* use for setting genetic variance prior */
    bin.sum2pq = (double*)calloc(model.nbin * model.npop, sizeof(double));
    sum2pq_group(model.npop, M.nmrk, model.nbin, model.nsnps, M.frq, bin.sum2pq);

    /* get initial variance components according to phenotypic variance */
    double* pvar = (double*)calloc(M.npop, sizeof(double));
    var(model.yc, model.nind, M.npop, phe.miss, 0, pvar);
    priorMat(M.npop, pvar, model.h2, model.rg, model.Va, model.nbin, bin.sum2pq);
    // copy_d(model.Va, 1, model.Va + n_pp, 1, n_pp, model.nbin - 1); /* same in all bins */
    priorMat(M.npop, pvar, 1.0 - model.h2, 0.0, model.Re, 1, NULL);
  }

  // /* genetic variance prior */
  // double sum2pq = sumd(bin.sum2pq, model.nbin, 1);
  // for (size_t i = 0; i < model.nbin; i++) {
  //   shrink(model.Va + i * n_pp, n_pp, 1, bin.sum2pq[i] / sum2pq, model.Va + i * n_pp, 1);
  // }

  /* inversion of random variance matrix (precision matrix) */
  // output_ddataframe(model.Va, model.nbin, n_pp, "var_prior.txt", 1.0);  // debug
  model.VaInv = (double*)calloc(n_pp * model.nbin, sizeof(double));
  for (size_t i = 0; i < model.nbin; i++) {
    solve(M.npop, model.Va + i * n_pp, model.VaInv + i * n_pp);
  }
  model.Re_inv = (double*)calloc(n_pp, sizeof(double));
  solve(M.npop, model.Re, model.Re_inv);

  /* scale parameters of additive effect variance */
  model.priVaS = (double*)calloc(n_pp * model.nbin, sizeof(double));
  for (size_t i = 0; i < model.nbin; i++) {
    shrink(model.Va + i * n_pp, n_pp, 1, model.df_ve - M.npop - 1, model.priVaS + i * n_pp, 1);
  }

  /* scale parameters of residual effect variance */
  model.priReS = (double*)calloc(n_pp, sizeof(double));
  shrink(model.Re, n_pp, 1, model.df_va - M.npop - 1, model.priReS, 1);

  /* Dot multiplication of column vectors of fix effect incidence matrix */
  model.effsqs = (double*)calloc(model.fixLa + M.nmrk, sizeof(double));
  vec_square_sum(model.X, npind, model.fixLa, model.effsqs, 1);
  vec_square_sum(model.Mr, M.nmrk, model.nind, model.effsqs + model.fixLa, 0);

  /* Whether to sample the scale parameter of the random effects variance distribution */
  char* prioSaPdf[3] = {"", "wishart"};                /* Currently supported distribution for genetic variance */
  char* prioSePdf[3] = {"", "wishart", "gamma"};       /* Currently supported distribution for residual variance */
  int SamVaSigma = in_char(para.prioSa, prioSaPdf, 3); /* estimate scale parameter of genetic variance */
  int SamReSigma = in_char(para.prioSe, prioSePdf, 3); /* estimate scale parameter of residual variance */

  /* Report parameter information */
  printf("\nInput genotype file  : %s\n", geno.filename == NULL ? geno.prefix : geno.filename);
  printf("Input phenotype file : %s\n", phe.filename);
  printf("output directory     : %s\n", para.outdir);
  printf("No. of individuals   : %d\n", model.nind);
  printf("No. of loci          : %d\n", M.nmrk);
  printf("No. of regions       : %d\n", model.nbin);
  printf("Missing phenotype    : %lf\n", phe.miss);
  printf_array("Prior Vara           : ", model.Va, n_pp, 4);
  printf_array("Prior Vare           : ", model.Re, n_pp, 4);
  printf("No. of cycles        : %d\n", para.iterAll);
  printf("Burnin               : %d\n", para.burnin);
  printf("Thinning rate        : %d\n", para.thin);
  printf("Seed                 : %ld\n", para.seed);

  /* Variable declaration and definition */
  int fixLa = model.fixLa;                           // check const MKL_INT
  int ntheta = fixLa + model.nRan + M.nmrk * M.npop; /* Sum of all effect levels */
  model.Se = (double*)calloc(n_pp, sizeof(double));
  model.Sa = (double*)calloc(n_pp, sizeof(double));
  model.Lhs = (double*)calloc(pow(fixLa, 2), sizeof(double));
  model.Rhs = (double*)calloc(fixLa, sizeof(double));
  model.theta = (double*)calloc(ntheta, sizeof(double)); /* effects */
  model.tmp_p = (double*)calloc(M.npop, sizeof(double)); /* Auxiliary variable */
  int index;
  int snp_start = 0; /* Starting SNP tag index of each interval */
  int snp_end = 0;   /* End SNP tag index for each interval */
  int elasp = 0;     /* Time spent in each para.report round */
  int eta = 0;
  double alpha_gik = 0.0;
  double beta_gik = 0.0;
  double* E_alphaij = (double*)calloc(M.npop, sizeof(double));
  double* ycAlphaj = (double*)calloc(M.npop, sizeof(double));
  double* V_alphaij = (double*)calloc(n_pp, sizeof(double));
  double* g = (double*)calloc(model.nbin * n_pp, sizeof(double));
  double* theta_mean = (double*)calloc(ntheta, sizeof(double));
  double* var_mean = (double*)calloc((model.nbin + 1) * n_pp, sizeof(double));
  double* precAlphaInv = (double*)calloc(model.nbin * n_pp, sizeof(double));
  double* ReInvYcorr = (double*)calloc(model.nbin * M.npop, sizeof(double));
  double* wishartIS = (double*)calloc(n_pp, sizeof(double));
  double* wishartScale = (double*)calloc(n_pp, sizeof(double));
  double* cover_judge = (double*)calloc((para.iterAll / para.report_sep) * para.chain_var, sizeof(double));
  double* alpha_last = (double*)calloc(M.npop, sizeof(double));
  int judge_i = 0, errcode = 0;
  int varnum = model.nbin * n_pp;

  /* Missing phenotype sampling parameters */
  struct RseSample resmiss;
  resmiss.nmiss = M.npop - 1;
  resmiss.sample = (double*)calloc(model.nind * resmiss.nmiss, sizeof(double));
  resmiss.tmp = (double*)calloc(resmiss.nmiss, sizeof(double));
  resmiss.E_ym = (double*)calloc(resmiss.nmiss, sizeof(double));
  resmiss.V_ym = (double*)calloc(resmiss.nmiss * resmiss.nmiss, sizeof(double));
  resmiss.R0_mo = (double*)calloc(resmiss.nmiss, sizeof(double));
  resmiss.R0_mm = (double*)calloc(resmiss.nmiss * resmiss.nmiss, sizeof(double));
  resmiss.obs = (int*)calloc(n_pp, sizeof(int));
  copy_i(model.obs, 1, resmiss.obs, 1, n_pp, 1);

  /* Refresh the cache and output the log informations */
  fflush(stdout);

  // double* debug;  // debug
  // char buffer[500];  // debug 2
  // double *LRhs = (double*) calloc(model.fixLa * (model.fixLa + 1), sizeof(double));  // debug 2
  // double *SSE = (double*) calloc(n_pp * (model.nbin + 1), sizeof(double));  // debug 2
  // int nsample = 100000; // debug 2
  // char randn[200] = "samples10w.txt";
  // model.randn = (double*)calloc(nsample, sizeof(double)); // debug 2
  // model.randidx = 0; // debug 2
  // read_dataframe_double(randn, nsample, 1, model.randn, 0, &nsample, 0); // debug 2

  /* Start MCMC sampling */
  printf("\n### Begin Gibbs sampling...\n");
  long start_all = getTimeUsec();
  long start_t = start_all;
  for (int iter = 1; iter <= para.iterAll; iter++) {
    /* reckon by time */
    // start_t = getTimeUsec();

    /* ----------- imputation of Missing residuals (phenotypes) (Gianola and Fernando 2020) ----------- */
    // sprintf(buffer, "yc%d.txt", iter);                          // debug 2
    // output_ddataframe(model.yc, M.npop * model.nind, 1, buffer, 1.0);  // debug 2
    // read_dataframe_double("ycJ.txt", M.npop * model.nind, 1, model.yc, 0, &iter, 0);  // debug
    if (model.npop > 1) sampleMissRes(&model, &resmiss, stream);
    // read_dataframe_double("ycJ.txt", M.npop * model.nind, 1, model.yc, 0, &iter, 0);  // debug
    // sprintf(buffer, "yc%d.txt", iter);                          // debug
    // output_ddataframe(model.yc, M.npop * model.nind, 1, buffer, 1.0);  // debug

    /* ----------- Location Parameters ----------- */
    yc_local_update(model.yc, model.X, model.theta, npind, fixLa, 1.0);
    /* update left and rigth hand side in MME */
    left_right_hand_update(&model, iter);
    // copy_d(model.Lhs, 1, LRhs, 1, model.fixLa * model.fixLa, 1);                 // debug 2
    // copy_d(model.Rhs, 1, LRhs + model.fixLa * model.fixLa, 1, model.fixLa, 1);   // debug 2
    // sprintf(buffer, "LRhs%d.txt", iter);                                       // debug 2
    // output_ddataframe(LRhs, model.fixLa * (model.fixLa + 1), 1, buffer, 1.0);  // debug 2
    fix_sample(&model, stream);
    // read_dataframe_double("solJ.txt", fixLa, 1, model.theta, 0, &iter, 0);  // debug
    // sprintf(buffer, "sol%d.txt", iter);                                        // debug 2
    // output_ddataframe(model.theta, fixLa, 1, buffer, 1.0);                     // debug 2
    yc_local_update(model.yc, model.X, model.theta, npind, fixLa, -1.0);
    // sprintf(buffer, "yc%d.txt", iter);                          // debug
    // output_ddataframe(model.yc, M.npop * model.nind, 1, "yc.txt", 1.0);  // debug

    /* ----------- Marker effect sampling ----------- */
    snp_start = 0;
    for (size_t bin_i = 0; bin_i < model.nbin; bin_i++) {
      if (bin_i > 0) {
        snp_start += model.nsnps[bin_i - 1];
      }
      snp_end = snp_start + model.nsnps[bin_i];

      for (size_t snp_mj = snp_start; snp_mj < snp_end; snp_mj++) {
        /* Update correction phenotype */
        // yc_alpha_update(&model, snp_mj, 1.0);
        copy_d(model.theta + fixLa + snp_mj * M.npop, 1, alpha_last, 1, M.npop, 1);

        /* (\mathbf{m}_{ij}^{'}\bigotimes \mathbf{I}_p)\mathbf{y}^{\S} */
        memset(ycAlphaj, 0, M.npop * sizeof(double));
        for (size_t pop_i = 0; pop_i < M.npop; pop_i++) {
          ycAlphaj[pop_i] = cblas_ddot(model.nind, model.Mr + snp_mj * model.nind, 1, model.yc + pop_i * model.nind, 1);
          ycAlphaj[pop_i] += model.effsqs[fixLa + snp_mj] * model.theta[fixLa + snp_mj * M.npop + pop_i];
        }

        /* R_0^{-1}(\mathbf{m}_{ij}^{'}\bigotimes \mathbf{I}_p)\mathbf{y}^{\S} */
        mat_vec_product(M.npop, M.npop, model.Re_inv, M.npop, ycAlphaj, 1, ReInvYcorr);

        /* variance of marker snp_mj */
        add(n_pp, model.VaInv + bin_i * n_pp, model.Re_inv, model.effsqs[snp_mj + fixLa], precAlphaInv);
        solve(M.npop, precAlphaInv, V_alphaij);

        /* V^{-1}R_0^{-1}(\mathbf{m}_{ij}^{'}\bigotimes \mathbf{I}_p)\mathbf{y}^{\S} */
        mat_vec_product(M.npop, M.npop, V_alphaij, M.npop, ReInvYcorr, 1, E_alphaij);

        /* Sampling of the j mark in bin i */
        alpha_sample(&model, snp_mj, bin_i, precAlphaInv, ycAlphaj, stream);
      }

      // debug = model.theta + fixLa + snp_start * M.npop;  // debug
      // sprintf(buffer, "alpha%d%zu.txt", bin_i + 1);                                 // Debug
      // read_dataframe_double("ycJ.txt", M.npop * model.nind, 1, model.yc, 0, &iter, 0);  // Debug
      // read_dataframe_double(buffer, model.nsnps[bin_i], M.npop, debug, 0, &iter, 0);           // Debug
      // sprintf(buffer, "alpha%d%zu.txt", iter, bin_i + 1);                                  // debug
      // output_ddataframe(debug, model.nsnps[bin_i], M.npop, buffer, 1.0);       // debug
      // sprintf(buffer, "yc%d.txt", iter);                          // debug
      // output_ddataframe(model.yc, M.npop * model.nind, 1, buffer, 1.0);  // debug

      /* ----------- Marker effect variance sampling ----------- */
      /* Dot product of marker effect among breeds \sum_{j=1}^{m_i}\alpha_{ij}\alpha_{ij'} */
      // memset(model.Sa, 0, n_pp * sizeof(double));
      mat_cblas_ddot(model.theta + fixLa + snp_start * M.npop, model.nsnps[bin_i], M.npop, model.Sa);
      add(n_pp, model.Sa, model.priVaS + bin_i * n_pp, 1.0, model.Sa);
      // copy_d(model.Sa, 1, SSE + bin_i * n_pp, 1, n_pp, 1);  // debug 2

      /* Sampling the variance covariance matrix of bin i */
      invWishart(model.Sa, model.df_va + model.nsnps[bin_i], M.npop, model.Va + bin_i * n_pp, stream);
      // vdRngInvWishart(model.Sa, model.df_va + model.nsnps[bin_i], M.npop, model.Va + bin_i * n_pp, stream);
      // sprintf(buffer, "Va%zuJ.txt", bin_i + 1);  // debug
      // read_dataframe_double(buffer, M.npop, M.npop, model.Va + bin_i * n_pp, 0, &iter, 0);  // debug
      // sprintf(buffer, "Va%d%zu.txt", iter, bin_i + 1);  // debug
      // output_ddataframe(model.Va + bin_i * n_pp, M.npop, M.npop, buffer, 1.0);  // debug

      /* inversion */
      solve(M.npop, model.Va + bin_i * n_pp, model.VaInv + bin_i * n_pp);
    }
    // read_dataframe_double("ycJ.txt", M.npop * model.nind, 1, model.yc, 0, &iter, 0);  // debug
    // sprintf(buffer, "yc%d.txt", iter);                                                    // debug
    // output_ddataframe(model.yc, M.npop * model.nind, 1, buffer, 1.0);                 // debug
    // sprintf(buffer, "alpha%d.txt", iter);                                                    // debug 2
    // output_ddataframe(model.theta + model.fixLa, M.nmrk, M.npop, buffer, 1.0);       // debug 2
    // sprintf(buffer, "Va%d.txt", iter);                                                       // debug 2
    // output_ddataframe(model.Va, model.nbin * M.npop, M.npop, buffer, 1.0);           // debug 2

    /* ----------- Residual effect covariance ----------- */
    /* \sum\boldsymbol e_l\boldsymbol e_l\prime */
    memset(model.Se, 0, n_pp * sizeof(double));
    for (size_t i = 0; i < M.npop; i++) {
      for (size_t j = 0; j <= i; j++) {
        /* constrain to covariance */
        if (para.constrain && i != j) continue;

        /* cblas_ddot There is a difference in the order of 1e-12 between manual loop addition and manual loop addition */
        model.Se[i * M.npop + j] = cblas_ddot(model.nind, model.yc + i * model.nind, 1, model.yc + j * model.nind, 1);
        model.Se[j * M.npop + i] = model.Se[i * M.npop + j];
      }
    }

    // Add the priori
    add(n_pp, model.Se, model.priReS, 1.0, model.Se);
    // copy_d(model.Se, 1, SSE + (model.nbin) * n_pp, 1, n_pp, 1);         // debug 2
    // sprintf(buffer, "SSE%d.txt", iter);                               // debug 2
    // output_ddataframe(SSE, (model.nbin + 1) * n_pp, 1, buffer, 1.0);  // debug 2

    /* inv-wishart */
    if (para.constrain) {
      for (size_t i = 0; i < M.npop; i++) {
        index = vdRngChiSquare(METHOD_C, stream, 1, model.Re + i * M.npop + i, model.nind + model.df_ve);
        model.Re[i * M.npop + i] = model.Se[i * M.npop + i] / model.Re[i * M.npop + i];
        if (index != 0) {
          printf("Error in vdRngChiSquare, status = %d", index);
        }
      }
    }
    else {
      invWishart(model.Se, model.df_ve + model.nind, M.npop, model.Re, stream);  // TODO iter>=2 Only slightly larger
    }

    /* vdRngInvWishart(Se, model.df_ve + model.nind, M.npop, model.Re, stream); */
    // read_dataframe_double("RJ.txt", M.npop, M.npop, model.Re, 0, &iter, 0);  // debug
    // sprintf(buffer, "R%d.txt", iter);                                  // debug 2
    // output_ddataframe(model.Re, M.npop, M.npop, buffer, 1.0);  // debug 2

    /* inversion */
    solve(M.npop, model.Re, model.Re_inv);

    /* ----------- scale parameter of residual effect covariance ----------- */
    if (SamReSigma) {
      /* scale parameter */
      add_diag(M.npop, model.Re_inv, wishartIS, 1.0);
      inverse(wishartIS, M.npop, wishartScale);
      /* sample */
      wishart(wishartScale, (model.df_ve + M.npop), M.npop, model.priReS, stream);
    }

    /* ----------- A priori scaling parameter of additive effect covariance matrix ----------- */
    // memset(wishartIS, 0, n_pp * sizeof(double));
    if (SamVaSigma) {
      if (strcmp(para.prioSa, "wishart") == 0) {
        memset(wishartIS, 0, n_pp * sizeof(double));
        /* Sum of SNP effect variances */
        for (size_t pop_i = 0; pop_i < M.npop; pop_i++) {
          /* a identity matrix as prior */
          wishartIS[pop_i * M.npop + pop_i] = 1.0;
          for (size_t pop_j = 0; pop_j < (pop_i + 1); pop_j++) {
            wishartIS[pop_i * M.npop + pop_j] +=
                cblas_ddot(model.nbin, model.nsnps, 1, model.VaInv + pop_i * M.npop + pop_j, n_pp);
            /*for (size_t snp_mj = snp_start; snp_mj < snp_end; snp_mj++)
            {
                wishartIS[pop_i * M.npop + pop_j] += model.theta[model.fixL + snp_start * M.npop + pop_i] *
            model.theta[model.fixL + snp_start * M.npop + pop_j];
            }*/
            wishartIS[pop_j * M.npop + pop_i] = wishartIS[pop_i * M.npop + pop_j];
          }
        }
        solve(M.npop, wishartIS, wishartScale);
        wishart(wishartScale, model.df_va * M.nmrk + M.npop, M.npop, model.priVaS, stream);
      }
      else if (strcmp(para.prioSa, "gamma") == 0) {
        alpha_gik = (model.df_va + M.npop) / 2.0; /* \alpha */
        for (int bin_i = 0; bin_i < model.nbin; bin_i++) {
          for (size_t pop_i = 0; pop_i < M.npop; pop_i++) {
            /* \beta */
            beta_gik = model.df_va * (1 / model.VaInv[bin_i * n_pp + pop_i * M.npop + pop_i]) + 1 / (para.A * para.A);
            /* Sampling in gamma distribution */
            errcode =
                vdRngGamma(METHOD_G, stream, 1, g + n_pp * bin_i + pop_i * M.npop + pop_i, alpha_gik, 0, beta_gik);

            if (errcode != VSL_STATUS_OK) {
              printf("error in vdRngGamma: %d\n", errcode);
              exit(errcode);
            }
          }
        }
      }
      else {
        printf("parameter prioSa can only be wishart or gamma.\n");
        exit(ERROR_PARA);
      }
    }

    /* check whether the specified number of cycles has been reached */
    if ((iter > para.burnin) && (iter % para.thin == 0)) {
      index = (iter - para.burnin) / para.thin;

      /* Total effect value */
      for (size_t i = 0; i < ntheta; i++) {
        theta_mean[i] += (model.theta[i] - theta_mean[i]) / (double)index;
      }

      /* Additive effect variance */
      for (size_t i = 0; i < varnum; i++) {
        var_mean[i] += (model.Va[i] - var_mean[i]) / (double)index;
      }

      /* Residual effect variance */
      for (size_t i = 0; i < n_pp; i++) {
        var_mean[i + varnum] += (model.Re[i] - var_mean[i + varnum]) / (double)index;
      }
    }

    /* Reporting at regular intervals */
    if ((iter % para.report_sep) == 0) {
      /* Output the specified number of variables in the MCMC chain */
      if (para.chain_var > 0) {
        for (size_t var_i = 0; var_i < para.chain_var; var_i++) {
          cover_judge[judge_i * para.chain_var + var_i] = model.theta[var_i];
        }
        judge_i++;
      }

      index = (getTimeUsec() - start_t) / 1000000;
      elasp = (getTimeUsec() - start_all) / 1000000 / 60;
      eta = ((para.iterAll - iter) / para.report_sep * index) / 60;
      printf(
          "Iter %d: %d s\telasp: %d min\teta: %d min\tFirst 2 fixed and 2 random: %.3f "
          "%.3f %.3f %.3f\n",
          iter, index, elasp, eta, model.theta[0], model.theta[model.fixL], model.theta[fixLa], model.theta[fixLa + 1]);

      /* Refresh the cache and output the log informations */
      fflush(stdout);
      start_t = getTimeUsec();
    }
  }

  printf("\nMCMC finished.\n\n");

  /* output the posterior mean of the fixed effects */
  sprintf(filename, "%s/%s", para.outdir, "sol.txt");
  output_ddataframe(theta_mean, model.fixL, M.npop, filename, 1.0);

  /* output the posterior mean of the marker effects */
  sprintf(filename, "%s/%s", para.outdir, para.effOutf);
  output_ddataframe(theta_mean + fixLa, M.nmrk, M.npop, filename, 1.0);

  /* output the posterior mean of the variance */
  sprintf(filename, "%s/%s", para.outdir, para.varOutf);
  output_ddataframe(var_mean, model.nbin + 1, n_pp, filename, 1.0);

  /* Output sample changes */
  if (para.chain_var > 0) {
    sprintf(filename, "%s/%s", para.outdir, para.mcmcOut);
    output_ddataframe(cover_judge, judge_i, para.chain_var, filename, 1.0);
  }

  /* Calculate estimated breeding value */
  out.gebv = (double*)calloc(M.nind * M.npop, sizeof(double));
  // cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, M.nind, M.npop, M.nmrk, 1.0, geno.M, M.nmrk,
  //             model.theta + fixLa, M.npop, 1.0, out.gebv, M.npop);
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, M.nind, M.npop, M.nmrk, 1.0, M.M, M.nmrk, theta_mean + fixLa,
              M.npop, 1.0, out.gebv, M.npop);

  /* output the gebv */
  for (size_t i = 0; i < M.npop; i++) {
    sprintf(filename, "%s/%s_y%zu.txt", para.outdir, para.gebvOutf, i + 1);
    output_ddataframe_names(out.gebv + i, M.nind, 1, M.npop, filename, "ID\tEBV", geno.iid);
  }

  /* free memory */
  free(ycAlphaj);
  free(theta_mean);
  free(var_mean);
  free(precAlphaInv);
  free(ReInvYcorr);
  free(wishartIS);
  free(wishartScale);
  free(model.Mr);
  free(model.Re);
  free(g);
  free(E_alphaij);
  free(V_alphaij);
  free(cover_judge);

  time(&t);
  ctime_r(&t, time_buf);
  printf("Run ended at %s", time_buf);

  /* close the log file */
  if (strlen(para.logf) != 0) {
    fclose(stdout);
  }

  return 0;
}

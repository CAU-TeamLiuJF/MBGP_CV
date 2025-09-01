/**
 * @file multi_pop_fun.h
 * @brief Declare functions for model evaluation
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

#include <argp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "mkl.h"
#include "mkl_pblas.h"

#define BRNG VSL_BRNG_MCG31
#define METHOD_N VSL_RNG_METHOD_GAUSSIAN_ICDF
#define METHOD_G VSL_RNG_METHOD_GAMMA_GNORM
#define METHOD_C VSL_RNG_METHOD_CHISQUARE_CHI2GAMMA
#define STOR_N VSL_MATRIX_STORAGE_FULL
#define MAX_ID_LEN 30      /* Maximum length of individual ID */
#define MAX_POPS 10        /* Maximum number of populations */
#define MISS_GENE_CODE 3   /* Missing code in allele content matrix */

/* -------------- structure definition -------------- */

/**
 * @brief information for allele content matrix
 */
typedef struct
{
  int nbin;       /* number of bins */
  int nsnp;       /* Number of SNPs when grouping with a fixed number of SNPs */
  double *nsnps;  /* number of SNPs in each bin */
  int *nbins;     /* number of bins in each chromosome */
  int *chr;       /* Which chromosome is bin located on */
  int *start;     /* Physical location where bin starts */
  int *end;       /* Physical location where bin ends */
  int *mean;      /* The mean value of r2 at the breakpoint */
  double *sum2pq; /* \sum{2pq} in each bin */
  double thresh;  /* The r2 threshold value used for bins division based on LD [0.1] */
  int win;        /* R2 will be computed between SNPs that are at most this far apart [100] */
  int min;        /* minimum number of SNPs in each bin [50] */
  char *type;     /* Method for dividing snp into different groups, LD/fix [LD] */
} Bin;

/**
 * @brief Parameters related to model operation
 * @details
 */
struct Model
{
  int nind_gen;   /* number of individuals with genotype information */
  int nind_ebv;   /* Number of individuals with only genetic information but no phenotypic information */
  int nind;       /* Number of individuals with both genetic information and phenotypic information */
  int nmrk;       /* Number of markers */
  double df_va;   /* Degree of freedom of prior distribution of genetic variance */
  double df_ve;   /* Degree of freedom of prior distribution of residual variance */
  double h2;      /* Heritability [0.5] */
  double rg;      /* Genetic correlation */
  int npop;       /* Number of populations */
  int *npops;     /* Number of individuals in each population */
  int nfix;       /* Number of fixed effects */
  int *fixLs;     /* Number of levels of each fixed effect */
  int *fix_nobs;  /* Number of observations of each levels in all fixed effect */
  int fixL;       /* Sum of levels of all fixed effects within breed */
  int fixLa;      /* Sum of levels of all fixed effects among breeds */
  int ran_level;  /* Sum of levels of all fixed effects */
  double *effsqs; /* Sum of squares of effect values within a single level (marker) */
  int nRan;       /* Number of random effects */
  int nbin;       /* Number of genomic intervals */
  double *Va;     /* Additive genetic variance */
  double *VaInv;  /* Inverse of additive genetic variances */
  double *nsnps;  /* Number of SNPs in each genomic interval */
  double *X;      /* incidence matrix, PS: MKL call needs float / double type */
  double *Mr;     /* Additive effect gene content matrix with 0/1/2 of reference population */
  int Mr_store;   /* storage format in reference group, same as above */
  double *M_ebv;  /* gene content matrix of individuals without phenotype */
  double *yc;     /* corrected phenotype */
  int *obs;       /* Column of non-missing phenotype */
  int *ntype_obs; /* Number of individuals with specific non-missing phenotype types */
  double *Re;     /* Residual effect variance matrix */
  double *Re_inv; /* inversion of residual effect variance matrix */
  double *R;      /* Individual residual correlation matrix */
  double *priVaS; /* scale factor for genetic variance */
  double *priReS; /* scale factor for residual variance */
  double *Sa;     /* genetic effect vector inner product */
  double *Se;     /* residual effect vector inner product */
  double *Lhs;    /* Left-hand side of MME */
  double *Rhs;    /* Right-hand side of MME */
  double *theta;  /* parameters */
  double *tmp_pp; /* temporary array */
  double *tmp_p;  /* temporary array */
  double *X1pX1;  /* \sum{2pq} */
  int ld;         /* leading dimension of y and yc */
  int inc;        /* increment of value in y and yc */
  struct Genotype *geno;   /* Structure for storing genotype parameters */

  // double *randn;   // debug 2
  // int randidx;     // debug 2
};

/**
 * @brief Parameters of parse_opt
 */
struct arguments
{
  int thread;          /* number of threads for the OpenMP run-time library [1] */
  int require_para;    /* Number of parameters that must be provided [3] */
  int A;               /* A sufficiently large value, which is used when priosa is gamma [100000] */
  int nsnp_bin;        /* if binf file is not provided, bin is defined according to a fixed number of markers [100] */
  int iterAll;         /* length of MCMC chain [1000] */
  int burnin;          /* burnin steps [100] */
  int thin;            /* thinning rate [10] */
  int report_sep;      /* MCMC sampling report interval [100] */
  int chain_var;       /* output first c prediction in MCMC chain [NULL] */
  int constrain;       /* constrain the covariance between the residuals to zero [FALSE] */
  long seed;           /* random number seed [81239] */
  char *outdir;        /* Output file directory [NULL] */
  int binPos;          /* column containing the number of SNPs in each bins at binf [1] */
  char *raw;           /* plink raw file of individuals [NULL] */
  char *prefix;        /* plink binary file prefix of individuals [NULL] */
  char *binf;          /* file contains number of SNPs per interval [NULL] */
  char *varf;          /* SNP variance prior [varp * 0.5] */
  char *effOutf;       /* posteriori mean of location effects output file name [theta.txt] */
  char *varOutf;       /* posteriori mean of dispersion matrix output file name [var.txt] */
  char *gebvOutf;      /* breeding value output file [gebv.txt] */
  char *mcmcOut;       /* output file names of prediction in MCMC chain [MCMC_process.txt] */
  char *chr;           /* tmp */
  char *out;           /* tmp */
  char *logf;          /* redirect standard output to this file */
  char *prioSa;        /* distribution of scale of additive effect prior [NULL/wishart/gamma] */
  char *prioSe;        /* distribution of scale of residual effect prior [NULL/wishart] */
  double maf;          /* filters out all SNPs with MAF below the provided threshold [0.01] */
  double geno;         /* filters out all SNPs with missing call rates exceeding the provided value [0.1] */
  double mind;         /* filters out all samples with missing call rates exceeding the provided value [0.1] */
  struct Pheno *phe;   /* Structure for storing phenotype parameters */
  struct Model *model; /* Structure for storing model parameters */
  struct Genotype *gene;   /* Structure for storing genotype parameters */
  Bin *bin;                /* bin information */
};

/**
 * @brief parameters for obtaining information from phenotype file
 */
struct Pheno
{
  char *filename; /* filename of phenotype file */
  char *phePos;   /* phenotypes position (separated by comma) at phenotype file, e.g. 4,5 [NULL] */
  char *ranPos;   /* environmental random effect position at phenotype file [NULL] */
  double miss;    /* reals below this value are regarded as missing [-99] */
  int idPos;      /* id position at phenotype file [1] */
  int cols;       /* Number of columns in the file */
  int rows;       /* Number of rows in the file */
  int nind;       /* Number of individuals */
  double *y;      /* phenotypes */
  char *fixPos;   /* fixed effect position (separated by comma) at phenotype file, e.g. 2,3 [NULL] */
  int *fixCols;   /* Line number of fixed effect */
  int *triCols;   /* Line number of phenotype */
  int ntriCol;    /* number of lines in phenotype */
  int header;     /* whether the phenotype file has header */
  int nindPhe;    /* number of individuals with phenotypes */
  char **ID;      /* Individual ID in phenotype file */
  char **pheID;   /* Individual ID with phenotype */
  int type;       /* Phenotype storage format, 1 for matrix (nind x npop) and 0 for vector (npops x 1) */
  int *pheIDInx;  /* index of Individual ID with phenotype */
};

/**
 * @brief parameters for obtaining information from fam file
 */
typedef struct
{
  char *fid;    /* Family ID ('FID') */
  char *iid;    /* Within-family ID ('IID'; cannot be '0') */
  char *father; /* Within-family ID of father ('0' if mother isn't in dataset) */
  char *mother; /* Within-family ID of mother ('0' if mother isn't in dataset) */
  int sex;                 /* Sex code ('1' = male, '2' = female, '0' = unknown) */
  double phe;     /* Phenotype value ('1' = control, '2' = case, '-9'/'0'/non-numeric = missing data if case/control) */
  int nmiss;      /* Number of missing SNPs in this individual */
} Fam;

/**
 * @brief parameters for obtaining information from bim file
 */
typedef struct
{
  char *chr;   /* Chromosome code (either an integer, or 'X'/'Y'/'XY'/'MT'; '0' indicates unknown) or name */
  char *id;    /* Variant identifier */
  double cm;  /* Position in morgans or centimorgans (safe to use dummy value of '0') */
  int pos;    /* Base-pair coordinate (1-based; limited to 231-2) */
  char A1[2];    /* Allele 1 (corresponding to clear bits in .bed; usually minor) only snp */
  char A2[2];    /* Allele 2 (corresponding to set bits in .bed; usually major) */
  // int nmiss;  /* Number of missing SNPs in this SNP */
} Bim;

/**
 * @brief information for allele content matrix
 */
typedef struct
{
  double *M;      /* Additive effect gene content matrix with 0/1/2 of reference and validation population */
  int format;     /* storage format of gene content matrix, 1 for nmrk x nind and 0 for nind x nmrk [0] */
  int ld;         /* leading dimension of matrix M */
  int inc;        /* increment of value in matrix M */
  int nind;       /* number of individuals */
  int nmrk;       /* number of SNPs */
  int center;     /* Whether to centralize the gene content matrix [1] */
  int nchr;       /* Number of chromosomes */
  int *chr_pos;   /* Chromosome index of SNPs */
  int npop;       /* number of breeds */
  char **chr;     /* unique chromosome identifiers */
  int *nsnps;     /* Number of SNPs in each chromosome */
  int *snp_index; /* The snp index (0/1) in this matrix relative to the complete matrix (from the genotype file) */
  int *ind_index; /* The ind index (0/1) in this matrix relative to the complete matrix (from the genotype file) */
  double *frq;    /* gene frequency */
  Fam *fam;       /* fam file information of Ma */
  Bim *bim;       /* bim file information of Ma */
  int *ninds;      /* number of individuals in each breed */
  char **fids;    /* identifiers of breeds */
  // parameters for storage format
  int row;        /* number of lines in matrix M */
  int col;        /* number of columns in matrix M */
  int *row_index; /* row index */
  int *col_index; /* column index */
  int *mrk_miss;  /* number of missing samples in each SNP and populations */
} Gmat;

/**
 * @brief parameters for obtaining information from genotype file
 */
struct Genotype
{
  char *filename;  /* filename of plink raw format file */
  char *prefix;    /* plink binary file prefix */
  char *mapf;      /* plink map(or bim) file */
  char **iid;      /* individuals IDs, use to output GEBV */
  // int nind_org;    /* number of individuals */
  // int nind;        /* number of individuals */
  // int nmrk;        /* number of markers */
  // int nmrk_org;    /* number of markers */
  // int npop;        /* number of breeds */
  // int *ninds;      /* number of breeds */
  // int *pos;        /* Physical location of the markers */
  // char **chr;      /* chromosome identifier of each SNP */
  // char **chr_uniq; /* unique chromosome identifiers */
  // int nchr;        /* Number of chromosomes */
  // int *chr_pos;    /* Chromosome index of SNPs */
  // int *nsnp_chr;   /* Number of SNPs in each chromosome */
  // double *frq;     /* gene frequency */
  // double *sum2pq;  /* \sum{2pq} */
  // Fam *fam;        /* fam file information of Ma */
  // Bim *bim;        /* bim file information of Ma */
  Gmat Ma;         /* full gene content matrix from plink file */
  Gmat Mq;         /* gene content matrix after quality control */
  Gmat Mr;         /* gene content matrix of reference population */
  // int M_store;     /* storage format of gene content matrix, 1 for nmrk x nind and 0 for nind x nmrk [0] */
  // double *M;       /* Additive effect gene content matrix with 0/1/2 of reference and validation population */
  // int center;      /* Whether to centralize the gene content matrix [1] */
  /* quality control */
};

/**
 * @brief Parameters required to output the result
 */
struct Output
{
  double *gebv;      /* Genome breeding value */
  char **genoID;     /* Genotype individual ID */
  char **pheID;      /* Phenotype file individual ID */
  char **popN;       /* breed/population ID */
  int nind_phe;      /* Number of individuals in the phenotype file */
  char **snpid;      /* snp ID */
  char **fix_labels; /* Effect level labels */
};

/**
 * @brief Parameters required to output the result
 */
struct RseSample
{
  double *sample; /*  */
  double *tmp;    /*  */
  double *E_ym;   /*  */
  double *V_ym;   /*  */
  double *R0_mo;  /*  */
  double *R0_mm;  /*  */
  int *obs;    /*  */
  double R0_oo;  /*  */
  int nmiss;      /*  */
};

// -------------- declarations of functions -------------- */

void inverse(double *in, int n, double *a);
void wishart(double *V, double df, int p, double *C, VSLStreamStatePtr stream);
void invWishart(double *V, double df, int p, double *C, VSLStreamStatePtr stream);
void solve(int n, double *A, double *B);
void vdRngInvWishart(double *S, double df, int n, double *sample, VSLStreamStatePtr stream);
void vdRngNor(int num, int d, double *mu, double *Sigma, double *sample, VSLStreamStatePtr stream);
long getTimeUsec();
void priorMat(int n, double *pvar, double per, double cor, double *out, int nbin, double *sum2pq);
void gtcenter(double *M, int nmrk, int nind, int trans);
// void sampleMissRes(struct Model *model, VSLStreamStatePtr stream);
void sampleMissRes(struct Model* model, struct RseSample* resmp, VSLStreamStatePtr stream);
void make_incidence(struct Pheno *phe, struct Model *model, struct Output *output);
void make_incidence_old(struct Pheno *phe, struct Model *model, struct Output *output);
void obs_type(int npop, int *ninds, double *y, double miss, int type, int *obs_index);
void make_bins(struct Genotype *geno, Gmat *M, Bin *bin);
void make_bins_fix(Gmat *M, Bin *bin);
void make_bins_LD(Gmat *M, Bin *bin);
void one_col_to_multi(double* vec, int npop, int* npops, double miss, double* multi, int type);
void single_phe(struct Pheno *phe, struct Model *model, double miss);
void load_phenotype(struct Pheno *pheno, struct Model *model, char *phe_pos, double miss);
void multi_phe_to_vec(struct Pheno *phe, struct Model *model, double miss);
void single_phe_to_vec(struct Pheno *phe, struct Model *model, double miss);
void left_right_hand_update(struct Model* model, int iter);
// void fix_sample(double *Lhs, double *Rhs, double *theta, int n, VSLStreamStatePtr stream);
void fix_sample(struct Model* model, VSLStreamStatePtr stream);
void yc_local_update(double *yc, double *X, double *sol, const MKL_INT npind, const MKL_INT fixLa, const double sign);
void yc_alpha_update(struct Model *model, double *snp_last, int snp_pos);
void mat_cblas_ddot(double *mat, int row, int col, double *S);
void gene_frq(Gmat *M);
void alpha_sample(struct Model *model, int snp_pos, int bin_i, double *C, double *wRinv, VSLStreamStatePtr stream);
void get_ref_M(int nref, char **refID, char **ID, int nmrk, int nind, double *M, int type_M, double *Mref, int type_Mr);
void center(double *M, int nind, int nmrk, int type);
void read_fam(char *prefix, Gmat *M, char ***iid);
void read_bim(char *prefix, char *mapf, Gmat *M);
void read_bed(char *prefix, Gmat *M);
void genotype_qc(Gmat *org, Gmat *new, double maf, double nind, double geno);
void plink_binary_parse(char *prefix, struct Genotype *gene);
void plink_raw_parse(char *rawf, Gmat *M, char ***iid);
void sum2pq_group(int npop, int nmrk, int nbin, double *nsnps, double *frq, double *sum2pq);
void plink_parse(char *rawf, char* prefix, struct Genotype *gene);
double gene_content_correlation(double *va, double frqa, int inca, double *vb, double frqb, int incb, int n);
void gmat_format(Gmat *mat, int aim);
void copy_gmat(Gmat *to, Gmat *from);
void init_gmat(Gmat *M);
void remove_individual(Gmat *M, int popi, int location);
void remove_SNP(Gmat *M, int chri, int location);
void compute_r2(Gmat *M, int chr, int start, int nmrk, int win, double *mean);
void read_allele(char *to, char *from);
void fid_order(Gmat *M);

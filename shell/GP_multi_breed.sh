#!/usr/bin/bash


########################################################################################################################
## Version: 1.3.1
## Author:    Liweining liwn@jaas.ac.cn
## Orcid:     0000-0002-0578-3812
## Institute: College of Animal Science and Technology, China Agricul-tural University, Haidian, 100193, Beijing, China
## Date:      2024-07-05
##
## Function：
## To estimate the accuracy of the joint prediction, it is necessary to first complete the within-breed genomic prediction
##  (i.e. there are val */rep * folders in different breed (population) directories) before running this script
##
## Usage: ./GP_multi_breed.sh --pops "breedA breedB" ...(Please refer to --help for detailed parameters)
##
## Dependent software/environment:
##  1. R
##  2. plink/1.9
##  3. gmatrix
##  4. mbBayesABLD
##  5. Other R languages and Bash scripts
##
## License:
##  This script is licensed under the GPL-3.0 License.
##  See https://www.gnu.org/licenses/gpl-3.0.en.html for details.
########################################################################################################################

###################  Parameter processing  ###################
##############################################################
## NOTE: This requires GNU getopt.  On Mac OS X and FreeBSD, you have to install this
TEMP=$(getopt -o h4 --long ref:,rg_local,keep_all,nbin:,binf:,summs:,dirPre:,noCov,debug,traceplot,bfileW:,bfileM:,pops:,pedf:,priorVar:,DIR:,vidf:,method:,rg:,re:,type:,GmatM:,gmat:,gidf:,invA:,all_eff:,ran_eff:,tbvf:,tbv_col:,phereal:,add_rf:,fold:,rep:,miss:,bincol:,num_int:,win:,nsnp_win:,r2_merge:,bin:,LD_maxD:,r2:,inter:,fix_snp:,seed:,iter:,burnin:,thin:,report_sep:,code:,thread:,alpha:,out:,nchr:,updatepri:,prefix:,VarA,dmu4,overwrite,suffix,dense \
  -n 'javawrap' -- "$@")
if [ $? != 0 ]; then
  echo "Terminating..." >&2
  exit 1
fi
eval set -- "$TEMP"

## Parse arguments
while true; do
  case "$1" in
  --pops )      pops="$2";       shift 2 ;; ## Breed A identifiers, such as 'YY DD'
  --bfileW )    bfileW="$2";     shift 2 ;; ## Genotype files for each breed, like './data/YY ./data/DD'
  --bfileM )    bfileM="$2";     shift 2 ;; ## Prefix for merged PLINK binary files
  --binf )      binf="$2";       shift 2 ;; ## genome block (bin) file
  --bincol )    bincol="$2";     shift 2 ;; ## Column in the block file indicating the number of markers within block [1]
  --pedf )      pedf="$2";       shift 2 ;; ## Pedigree file, can be a single merged pedigree or a space-separated list files for each breed
  --summs )     summs="$2";      shift 2 ;; ## Prefix for GEMMA results files
  --priorVar )  priorVar="$2";   shift 2 ;; ## Initial values for variance components, can be null/predict/pheno/A/B/... [pheno]
  --h2 )        h2="$2";         shift 2 ;; ## Heritability used to generate initial variance components, like '0.5' or '0.3 0.1' [0.5]
  --priorRg )   rg="$2";         shift 2 ;; ## Prior for additive genetic correlation between Breed A and B (two-trait model) [0.001]
  --priorRe )   re="$2";         shift 2 ;; ## Prior for residual correlation between breed A and B (two-trait model) [0.001]
  --DIR )       DIR="$2";        shift 2 ;; ## Parameter card prefix [type]
  --method )    method="$2";     shift 2 ;; ## Genomic prediction model, options are BLUP/GBLUP/ssGBLUP/bayesR/mbBayesAB [GBLUP]
  --type )      type="$2";       shift 2 ;; ## Joint prediction model, options are single/multi [single]
  --bayes )     bayes="$2";      shift 2 ;; ## Bayesian method for multi-trait prediction, options are mbBayesAB or bayesR [mbBayesAB]
  --GmatM )     GmatM="$2";      shift 2 ;; ## Method to construct genotype matrix [multi/single]
  --gmat )      gmat="$2";       shift 2 ;; ## User-provided relationship matrix or inverse matrix (id id value)
  --gidf )      gidf="$2";       shift 2 ;; ## Genotype individual IDs, consistent with the IDs in the user-specified G matrix file
  --invA )      invA="$2";       shift 2 ;; ## A-inverse construction method (1/2/3/4/6), refer to DMU manual for details [1]
  --all_eff  )  all_eff="$2";    shift 2 ;; ## Third line in $MODEL in DIR (all effects), the first three are ignored, e.g., "2 1", separated by '-' if different between breeds
  --ran_eff  )  ran_eff="$2";    shift 2 ;; ## Fourth line in $MODEL in DIR (random effects), the first one is ignored, e.g., "1", separated by '-' if different between breeds
  --tbvf )      tbvf="$2";       shift 2 ;; ## File containing TBVss
  --tbv_col )   tbv_col="$2";    shift 2 ;; ## Column number of TBVs in the phenotype file, should be "same" if the TBVs is the same as the corresponding phenotype in the phenotype file
  --phereal )   phereal="$2";    shift 2 ;; ## Column number of the real phenotype in the phenotype file
  --add_rf )    add_rf="$2";     shift 2 ;; ## Group in which additive effects are located
  --fold )      fold="$2";       shift 2 ;; ## Cross-validation fold number
  --rep )       rep="$2";        shift 2 ;; ## Number of repetitions
  --miss )      miss="$2";       shift 2 ;; ## Missing phenotype identifier [-99]
  --num_int )   num_int="$2";    shift 2 ;; ## Number of integer variables in the phenotype file [obtained from breed A's phenotype file]
  --nsnp_win )  nsnp_win="$2";   shift 2 ;; ## Number of SNPs per window before merging
  --win )       win="$2";        shift 2 ;; ## Window size when calculating the mean R2 for specified SNPs [50]
  --bin )       bin="$2";        shift 2 ;; ## Whether to merge adjacent windows, options are fix/frq/ld/lava/cubic [lava]
  --maf )       maf="$2";        shift 2 ;; ## Filter SNPs based on the specified minor allele frequency threshold when bin is cubic [-0.01]
  --nbin )      nbin="$2";       shift 2 ;; ## Number of blocks, constrained to be close to nbin when bin is fixed
  --ref )       ref="$2";        shift 2 ;; ## Reference panel for dividing blocks, options are M/1[index of breed]/2/... [M]
  --r2_merge )  r2_merge="$2";   shift 2 ;; ## LD threshold for merging adjacent windows
  --LD_maxD )   LD_maxD="$2";    shift 2 ;; ## PLINK parameter, maximum distance in kb for SNP pairs when calculating LD
  --r2 )        r2="$2";         shift 2 ;; ## PLINK parameter, R2 threshold in the output file when calculating LD
  --inter )     inter="$2";      shift 2 ;; ## PLINK parameter, LD is not calculated for marker pairs with more than inter intervening loci
  --seed )      seed="$2";       shift 2 ;; ## Random seed for MCMC sampling
  --iter )      iter="$2";       shift 2 ;; ## Total number of MCMC iterations
  --burnin )    burnin="$2";     shift 2 ;; ## Number of MCMC burn-in iterations
  --thin )      thin="$2";       shift 2 ;; ## MCMC sampling interval
  --report )    report_sep="$2"; shift 2 ;; ## MCMC report interval
  --code )      code="$2";       shift 2 ;; ## Script path
  --thread )    thread="$2";     shift 2 ;; ## Number of parallel DMU tasks
  --alpha )     alpha="$2";      shift 2 ;; ## ssGBLUP relationship matrix alpha
  --vidf )      vidf="$2";       shift 2 ;; ## Filename of the validation breed individual ID file
  --out )       out="$2";        shift 2 ;; ## Prefix for the output accuracy filename
  --dirPre )    dirPre="$2";     shift 2 ;; ## Prefix added to the output folder for EBV
  --updatepri ) updatepri="$2";  shift 2 ;; ## Update prior variance-covariance matrix at the specified round [0]
  --nchr )      nchr="$2" ;      shift 2 ;; ## Number of chromosomes [30]
  --prefix )    prefix="$2";     shift 2 ;; ## Add the specified prefix to the single/multi folder
  --suffix )    suffix=true;     shift   ;; ## Add breed name suffix to single/multi folder
  --dense )     dense=true;      shift   ;; ## Set method in DMU's ANALYSE to 31 to utilize multi-threading for variance component calculations
  --debug )     debug=true;      shift   ;; ## Skip time-consuming steps like DMU, gmatrix, bayes, etc.
  --rg_local )  rg_local=true;   shift   ;; ## Different covariance for each genomic region, for local LD, frequency correlation coefficients
  --noCov )     noCov=true;      shift   ;; ## Constrain residual effects between traits to 0
  --keep_all )  keep_all=true;   shift   ;; ## Keep all Bayes result files
  --traceplot ) traceplot=true;  shift   ;; ## Trace plot after MCMC burn-in period
  --overwrite ) overwrite=true;  shift   ;; ## Overwrite if preparation files (phenotype, ID files, etc.) already exist
  -4 | --dmu4 ) dmu4=true;       shift   ;; ## Use dmu4 for prediction, variance components will not be estimated
  -h | --help)  grep ";; ##" $0 | grep -v help && exit 1 ;;
  -- ) shift; break ;;
  * ) break ;;
  esac
done

## Working directory
workdir=$(pwd)

## Log directory
logp=${workdir}/log
mkdir -p ${logp}

## Script directory
if [[ ${code} ]]; then
  [[ ! -d ${code} ]] && echo "${code} not exists! " && exit 5
else
  script_path=$(dirname "$(readlink -f "$0")")
  code=$(dirname "$script_path")
fi

## Scripts
pheMerge=${code}/R/multibreed_pheno.R
accur_cal=${code}/R/accuracy_bias_calculation.R
fix_frq_ld_bolck=${code}/R/fix_frq_ld_bolck.R
lava_cubic_bolck=${code}/shell/lava_cubic_bolck.sh
job_pool=${code}/shell/parallel_jobs_func.sh
func=${code}/shell/function.sh
multiG=${code}/R/multibreed_relationship_matrix.R
variance_prior=${code}/R/variance_prior_setting.R
MCMC_polt=${code}/R/MCMC_chain_plot.R
local_rg=${code}/R/local_block_rg.R
bayesR=${code}/shell/BayesR.sh

## Load custom functions
[[ ! -s ${func} ]] && echo "Error: ${func} not found! " && exit 5
source ${job_pool} ## Parallel computing
source ${func}

## Add program paths to environment variables
export PATH=${code}/bin:$PATH

## Check if required programs are in the environment path and executable
check_command plink gmatrix mbBayesABLD LD_mean_r2 run_dmu4 run_dmuai

## Check if required script files exist and are executable
check_command $pheMerge $accur_cal $fix_frq_ld_bolck $lava_cubic_bolck $job_pool $func
check_command $multiG $variance_prior $MCMC_polt $local_rg

## Prevent warnings when executing R scripts ("ignoring environment value of R_HOME")
unset R_HOME

## Directory suffix
if [ "$suffix" = true ]; then
  suffix="_$(echo "${pops}" | tr ' ' '_')"
else
  suffix=""
fi

## Directory to store results
tpath=${workdir}/${prefix}${type}${suffix}
mkdir -p ${tpath}

## Default parameters
nchr=${nchr:=30}
updatepri=${updatepri:=0}
LD_maxD=${LD_maxD:=10000}
r2=${r2:=0}
inter=${inter:=99999}
nsnp_win=${nsnp_win:=100}
win=${win:=50}
bin=${bin:=lava}
ref=${ref:=M}
bincol=${bincol:=1}
fold=${fold:=1}
rep=${rep:=1}
miss=${miss:=-99}
type=${type:=single}
phereal=${phereal:=1}
h2=${h2:=0.5}
rg=${rg:=0.001}
re=${re:=0.001}
iter=${iter:=30000}
burnin=${burnin:=20000}
thin=${thin:=10}
DIR=${DIR:=${type}}
report_sep=${report_sep:=100}
GmatM=${GmatM:="single"}
invA=${invA:=1}
alpha=${alpha:=0.05}
method=${method:=GBLUP}
maf=${maf:=-0.01}
vidf=${vidf:=val.id}
r2_merge=${r2_merge:=0.1}
add_rf=${add_rf:=1}
geno=${geno:=0.2}
mind=${mind:=0.2}
code=${code:=${HOME}/liwn/code}
[[ ${tbvf} && ! ${tbv_col} ]] && tbv_col=2 ## If the TBVs file exists, assume the TBVs is in the 2rd column by default, with IDs in the 1st column.

## Parameter initialization
[[ ${noCov} ]] && noCov=" --nocov "
# [[ ${all_samps} ]] && all_samps=" --all_samp_out "
[[ ${keep_all} ]] && keep_all=" --keep_all "
if [[ $(echo ${h2} | awk '{print NF}') -gt 1 ]]; then
  h2B=$(echo ${h2} | awk '{print $2}')
  h2=$(echo ${h2} | awk '{print $1}')
else
  h2B=${h2}
fi

## Random seed
if [[ ! ${seed} ]]; then
  seed=$RANDOM
  echo "$seed" >MCMC_seed.txt
fi

## Accuracy output file name
if [[ ! ${out} ]]; then
  out=accur_${method}
fi

## Set the number of parallel jobs
if [[ ! ${thread} ]]; then
  if [[ ${SLURM_CPUS_ON_NODE} ]]; then thread=${SLURM_CPUS_ON_NODE}; else thread=$((rep * fold)); fi
fi

## Initialize job pool
[[ ! ${debug} ]] && job_pool_init ${thread} 0

##################  Parse command line parameters  ##############
#################################################################
np=$(echo ${pops} | awk '{print NF}')
IFS=" " read -r -a popN <<<"$pops"
IFS=" " read -r -a bfileWa <<<"$bfileW"

## DMU multi-threading
if [[ ${method} == "GBLUP" && ${dense} ]]; then
  method_dmu="31"
else
  method_dmu="1"
fi

###############  Check if within-breed prediction is complete  ###########
##########################################################################
## Cross-validation and repetition
rep=$(find ${workdir}/${popN[0]}/val1/rep* -type d | wc -l)
fold=$(find ${workdir}/${popN[0]}/val*/rep1 -type d | wc -l)

## Check if files exist in each subset
for r in $(seq 1 ${rep}); do
  for f in $(seq 1 ${fold}); do
    for b in "${popN[@]}"; do
      for file in pheno.txt ${vidf}; do
        [[ ! -s ${workdir}/${b}/val${f}/rep${r}/${file} ]] &&
          echo "${workdir}/${b}/val${f}/rep${r}/${file} not found! " &&
          exit 1

        ## Get the number of integer and real variables in the phenotype file
        if [[ ! ${num_int} ]]; then
          phef_within=${workdir}/${b}/val${f}/rep${r}/pheno.txt
          ncol=$(awk 'END{print NF}' ${phef_within})
          for i in $(seq 1 ${ncol}); do
            dot=$(awk -vl=${i} '{print $l}' ${phef_within} | grep -c "\.")
            [[ ${dot} -gt 0 ]] && num_int=$((i - 1)) && break
          done
        fi
        [[ ! ${num_real} ]] && num_real=$(($(awk 'END{print NF}' ${phef_within}) - num_int))
      done
    done
  done
done

##################  Retrieve Model Effect Settings Within Breeds  #############
###############################################################################
## Retrieve effect settings from one of the breed subsets (effects are the same across breeds, use the first breed as reference)
firstB=$(echo ${pops} | cut -d ' ' -f 1)
within_DIR=$(find ${workdir}/${firstB}/val1/rep1 -name "*.DIR" | head -n 1)
[[ ! -s ${within_DIR} ]] && echo "${within_DIR} not found! " && exit 1
model_line=$(grep -n "MODEL" ${within_DIR} | head -n 1 | cut -d ':' -f 1)
all_eff=$(sed -n "$((model_line + 3))p" ${within_DIR} | cut -d ' ' -f '4-')
ran_eff=$(sed -n "$((model_line + 4))p" ${within_DIR} | cut -d ' ' -f '2-')
## Number of Effects
nA=$(echo ${all_eff} | awk '{print NF}')
nR=$(echo ${ran_eff} | awk '{print NF}')
## Add fixed effects for breeds
((nA++))
all_eff="$((num_int + 1)) ${all_eff}"

################  Genotype Files  #############
###############################################
## Paths to genotype files within breeds
# IFS=" " read -r -a bfiles <<<"${bfileM} ${bfileW[*]}"
if [[ ! ${debug} ]]; then
  if [[ ! ${bfileW} ]]; then
    IF="\n" mapfile -t bfiles < <(printf "%s\n" "${popN[@]}" | xargs -I {} echo "${workdir}/{}/{}")
  else
    IFS=" " read -r -a bfiles <<<"${bfileW}"
  fi
  ## Check if fam/bim/bed files exist
  check_plink "${bfiles[@]}" ${nchr}
  ## Generate merged plink files for all breeds
  if [[ ! -s ${gmat} ]]; then
    if [[ ${bfileM} ]]; then
      check_plink "${bfileM}" ${nchr}
    else
      ## Merge plink files
      bfileM=${tpath}/merge
      merge_plink "${bfiles[*]}" ${bfileM}
    fi
  fi
fi

##################  Reference Population Phenotype Files  ##############
########################################################################
$pheMerge \
  --pops "${pops}" \
  --phef "${workdir}/#breed#/val#val#/rep#rep#/pheno.txt" \
  --rep ${rep} \
  --fold ${fold} \
  --nInt ${num_int} \
  --type ${type} \
  --method ${method} \
  --overwri ${overwrite} \
  --pheCol ${phereal} \
  --fixPop \
  --out "${tpath}/val#val#/rep#rep#/pheno.txt"

#####################  Pedigree Files  #####################
############################################################
if [[ ${pedf} ]]; then
  for pedi in ${pedf}; do
    if [[ ! -s ${pedi} ]]; then
      echo "${pedi} not found! " && exit 1
    else
      cat ${pedi} >pedi_merge.txt
    fi
  done
  echo "Number of individuals in pedigree: $(wc -l <pedi_merge.txt)"
else
  [[ -s pedi_merge.txt ]] && rm pedi_merge.txt
fi

####################  Determining Genomic Regions  ##################
#####################################################################
if [[ ${method} == 'mbBayesAB' ]]; then
  cd ${tpath} || exit

  if [[ ! -s ${binf} || ${overwrite} ]]; then
    ## Filename prepare
    if [[ ${bin} == "fix" ]]; then
    bin_prefix=${bin}_${nsnp_win}
    elif [[ ${bin} == "lava" ]]; then
    bin_prefix=${bin}_${ref}_${nsnp_win}
    elif [[ ${bin} == "cubic" ]]; then
    bin_prefix=${bin}_${ref}_${win}
    elif [[ ${bin} == "frq" || ${bin} == "ld" ]]; then
    bin_prefix=${bin}_${r2_merge}_${nsnp_win}
    else
    echo "${bin} can only be fix, frq, ld or lava! "
    exit 1
    fi
    binf=${tpath}/${bin_prefix}.txt

    ## Partition blocks
    echo 'Generating genome regions file...'
    if [[ ${bin} == "fix" ]]; then
      if [[ ${nbin} ]]; then
        nsnp=$(wc -l <${bfileM}.bim)
        ## Ensure the final number of partitions is close to nbin
        nsnp_win=$((nsnp / nbin))
        echo "Number of SNPs per window set to: ${nsnp_win}"
      fi

      [[ ! ${debug} ]] && \
        $fix_frq_ld_bolck \
          --win ${nsnp_win} \
          --map ${bfileM}.bim \
          --bin_merge ${bin} \
          --out ${binf}
    elif [[ ${bin} == "frq" || ${bin} == "ld" ]]; then
      # ## Check if genotype files for each breed are provided
      # [[ ! "${bfileW[*]}" ]] && echo "Required parameter 'bfileW' is missing! " && exit 1

      ## Common markers among breeds
      awk '{print $2}' ${bfileM}.bim >SNP_share_id.txt

      ## LD calculation
      for i in $(seq 0 $((np - 1))); do
        plink --bfile ${bfiles[i]} \
          --chr-set ${nchr} \
          --extract SNP_share_id.txt \
          --r2 \
          --freq \
          --ld-window-kb ${LD_maxD} \
          --ld-window ${inter} \
          --ld-window-r2 ${r2} \
          --out ${popN[i]} >>${logp}/plink_ld_frq.log
      done

      ## Define blocks as required
      [[ ! ${debug} ]] && \
        $fix_frq_ld_bolck \
          --prefixs "${popN[*]}" \
          --bin_merge ${bin} \
          --win ${nsnp_win} \
          --seg ${r2_merge} \
          --map ${bfileM}.bim \
          --out ${binf}~

      ## Extract a single-column file indicating the number of SNPs in each block
      awk '{print $3}' ${binf}~ >${binf}
    elif [[ ${bin} == "lava" || ${bin} == "cubic" ]]; then
      ## Reference panel (population genomic data) used for partitioning
      if [[ ${ref} == "M" ]]; then
        ## Use the merged panel for all breeds
        bfile_block=bfile${ref}  
        bfile_block=${!bfile_block}
      else
        ## Use the panel of the reference breed, where the order matches fid in pops
        bfile_block=${bfiles[((ref - 1))]}
      fi

      ## Generate partitioning file
      [[ ! ${debug} ]] && \
        $lava_cubic_bolck \
          --bfile ${bfile_block} \
          --win ${win} \
          --maf ${maf} \
          --type ${bin} \
          --minSize ${nsnp_win} \
          --out ${binf}~ >${logp}/${bin}_${ref}_block_${SLURM_JOB_ID}.log

      ## Extract a single-column file indicating the number of SNPs in each bin
      sed '1d' ${binf}~ | awk '{print $5}' >${binf}
    fi

    ## Check for errors
    if [[ ! -s ${binf} ]]; then
      echo 'error in creating bins file! '
      exit 1
    fi
  else
    # binf_base=$(basename "${binf}")
    # bin_prefix="${binf_base%.*}"
    bin_prefix=self_bin
    if [[ ${bincol} != 1 ]]; then
      awk -v col=${bincol} '{print $col}' ${binf} > bin_col1.txt
      binf=$(pwd)/bin_col1.txt
    fi
  fi
fi

#####################  Local Genetic Correlation  #####################
#######################################################################
if [[ ${rg_local} && ${bin} == "lava" ]]; then
  summA=
  summB=
  if [[ -s ${summA}.assoc.txt && -s ${summB}.assoc.txt ]]; then
    echo "calculating local genetic correlations..."
    ## Software needs modification
    $local_rg \
      --pops ${pops} \
      --summ1 ${summA}.assoc.txt \
      --summ2 ${summB}.assoc.txt \
      --bfile ${!bfile_block} \
      --block ${binf}~ \
      --out ${binf}~ &>${logp}/${bin}_local_rg_${SLURM_JOB_ID}.log
  else
    echo "Error: ${summA}.assoc.txt or ${summB}.assoc.txt not found! "
    exit 1
  fi
fi

#####################  Variance Components Preparation  ####################
############################################################################
if [[ ${priorVar} ]]; then
  ## Variance components file name
  parfA=$(basename "$(find ${workdir}/${A} -name "*PAROUT" | head -n 1)")
  parfB=$(basename "$(find ${workdir}/${B} -name "*PAROUT" | head -n 1)")

  ## Different covariance priors for each block
  [[ ${rg_local} ]] && rg_local=" --rg_local ${binf}~ "

  if [[ ${priorVar} == "pheno" ]]; then
    ## Calculate genetic and residual variances based on phenotypic variance
    varf_para1="${workdir}/${A}/val#val#/rep#rep#/pheno.txt"
    varf_para2="${workdir}/${B}/val#val#/rep#rep#/pheno.txt"
  elif [[ ${priorVar} == "predict" ]]; then
    ## Calculate multi-breed genetic and residual variances based on within-breed prediction variances
    varf_para1="${workdir}/${A}/val#val#/rep#rep#/${parfA}"
    varf_para2="${workdir}/${B}/val#val#/rep#rep#/${parfB}"
  elif [[ ${priorVar} == "A" ]]; then
    ## Use variance components estimated within breed A
    varf_para1="${workdir}/${A}/val#val#/rep#rep#/${parfA}"
    varf_para2="null"
  elif [[ ${priorVar} == "B" ]]; then
    ## Use variance components estimated within breed B
    varf_para1="${workdir}/${B}/val#val#/rep#rep#/${parfB}"
    varf_para2="null"
  fi

  ## Generate variance components (initial values) file
  $variance_prior \
    --filea ${varf_para1} \
    --fileb ${varf_para2} \
    --rep ${rep} \
    --fold ${fold} \
    --h2 ${h2} \
    --h2B ${h2B} \
    --rg ${rg} \
    --re ${re} \
    --type ${type} \
    --method ${method} \
    --var ${priorVar} \
    --norec \
    --add_rf ${add_rf} \
    --overwri ${overwrite} \
    --pcol $((num_int + phereal)) \
    ${rg_local} \
    --out ${tpath}/val#val#/rep#rep#/${type}_${dirPre}${bin}_prior.txt
fi

#####################  Relationship Matrix Preparation  #####################
#############################################################################
if [[ ${method} =~ 'GBLUP' ]]; then
  ## Generate the Genotype Relationship Matrix (G-matrix)
  if [[ -s ${gmat} ]]; then
    if [[ ${method} == "ssGBLUP" && ! -s ${gidf} ]]; then
      echo "${gidf} not found! " && exit 1
    else
      echo "Use user provided relationship matrix: ${gmat}"
    fi
  elif [[ ${GmatM} == 'multi' ]]; then
    ## Multi-breed relationship matrix
    [[ ! "${bfileW[*]}" ]] && echo "Required parameter 'bfileW' is missing! " && exit 1

    :> all_breed.id
    for i in $(seq 0 $((np - 1))); do
      awk '{print $2}' ${bfileWa[i]}.fam >${bfileW[i]}.ids
      cat ${bfileW[i]}.ids >>all_breed.id
      idf="${bfileW[i]}.ids ${idf}"
    done

    ## Format conversion, remove loci with missing values (loci present in only one breed)
    plink \
      --bfile ${bfileM} \
      --chr-set ${nchr} \
      --geno 0 \
      --recode A \
      --out merge >>${logp}/plink_recodeA.log

    ## Multi-breed relationship matrix (need to be modified)
    $multiG \
      --rawf merge.raw \
      --idf ${idf} \
      --out merge
    mv merge.grm merge.agrm.id_fmt
  elif [[ ${GmatM} == 'single' ]]; then
    ## Generate the Genotype Relationship Matrix (G-matrix)
    if [[ ${method} == "GBLUP" ]]; then
      gmat_inv="--inv"
      [[ ! ${gmat} ]] && gmat=${tpath}/merge.agiv.id_fmt
    elif [[ ${method} == "ssGBLUP" ]]; then
      gmat_inv=""
      [[ ! ${gmat} ]] && gmat=${tpath}/merge.agrm.id_fmt
      [[ ! ${gidf} ]] && gidf=${tpath}/merge.id
    fi

    ## Construct the Genomic Relationship Matrix
    echo "Read the plink bed file and Calculate the additive G matrix..."
    [[ ! ${debug} ]] && gmatrix --bfile ${bfileM} --grm agrm --out ${tpath}/merge ${gmat_inv}
    echo "G matrix created."

    ## Rename
    if [[ ${method} == "GBLUP" ]]; then
      ## If gmat file name is specified, rename (and possibly move) the file
      [[ ${gmat} && ! -s ${gmat} ]] && \
        mv merge.agiv.id_fmt ${gmat} && \
        echo "gmat created and has been renamed to: ${gmat}"
    elif [[ ${method} == "ssGBLUP" ]]; then
      ## If gmat file name is specified, rename (and possibly move) the file
      [[ ${gmat} && ! -s ${gmat} ]] && \
        mv merge.agrm.id_fmt ${gmat} && \
        echo "gmat created and has been renamed to: ${gmat}"
      [[ ${gidf} && ! -s ${gidf} ]] && mv ${tpath}/merge.id ${gidf}
    fi
  else
    echo "Gmat can only be multi or single! "
    exit 1
  fi
fi

###################  DMU Parameter File  ####################
#############################################################
if [[ ${method} =~ 'GBLUP' ]]; then
  cd ${tpath} || exit

  ## ANALYSE
  [[ ${dmu4} ]] && ANALYSE="11 9 0 0" || ANALYSE="1 ${method_dmu} 0 0"

  ## MODEL
  if [[ ${type} == 'single' ]]; then
    num_real=1
    MODEL="1\n0\n1 0 ${nA} ${all_eff}\n${nR} ${ran_eff}\n0\n0"
  else
    ABSORB="0"
    MODELS="1 0 ${nA} ${all_eff}"
    RANDOMS="${nR} ${ran_eff}"
    REGRES="0"
    NOCOV=""
    nNOCOV=0
    for i in $(seq 2 ${np}); do
      ABSORB="${ABSORB}\n0"
      MODELS="${MODELS}\n${i} 0 ${nA} ${all_eff}"
      RANDOMS="${RANDOMS}\n${nR} ${ran_eff}" # RANDOM is a shell keyword, so renamed
      REGRES="${REGRES}\n0"
      for j in $(seq 1 $((i - 1))); do
        NOCOV="${NOCOV}\n${j} ${i}"
        ((nNOCOV++))
      done
    done

    num_real=${np}
    MODEL="${np}\\n${ABSORB}\\n${MODELS}\\n${RANDOMS}\\n${REGRES}\\n${nNOCOV}${NOCOV}"
  fi

  ## Write out the parameter file
  [[ -s ${DIR}.DIR ]] && echo "warn: ${DIR}.DIR has been overwritten! "
  {
    echo "\$COMMENT"
    echo "get EBV of individuals in validation population with reduced phenotypes"
    echo "\$ANALYSE ${ANALYSE}"
    echo "\$DATA  ASCII ($((num_int + 1)), ${num_real}, ${miss}) pheno.txt"
    echo -e "\$MODEL\n${MODEL}"
    echo "\$PRIOR %varf%"
  } >${DIR}.DIR

  ## Variance component structure
  if [[ ${method} == "GBLUP" ]]; then
    echo "\$VAR_STR ${add_rf} GREL ASCII ${gmat}" >>${DIR}.DIR
    add_sol=3
  elif [[ ${method} == "ssGBLUP" ]]; then
    [[ ! -s pedi_merge.txt ]] && echo "pedi_merge.txt not found! " && exit 1
    echo "\$VAR_STR ${add_rf} PGMIX ${invA} ASCII pedi_merge.txt ${gidf} ${tpath}/merge.agrm.id_fmt ${alpha} G-ADJUST" >>${DIR}.DIR
    add_sol=4
  elif [[ ${method} == "BLUP" ]]; then
    [[ ! -s pedi_merge.txt ]] && echo "pedi_merge.txt not found! " && exit 1
    echo "\$VAR_STR ${add_rf} PED ${invA} ASCII pedi_merge.txt" >>${DIR}.DIR
    add_sol=4
  else
    echo "method can only be BLUP/GBLUP/ssGBLUP! "
    exit 1
  fi

  ## Output effect values
  echo "\$SOLUTION" >>${DIR}.DIR
fi


#####################  Subset Processing and prediction  #####################
##############################################################################
for r in $(seq 1 ${rep}); do
  for f in $(seq 1 ${fold}); do
    ## Change working directory
    vali_path=${tpath}/val${f}/rep${r}
    cd ${vali_path} || exit

    ## Skip breeding value estimation if in debug mode
    [[ ${debug} ]] && continue

    ## (ss)GBLUP Model
    if [[ ${method} =~ 'GBLUP' ]]; then
      ## Copy DMU parameter file
      cp ${tpath}/${DIR}.DIR .

      ## Provide variance components (or initial values) as needed
      if [[ -s ${type}_${dirPre}${bin}_prior.txt ]]; then
        sed -i "s#%varf%#${type}_${dirPre}${bin}_prior.txt#" ${DIR}.DIR
      else
        sed -i '/$PRIOR/d' ${DIR}.DIR
      fi

      ## Breeding value estimation
      if [[ ${dmu4} ]]; then
        job_pool_run /usr/bin/time -v run_dmu4 ${DIR} ${vali_path}
      else
        job_pool_run /usr/bin/time -v run_dmuai ${DIR} ${vali_path}
      fi
    fi

    ## Fixed effects
    fix_eff=${all_eff%" ${ran_eff}"}

    ## MT-BayesABLD Model
    if [[ ${method} == 'mbBayesAB' ]]; then
      ## Get GEBVs
      job_pool_run /usr/bin/time -v mbBayesABLD \
        --bfile ${bfileM} \
        --phef ${vali_path}/pheno.txt \
        --fix "${fix_eff}" \
        --binf ${binf} \
        --iter ${iter} \
        --burnin ${burnin} \
        --thin ${thin} \
        --outdir ${vali_path} \
        --report ${report_sep} \
        --varOut var_${bin}_${r2_merge}_${nsnp_win}.txt \
        --effOut effect_${bin}_${r2_merge}_${nsnp_win}.txt \
        --gebvOut EBV_${bin} \
        --mcmcOut MCMC_process_${bin}.txt \
        --seed ${seed} \
        --logf ${bin}_${r2_merge}_${nsnp_win}_gibs_${SLURM_JOB_ID}.log \
        ${noCov}
    fi

    if [[ ${method} == 'bayesR' ]]; then
      if [[ ${type} == 'single' ]]; then
        ## Get GEBVs
        job_pool_run /usr/bin/time -v $bayesR \
          --proj "${vali_path}" \
          --bfile "${bfileM}" \
          --fix "${fix_eff}" \
          --iter ${iter} \
          --burnin ${burnin} \
          --seed ${seed} \
          --gebv_out "EBV_bayesR.txt"
      else
        echo "Not currently supported! "
        exit 2
      fi
    fi
  done
done

## Wait for background jobs to finish
[[ ! ${debug} ]] && job_pool_wait
## Release cpus
[[ ! ${debug} ]] && job_pool_shutdown

####################  Accuracy Calculation  ##################
##############################################################
## Change working directory
cd ${tpath} || exit

## Accuracy calculation parameters
if [[ ${tbv_col} == "same" ]]; then
  option=" --tbv_col $((num_int + phereal))"
elif [[ ${tbv_col} ]]; then
  option=" --tbv_col ${tbv_col}"
fi
[[ ${tbvf} ]] && option="${option} --tbvf ${tbvf} "

## Other parameters
if [[ ${method} =~ 'BLUP' ]]; then
  option="${option} --add_sol ${add_sol} --dir_val ${tpath}/val#val#/rep#rep#/${DIR}"
  if [[ ${type} == 'single' ]]; then
    ## ST-GBLUP model
    ebv_col=1
  else
    ## MT-GBLUP model
    ebv_col=%i%
  fi
else
  ## GEBVs filename
  if [[ ${method} == "bayesR" ]]; then
    gebvf=EBV_bayesR.txt
  else
    gebvf=EBV_${bin}_y%i%.txt
  fi

  ## MT-Bayes model
  option="${option} --ebvf ${tpath}/val#val#/rep#rep#/${gebvf}"
  ebv_col=2
fi

## Calculate accuracy
for i in $(seq 0 $((np - 1))); do
  ip1=$((i+1))
  ## Output file name
  if [[ ${method} == 'mbBayesAB' ]]; then
    outf=${out}_${dirPre}_${bin}_${popN[${i}]}.txt
    [[ ${bin} == 'lava' ]] && outf=${out}_${ref}_${dirPre}_${bin}_${popN[${i}]}.txt
  else
    outf=${out}_${popN[${i}]}.txt
  fi

  ## Replace possible double underscores
  outf=${outf/__/_}

  ## Calculate accuracy
  $accur_cal \
    --ebv_col ${ebv_col/\%i\%/${ip1}} \
    --val_idf ${workdir}/${popN[${i}]}/val#val#/rep#rep#/${vidf} \
    --fid ${popN[${i}]} \
    --rep ${rep} \
    --fold ${fold} \
    ${option/\%i\%/${ip1}} \
    --out ${tpath}/${outf}
done

#################  BayesABLD Model Convergence Plotting  #################
##########################################################################
## BayesABLD Model Convergence Plotting
if [[ ${type} == 'multi' ]]; then
  if [[ ${traceplot} ]]; then
    for r in $(seq 1 ${rep}); do
      ## Plot Bayes MCMC Chains
      for f in $(seq 1 ${fold}); do
        $MCMC_polt \
          --files ${vali_path}/${bin}_MCMC_process.txt \
          --start ${burnin} \
          --end ${iter} \
          --thin ${report_sep} \
          --names "μ1 μ2 alpha11 alpha12" \
          --out "${vali_path}/${bin}_MCMC_process"
      done
    done
  fi
fi

###################  Delete Intermediate Files  #####################
#####################################################################
## Genotype files
# [[ ${bfileM} && ${bfileM} =~ rmMiss ]] && rm ${bfileM}.*

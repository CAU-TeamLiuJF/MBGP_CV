#!/usr/bin/bash

########################################################################################################################
## Version: 1.3.0
## Author:    Liweining liwn@cau.edu.cn
## Orcid:     0000-0002-0578-3812
## Institute: College of Animal Science and Technology, China Agricul-tural University, Haidian, 100193, Beijing, China
## Date:      2024-08-20
##
## Function：
## Used for accuracy calculation based on BLUP/Bayes cross validation
##
##
## ./single_breed.sh --phef /path/to/your/phenotype ...(Please refer to --help for detailed parameters)
##
## Dependent on software/environment:
##  1. R
##  2. plink/1.9
##  3. gmatrix
##  4. mbBayesABLD
##  5. Other R languages and Bash scripts
##
##
## License:
##  This script is licensed under the GPL-3.0 License.
##  See https://www.gnu.org/licenses/gpl-3.0.en.html for details.
########################################################################################################################


###################  Parameter processing  #####################
####################################################
## NOTE: This requires GNU getopt.  On Mac OS X and FreeBSD, you have to install this
## Parse arguments
TEMP=$(getopt -o 4hp:b:v:m:o: --long rmNeg,intercept,debug,label:,phef:,bfile:,varf:,pedf:,gmat:,DIR:,rep:,fold:,add_rf:,gen:,year:,iyse:,tbvf:,tbv_col:,phereal:,binf:,seed:,all_eff:,iter:,burnin:,bin:,thin:,report:,ran_eff:,method:,thread:,miss:,gidf:,add_sol:,invA:,num_int:,code:,alpha:,out:,overWri,valphe,dmu4,dense,help \
  -n 'javawrap' -- "$@")
if [ $? != 0 ]; then
  echo "Terminating..." >&2
  exit 1
fi
eval set -- "$TEMP"
##  parse parameter
while true; do
  case "$1" in
  -p |--phef)    phef="$2";    shift 2 ;; ## Phenotype file
  -b |--bfile)   bfile="$2";   shift 2 ;; ## plink binary file prefix
      --label)   label="$2";   shift 2 ;; ## Breed label
  -v |--varf)    varf="$2";    shift 2 ;; ## Variance component file (if not provided, estimated from pedigree/SNP information)
      --pedf)    pedf="$2";    shift 2 ;; ## Pedigree file (if not provided, run GBLUP)
  -m |--gmat)    gmat="$2";    shift 2 ;; ## User-provided relationship matrix or inverse matrix (id id value)
      --gidf)    gidf="$2";    shift 2 ;; ## Genotyped individual IDs, consistent with the IDs in the user-specified G-matrix file
      --DIR)     DIR="$2";     shift 2 ;; ## Parameter file prefix [within]
      --rep)     rep="$2";     shift 2 ;; ## Cross-validation repeat number [1]
      --fold)    fold="$2";    shift 2 ;; ## Number of cross-validation folds [1]
      --add_rf)  add_rf="$2";  shift 2 ;; ## Additive effect group, default is the 1st random effect group [1]
      --gen)     gen="$2";     shift 2 ;; ## Individuals of the last gen generations with phenotypes are the validation group [ALL]
      --year)    year="$2";    shift 2 ;; ## Individuals born after this year with phenotypes are the validation group
      --iyse)    iyse="$2";    shift 2 ;; ## Required when providing year, format as idcol:yearcol:years_start:years_end
      --tbvf)    tbvf="$2";    shift 2 ;; ## TBV file; if tbv_col is provided but tbvf is missing, tbvf is set to the phenotype file
      --tbv_col) tbv_col="$2"; shift 2 ;; ## Column of TBVs in the phenotype file; if TBV is the same as phenotype, set to "same"
      --phereal) phereal="$2"; shift 2 ;; ## Column position of the phenotype in the phenotype file [1]
      --all_eff) all_eff="$2"; shift 2 ;; ## All effects columns in $MODEL in the 3rd line of DIR; the first 3 digits are not needed,
                                          ## only the effect columns are required, e.g., "2 3 1"
      --ran_eff) ran_eff="$2"; shift 2 ;; ## All random effects group in $MODEL in the 4th line of DIR; the first digit is not needed, 
                                          ## only the random effect groups are required, e.g., "1" [1]
      --method)  method="$2";  shift 2 ;; ## Method for breeding value estimation, can be BLUP/GBLUP/ssGBLUP/mbBayesAB/BayesR [GBLUP]
      --bin )    bin="$2";     shift 2 ;; ## Whether to merge adjacent windows, fix/frq/ld/ind/cubic [ind]
      --binf )   binf="$2";    shift 2 ;; ## Regions file
      --seed )   seed="$2";    shift 2 ;; ## Random seed for MCMC sampling and validation group partitioning [40296]
      --dense )  dense=true;   shift   ;; ## Set method in DMU's ANALYSE to 31, use multi-threading to calculate variance components
      --iter )   iter="$2";    shift 2 ;; ## Total number of MCMC iterations
      --burnin ) burnin="$2";  shift 2 ;; ## MCMC burn-in iterations
      --thin )   thin="$2";    shift 2 ;; ## MCMC sampling interval
      --report ) report="$2";  shift 2 ;; ## MCMC reporting interval
      --thread)  thread="$2";  shift 2 ;; ## Number of parallel dmu tasks
      --miss)    miss="$2";    shift 2 ;; ## Missing phenotype identifier [-99]
      --invA)    invA="$2";    shift 2 ;; ## A-inverse construction method (1/2/3/4/6), 1 considers inbreeding, 2 does not. (see DMU manual) [1]
      --num_int) num_int="$2"; shift 2 ;; ## Number of integer columns [calculated based on the file]
      --code)    code="$2";    shift 2 ;; ## Path to code [/public/home/liujf/liwn/code]
      --alpha)   alpha="$2";   shift 2 ;; ## G-matrix correction coefficient (whether to consider inbreeding) [0.05]
      --intercept) mean=true;  shift   ;; ## Add a population mean column after the last integer column
      --rmNeg)    rmNeg=true;  shift   ;; ## Remove negative values after the last integer column
      --valphe)   valphe=true; shift   ;; ## Individuals in the validation group only need to have phenotypes, not genotypes
      --overWri) overWri=true; shift   ;; ## overWri if the file to be generated already exists
      --debug)    debug=true;  shift   ;; ## Do not run DMU and gmatrix
  -o | --out)     out="$2";    shift 2 ;; ## Output corrected phenotype file name
  -4 | --dmu4)    dmu4=true;   shift   ;; ## Use DMU4 model to estimate breeding values
  -h | --help)    grep ";; ##" $0 | grep -v help && exit 1 ;;
  -- ) shift; break ;;
  * ) break ;;
  esac
done

## Avoid warnings when running R scripts ("ignoring environment value of R_HOME")
unset R_HOME

## Default parameters
seed=${seed:=40296}                           ## Random seed
nchr=${nchr:=30}                              ## Number of chromosomes
DIR=${DIR:=within}                            ## Parameter card name
iter=${iter:=30000}                           ## Total number of MCMC iterations
burnin=${burnin:=20000}                       ## MCMC burn-in iterations
thin=${thin:=10}                              ## MCMC sampling interval
report=${report:=100}                         ## MCMC reporting interval
alpha=${alpha:=0.05}                          ## G-matrix parameter
thread=${thread:=1}                           ## Number of parallel dmu tasks
phereal=${phereal:=1}                         ## Phenotype column
add_rf=${add_rf:=1}                           ## Additive random effect group
ran_eff=${ran_eff:=1}                         ## Random effects
rep=${rep:=1}                                 ## Repeats
fold=${fold:=1}                               ## Number of validation groups
miss=${miss:=-99}                             ## Missing phenotype identifier
method=${method:=GBLUP}                       ## Breeding value estimation method
out=${out:=accuracy.txt}                      ## Accuracy output file name
invA=${invA:=2}                               ## A-inverse construction method, default considers inbreeding
[[ ${tbv_col} && ! ${tbvf} ]] && tbvf=${phef} ## File containing TBVs
[[ ${tbvf} && ! ${tbv_col} ]] && tbv_col=2    ## TBV file exists, default TBV in the second column, id in the first column

## Main folder
workdir=$(pwd)

## Log folder
logp=${workdir}/log
mkdir -p ${logp}

## Script folder
if [[ ${code} ]]; then
  [[ ! -d ${code} ]] && echo "${code} not exists! " && exit 5
else
  script_path=$(dirname "$(readlink -f "$0")")
  code=$(dirname "$script_path")
fi

## Scripts
phe_group=${code}/R/validation_population_define.R
accur_cal=${code}/R/accuracy_bias_calculation.R
keep_phe_gid=${code}/R/keep_pheno_geno_individuals.R
fix_file=${code}/R/bayesR_file.R  ## Genotype file and fixed Structural matrix for bayesR
bayesR=${code}/shell/BayesR.sh
job_pool=${code}/shell/parallel_jobs_func.sh
func=${code}/shell/function.sh

## Load custom functions
[[ ! -s ${func} ]] && echo "Error: ${func} not found! " && exit 5
source ${func}
source ${job_pool}

## Add program path to environment variables
export PATH=${code}/bin:$PATH

## Check if the necessary parameters are provided
if [[ ! -s ${phef} ]]; then
  echo "phenotype file ${phef} not found! "
  exit 1
elif [[ ! -s ${bfile}.fam && ! -s ${bfile}.map && ! -s ${gidf} && ! -s ${pedf} ]]; then
  echo "plink file ${bfile}.fam, pedigree file ${pedf} or genotyped individuals id file ${gidf} not found! "
  exit 1
elif [[ -s ${gmat} && ${method} == "ssGBLUP" && ! -s ${gidf} ]]; then
  echo "genotyped individuals id file ${gidf} not found! "
  exit 1
elif [[ ! ${label} ]]; then
  echo "breed label parameter --label not found! "
  exit 1
fi

## Check if the required programs are available in the environment path and executable
check_command plink gmatrix mbBayesABLD LD_mean_r2 run_dmu4 run_dmuai

## Check if the required script files exist and have execution permissions
check_command $phe_group $accur_cal $keep_phe_gid $job_pool $func $bayesR


##################  Parse command-line arguments  ##############
################################################################

## DMU multi-threading
if [[ ${method} == "GBLUP" || ${dense} ]]; then
  method_dmu="31"
else
  method_dmu="1"
fi

# ## Job threads
# [[ ${SLURM_CPUS_ON_NODE} ]] && thread=${SLURM_CPUS_ON_NODE}

## Random seed
if [[ ! ${seed} ]]; then
  seed=$RANDOM
  echo "$seed" >MCMC_seed.txt
fi

## Determine if validation group partitioning and other steps have been completed
val_phe="${workdir}/val1/rep1/pheno.txt"
if [[ ! -s ${val_phe} || ${overWri} ]]; then
  overWri=true
else
  unset overWri
fi

## Initialize the job pool
job_pool_init ${thread} 0

## Copy the phenotype file to the working directory
cp ${phef} ${workdir}/pheno_within.txt
phef=${workdir}/pheno_within.txt

## Number of integer and real variable columns in the phenotype file
if [[ ! ${num_int} ]]; then
  ncol=$(awk 'END{print NF}' ${phef})
  for i in $(seq 1 ${ncol}); do
    dot=$(awk -vl=${i} '{print $l}' ${phef} | grep -c "\.")
    [[ ${dot} -gt 0 ]] && num_int=$((i - 1)) && break
  done
fi
num_real=$(($(awk 'END{print NF}' ${phef}) - num_int))

## Number of effects
nA=$(echo ${all_eff} | awk '{print NF}')
nR=$(echo ${ran_eff} | awk '{print NF}')

###################  Genotype File  ####################
########################################################
if [[ ${bfile} && ! -s ${bfile}.bim ]]; then
  if [[ ! -s ${bfile}.map ]]; then
    echo "plink file ${bfile}.map not found! "
    exit 1
  else
    ## Extract markers on the specified chromosomes
    plink \
      --file ${bfile} \
      --chr-set ${nchr} \
      --make-bed --out ${label} >${logp}/plink_single_copy.log
    bfile=$(pwd)/${label}
  fi
elif [[ -s ${bfile}.bim ]]; then
  plink \
    --bfile ${bfile} \
    --chr-set ${nchr} \
    --make-bed --out ${label} >${logp}/plink_single_copy.log
  bfile=$(pwd)/${label}
fi

#####################  Phenotype File  ####################
###########################################################
phe_col=$((phereal + num_int))
## Select individuals with both genotype and phenotype as the reference group
if [[ -s ${bfile}.fam ]]; then
  $keep_phe_gid \
    --famf "${bfile}.fam" \
    --phef ${phef} \
    --rm single \
    --num_int ${num_int} \
    --phec ${phe_col} \
    --rmOut gid_miss_phe.txt
fi
## Remove genotyped individuals without phenotype information
if [[ -s gid_miss_phe.txt ]]; then
  echo "remove $(wc -l <gid_miss_phe.txt) individuals without phenotype"
  plink \
    --bfile ${bfile} \
    --chr-set ${nchr} \
    --remove gid_miss_phe.txt \
    --make-bed --out ${bfile} >${logp}/plink_single_rm_phemiss.log
fi
## Add a column of intercept (all "1") after integer columns in the phenotype file
if [[ ${mean} ]]; then
  ## Parameter card information
  ((nA++))
  ((num_int++))
  all_eff="${num_int} ${all_eff}"

  ## Phenotype file
  awk -v column="${num_int}" -v value="1" '
    BEGIN {
        FS = OFS = " ";
    }
    {
        for ( i = NF + 1; i > column; i-- ) {
            $i = $(i-1);
        }
        $i = value;
        print $0;
    }
    ' ${phef} >${phef}.tmp
  echo "add populations mean column (${num_int}) in the phenotype file."
  mv ${phef}.tmp ${phef}
fi
## Check if there are any individuals left in the phenotype file
if [[ ! -s ${phef} ]]; then
  echo "no individuals in phenotype file ${phef}! "
  [[ ${workdir} =~ "mbGS" ]] && rm -r ${workdir}
  exit 1
fi

###################  DMU Parameter Card Template  ####################
######################################################################
if [[ ${method} =~ "BLUP" ]]; then
  ## ANALYSE
  if [[ ${dmu4} ]]; then
    ANALYSE="11 9 0 0"
  else
    ANALYSE="1 ${method_dmu} 0 0"
  fi

  [[ -s ${DIR}.DIR ]] && echo "warn: ${DIR}.DIR will be overwritten! "
  {
    echo "\$COMMENT"
    echo "get EBV of individuals in validation population with reduced phenotypes"
    echo "\$ANALYSE ${ANALYSE}"
    echo "\$DATA  ASCII (${num_int}, ${num_real}, ${miss}) pheno.txt"
    echo -e "\$MODEL\n1\n0\n${phereal} 0 ${nA} ${all_eff}\n${nR} ${ran_eff}\n0\n0"
    echo '$VAR_STR %VAR_STR%'
    echo '$PRIOR %PRIOR%'
    echo '$SOLUTION'
  } >${DIR}.DIR
fi

#####################  Variance Structure  #####################
################################################################
if [[ ${method} == "GBLUP" ]]; then
  ## Variance-covariance structure file
  if [[ ${gmat} ]]; then
    if [[ ! -s ${gmat} ]]; then
      echo "${gmat} not found! "
      exit 1
    else
      echo "Use the user specified genome relationship matrix: ${gmat}"
    fi
  else
    echo "Read the plink bed file and Calculate the additive G matrix..."
    [[ ! ${debug} ]] && gmatrix --bfile ${bfile} --grm agrm --out full --inv
    echo "G matrix created."
    gmat=${workdir}/full.agiv.id_fmt
  fi
  ## Parameter card
  sed -i "s#%VAR_STR%#${add_rf} GREL ASCII ${gmat}#g" ${DIR}.DIR
  ## Additive effect code in the SOL results file
  add_sol=3
elif [[ ${method} == "BLUP" ]]; then
  ## BLUP
  sed -i "s#%VAR_STR%#${add_rf} PED ${invA} ASCII ${pedf}#g" ${DIR}.DIR
  add_sol=4
elif [[ ${method} == "ssGBLUP" ]]; then
  if [[ ${gmat} ]]; then
    [[ ! -s ${gmat} ]] && echo "${gmat} not found! " && exit 1
    [[ ! -s ${gidf} ]] && echo "${gidf} not found! " && exit 1
  else
    echo "Read the plink bed file and Calculate the additive G matrix..."
    [[ ! ${debug} ]] && gmatrix --bfile ${bfile} --grm agrm --out full
    echo "G matrix created."
    gmat=${workdir}/full.agrm.id_fmt
    gidf=${workdir}/full.id
  fi
  sed -i "s#%VAR_STR%#${add_rf} PGMIX ${invA} ASCII ${pedf} ${gidf} ${gmat} ${alpha} G-ADJUST#g" ${DIR}.DIR
  add_sol=4
fi

## Specify variance components (initial values)
if [[ ${method} =~ "BLUP" ]]; then
  if [[ -s ${varf} ]]; then
    sed -i "s#%PRIOR%#${varf}#g" ${DIR}.DIR
  else
    if [[ ${dmu4} ]]; then
      echo "${varf} not found! "
      exit 1
    else
      [[ ${varf} ]] && echo "${varf} not found! "
      sed -i '/$PRIOR.*/d' ${DIR}.DIR
    fi
  fi
fi

###################  Validation Group Division  ####################
####################################################################
## Parameter settings
if [[ ${overWri} ]]; then
  option=' '
  [[ ${gen} ]] && option="${option} --gen ${gen} --pedf ${pedf}"               ## Divide by generations
  [[ ${year} && ${iyse} ]] && option="${option} --year ${year} --iyse ${iyse}" ## Divide by birth year
  [[ ! ${valphe} && -s ${bfile}.fam ]] && option="${option} --fam ${bfile}.fam"

  ## Divide into reference group and validation group
  $phe_group \
    --phef "${phef}" \
    --nonmiss "${all_eff}" \
    --rep ${rep} \
    --fold ${fold} \
    --seed ${seed} \
    --outvid "val.id" \
    --outdir "${workdir}/val#val#/rep#rep#" \
    --rminvail \
    --keepid "${workdir}/keep_fid_id.txt" \
    --pheCol ${phe_col} \
    ${option}

  ## Keep only the individuals that can serve as the reference group (with phenotype, genotype, and non-missing fixed effects)
  if [[ -s ${workdir}/keep_fid_id.txt ]]; then
    plink --bfile ${bfile} \
      --chr-set ${nchr} \
      --keep ${workdir}/keep_fid_id.txt \
      --make-bed --out ${bfile} >${logp}/plink_full_set.log

    rm ${bfile}.fam~ ${bfile}.bim~ ${bfile}.bed~
    echo "keep $(wc -l <${workdir}/keep_fid_id.txt) individuals in geneotype file."
  fi
fi

###################  Estimation of Breeding Value for Validation Group  #####################
#############################################################################################
for r in $(seq 1 ${rep}); do
  for f in $(seq 1 ${fold}); do
    ## Replace missing values
    sed -i "s/NA/${miss}/Ig" ${workdir}/val${f}/rep${r}/pheno.txt
    check_alphabet ${workdir}/val${f}/rep${r}/pheno.txt

    ## Fixed effects
    fix_eff=${all_eff%" ${ran_eff}"}

    ## mbBayesAB model
    if [[ ${method} == "mbBayesAB" ]]; then
      ## If binf file is not provided, run BayesA
      if [[ ! -s ${binf} ]]; then
        binf=${workdir}/val${f}/rep${r}/fixed_bins_1.txt
        awk '{print "1"}' ${bfile}.bim >${binf}
        echo "Fit bayesA with only one snp in each bins"
      else
        ## Divide types
        case ${binf} in
          *fix*) bin="fix" ;;
          *lava*) bin="lava" ;;
          *cubic*) bin="cubic" ;;
          *ld*) bin="ld" ;;
          *frq*) bin="frq" ;;
          *) bin="unknown" ;;
        esac
      fi

      ## Run Bayes model
      job_pool_run mbBayesABLD \
        --bfile ${bfile} \
        --phef ${workdir}/val${f}/rep${r}/pheno.txt \
        --fix "${fix_eff}" \
        --phe ${phe_col} \
        --binf ${binf} \
        --iter ${iter} \
        --burnin ${burnin} \
        --thin ${thin} \
        --outdir ${workdir}/val${f}/rep${r} \
        --report ${report} \
        --varOut var_${bin}.txt \
        --effOut effect_${bin}.txt \
        --gebvOut EBV_${bin} \
        --mcmcOut MCMC_process_${bin}.txt \
        --seed ${seed} \
        --logf ${bin}_gibs_${SLURM_JOB_ID}.log
      gebvf=EBV_${bin}_y1.txt
    fi

    ## mbBayesAB model
    if [[ ${method} == "bayesR" ]]; then
      job_pool_run $bayesR \
        --proj "${workdir}/val${f}/rep${r}" \
        --bfile "${bfile}" \
        --phe_col ${phe_col} \
        --fix "${fix_eff}" \
        --iter ${iter} \
        --burnin ${burnin} \
        --seed ${seed} \
        --gebv_out "EBV_bayesR.txt"
      gebvf=EBV_bayesR.txt
    fi

    if [[ ${method}  =~ "BLUP" ]]; then
      ## Parameter card
      cp ${workdir}/${DIR}.DIR ${workdir}/val${f}/rep${r}/${DIR}.DIR

      ## Run BLUP model
      if [[ ${dmu4} ]]; then
        [[ ! ${debug} ]] && job_pool_run run_dmu4 ${DIR} ${workdir}/val${f}/rep${r}
      else
        [[ ! ${debug} ]] && job_pool_run run_dmuai ${DIR} ${workdir}/val${f}/rep${r}
      fi
    fi
  done
done

## Wait for background processes to finish
job_pool_wait
## Release threads
job_pool_shutdown

####################  Accuracy Calculation  ##################
##############################################################
## Parameters
[[ ${rmNeg} ]] && rmNeg=" --rmNeg " ## Whether to remove negative values (abnormalities) from accuracy calculation results
if [[ ${tbv_col} == "same" ]]; then
  ## Use the original phenotype column in the file as the adjusted phenotype
  option=" --tbv_col ${phe_col}"
elif [[ ${tbv_col} ]]; then
  option=" --tbv_col ${tbv_col}"
fi
[[ ${tbvf} ]] && option="${option} --tbvf ${tbvf} "
if [[ ${method} =~ 'BLUP' ]]; then
  ## BLUP model
  option="${option} --add_sol ${add_sol} --dir_val ${workdir}/val#val#/rep#rep#/${DIR}"
  ebv_col=1
else
  ## BayesABLD model
  ebv_col=2
  option="${option} --ebvf ${workdir}/val#val#/rep#rep#/${gebvf}"
fi
## Calculate accuracy
$accur_cal \
  --ebv_col ${ebv_col} \
  --val_idf ${workdir}/val#val#/rep#rep#/val.id \
  --rep ${rep} \
  --fold ${fold} \
  ${rmNeg} \
  ${option} \
  --out ${workdir}/${out}

###################  Delete Intermediate Files  #####################
#####################################################################
## Parameter card
[[ -s ${workdir}/${DIR}.DIR ]] && rm ${workdir}/${DIR}.DIR

[[ $? -ne 0 ]] && echo "Accuracy calculation completed, file output to: ${out}"

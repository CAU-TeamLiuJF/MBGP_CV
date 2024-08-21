#!/usr/bin/bash

########################################################################################################################
## Version:   1.2.0
## Author:    Liweining liwn@cau.edu.cn
## Orcid:     0000-0002-0578-3812
## Institute: College of Animal Science and Technology, China Agricul-tural University, Haidian, 100193, Beijing, China
## Date:      2024-08-20
##
## Function:
## Obtain phenotypes corrected for all fixed and non genetic effects (Christensen et al., 2012)
##
## Usage: ./dmu_get_pheno_adj.sh --phef "pheno.txt" ...(Please refer to --help for detailed parameters)
## 
## Dependent software/environment:
##  1. R/4.1.0
##  2. plink/1.9
##  3. gmatrix
##  4. mbBayesABLD
##  5. Other R languages and Bash scripts
## 
## License:
##  This script is licensed under the GPL-3.0 License.
##  See https://www.gnu.org/licenses/gpl-3.0.en.html for details.
########################################################################################################################

################  Parameter processing  ##################
##########################################################
## NOTE: This requires GNU getopt.  On Mac OS X and FreeBSD, you have to install this
## Command-line parameters
TEMP=$(getopt -o 4vh\? --long append,ped_var,dmu4,help,phef:,bfile:,gmat:,gidf:,pedf:,all_eff:,ran_eff:,add_rf:,invA:,varf:,phereal:,miss:,num_int:,code:,DIR:,alpha:,out:,nchr:,intercept,debug \
    -n 'javawrap' -- "$@")
if [ $? != 0 ]; then
    echo "Open script $0 to view instructions"
    echo "Terminating..." >&2
    exit 1
fi
eval set -- "$TEMP"

## Parse parameters
while true; do
  case "$1" in
    --phef )     phef="$2";       shift 2 ;;  ## Phenotype file
    --bfile )    bfile="$2";      shift 2 ;;  ## plink binary file prefix
    --pedf )     pedf="$2";       shift 2 ;;  ## Pedigree file (GBLUP runs if not provided)
    --gmat )     gmat="$2";       shift 2 ;;  ## User-provided relationship matrix or inverse matrix (id id value)
    --gidf )     gidf="$2";       shift 2 ;;  ## Genotype individual ID, consistent with the individual ID in the user-specified G matrix file
    --all_eff )  all_eff="$2";    shift 2 ;;  ## Row 3 in $MODEL within DIR (all effects), first 3 positions are not needed, only the columns
                                              ## where all effects are located, e.g., "2 3 1"
    --ran_eff )  ran_eff="$2";    shift 2 ;;  ## Row 4 in $MODEL within DIR (random effects group), first position is not needed, only the group
                                              ## where all random effects are located, e.g., "1"
    --add_rf )   add_rf="$2";     shift 2 ;;  ## Group where additive effects are located
    --add_sol )  add_sol="$2";    shift 2 ;;  ## Additive effects identifier in the first column of the SOL file
    --invA )     invA="$2";       shift 2 ;;  ## Method for constructing the inverse of A (1/2/3/4/6), 1 considers inbreeding, 
                                              ## 2 does not consider inbreeding, others see DMU manual
    --varf )     varf="$2";       shift 2 ;;  ## Variance components file (use pedigree/SNP information for estimation if not provided)
    --phereal )  phereal="$2";    shift 2 ;;  ## Position of the real-number column in the phenotype file
    --miss )     miss="$2";       shift 2 ;;  ## Missing phenotype indicator
    --num_int )  num_int="$2";    shift 2 ;;  ## Number of integer columns
    --code )     code="$2";       shift 2 ;;  ## Path to the code
    --nchr )     nchr="$2";       shift 2 ;;  ## Number of chromosomes used in the species [30]
    --DIR )      DIR="$2";        shift 2 ;;  ## Parameter card file prefix
    --alpha )    alpha="$2";      shift 2 ;;  ## G matrix correction coefficient (whether to consider inbreeding)
    --out )      out="$2";        shift 2 ;;  ## Output corrected phenotype file name
    --append )   append=true;     shift   ;;  ## Append to the result file instead of overwriting
    --debug)     debug=true;      shift   ;;  ## Add a population mean column after the last integer column
    --intercept) mean=true;       shift   ;;  ## Add a population mean column after the last integer column
  -v | --ped_var )  ped_var=true; shift   ;;  ## Estimate variance components using pedigree information
  -4 | --dmu4 )     dmu4=true;    shift   ;;  ## Use DMU4 model to estimate breeding values
  -h | --help | -\? )  echo "Open script $0 to view instructions" && exit 1 ;;
  -- ) shift; break ;;
  * ) break ;;
  esac
done

## Directory of the script
if [[ ${code} ]]; then
  [[ ! -d ${code} ]] && echo "${code} not exists! " && exit 5
else
  script_path=$(dirname "$(readlink -f "$0")")
  code=$(dirname "$script_path")
fi

## Add the program path to the environment variable
export PATH=${code}/bin:$PATH

## Default parameter values
ran_eff=${ran_eff:="1"}        ## Default is only one effect for additive genetics
all_eff=${all_eff:="2 1"}      ## Default is only one fixed effect for the mean, located in the second column of the phenotype (all set to 1)
phereal=${phereal:="1"}        ## Phenotype column
add_rf=${add_rf:="1"}          ## Group where additive random effects are located
miss=${miss:="-99"}            ## Missing phenotype indicator
invA=${invA:="1"}              ## Method for constructing the inverse of A (whether to consider inbreeding)
alpha=${alpha:="0.05"}         ## G matrix correction coefficient (whether to consider inbreeding)
DIR=${DIR:="phe_adj"}          ## Parameter card file prefix
out=${out:="phe_adj.txt"}      ## Output file name
nchr=${nchr:="30"}             ## Number of chromosomes
append=${append:=false}        ## Append results to the existing file instead of overwriting
varf=${varf:=}                 ## Prevent vscode from reporting errors

## Avoid warnings when executing the R script ("ignoring environment value of R_HOME")
unset R_HOME

## Check if the necessary parameters are provided
if [[ ! -s ${phef} ]]; then
  echo "phenotype file ${phef} not found! "
  exit 1
elif [[ ! -s ${bfile}.fam && ! -s ${gidf} && ! -s ${pedf} ]]; then
  echo "plink file ${bfile}.fam, pedigree file ${pedf} or genotyped individuals id file ${gidf} not found! "
  exit 1
elif [[ -s ${gmat} && ! -s ${gidf} ]]; then
  echo "genotyped individuals id file ${gidf} not found! "
  exit 1
fi

## Main directory
workdir=$(pwd)

## Log directory
logp=${workdir}/log
mkdir -p ${logp}

## Dependency scripts
keep_phe_gid=${code}/R/keep_pheno_geno_individuals.R
miss_phe=${code}/R/pheno_miss_remove.R
phe_adj=${code}/R/adj_pheno_cal.R
func=${code}/shell/function.sh

## Load custom functions
[[ ! -s ${func} ]] && echo "Error: ${func} not found! " && exit 5
source ${func}

## Number of integer and real columns in the phenotype file
ncol=$(awk 'END{print NF}' ${phef})
for i in $(seq 1 ${ncol}); do
  dot=$(awk -vl=${i} '{print $l}' ${phef} | grep -c "\.")
  [[ ${dot} -gt 0 ]] && num_int=$((i - 1)) && break
done
num_real=$(($(awk 'END{print NF}' ${phef}) - num_int))

## Number of effects
nA=$(echo ${all_eff} | awk '{print NF}')
nR=$(echo ${ran_eff} | awk '{print NF}')


###################  Phenotype File Processing  #####################
#####################################################################
## Remove individuals with missing phenotypes
echo "removing individuals missing phenotypes in the phenotype file..."
[[ -s ${bfile}.fam ]] && option="--map ${bfile}.fam"
$miss_phe \
  --file ${phef} \
  --col $((phereal + num_int)) \
  --miss ${miss} \
  --missid ${workdir}/miss_phe.id \
  --out "${workdir}/pheno_adj.txt" \
  ${option}

phef=${workdir}/pheno_adj.txt
if [[ ${bfile} ]]; then
  ## Remove individuals with missing phenotypes from the genotype file
  if [[ -s ${workdir}/miss_phe.id ]]; then
    n_miss_phe=$(cat ${workdir}/miss_phe.id | wc -l)
    echo "remove ${n_miss_phe} individuals with the missing value in plink files"

    ## Rename and back up the original files
    plink \
      --bfile ${bfile} \
      --chr-set ${nchr} \
      --make-bed --out ${bfile}.org >${logp}/plink_rename.log

    ## Remove individuals with missing phenotypes
    plink \
      --bfile ${bfile}.org \
      --chr-set ${nchr} \
      --remove ${workdir}/miss_phe.id \
      --make-bed --out ${bfile} >>${logp}/plink_rm_miss_phe.log
  fi
fi

## Set the intercept column (set entire column to "1")
if [[ ${mean} ]]; then
  ## Effect parameters
  ((nA++))
  ((num_int++))
  all_eff="${num_int} ${all_eff}"

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
  echo "add populations mean column ${num_int} in the phenotype file."
  mv ${phef}.tmp ${phef}
fi


###################  dmu parameter card template  ####################
######################################################################
[[ -s ${DIR}.DIR ]] && echo "warn: ${DIR}.DIR will be overwrited! "
{
  echo "\$COMMENT"
  echo "creating phenotypes corrected for fixed effects and non-genetic random effects"
  echo "\$ANALYSE %ANALYSE%"
  echo "\$DATA  ASCII (${num_int}, ${num_real}, ${miss}) %phef%"
  echo -e "\$MODEL\n1\n0\n${phereal} 0 ${nA} ${all_eff}\n${nR} ${ran_eff}\n0\n0"
  echo -e "\$VAR_STR %VAR_STR%"
  echo -e "\$PRIOR %PRIOR%"
  echo -e "\$RESIDUALS ASCII"
  echo -e "\$SOLUTION"
} >${DIR}.DIR


##################  File validity check  ####################
#############################################################
## Phenotype file (check for non-numeric characters)
sed -i "s/na/${miss}/Ig" ${phef}
check_alphabet ${phef}
[[ -s ${pedf} ]] && check_alphabet ${pedf}
[[ -s ${gmat} ]] && check_alphabet ${gmat}
[[ -s ${bfile}.fam ]] && check_alphabet ${bfile}.fam 2


###################  Variance component estimation  #####################
#########################################################################
## Check whether variance components need to be estimated (using pedigree)
if [[ ${ped_var} ]]; then
  ## Check files
  [[ ! -s ${pedf} ]] && echo "${pedf} not found! " && exit 1

  ## Replace parameter card information
  sed 's#%ANALYSE%#1 1 0 0#g' ${DIR}.DIR >ped_var.DIR
  sed -i '/\$RESIDUALS.*/d' ped_var.DIR
  sed -i '/\$SOLUTION.*/d' ped_var.DIR
  sed -i "s#%VAR_STR%#${add_rf} PED ${invA} ASCII ${pedf}#g" ped_var.DIR
  ## Variance components
  if [[ -s ${varf} ]]; then
    sed -i "s#%PRIOR%#${varf}#g" ped_var.DIR
  else
    sed -i '/\$PRIOR.*/d' ped_var.DIR
  fi

  ## Estimate variance components using AIREML
  run_dmuai ped_var

  ## Variance component file
  sed -i 's#%PRIOR%#ped_var.PAROUT#g' ${DIR}.DIR
elif [[ -s ${varf} ]]; then
  sed -i 's#%PRIOR%#\${varf}#g' ${DIR}.DIR
else
  sed -i '/\$PRIOR.*/d' ${DIR}.DIR
fi


###################  Variance component structure  ####################
#######################################################################
if [[ -s ${bfile}.fam || -s ${gmat} ]] && [[ ! -s ${pedf} ]]; then
  ## GBLUP
  method=GBLUP
  gmat=${gmat:=full.agiv.id_fmt} ## 3-column format G inverse matrix file (id id value)
  sed -i "s#%VAR_STR%#${add_rf} GREL ASCII ${gmat}#g" ${DIR}.DIR

  ## Keep only phenotypes with genotyped individuals
  $keep_phe_gid \
    --famf "${bfile}.fam" \
    --phef ${phef} \
    --out ${phef}_gid
  phef=${phef}_gid
  add_sol=3
elif [[ ! -s ${bfile}.fam && -s ${pedf} ]]; then
  ## BLUP
  method=BLUP
  sed -i "s#%VAR_STR%#${add_rf} PED ${invA} ASCII ${pedf}#g" ${DIR}.DIR
  add_sol=4
else
  method=ssGBLUP
  gmat=${gmat:=full.agrm.id_fmt} ## 3-column format G matrix file (id id value)
  gidf=${gidf:=full.id}          ## 1-column id
  sed -i "s#%VAR_STR%#${add_rf} PGMIX ${invA} ASCII ${pedf} ${gidf} ${gmat} ${alpha} G-ADJUST#g" ${DIR}.DIR
  add_sol=4
fi


###################  Relationship matrix construction  #####################
############################################################################
## Construct G inverse matrix (also output genotype individual ids) using gmatrix software
if [[ ! -s ${gmat} && ${bfile} ]]; then
  [[ ${method} == "GBLUP" ]] && inv=" --inv" || inv=""
  echo "Read the plink bed file and Calculate the additive G matrix..."
  
  ## Generate G matrix using gmatrix software
  gmatrix \
    --bfile ${bfile} \
    --grm agrm \
    --out full ${inv} >${logp}/gmatrix.log

  if [[ $? -ne 0 ]]; then
    echo "G matrix calculation error! "
    exit 1
  else
    echo "G matrix created."
  fi
fi


#################  Breeding value estimation  ###################
#################################################################
echo "estimating breeding values using the ${method} model..."
sed -i "s#%phef%#${phef}#g" ${DIR}.DIR
if [[ ! ${debug} ]]; then
  if [[ ${dmu4} ]]; then
    [[ ! -s ${varf} ]] && echo "${varf} not found! " && exit 1
    ## Run dmu4
    sed -i 's#%ANALYSE%#11 9 0 0#g' ${DIR}.DIR
    run_dmu4 ${DIR}
    [[ $? -ne 0 ]] && echo "error in dmu4! " && exit 1
  else
    sed -i 's#%ANALYSE%#1 31 0 0#g' ${DIR}.DIR
    run_dmuai ${DIR}
    [[ $? -ne 0 ]] && echo "error in dmuai! " && exit 1
  fi
fi


#################  Calculate adjusted phenotype  #################
##################################################################
echo "calculating adjusted phenotype..."
$phe_adj \
  --DIR ${DIR} \
  --phe ${phef} \
  --add_sol ${add_sol} \
  --out ${out} \
  --append ${append}

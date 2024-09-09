#!/bin/bash

########################################################################################################################
## Version: 1.3.0
## Author:    Liweining liwn@jaas.ac.cn
## Orcid:     0000-0002-0578-3812
## Institute: College of Animal Science and Technology, China Agricul-tural University, Haidian, 100193, Beijing, China
## Date:      2024-08-20
##
## Function：
## Run BayesR
##
##
## Usage:
## ./BayesR.sh --proj "/path/to/project" ...(Please refer to --help for detailed parameters)
##
##
## License:
##  This script is licensed under the GPL-3.0 License.
##  See https://www.gnu.org/licenses/gpl-3.0.en.html for details.
########################################################################################################################


###################  Parameter processing  #####################
####################################################
## NOTE: This requires GNU getopt.  On Mac OS X and FreeBSD, you have to install this
## Parameters
TEMP=$(getopt -o h --long proj:,bfile:,fix:,phef:,phe_col:,code:,iter:,burnin:,nthread:,seed:,gebv_out:,help \
              -n 'javawrap' -- "$@")
if [ $? != 0 ] ; then echo "Terminating..." >&2 ; exit 1 ; fi
eval set -- "$TEMP"
## Parse parameters
while true; do
  case "$1" in
    --proj )     proj="$2";     shift 2 ;; ## Project directory [required]
    --bfile )    bfile="$2";    shift 2 ;; ## plink genotype file [required]
    --fix )      fix="$2";      shift 2 ;; ## fix effects [2]
    --phef )     phef="$2";     shift 2 ;; ## genotype file [pheno.txt]
    --phe_col )  phe_col="$2";  shift 2 ;; ## phenotype column [last column]
    --code )     code="$2";     shift 2 ;; ## Scripts directory, e.g., /BIGDATA2/cau_jfliu_2/liwn/code [NULL]
    --miss )     miss="$2";     shift 2 ;; ## Missing phenotype identifier [-99]
    --iter )     iter="$2";     shift 2 ;; ## MCMC iteration [30000]
    --burnin )   burnin="$2";   shift 2 ;; ## MCMC burn-in [20000]
    --nthread )  nthread="$2";  shift 2 ;; ## Number of cores [1]
    --seed )     seed="$2";     shift 2 ;; ## Random number seed [8123]
    --gebv_out ) gebv_out="$2"; shift 2 ;; ## GEBV output file name [bayesR_gebv.txt]
  -h | --help )   grep ";; ##" $0 | grep -v help && exit 1 ;;
  -- ) shift; break ;;
  * ) shift; break ;;
  esac
done

## Check the scripts directory
if [[ ${code} ]]; then
  [[ ! -d ${code} ]] && echo "${code} not exists! " && exit 5
else
  script_path=$(dirname "$(readlink -f "$0")")
  code=$(dirname "$script_path")
fi

## Check if required parameters are provided
if [[ ! -d ${proj} ]]; then
  echo "${proj} not found! "
  exit 1
else
  cd ${proj} ||exit 1
fi

## default parameters
phef=${phef:=pheno.txt}
miss=${miss:=-99}
fix=${fix:=2}
nthread=${nthread:=1}
seed=${seed:=8123}
iter=${iter:=300}
burnin=${burnin:=200}
gebv_out=${gebv_out:=bayesR_gebv.txt}

## scripts
func=${code}/shell/function.sh
fix_file=${code}/R/bayesR_file.R ## genotype file and fixed Structural matrix

## Load self-defined functions
[[ ! -s ${func} ]] && echo "Error: ${func} not found! " && exit 5
source ${func}

## Add the program path to the environment variable
export PATH=${code}/bin:$PATH

## Check if the required program can be retrieved and executed in the environment variables
check_command bayesRv2

## Check if the required script files exist and have execution permissions
check_command $fix_file

## phenotype column
if [[ ! ${phe_col} ]]; then
  phe_col=$(awk 'NR==1{print NF}' ${phef})
fi

## get structure matrix of fixed effect
$fix_file \
  --bfile "${bfile}" \
  --phe "${phef}" \
  --miss "${miss}" \
  --fix "${fix}" \
  --phe_col "${phe_col}" \
  --out_dir "${proj}" \
  --gt_out "genotype" \
  --fix_out fix_BayesR.txt
[[ -s fix_BayesR.txt ]] && fix_effect=" -covar fix_BayesR.txt "

## run bayesR
echo "Runing bayesR..."
bayesRv2 \
  -bfile "genotype" \
  -numit ${iter} \
  -burnin ${burnin} \
  -out bayesR \
  -seed ${seed} \
  -nthreads ${nthread} \
  ${fix_effect}

## get gebv
bayesRv2 \
  -bfile "genotype" \
  -predict \
  -model bayesR.model \
  -freq bayesR.frq \
  -param bayesR.param \
  -out bayesR_val

## paste
awk '{print $2}' genotype.fam >id_tmp.txt
paste -d" " id_tmp.txt bayesR_val.gv >${gebv_out}
rm id_tmp.txt
[[ $? -eq 0 ]] && echo "GEBV result has been output to: ${gebv_out}"

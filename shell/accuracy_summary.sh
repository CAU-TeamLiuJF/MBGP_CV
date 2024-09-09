#!/bin/bash

########################################################################################################################
## Version: 1.3.0
## Author: Weining Li liwn@jaas.ac.cn
## Date: 2023-07-05
## 
## Statistical results of cross-validation accuracy under different conditions
## 
## Usage: ./varcomp_summary.sh --help
## 
## License:
##  This script is licensed under the GPL-3.0 License.
##  See https://www.gnu.org/licenses/gpl-3.0.en.html for details.
########################################################################################################################


###################  Parameter processing  #####################
################################################################
## NOTE: This requires GNU getopt.  On Mac OS X and FreeBSD, you have to install this
## Parameters
TEMP=$(getopt -o h --long proj:,breeds:,rep:,dist:,cor:,traits:,bin:,dirPre:,out:,help \
              -n 'javawrap' -- "$@")
if [ $? != 0 ] ; then echo "Terminating..." >&2 ; exit 1 ; fi
eval set -- "$TEMP"
## Parse parameters
while true; do
  case "$1" in
    --proj )     proj="$2";     shift 2 ;; ## Project directory [required]
    --breeds )   breeds="$2";   shift 2 ;; ## Label of pulation/breed, e.g., 'YY DD' [required]
    --traits )   traits="$2";   shift 2 ;; ## Trait names, e.g., "DF DPM" [""]
    --rep )      rep="$2";      shift 2 ;; ## The number of repeats [optional]
    --dist )     dist="$2";     shift 2 ;; ## Distribution of additive genetic correlation [""]
    --cor )      cor="$2";      shift 2 ;; ## Additive genetic correlation size [""]
    --dirPre )   dirPre="$2";   shift 2 ;; ## Extra prefix for EBV folder [""]
    --bin )      bins="$2";     shift 2 ;; ## Block partitioning methods for multi-breed genomic prediction, fix/frq/ld/lava/cubic ["fix lava cubic"]
    --out )      out="$2";      shift 2 ;; ## Output filename for accuracy [accuracy_$date.txt]
  -h | --help)    grep ";; ##" $0 | grep -v help && exit 1 ;;
  -- ) shift; break ;;
  * ) shift; break ;;
  esac
done

## Check if required parameters are provided
if [[ ! -d ${proj} ]]; then
  echo "${proj} not found! "
  exit 1
elif [[ ! ${breeds} ]]; then
  echo "para --breeds is reduired! "
  exit 1
fi

## Date
today=$(date +%Y%m%d)

## Suppress warnings when executing R scripts("ignoring environment value of R_HOME")
unset R_HOME

## Default parameters
out=${out:=${proj}/accuracy_${today}.txt}
bins=${bins:="fix lava cubic"}
dirPre=${dirPre:=""}
traits=${traits:="/"}
rep=${rep:="/"}
dist=${dist:="/"}
cor=${cor:="/"}

## Parse parameters
read -ra breeds_array <<<"$breeds"
read -ra bins_array <<<"$bins"
read -ra traits_array <<<"$traits"
read -ra reps_array <<<"$rep"
read -ra dists_array <<<"$dist"
read -ra cors_array <<<"$cor"

############## Accuracy results statistics (different combinations) ##########
##############################################################################
echo "Simulation_repeats correlation_distribution correlation prediction_type model reference_population breed trait number_of_repeats number_of_cross-validation_folds accuracy unbiasedness rank_correlation validation_population_size" >${out}
for t in "${traits_array[@]}"; do
  for r in "${reps_array[@]}"; do
    for d in "${dists_array[@]}"; do
      for c in "${cors_array[@]}"; do
        for b in "${breeds_array[@]}"; do
          path=${proj}/${t}

          ## Path settings for simulated scenarios
          [[ ${r} != "/" ]] && path=${path}/rep${r}
          [[ ${d} != "/" ]] && path=${path}/${d}
          [[ ${c} != "/" ]] && path=${path}/cor${c}

          ## Handle cases with different types of "/" in the path
          path=$(echo "$path" | sed 's#/\{2,\}#/#g; s#/$##')

          ## Check if the directory exists
          [[ ! -d ${path} ]] && continue

          ## Cross-validation parameters
          [[ ! -d ${path}/${b}/val1 ]] && continue
          rep=$(find ${path}/${b}/val1 -name "rep*" -type d | wc -l)
          fold=$(find ${path}/${b}/val* -name "rep1" -type d | wc -l)

          ## within population
          wf=${path}/${b}/accur_GBLUP.txt
          bf=${path}/${b}/accur_mbBayesAB.txt
          # if [[ -s ${wf} ]]; then
          #   sed '1d' ${wf} | awk '{print "within","'${b}'","'${b}'","'${t}'", $0}' >>${out}
          # # else
          #   # echo "${accf} not found! "
          # fi

          ## single-trait joint prediction models
          accf=$(find ${path} -path "*/singl*/*" -name "accur_*${b}.txt" 2>/dev/null)
          if [[ ${accf} ]]; then
            # accf=${path}/blend/accur_GBLUP_${b}.txt
            for f in ${accf}; do
              if [[ -s ${f} ]]; then
                ## reference populations
                ref=$(dirname ${f})
                ref=$(basename ${ref})
                ref=${ref/single_/}

                ## model
                model=$(basename ${f})
                model=${model/accur_/}
                model=${model/_${b}.txt/}

                ## print to file
                {
                  awk '{print "'${r}'","'${d}'","'${c}'","single","'${model}'","'${ref}'","'${b}'","'${t}'",$0}' ${f}
                  [[ -s ${wf} ]] && \
                    awk '{print "'${r}'","'${d}'","'${c}'","within","w-GBLUP","'${ref}'","'${b}'","'${t}'", $0}' ${wf}
                  [[ -s ${bf} ]] && \
                    awk '{print "'${r}'","'${d}'","'${c}'","within","w-mbBayesAB","'${ref}'","'${b}'","'${t}'", $0}' ${bf}
                } >>${out}
              # else
                # echo "${accf} not found! "
              fi
            done
          fi

          ## multi-trait joint prediction models
          accf=$(find ${path} -path "*/mult*/*" -name "accur_*${b}.txt" 2>/dev/null)
          # accf=$(find ${path}/unio* -name "accur_*${b}*txt" 2>/dev/null)
          if [[ ${accf} ]]; then
            for f in ${accf}; do
              if [[ -s ${f} ]]; then
                ## reference populations
                ref=$(dirname ${f})
                ref=$(basename ${ref})
                ref=${ref/multi_/}

                ## model
                model=$(basename ${f})
                model=${model/accur_/}
                model=${model/_${b}.txt/}

                ## print to file
                {
                  awk '{print "'${r}'","'${d}'","'${c}'","multi","'${model}'","'${ref}'","'${b}'","'${t}'",$0}' ${f}
                  [[ -s ${wf} ]] && \
                    awk '{print "'${r}'","'${d}'","'${c}'","within","w-GBLUP","'${ref}'","'${b}'","'${t}'", $0}' ${wf}
                  [[ -s ${bf} ]] && \
                    awk '{print "'${r}'","'${d}'","'${c}'","within","w-mbBayesAB","'${ref}'","'${b}'","'${t}'", $0}' ${bf}
                } >>${out}
              # else
              #   echo "${accf} not found! "
              fi
            done
          fi
        done
      done
    done
  done
done

## Remove the first lines
sed -i '/rep/d' ${out}
sed -i '/mean/d' ${out}
## Replace spaces
sed -i 's/ /\t/g' ${out}

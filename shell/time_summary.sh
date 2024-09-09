#!/bin/bash

########################################################################################################################
## Version: 1.0.0
## Author:    Liweining liwn@jaas.ac.cn
## Orcid:     0000-0002-0578-3812
## Institute: College of Animal Science and Technology, China Agricul-tural University, Haidian, 100193, Beijing, China
## Date:      2024-09-02
##
## Function：
## Statistical variance component results
##
##
## Usage:
## ./varcomp_summary.sh --proj "/path/to/project" ...(Please refer to --help for detailed parameters)
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
TEMP=$(getopt -o h --long code:,proj:,breeds:,rep:,dist:,cor:,traits:,bin:,dirPre:,out:,help \
              -n 'javawrap' -- "$@")
if [ $? != 0 ] ; then echo "Terminating..." >&2 ; exit 1 ; fi
eval set -- "$TEMP"
##  parse parameter
while true; do
  case "$1" in
    --proj )     proj="$2";     shift 2 ;; ## Project directory [required parameter]
    --breeds )   breeds="$2";   shift 2 ;; ## Breed identifiers, e.g., 'YY DD' [required parameter]
    --traits )   traits="$2";   shift 2 ;; ## Trait names, e.g., "DF DPM" ["/"]
    --rep )      rep="$2";      shift 2 ;; ## Replication number ["/"]
    --dist )     dist="$2";     shift 2 ;; ## Distribution of additive genetic relationships ["/"]
    --cor )      cor="$2";      shift 2 ;; ## Magnitude of additive genetic correlation ["/"]
    --dirPre )   dirPre="$2";   shift 2 ;; ## Prefix added to the EBV output folder [""]
    --bin )      bins="$2";     shift 2 ;; ## Method of region division for multi-breed prediction, fix/frq/ld/lava/cubic ["fix lava cubic"]
    --code )     code="$2";     shift 2 ;; ## Directory where the script file is located, e.g., /BIGDATA2/cau_jfliu_2/liwn/code [NULL]
    --out )      out="$2";      shift 2 ;; ## Output filename for accuracy [accuracy_$date.txt]
  -h | --help)    grep ";; ##" $0 | grep -v help && exit 1 ;;
  -- ) shift; break ;;
  * ) shift; break ;;
  esac
done

## Check the required parameters
if [[ ! -d ${proj} ]]; then
  echo "${proj} not found! "
  exit 1
elif [[ ! ${breeds} ]]; then
  echo "para --breeds is reduired! "
  exit 1
fi

## Date
today=$(date +%Y%m%d)

## Default parameters
out=${out:=${proj}/time_${today}.txt}
bins=${bins:="fix lava cubic"}
dirPre=${dirPre:=""}
traits=${traits:="/"}
rep=${rep:="/"}
dist=${dist:="/"}
cor=${cor:="/"}

## Avoid warnings when running R scripts ("ignoring environment value of R_HOME")
unset R_HOME

## DMU parameter card name
DIR_full=phe_adj_BLUP
DIR_within=within
DIR_single=single
DIR_multi=multi
bayesR_prefix=bayesR

##  parse parameter
read -ra breeds_array <<<"$breeds"
read -ra bins_array <<<"$bins"
read -ra traits_array <<<"$traits"
read -ra reps_array <<<"$rep"
read -ra dists_array <<<"$dist"
read -ra cors_array <<<"$cor"

##############  Variance component sorting  ##############
##########################################################
for p in "${traits_array[@]}"; do # p=${traits_array[0]};re=${reps_array[0]};d=${dists_array[0]};c=${cors_array[0]};b=${breeds_array[0]}
  for re in "${reps_array[@]}"; do
    for d in "${dists_array[@]}"; do
      for c in "${cors_array[@]}"; do
        path=${proj}/${p}

        ## Set the path under the simulation scenario
        [[ ${re} != "/" ]] && path=${path}/rep${re}
        [[ ${d} != "/" ]] && path=${path}/${d}
        [[ ${c} != "/" ]] && path=${path}/cor${c}

        ## Handle cases where there are multiple slashes in the path
        path=$(echo "$path" | sed 's#/\{2,\}#/#g; s#/$##')

        ## Check if the directory exists
        [[ ! -d ${path} ]] && continue

        ## Time consumption within breed
        for b in "${breeds_array[@]}"; do
          ## Cross-validation parameters
          rep=$(find ${path}/${b}/val1 -name "rep*" -type d | wc -l)
          fold=$(find ${path}/${b}/val* -name "rep1" -type d | wc -l)

          ## Time consumption for each subset
          for r in $(seq 1 ${rep}); do # r=1;f=1
            for f in $(seq 1 ${fold}); do
                ## within prediction --- GBLUP
                lst=${path}/${b}/val${f}/rep${r}/${DIR_within}.lst
                if [[ -s ${lst} ]]; then
                  start=$(grep 'start:' ${lst} | awk '{print $2,$3}')
                  end=$(grep 'end:' ${lst} | awk '{print $2,$3}')

                  ## Convert time to seconds
                  start_seconds=$(date -d "$start" +%s)
                  end_seconds=$(date -d "$end" +%s)

                  ## Calculate time interval (seconds)
                  interval_sec=$((end_seconds - start_seconds))
                  time_wG=$(echo "scale=1; $interval_sec / 60" | bc)

                  echo "${re} ${d} ${c} w-GBLUP ${b} ${p} ${r} ${f} ${time_wG}" >>${out}
                fi

                ## within prediction --- Bayes
                mapfile -t logf < <(find ${path}/${b}/val${f}/rep${r} -name "*gibs*.log" 2>/dev/null)
                [[ ! -s ${logf[0]} ]] && continue

                start=$(grep 'Run started' ${logf[0]} | awk '{print $5,$6,$7}')
                end=$(grep 'Run ended' ${logf[0]} | awk '{print $5,$6,$7}')

                ## Check if the program has ended normally
                [[ ! ${start} || ! ${end} ]] && \
                  echo "warn: ${logf[0]} Running abnormally, please check! " && \
                  continue

                ## Convert time to seconds
                start_seconds=$(date -d "$start" +%s)
                end_seconds=$(date -d "$end" +%s)

                ## Calculate time interval (seconds)
                interval_sec=$((end_seconds - start_seconds))
                time_wR=$(echo "scale=1; $interval_sec / 60" | bc)

                echo "${re} ${d} ${c} w-fix ${b} ${p} ${r} ${f} ${time_wR}" >>${out}
            done
          done
        done

        ## Breed combinations
        mapfile -t single < <(find ${path} -name "singl*" -type d 2>/dev/null)

        for t in "${single[@]}"; do
          dir=$(basename ${t})
          types=${dir/single_/}
          IFS='_' read -r -a breeds_sub <<<"$types"
          nb=${#breeds_sub[@]}

          ## Types
          if [[ ${types} == "single" ]]; then
            types=""
            comb=${breeds// /_}
          else
            types="_${types}"
            comb=$(IFS="_"; echo "${breeds_sub[*]}")
          fi

          ## each model
          for r in $(seq 1 ${rep}); do
            for f in $(seq 1 ${fold}); do
              ## Single-trait GBLUP joint prediction
              lst=${path}/single${types}/val${f}/rep${r}/${DIR_single}.lst
              if [[ -s ${lst} ]]; then
                  start=$(grep 'start:' ${lst} | awk '{print $2,$3}')
                  end=$(grep 'end:' ${lst} | awk '{print $2,$3}')

                  ## Check if the program has ended normally
                  if [[ ! ${start} || ! ${end} ]]; then
                    echo "warn: ${lst} Running abnormally, please check! "
                  else
                    ## Convert time to seconds
                    start_seconds=$(date -d "$start" +%s)
                    end_seconds=$(date -d "$end" +%s)

                    ## Calculate time interval (seconds)
                    interval_sec=$((end_seconds - start_seconds))
                    time_sG=$(echo "scale=1; $interval_sec / 60" | bc)
                    echo "${re} ${d} ${c} m-STGBUP ${comb} ${p} ${r} ${f} ${time_sG}" >>${out}
                  fi
              else
                  echo "${lst} not found! "
              fi

              ## Single-trait bayesR joint prediction
              logfBR=${path}/single${types}/val${f}/rep${r}/${bayesR_prefix}.log
              if [[ -s ${logfBR} ]]; then
                  start=$(grep 'Run started' ${logfBR} | awk '{print $4,$5}')
                  end=$(grep 'Run ended' ${logfBR} | awk '{print $4,$5}')

                  ## Check if the program has ended normally
                  if [[ ! ${start} || ! ${end} ]]; then
                    echo "warn: ${logfBR} Running abnormally, please check! "
                  else
                    ## Convert time to seconds
                    start_seconds=$(date -d "$start" +%s)
                    end_seconds=$(date -d "$end" +%s)

                    ## Calculate time interval (seconds)
                    interval_sec=$((end_seconds - start_seconds))
                    time_sR=$(echo "scale=1; $interval_sec / 60" | bc)

                    echo "${re} ${d} ${c} m-bayesR ${comb} ${p} ${r} ${f} ${time_sR}" >>${out}
                    fi
              else
                  echo "${logfBR} not found! "
              fi

              ## multi-breed GBLUP joint prediction
              lst=${path}/multi${types}/val${f}/rep${r}/${DIR_multi}.lst
              if [[ -s ${lst} ]]; then
                  start=$(grep 'start:' ${lst} | awk '{print $2,$3}')
                  end=$(grep 'end:' ${lst} | awk '{print $2,$3}')

                  ## Convert time to seconds
                  start_seconds=$(date -d "$start" +%s)
                  end_seconds=$(date -d "$end" +%s)

                  ## Check if the program has ended normally
                  if [[ ! ${start} || ! ${end} ]]; then
                    echo "warn: ${lst} Running abnormally, please check! "
                  else
                    ## Calculate time interval (seconds)
                    interval_sec=$((end_seconds - start_seconds))
                    time_mG=$(echo "scale=1; $interval_sec / 60" | bc)
                    echo "${re} ${d} ${c} m-MTGBUP ${comb} ${p} ${r} ${f} ${time_mG}" >>${out}
                  fi
              else
                  echo "${lst} not found! "
              fi

              for bin in "${bins_array[@]}"; do
                multi_path=${path}/multi${types}/val${f}/rep${r}

                [[ ! -d ${multi_path} ]] && continue

                mapfile -t "logf" < <(find ${multi_path} -name "${bin}*gibs*log" 2>/dev/null)

                [[ ${#logf[@]} -gt 1 ]] && echo "warn: multiple ${bin} log file not found in ${multi_path}"
                [[ ! ${logf[0]} ]] && echo "warn: ${bin} log file not found in ${multi_path}" && continue

                start=$(grep 'Run started' ${logf[0]} | awk '{print $5,$6,$7}')
                end=$(grep 'Run ended' ${logf[0]} | awk '{print $5,$6,$7}')

                ## Check if the program has ended normally
                [[ ! ${start} || ! ${end} ]] && \
                  echo "warn: ${logf[0]} Running abnormally, please check! " && \
                  continue

                ## Convert time to seconds
                start_seconds=$(date -d "$start" +%s)
                end_seconds=$(date -d "$end" +%s)

                ## Calculate time interval (seconds)
                interval_sec=$((end_seconds - start_seconds))
                time_mR=$(echo "scale=1; $interval_sec / 60" | bc)

                echo "${re} ${d} ${c} mbBayesAB-${bin} ${comb} ${p} ${r} ${f} ${time_mR}" >>${out}
              done
            done
          done
        done
      done
    done
  done
done

## Remove rows with NA values
sed -i '/ $/d' ${out}
sed -i 's/ /\t/g' ${out}

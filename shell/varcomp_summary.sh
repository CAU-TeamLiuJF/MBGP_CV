#!/bin/bash

########################################################################################################################
## Version: 1.3.0
## Author:    Liweining liwn@cau.edu.cn
## Orcid:     0000-0002-0578-3812
## Institute: College of Animal Science and Technology, China Agricul-tural University, Haidian, 100193, Beijing, China
## Date:      2024-08-20
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
out=${out:=${proj}/varcomp_${today}.txt}
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
for p in "${traits_array[@]}"; do
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

        ## Heritability within breed
        for b in "${breeds_array[@]}"; do
          ## Cross-validation parameters
          rep=$(find ${path}/${b}/val1 -name "rep*" -type d | wc -l)
          fold=$(find ${path}/${b}/val* -name "rep1" -type d | wc -l)

          ## Variance components for each subset
          for r in $(seq 1 ${rep}); do
            for f in $(seq 1 ${fold}); do
                ## within prediction --- GBLUP
                lst=${path}/${b}/val${f}/rep${r}/${DIR_within}.lst
                if [[ -s ${lst} ]]; then
                  h2_within=$(grep -A2 'Trait  correlation' ${lst} | tail -n 1 | awk '{print $2}')
                  echo "${re} ${d} ${c} w-GBLUP ${b} ${b} ${p} ${r} ${f} h2 ${h2_within}" >>${out}
                fi

                ## within prediction --- Bayes
                ebvf=${path}/${b}/val${f}/rep${r}/EBV_fix_y1.txt
                varf=${path}/${b}/val${f}/rep${r}/var_fix.txt
                [[ ! -s ${ebvf} || ! -s ${varf} ]] && continue
                varg=$(sed '1d' ${ebvf} | awk '{sum+=$2; sumsq+=($2)^2} END {print (sumsq/NR-(sum/NR)^2)}')
                vare=$(tail -n 1 ${varf})
                h2_Bayes=$(echo "scale=4; $varg / ($varg + $vare)" | bc)
                echo "${re} ${d} ${c} w-fix ${b} ${b} ${p} ${r} ${f} h2 ${h2_Bayes}" >>${out}
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

          ## Each breed in the combination
          for i in "${!breeds_sub[@]}"; do
            ## Subsets
            for r in $(seq 1 ${rep}); do
              for f in $(seq 1 ${fold}); do
                ## Single-trait GBLUP joint prediction
                lst=${path}/single${types}/val${f}/rep${r}/${DIR_single}.lst
                if [[ -s ${lst} ]]; then
                  h2_sG=$(grep -A2 'Trait  correlation' ${lst} | tail -n 1 | awk '{print $2}')
                else
                  h2_sG=""
                  echo "${lst} not found! "
                fi

                ## Single-trait bayesR joint prediction
                model=${path}/single${types}/val${f}/rep${r}/${bayesR_prefix}.model
                if [[ -s ${model} ]]; then
                  Va=$(grep 'Va' ${model} | awk '{printf("%f",$2)}')
                  Ve=$(grep 'Ve' ${model} | awk '{printf("%f",$2)}')
                  h2_sR=$(echo "scale=4; $Va / ($Va + $Ve)" | bc)
                else
                  h2_sR=""
                  echo "${model} not found! "
                fi

                ## multi-breed GBLUP joint prediction
                lst=${path}/multi${types}/val${f}/rep${r}/${DIR_multi}.lst
                if [[ -s ${lst} ]]; then
                  message="Correlation matrix for random"
                  h2_mG=$(grep -A $((i + 2)) 'Trait  correlation' ${lst} | tail -n 1 | awk '{print $2}')

                  for j in $(seq $((i + 1)) $((nb - 1))); do
                    rg_union=$(grep -A $((j + 2)) "${message}" ${lst} | tail -n 1 | awk -v col=$((i + 2)) '{print $col}')
                    echo "${re} ${d} ${c} m-MTGBUP ${comb} ${breeds_sub[i]}_${breeds_sub[j]} ${p} ${r} ${f} rg ${rg_union}" >>${out}
                  done
                else
                  h2_mG=""
                  echo "${lst} not found! "
                fi

                for bin in "${bins_array[@]}"; do
                  # for dirPre in /; do
                    multi_path=${path}/multi${types}/val${f}/rep${r}

                    ebvfi=$(find ${multi_path} -name "EBV_${bin}_y$((i + 1)).txt" 2>/dev/null)
                    varfi=$(find ${multi_path} -name "var_${bin}*.txt" 2>/dev/null)
                    [[ ! -s ${ebvfi} || ! -s ${varfi} ]] && continue

                    ## Heritability
                    varg=$(sed '1d' ${ebvfi} | awk '{sum+=$2; sumsq+=($2)^2} END {print (sumsq/NR-(sum/NR)^2)}')
                    vare=$(tail -n 1 ${varfi} | awk -v col=$((i * nb + i + 1)) '{printf("%f",$col)}')
                    h2_multi=$(echo "scale=4; $varg / ($varg + $vare)" | bc)
                    echo "${re} ${d} ${c} ${bin} ${comb} ${breeds_sub[i]} ${p} ${r} ${f} h2 ${h2_multi}" >>${out}

                    ## Genetic correlation
                    for j in $(seq $((i + 1)) $((nb - 1))); do
                      ebvfj=$(find ${multi_path} -name "EBV_${bin}_y$((j + 1)).txt" 2>/dev/null)
                      [[ ! -s ${ebvfj} ]] && continue

                      # Calculate values in column 2 of file A minus the mean
                      meanA=$(awk 'NR>1{sum+=$2}END{print sum/(NR-1)}' ${ebvfi})
                      awk -v meanA="$meanA" 'NR>1{print $2-meanA}' ${ebvfi} >A.tmp
                      # Calculate values in column 2 of file B minus the mean
                      meanB=$(awk 'NR>1{sum+=$2}END{print sum/(NR-1)}' ${ebvfj})
                      awk -v meanB="$meanB" 'NR>1{print $2-meanB}' ${ebvfj} >B.tmp
                      # Calculate covariance
                      cov=$(paste A.tmp B.tmp | awk '{sum+=($1*$2)}END{print sum/(NR - 1)}')
                      var2=$(sed '1d' ${ebvfj} | awk '{sum+=$2; sumsq+=($2)^2} END {print (sumsq/NR-(sum/NR)^2)}')

                      ## Calculate genetic correlation
                      if [[ $(echo "$varg < 0" | bc -l) -eq 1 || $(echo "$var2 < 0" | bc -l) -eq 1 ]]; then
                        echo "varg=$varg var2=$var2"
                        rg_multi=0
                      else
                        rg_multi=$(echo "scale=4; $cov / sqrt($varg * $var2)" | bc | xargs printf "%.4f")
                      fi

                      echo "${re} ${d} ${c} ${bin} ${comb} ${breeds_sub[i]}_${breeds_sub[j]} ${p} ${r} ${f} rg ${rg_multi}" >>${out}
                    done
                  # done
                done

                ## Write to file
                {
                  echo "${re} ${d} ${c} bayesR ${comb} ${breeds_sub[i]} ${p} ${r} ${f} h2 ${h2_sR}"
                  echo "${re} ${d} ${c} m-STGBUP ${comb} ${breeds_sub[i]} ${p} ${r} ${f} h2 ${h2_sG}"
                  echo "${re} ${d} ${c} m-MTGBUP ${comb} ${breeds_sub[i]} ${p} ${r} ${f} h2 ${h2_mG}"
                } >>${out}
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

[[ -s A.tmp ]] && rm A.tmp B.tmp

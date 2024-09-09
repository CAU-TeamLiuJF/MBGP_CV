#!/usr/bin/bash

########################################################################################################################
## Version: 1.3.0
## Author:    Liweining liwn@jaas.ac.cn
## Orcid:     0000-0002-0578-3812
## Institute: College of Animal Science and Technology, China Agricul-tural University, Haidian, 100193, Beijing, China
## Date:      2024-08-20
##
## Function：
## Generate the parameter file *. prm required for QMSim software based on the provided parameters
##
##
## Usage:
## ./QMSim_genome_process.sh --proj "/path/to/project" ...(Please refer to --help for detailed parameters)
##
##
## License:
##  This script is licensed under the GPL-3.0 License.
##  See https://www.gnu.org/licenses/gpl-3.0.en.html for details.
########################################################################################################################


###################  Parameter processing  #####################
################################################################
## NOTE: This requires GNU getopt.  On Mac OS X and FreeBSD, you have to install this
## Parse arguments
TEMP=$(getopt -o h --long proj:,breeds:,rep:,nsnp:,geno_gen:,geno_sel:,binDiv:,binThr:,maf:,nginds:,last_litters:,last_females:,code:,out:,help \
              -n 'javawrap' -- "$@")
if [ $? != 0 ] ; then echo "Terminating..." >&2 ; exit 1 ; fi
eval set -- "$TEMP"
##  parse parameter
while true; do
  case "$1" in
    --proj )         proj="$2";     shift 2 ;; ## Project directory [required parameter]
    --breeds )       breeds="$2";   shift 2 ;; ## Population/Breed identifiers, e.g., 'A B C' ["A B"]
    --rep )          rep="$2";      shift 2 ;; ## Which replicate of QMSim output to process, e.g., 1 for lm_mrk_001.txt [1]
    --nsnp )         nsnp="$2";     shift 2 ;; ## Number of SNPs to select [50000]
    --geno_gen )     geno_gen="$2"; shift 2 ;; ## Generation(s) to output genotyped individuals [8-10]
    --geno_sel )     geno_sel="$2"; shift 2 ;; ## From which generations to select genotyped individuals [geno_gen]
    --binDiv )       binDiv="$2";   shift 2 ;; ## Criterion for dividing regions when sampling SNPs, pos/frq [pos]
    --binThr )       binThr="$2";   shift 2 ;; ## length for dividing regions when sampling SNPs, in physical position cM or frequency step [10]
    --maf )          maf="$2";      shift 2 ;; ## Minimum allowable allele frequency when sampling SNPs [0.01]
    --nginds )       nginds="$2";   shift 2 ;; ## Number of genotyped individuals to select per breed ["600 600 ..."]
    --nchr )         nchr="$2" ;    shift 2 ;; ## Number of chromosomes [30]
    --last_litters ) litters="$2";  shift 2 ;; ## Number of individuals per litter during the LD stabilization phase for each breed ["10 10 ..."]
    --last_females ) females="$2";  shift 2 ;; ## Number of females in the population in the last phase for each breed ["500 500 ..."]
    --code )         code="$2";     shift 2 ;; ## Directory where script files are located, e.g., /BIGDATA2/cau_jfliu_2/liwn/code [NULL]
    --out )          out="$2";      shift 2 ;; ## Prefix for the final output genotype files [merge]
  -h | --help)    grep ";; ##" $0 | grep -v help && exit 1 ;;
  -- ) shift; break ;;
  * ) shift; break ;;
  esac
done

## Check the required parameters
if [[ ! -d ${proj} ]]; then
  echo "project path ${proj} not found! "
  exit 1
fi

## Avoid warnings when running R scripts ("ignoring environment value of R_HOME")
unset R_HOME

## Directory of the script
if [[ ${code} ]]; then
  [[ ! -d ${code} ]] && echo "${code} not exists! " && exit 5
else
  script_path=$(dirname "$(readlink -f "$0")")
  code=$(dirname "$script_path")
fi

## Path/Scripts
func=${code}/shell/function.sh
mrk_sel=${code}/R/QMSim_mrk_select.R
gind_sel=${code}/R/geno_individuals_select.R
PCA_cal=${code}/shell/pca_multi_pop.sh
geno_dist=${code}/shell/distance_multi_pop.sh
PCA_plot=${code}/R/PCA_plot.R
LD_cor=${code}/R/LD_decay_plot.R
# block_LD=${code}/R/block_LD_cor.R
corr_cal=${code}/R/columns_correlation.R

## Add the program path to the environment variable
export PATH=${code}/bin:$PATH

## Load custom functions
[[ ! -s ${func} ]] && echo "Error: ${func} not found! " && exit 5
source ${func}

## Check if the required program can be retrieved and executed in the environment variables
check_command plink QMSim_selected

## Check if the required script files exist and have execution permissions
check_command $gind_sel $PCA_cal $geno_dist $PCA_plot $LD_cor $mrk_sel $corr_cal

## Default parameters
out=${out:=merge}
windows=${windows:="10000"}
r2=${r2:="0"}
inter=${inter:="99999"}
nsnp=${nsnp:="50000"}
geno_gen=${geno_gen:="8-10"}
maf=${maf:="0.01"}
binDiv=${binDiv:="pos"}
binThr=${binThr:="10"}
rep=${rep:="1"}
breeds=${breeds:="A B"}
geno_sel=${geno_sel:="${geno_gen}"}

## Obtain the number of breeds
read -ra breeds <<<"$breeds"
np=${#breeds[@]}

## Default parameters determined by the number of breeds
nginds=${nginds:=$(printf "%${np}s" | sed "s/ /600 /g" | sed 's/ *$//')}
litters=${litters:=$(printf "%${np}s" | sed "s/ /10 /g" | sed 's/ *$//')}
females=${females:=$(printf "%${np}s" | sed "s/ /200 /g" | sed 's/ *$//')}
# echo "line113: nginds=${nginds}"
# echo "line114: females=${females}" && exit

##  parse parameter
read -ra nginds <<<"$nginds"
read -ra litters <<<"$litters"
read -ra females <<<"$females"
qrep=$(printf "%03d" "$rep")
mapfile -t gid_gen < <(seq ${geno_gen:0:1} ${geno_gen:2:3})

## Working directory
cd ${proj} || exit

## Retrieve the first 6 bits of the first random number from QMSim as the seed for the subsequent random process
seed=$(sed -n '2p' seed | awk '{print $2}')
seed=${seed:0:6}
echo ${seed} >random.seed

## Select markers
$mrk_sel \
  --popN "${breeds[*]}" \
  --freqf "%pop%_freq_mrk_${qrep}.txt" \
  --mapf "lm_mrk_${qrep}.txt" \
  --nsel ${nsnp} \
  --binDiv ${binDiv} \
  --binThr ${binThr} \
  --maf ${maf} \
  --out mrk_sel_index.txt \
  --outmapf ${breeds[0]}.map \
  --seed ${seed}

## Select individuals and extract marker information
for j in $(seq 0 $((np - 1))); do
  ## pedigree
  sed '1d' ${breeds[${j}]}_data_${qrep}.txt | awk '{print $1,$2,$3}' >${breeds[${j}]}_pedi.txt

  ## Select genotype individuals
  $gind_sel \
    --dataf ${breeds[${j}]}_data_${qrep}.txt \
    --gen_all ${geno_gen} \
    --gen_sel ${geno_sel} \
    --nsel ${nginds[${j}]} \
    --outIndex \
    --seed ${seed} \
    --out ${breeds[${j}]}_Ind_sel_index.txt

  ## map file
  if [[ ${j} == "0" ]]; then
    ## Convert the location information of Scientific notation that may exist in the map file to an integer
    awk '{printf "%s %s %s %d\n",$1,$2,$3,int($4)}' ${breeds[${j}]}.map >tmp.map
    mv tmp.map ${breeds[${j}]}.map
  else
    cp ${breeds[0]}.map ${breeds[${j}]}.map
  fi

  ## Number of genotype individuals
  ngid=$((${#gid_gen[@]} * ${litters[${j}]} * ${females[${j}]}))

  ## ped file
  QMSim_selected \
    --indexf mrk_sel_index.txt \
    --indIndexf ${breeds[${j}]}_Ind_sel_index.txt \
    --mrkf ${breeds[${j}]}_mrk_${qrep}.txt \
    --nInd ${ngid} \
    --fid ${breeds[${j}]} \
    --out ${breeds[${j}]}.ped &
done
wait

## Quality control
for b in "${breeds[@]}"; do
  plink --file ${b} --maf 0.05 --make-bed --out ${b}q
done

## Screen out marker loci that have passed quality control in all breeds
awk '{print $2}' "${breeds[@]/%/q.bim}" | sort | uniq -c | awk -v n=${np} '$1==n {print $2}' >common.snp
for b in "${breeds[@]}"; do
  plink --bfile ${b}q --extract common.snp --make-bed --out ${b}m
done

## pca calculate
qc_files=("${breeds[@]/%/m}")
$PCA_cal --pre_list "${qc_files[*]}" --fids "${breeds[*]}" --fid

## pca plot
$PCA_plot \
  --eigv "$(printf '%s_' "${qc_files[@]}")pca.txt" \
  --out "$(printf '%s_' "${breeds[@]}")pca"

## LD calculate
for b in "${breeds[@]}"; do
  plink \
    --bfile ${b}m \
    --freq \
    --r2 \
    --ld-window-kb ${windows} \
    --ld-window ${inter} \
    --ld-window-r2 ${r2} \
    --out ${b}m
done

## LD result statistics and plotting
$LD_cor \
  --files "$(printf '%sm ' "${breeds[@]}")" \
  --popN "${breeds[*]}" \
  --bin1 50 \
  --breaks 1000 \
  --bin2 100 \
  --max 5000 \
  --out g_"$(printf '%s_' "${breeds[@]}")"_5Mb

## The gene frequency and LD correlation between breeds
for t in frq ld; do
  :>${t}_cor.txt
  for bi in $(seq 0 $((np - 1))); do
    for bj in $(seq ${bi} $((np - 1))); do
      [[ ${bi} == "${bj}" ]] && continue
      cor=$($corr_cal --file1 ${breeds[${bi}]}m.${t} --file2 ${breeds[${bj}]}m.${t})
      echo "${breeds[${bi}]} ${breeds[${bj}]} ${cor}" >>${t}_cor.txt
    done
  done
  echo "correlation of ${t}:"
  cat ${t}_cor.txt
done

## Merge genotypes of all breeds
: >plink_merge_list.txt
for b in "${breeds[@]}"; do
  [[ ${b} == "${breeds[0]}" ]] && continue
  echo "${b}m" >>plink_merge_list.txt
done
plink \
  --bfile ${breeds[0]}m \
  --merge-list plink_merge_list.txt \
  --maf 0.05 \
  --make-bed \
  --out ${out}

## Genetic distance between breeds
$geno_dist --bfile ${out} --out ${out}.dist.summ

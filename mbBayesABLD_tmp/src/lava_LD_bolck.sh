#!/usr/bin/bash

########################################################################################################################
## Version: 1.3.0
## Author:    Liweining liwn@jaas.ac.cn
## Orcid:     0000-0002-0578-3812
## Institute: College of Animal Science and Technology, China Agricul-tural University, Haidian, 100193, Beijing, China
## Date:      2024-08-20
##
## Function：
## Using ldblock to divide the genome into approximately independent ld blocks based on ld information
##
##
## Usage: ./lava_LD_bolck.sh --bfile "/path/to/plink/binary/file/prefix" ...(Please refer to --help for detailed parameters)
##
## Dependent on software/environment:
##  1. ldblock
##  2. plink/1.9
##  3. Other R languages and Bash scripts
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
TEMP=$(getopt -o h --long bfile:,type:,spar:,diff:,maf:,win:,jobs:,minSize:,code:,out:,keepTmp,plot \
  -n 'javawrap' -- "$@")
if [ $? != 0 ]; then
  echo "Terminating..." >&2
  exit 1
fi
eval set -- "$TEMP"

##  parse parameter
while true; do
    case "$1" in
    --bfile )   bfile="$2";    shift 2 ;; ## PLINK_PREFIX
    --type )    type="$2";     shift 2 ;; ## Method for determining bins [lava/LD]
    --spar )    spar="$2";     shift 2 ;; ## Smoothness (0-1) when fitting curves [0.2]
    --diff )    diff="$2";     shift 2 ;; ##  [0.05]
    --maf )     maf="$2";      shift 2 ;; ## MAF filtering threshold [0.01]
    --geno )    geno="$2";     shift 2 ;; ## remove SNPs with missing call rates [0.2]
    --mind )    mind="$2";     shift 2 ;; ## remove individuals with missing call rates [0.2]
    --win )     win="$2";      shift 2 ;; ## Correlation window in number of SNPs [50]
    --jobs )    jobs="$2";     shift 2 ;; ## Number of parallel jobs [5]
    --minSize ) minSize="$2";  shift 2 ;; ## The minimum size in number of SNPs a block must have [50]
    --code )    code="$2";     shift 2 ;; ## Script path [/BIGDATA2/cau_jfliu_2/liwn/code]
    --out )     out="$2";      shift 2 ;; ## The prefix used for output files [ldblock]
    --nchr )    nchr="$2";     shift 2 ;; ## number of chromosomes [30]
    --keepTmp ) keepTmp=true;  shift   ;; ## Keep intermediate files [false]
    --plot )    plot="--plot"; shift   ;; ## plot mean r2 [FALSE]
    -h | --help ) grep " ;; ## " $0 && exit 1 ;;
    -- ) shift; break ;;
    * ) break ;;
    esac
done

## Folder containing the scripts
if [[ ${code} ]]; then
  [[ ! -d ${code} ]] && echo "${code} not exists! " && exit 5
else
  script_path=$(dirname "$(readlink -f "$0")")
  code=$(dirname "$script_path")
fi

## Software/scripts
RLD=${code}/src/LD_smoothing_block.R
job_pool=${code}/src/parallel_jobs_func.sh
func=${code}/src/function.sh
plink=${code}/bin/plink
ldblock=${code}/bin/ldblock
LD_mean_r2=${code}/bin/LD_mean_r2

## execute permission
chmod 0770 $plink $ldblock $LD_mean_r2

## Load custom functions
[[ ! -s ${func} ]] && echo "Error: ${func} not found! " && exit 5
source ${func}

## Avoid warnings when running R scripts ("ignoring environment value of R_HOME")
unset R_HOME

## Default parameters
out=${out:=ldblock.txt}
type=${type:=LD}
spar=${spar:=0.2}
diff=${diff:=0.05}
maf=${maf:=-0.01}
frq=${frq:=0}
geno=${geno:=1.0}
mind=${mind:=1.0}
win=${win:=50}
minSize=${minSize:=50}
jobs=${jobs:=5}
nchr=${nchr:=30}

## Check file format
if [[ ! -s ${bfile}.bim ]]; then
  if [[ -s ${bfile}.map ]]; then
    $plink --file ${bfile} --chr-set ${nchr} --make-bed --out ${bfile}
  else
    echo "plink file ${bfile}.map(.bim) not found! "
    exit 1
  fi
fi

## Number of parallel jobs
source ${job_pool}
job_pool_init ${jobs} 0

## Random seeds to prevent interference between different processes
seed=$RANDOM

## Chromosome information
chrs=$(awk '{print $1}' ${bfile}.bim | sort -n | uniq)
nchr=$(echo ${chrs} | tr " " "\n" | wc -l)

## Report
echo "Panel for bin definition: $bfile"
echo "Number of chromosomes in panel: $nchr"

## Divide each chromosome into blocks
for chr in ${chrs}; do
  ## Extract chromosome information
  $plink --bfile ${bfile} --chr-set ${nchr} --chr ${chr} --make-bed --out ld_block_tmp.${seed}.${chr} >/dev/null

  ## If the number of SNPs is less than the minimum number of markers in the interval, it will be skipped.
  ## If no results are detected, the chromosome will be separately divided into an interval in the result file
  [[ $(wc -l <ld_block_tmp.${seed}.${chr}.bim) -lt ${minSize} ]] && continue

  ## Genome partitioning
  if [[ ${type} == "lava" ]]; then
    $ldblock \
      ld_block_tmp.${seed}.${chr} \
      -frq ${frq} \
      -win ${win} \
      -min-size ${minSize} \
      -out ld_block_tmp.${seed}.${chr}
  else
    $LD_mean_r2 \
      --bfile ld_block_tmp.${seed}.${chr} \
      --maf ${maf} \
      --geno ${geno} \
      --mind ${mind} \
      --win ${win} \
      --min ${minSize} \
      --out ld_mean.${seed}
    [[ -f ld_mean.${seed}_chr${chr}.txt ]] && mv ld_mean.${seed}_chr${chr}.txt ld_block_tmp.${seed}.${chr}.breaks
  fi
done

## Waiting for the backend program to finish running
job_pool_wait
job_pool_shutdown

## Prepare block files
: >ld_block_tmp.${seed}.breaks
for chr in ${chrs}; do
  block_file=ld_block_tmp.${seed}.${chr}.breaks

  if [[ -s ${block_file} ]]; then
    if [[ ${type} == "lava" ]]; then
      nblock=$(($(wc -l <${block_file}) - 2))

      ## Prepare block files
      mapfile -t -O 1 start < <(awk -v lines="$(wc -l <${block_file})" 'NR>1 && NR<lines {print $6}' ${block_file})
      mapfile -t -O 1 stop < <(awk 'NR>2 {print $6-1}' ${block_file})
      mapfile -t -O 1 nsnp < <(awk 'NR > 2 {print $5 - prev} {prev = $5}' ${block_file})

      ## Output to file
      for line in $(seq ${nblock}); do
        echo "${chr} ${start[line]} ${stop[line]} ${nsnp[line]}" >>ld_block_tmp.${seed}.breaks
      done
    else
      ${RLD} \
        --r2 ${block_file} \
        --bim ld_block_tmp.${seed}.${chr}.bim \
        --spar ${spar} \
        --diff ${diff} \
        --out ld_block_tmp.${seed}.LD \
        ${plot}
      cat ld_block_tmp.${seed}.LD >>ld_block_tmp.${seed}.breaks
    fi
  else
    starti=$(head -n 1 ld_block_tmp.${seed}.${chr}.bim | awk '{print $4}')
    stopi=$(tail -n 1 ld_block_tmp.${seed}.${chr}.bim | awk '{print $4}')
    echo "${chr} ${starti} ${stopi} $(wc -l <ld_block_tmp.${seed}.${chr}.bim)" >>ld_block_tmp.${seed}.breaks
  fi
done

## Add sequence lines
awk '{print NR, $0}' ld_block_tmp.${seed}.breaks >${out}

## Check the number of markers
nsnp_bin=$(awk '{sum+=$5}END{print sum}' ${out})
if [[ ${nsnp_bin} -ne $(wc -l <${bfile}.bim) ]]; then
  echo "number of snp in bin file (${nsnp_bin}) are not equal to bfile ($(wc -l <${bfile}.bim))! "
  # rm ${out}
  rm ld_block_tmp.${seed}*
  exit 2
fi

## Insert title line
sed -i '1i LOC CHR START STOP nSNP' ${out}

## Delete intermediate files
[[ ! ${keepTmp} ]] && rm ld_block_tmp.${seed}*

## Report
echo "number of blocks: $(sed '1d' ${out} | wc -l)"
echo "blocks information file output to: ${out}"

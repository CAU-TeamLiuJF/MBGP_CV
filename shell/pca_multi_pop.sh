#!/usr/bin/bash

########################################################################################################################
## Version: 1.3.0
## Author:    Liweining liwn@cau.edu.cn
## Orcid:     0000-0002-0578-3812
## Institute: College of Animal Science and Technology, China Agricul-tural University, Haidian, 100193, Beijing, China
## Date:      2024-08-20
##
## Function：
## Perform PCA analysis and plot for any number of breeds provided
##
##
## Usage:
## ./pca_multi_pop.sh --pre_list "/path/to/plink/binary/files/prefix" ...(Please refer to --help for detailed parameters)
##
## Dependent on software/environment:
##  1. R
##  2. plink/1.9
##  3. Other R languages and Bash scripts
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
TEMP=$(getopt -o h --long pre_list:,fids:,out:,nchr:,help,fid,plot \
              -n 'javawrap' -- "$@")
if [ $? != 0 ] ; then echo "Terminating..." >&2 ; exit 1 ; fi
eval set -- "$TEMP"
##  parse parameter
while true; do
  case "$1" in
          --pre_list )   pre_list="$2";  shift 2 ;;  ## plink file prefix, e.g., "/public/home/popA /public/home/popB"
          --fids )       fids="$2";      shift 2 ;;  ## Breed (population) identifiers, e.g., "popA popB"
          --nchr )       nchr="$2" ;     shift 2 ;;  ## Number of chromosomes [30]
          --out )        out="$2" ;      shift 2 ;;  ## Output file name
          --fid )        fid=true;       shift   ;;  ## Breed (population) identifiers provided in the bfile (i.e., fid)
          --plot )       plot=true;      shift   ;;  ## Plot based on PCA results
    -h | --help )        grep " shift " $0 && exit 1 ;;
    -- ) shift; break ;;
    * ) shift; break ;;
  esac
done

## Directory of the script
if [[ ${code} ]]; then
  [[ ! -d ${code} ]] && echo "${code} not exists! " && exit 5
else
  script_path=$(dirname "$(readlink -f "$0")")
  code=$(dirname "$script_path")
fi

## Scripts
pca_plot=${code}/R/PCA_plot.R
func=${code}/shell/function.sh

## path to log information
logp=${code}/log

## Load custom functions
[[ ! -s ${func} ]] && echo "Error: ${func} not found! " && exit 5
source ${func}

## Check if the required programs are available in the environment and executable
check_command plink

## Check required parameters
[[ ! ${pre_list} ]] && echo "Open script $0 to view instructions" && exit 1

## Default parameters
nchr=${nchr:="30"}

## If only one file is provided, group by fid in the fam file
files_num=$(echo ${pre_list} | tr " " "\n" | wc -l)
if [[ ${files_num} -le 1 ]]; then
  unset fids
  if [[ ! ${fid} ]]; then
    echo "please provide the --fid option if the .ped(.fam) file contains fid"
    exit 1
  fi

  ## Check if the file exists
  check_plink "${pre_list}" ${nchr}

  ## Family ID in the plink genotype file
  fid_uniq=$(awk '{print $1}' ${pre_list}.fam | sort | uniq)

  ## Extract genotype information for different breeds
  for name in ${fid_uniq}; do
    ## Family ID
    awk ${name} >${name}_fid.txt

    ## Extract individuals with specified family ID
    plink \
      --bfile ${pre_list} \
      --keep-fam ${name}_fid.txt \
      --chr-set ${nchr} \
      --make-bed --out ${name} >${logp}/plink_pca_keep_fam.log
    rm ${name}_fid.txt
  done

  pre_list=("${fid_uniq}")
  # Name_list="${pre_list[*]}"
fi

## Get the identifier (family id) for each population
# pre_list=("${pre_list[@]}")
IFS=' ' read -ra pre_list <<< "${pre_list[@]}"
index=0
for prefix in "${pre_list[@]}"; do
  ((index++))

  ## Check if the file exists
  check_plink "${prefix}" ${nchr}

  ## Population id
  if [[ ${fids[*]} ]]; then
    fidi=$(echo "${fids[*]}" | cut -d ' ' -f ${index})
  else
    ## Extract the first row's iid and fid, and check if fid is represented by iid or all zeros. If so, assign a new fid (e.g., pop1)
    iid=$(head -n 1 ${prefix}.fam | awk '{print $2}')
    fidi=$(head -n 1 ${prefix}.fam | awk '{print $1}')
    if [[ ${fidi} == "${iid}" || ${fidi} == '0' ]]; then
      fidi=pop${index}
    fi
    # Name_list="${Name_list} ${fidi}"
  fi

  ## Output the mapping table of iid and fid
  if [[ ${index} -eq 1 ]]; then
    awk '{print $2,"'${fidi}'"}' ${prefix}.fam >iid_fid_pca.txt
  else
    awk '{print $2,"'${fidi}'"}' ${prefix}.fam >>iid_fid_pca.txt

    ## File list needed for merging populations
    if [[ ${index} -ge 2 ]]; then
      echo ${prefix}.bed ${prefix}.bim ${prefix}.fam >>merge_list.txt
    else
      [[ -f merge_list.txt ]] && rm merge_list.txt
    fi
  fi
done

## Output file name
[[ ! ${out} ]] && out=$(echo "${pre_list[@]}" | tr " " "_")_pca.txt

## Merge all population plink files
prefix1=$(echo "${pre_list[@]}" | cut -d ' ' -f 1)
num=$(echo "${pre_list[@]}" | wc -w)
plink \
  --bfile ${prefix1} \
  --merge-list merge_list.txt \
  --chr-set ${nchr} \
  --make-bed --out pop${num}_merge >${logp}/plink_pca_merge.log

## Multiple allele sites exist
if [[ $? -ne 0 ]]; then
  echo "please remove all offending variants (merged.missnp) in all populations"
  exit 2
fi

## PCA calculation
plink \
  --bfile pop${num}_merge \
  --chr-set ${nchr} \
  --pca 3 header \
  --out pop${num}_merge >${logp}/plink_pca.log

## Match population identifiers
if [[ ! ${fid} ]]; then
  sort -k 1n iid_fid_pca.txt -o iid_fid_pca.txt2
  sed '1d' pop${num}_merge.eigenvec | sort -k 2n > pop${num}_merge.eigenvec2
  join -1 2 -2 1 pop${num}_merge.eigenvec2 iid_fid_pca.txt2 >"${out}"
else
  cp pop${num}_merge.eigenvec "${out}" 
fi

## Plotting
if [[ ${plot} ]]; then
  $pca_plot --eigv ${out}
fi

## Delete intermediate files
rm pop${num}_merge.*
rm merge_list.*
rm iid_fid_pca.*

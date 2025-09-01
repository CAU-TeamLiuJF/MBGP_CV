#!/bin/bash

########################################################################################################################
## Version: 1.3.0
## Author:    Liweining liwn@jaas.ac.cn
## Orcid:     0000-0002-0578-3812
## Institute: College of Animal Science and Technology, China Agricul-tural University, Haidian, 100193, Beijing, China
## Date:      2024-08-20
##
## Function：
## Save custom functions
##
##
## Usage:
## ./function.sh
##
##
## License:
##  This script is licensed under the GPL-3.0 License.
##  See https://www.gnu.org/licenses/gpl-3.0.en.html for details.
########################################################################################################################


check_alphabet() {
  ## Check if the file contains any alphabetic characters (DMU software does not accept character-type values)
  [[ ! -s ${1} ]] && echo "${1} not found! " && exit 1
  if [[ ${2} ]]; then
    NotNum=$(awk -vl=${2} '{print $l}' ${1} | grep -c "[a-zA-Z]")
  else
    NotNum=$(grep -v "[0-9]e" <${1} | grep -c "[a-zA-Z]")
  fi

  if [[ ${NotNum} -gt 0 ]]; then
    echo "Non-numeric characters exist in ${1} file, please check!"
    exit 1
  fi
}

check_command() {
  # Check if the script file/software exists and has execution permission
  # Parameters:
  # $@：List of software names to check
  for cmd in "$@"; do
    ## Determine if it's a file and has execute permission
    if [[ -s ${cmd} ]]; then
      if [[ ! -x "$cmd" ]]; then
        # chmod 0777 "$cmd"
        echo "Error: \"$cmd\" does not have execute permission"
        exit 1
      else
        return 0
      fi
    fi

    ## Determine if it's a command and has execute permission
    if ! command -v "$cmd" &> /dev/null; then
      echo "Error: \"$cmd\" not found (not installed or not in PATH)"
      exit 1
    elif [[ ! -x $(command -v "$cmd") ]]; then
      # chmod 0777 "$cmd"
      echo "Error: \"$cmd\" does not have execute permission"
      exit 1
    fi
  done
}

check_plink() {
  ## Check if the plink binary files with the specified prefix exist
  ## If map and ped files exist, convert them to binary format
  ## Parameter 1: plink file name prefix, e.g., /path/to/plink_prefix
  ## Parameter 2: Number of chromosomes, default is 30

  ## Check if plink software is available
  check_command plink

  ## Set the number of chromosomes
  if [[ $2 ]]; then
    nchr=$2
  else
    nchr=30
  fi

  ## File format
  for prefix in $1; do
    if [[ ! -f ${prefix}.fam ]]; then
      if [[ -f ${prefix}.map ]]; then
        plink \
          --file ${prefix} \
          --chr-set ${nchr} \
          --make-bed \
          --out ${prefix}
      else
        echo "Error: ${prefix}.bed not found"
        # exit 1
      fi
    fi
  done
}

merge_plink() {
  ## Merge plink files
  ## Parameter 1: List of plink files, separated by spaces, e.g., "/path/to/plinkA_prefix /path/to/plinkB_prefix", should be enclosed in double quotes
  ## Parameter 2: Output file name prefix, e.g., /path/to/plink_merge_prefix
  ## Parameter 3: Number of chromosomes, default is 30
  ## Dependency: plink

  ## List of file names to merge
  local plink_list
  IF=" " read -r -a plink_list <<<"$1"

  ## Set the number of chromosomes
  if [[ $3 ]]; then
    nchr=$3
  else
    nchr=30
  fi

  ## Check the number of parameters
  nfiles=${#plink_list[@]}
  if [[ ${nfiles} -lt 2 ]]; then
    echo "Error: at least two plink files are needed to merge"
    return 1
  fi

  ## Temporary file to store the prefixes of the population files to be merged
  local seed_tmp=$RANDOM
  touch merge_${seed_tmp}.txt

  ## Check if plink files for each population exist
  check_plink "${plink_list[0]}" ${nchr}
  for i in $(seq 1 $((nfiles - 1))); do
    check_plink "${plink_list[${i}]}" ${nchr}
    echo ${plink_list[${i}]} >>merge_${seed_tmp}.txt
  done

  ## Check if the plink files with the specified output prefix already exist
  if [[ -f $2.bed ]]; then
    echo "Warning: $2.bed already exists and will be overwritten!"
  fi

  ## Merge files
  [[ -s merge_${seed_tmp}.txt ]] && \
    plink \
      --bfile ${plink_list[0]} \
      --chr-set ${nchr} \
      --merge-list merge_${seed_tmp}.txt \
      --make-bed \
      --out $2

  rm merge_${seed_tmp}.txt
}

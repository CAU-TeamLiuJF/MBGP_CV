#!/bin/bash
#SBATCH --job-name=QMSim

########################################################################################################################
## Version: 1.3.0
## Author:    Liweining liwn@cau.edu.cn
## Orcid:     0000-0002-0578-3812
## Institute: College of Animal Science and Technology, China Agricul-tural University, Haidian, 100193, Beijing, China
## Date:      2024-08-20
##
## Function：
## Generate the required parameter file *. prm for QMSim based on the specified parameters, and run QMSim for population simulation
##
##
## Usage:
## ./run_QMSim.sh --proj "/path/to/project" ...(Please refer to --help for detailed parameters)
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
TEMP=$(getopt -o h --long proj:,breeds:,thread:,seg_gens:,extentLDs:,last_males:,last_females:,founder_sel:,seg_sel:,last_sel:,last_litters:,geno_gen:,nchr:,nmloc:,nqloci:,QMSim_h2:,QMSim_qtlh2:,QMSim_rep:,bottleneck:,code:,prmpath:,sim_dir:,logf:,help \
              -n 'javawrap' -- "$@")
if [ $? != 0 ] ; then echo "Terminating..." >&2 ; exit 1 ; fi
eval set -- "$TEMP"
##  parse parameter
while true; do
  case "$1" in
    --proj )         proj="$2";         shift 2 ;; ## Project directory [required parameter]
    --breeds )       breeds="$2";       shift 2 ;; ## Population/Breed identifiers, e.g., 'YY DD' ["A B"]
    --thread )       thread="$2";       shift 2 ;; ## Number of threads for running QMSim [10]
    --seg_gens )     seg_gens="$2";     shift 2 ;; ## Number of generations after separation from the historical population for each breed ["40 10"]
    --extentLDs )    extentLDs="$2";    shift 2 ;; ## Number of generations in the last stable LD phase for each breed ["10 10"]
    --last_males )   last_males="$2";   shift 2 ;; ## Number of males in the population in the last phase for each breed ["100 10"]
    --last_females ) last_females="$2"; shift 2 ;; ## Number of females in the population in the last phase for each breed ["500 50"]
    --founder_sel )  founder_sel="$2";  shift 2 ;; ## Criteria for selecting individuals from the historical population for each breed ["tbv /h,tbv /l"]
    --seg_sel )      seg_sel="$2";      shift 2 ;; ## Criteria for selecting individuals during the generation selection phase for each breed ["phen /h,phen /l"]
    --last_sel )     last_sel="$2";     shift 2 ;; ## Criteria for selecting individuals during the LD stabilization phase for each breed ["rnd,rnd"]
    --last_litters ) last_litters="$2"; shift 2 ;; ## Number of individuals per litter during the LD stabilization phase for each breed ["10 10"]
    --geno_gen )     geno_gen="$2";     shift 2 ;; ## Generation(s) to output genotyped individuals [8-10]
    --nchr )         nchr="$2";         shift 2 ;; ## Number of chromosomes [18]
    --nmloc )        nmloc="$2";        shift 2 ;; ## Number of markers per chromosome [300000]
    --nqloci )       nqloci="$2";       shift 2 ;; ## Number of QTLs per chromosome [100]
    --QMSim_h2 )     QMSim_h2="$2";     shift 2 ;; ## Broad-sense heritability of the trait [0.3]
    --QMSim_qtlh2 )  QMSim_qtlh2="$2";  shift 2 ;; ## Narrow-sense heritability of the trait [0.3]
    --QMSim_rep )    QMSim_rep="$2";    shift 2 ;; ## Number of QMSim simulation replicates [1]
    --bottleneck )   bottleneck="$2";   shift 2 ;; ## Population size during the bottleneck phase in the historical population simulation [250]
    --logf )         logf="$2";         shift 2 ;; ## Log file name [QMSim_${sim_dir}.txt]
    --code )         code="$2";         shift 2 ;; ## Directory where script files are located, e.g., /BIGDATA2/cau_jfliu_2/liwn/code [automatically obtained]
    --prmpath )      prmpath="$2";      shift 2 ;; ## Path to the parameter template files [${code}/prm]
    --sim_dir )      sim_dir="$2";      shift 2 ;; ## Name of the output folder for simulation results; this folder should not already exist [rep1]
  -h | --help)    grep ";; ##" $0 | grep -v help && exit 1 ;;
  -- ) shift; break ;;
  * ) shift; break ;;
  esac
done

## Check the required parameters
if [[ ! -d ${proj} ]]; then
  echo "${proj} not found! "
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

## Add the program path to the environment variable
export PATH=${code}/bin:$PATH

## Path/Scripts
func=${code}/shell/function.sh
gind_sel=${code}/R/geno_individuals_select.R
PCA_cal=${code}/shell/pca_multi_pop.sh
geno_dist=${code}/shell/distance_multi_pop.sh
PCA_plot=${code}/R/PCA_plot.R
LD_cor=${code}/R/LD_decay_plot.R
mrk_sel=${code}/R/QMSim_mrk_select.R
block_LD=${code}/R/block_LD_cor.R

## Parameter card template file
prm_hist=templete_gloabal_history.prm
prm_sub=templete_subpopulation.prm
prm_geno=templete_genome_output.prm

## Load custom functions
[[ ! -s ${func} ]] && echo "Error: ${func} not found! " && exit 5
source ${func}

## Check if the required program can be retrieved and executed in the environment variables
check_command plink QMSim

## Check if the required script files exist and have execution permissions
check_command $gind_sel $PCA_cal $geno_dist $PCA_plot $LD_cor $mrk_sel $block_LD

## Default parameters
breeds=${breeds:="A B"}
sim_dir=${sim_dir:=rep1}
QMSim_rep=${QMSim_rep:="1"}
nchr=${nchr:="18"}
nmloc=${nmloc:="300000"}
geno_gen=${geno_gen:="8-10"}
nqloci=${nqloci:="100"}
QMSim_h2=${QMSim_h2:="0.3"}
QMSim_qtlh2=${QMSim_qtlh2:="${QMSim_h2}"}
bottleneck=${bottleneck:="250"}
thread=${thread:="10"}
prmpath=${prmpath:="${code}/prm"}
logf=${logf:="QMSim_${sim_dir}.txt"}

## Obtain the number of breeds
read -ra breeds <<<"$breeds"
np=${#breeds[@]}

## Default parameters determined by the number of breeds
founder_sel=${founder_sel:=$(printf "%${np}s" | sed "s/ /rnd,/g" | sed 's/,$//')}
seg_sel=${seg_sel:=${founder_sel}}
last_sel=${last_sel:=${founder_sel}}
last_litters=${last_litters:=$(printf "%${np}s" | sed "s/ /10 /g" | sed 's/ *$//')}
extentLDs=${extentLDs:=$(printf "%${np}s" | sed "s/ /10 /g" | sed 's/ *$//')}
last_males=${last_males:=$(printf "%${np}s" | sed "s/ /40 /g" | sed 's/ *$//')}
last_females=${last_females:=$(printf "%${np}s" | sed "s/ /200 /g" | sed 's/ *$//')}
seg_gens=${seg_gens:=$(printf "%${np}s" | sed "s/ /40 /g" | sed 's/ *$//')}

##  parse parameter
read -ra seg_gens <<<"$seg_gens"
read -ra last_males <<<"$last_males"
read -ra last_females <<<"$last_females"
read -ra extentLDs <<<"$extentLDs"
read -ra last_litters <<<"$last_litters"
IFS=, read -ra founder_sel <<<"$founder_sel"
IFS=, read -ra seg_sel <<<"$seg_sel"
IFS=, read -ra last_sel <<<"$last_sel"

mapfile -t gid_gen < <(seq ${geno_gen:0:1} ${geno_gen:2:3})
gen_all=${geno_gen:2:3}
rand=$RANDOM

## Check if the parameter card template file exists
if [[ ! -d ${prmpath} ]]; then
  echo "${prmpath} not exists! "
  exit 5
else
  for f in $prm_hist $prm_sub $prm_geno; do
    [[ ! -s ${prmpath}/$f ]] && echo "$f not found! " && exit 4
  done
fi

###################### Parameter card of simulated population  ######################
cd ${prmpath} || exit

## History population 
sed "s/%nthread%/${thread}/" $prm_hist >temp_${rand}.prm
sed -i "s/%rep%/${QMSim_rep}/" temp_${rand}.prm
sed -i "s/%h2%/${QMSim_h2}/" temp_${rand}.prm
sed -i "s/%qtlh2%/${QMSim_qtlh2}/" temp_${rand}.prm
sed -i "s/%bottleneck%/${bottleneck}/" temp_${rand}.prm

## Subpopulation
for i in $(seq 0 $((np - 1))); do
  sed "s/%pop%/${breeds[${i}]}/" $prm_sub >${breeds[${i}]}.prm
  sed -i "s#%founder_select%#${founder_sel[${i}]}#" ${breeds[${i}]}.prm
  sed -i "s#%seg_select%#${seg_sel[${i}]}#" ${breeds[${i}]}.prm
  sed -i "s#%last_select%#${last_sel[${i}]}#" ${breeds[${i}]}.prm
  sed -i "s/%md%/rnd/" ${breeds[${i}]}.prm
  sed -i "s/%accur%/${QMSim_qtlh2}/" ${breeds[${i}]}.prm
  sed -i "s/%seg_ng%/${seg_gens[${i}]}/" ${breeds[${i}]}.prm
  sed -i "s/%last_male%/${last_males[${i}]}/" ${breeds[${i}]}.prm
  sed -i "s/%last_female%/${last_females[${i}]}/" ${breeds[${i}]}.prm
  sed -i "s/%extentLD%/${extentLDs[${i}]}/" ${breeds[${i}]}.prm
  sed -i "s/%last_litter_size%/${last_litters[${i}]}/" ${breeds[${i}]}.prm
  sed -i "s/%geno_gen%/${gid_gen[*]}/g" ${breeds[${i}]}.prm
  sed -i "s/%freq_gen%/${gid_gen[*]}/g" ${breeds[${i}]}.prm

  ## Merge into the main parameter file
  cat ${breeds[${i}]}.prm >>temp_${rand}.prm
  rm ${breeds[${i}]}.prm
done

## Genome and output parameters
sed "s/%nmloci%/${nmloc}/" $prm_geno >genome.prm
sed -i "s/%nchr%/${nchr}/" genome.prm
sed -i "s/%nqloci%/${nqloci}/" genome.prm

## Final simulation parameter card
cat genome.prm >>temp_${rand}.prm
rm genome.prm

## Project directory
cd ${proj} || exit 5

## Output parameter card
sed "s#%out_dir%#${sim_dir}#" ${prmpath}/temp_${rand}.prm >QMSim.prm
rm ${prmpath}/temp_${rand}.prm

## QMSim simulation program running
QMSim QMSim.prm >> ${logf} 2>&1

## Change working path
rm QMSim.prm
cd ${proj}/${sim_dir} || exit

## Save QMSim simulation seeds
cp seed QMSim_${sim_dir}_seed.prv

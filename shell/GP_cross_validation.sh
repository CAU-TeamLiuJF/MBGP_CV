#!/usr/bin/bash
#SBATCH --job-name=GP_CV

########################################################################################################################
## Version:   1.2.0
## Author:    Liweining liwn@jaas.ac.cn
## Orcid:     0000-0002-0578-3812
## Institute: College of Animal Science and Technology, China Agricul-tural University, Haidian, 100193, Beijing, China
## Date:      2024-07-07
## 
## 
## Usage: GP_cross_validation.sh --help
## 
## License:
##  This script is licensed under the GPL-3.0 License.
##  See https://www.gnu.org/licenses/gpl-3.0.en.html for details.
########################################################################################################################



###################  Parameter processing  #####################
################################################################
## NOTE: This requires GNU getopt.  On Mac OS X and FreeBSD, you have to install this
## Parameters
TEMP=$(getopt -o h --long code:,proj:,type:,breeds:,thread:,traits:,trait:,h2s:,rg_sim:,rg_pri:,rg_dist:,means:,method:,phef:,pedf:,binf:,nqtl:,nbin_cor:,nsnp_cor:,nsnp_win:,prior:,bin:,bin_sim:,tbv_col:,all_eff:,ran_eff:,iter:,burnin:,ref:,dirPre:,nbin:,min:,bfile:,seed:,fold:,rep:,gen:,nsnp:,nsnp_sim:,out:,sim_dir:,nginds:,seg_gens:,extentLDs:,last_males:,last_females:,founder_sel:,seg_sel:,last_sel:,last_litters:,geno_gen:,maf:,binDiv:,binThr:,nchr:,nmloc:,nqloci:,QMSim_h2:,debug,suffix,dense,noCov,append,overlap,evenly,help \
              -n 'javawrap' -- "$@")
if [ $? != 0 ] ; then echo "Terminating..." >&2 ; exit 1 ; fi
eval set -- "$TEMP"
## Parse arguments
while true; do
  case "$1" in
    --proj )     proj="$2";     shift 2 ;; ## Project directory [required]
    --breeds )   breeds="$2";   shift 2 ;; ## Breed identifiers, e.g., 'YY DD' [required]
    --traits )   traits="$2";   shift 2 ;; ## All trait names in phenotype file, e.g., 'AGE BF' [required]
    --type )     type="$2";     shift 2 ;; ## Analysis type: bin/adj/gsim/geno/psim/within/single/multi/accur/var [required]
    --trait )    trait="$2";    shift 2 ;; ## Trait name to be analyzed, e.g., 'AGE' [traits]
    --bfile )    bfile="$2";    shift 2 ;; ## plink file prefix, e.g., "/public/home/merge" [NULL]
    --phef )     phef="$2";     shift 2 ;; ## Phenotype file [NULL]
    --pedf )     phef="$2";     shift 2 ;; ## Pedigree file [NULL]
    --rg_dist )  rg_dist="$2";  shift 2 ;; ## Distribution of additive genetic correlation size, uniform/normal/identical [identical]
    --rg_sim )   rg_sim="$2";   shift 2 ;; ## Additive genetic correlation size between breeds in trait simulation [0.2]
    --rg_pri )   rg_pri="$2";   shift 2 ;; ## Prior for additive genetic correlation between breeds in multi-trait model [NULL]
    --h2s )      h2s="$2";      shift 2 ;; ## Heritability size of traits for each breed in trait simulation [0.2]
    --means )    means="$2";    shift 2 ;; ## Population means for each breed in trait simulation [NULL]
    --iter )     iter="$2";     shift 2 ;; ## Total MCMC iterations [30000]
    --burnin )   burnin="$2";   shift 2 ;; ## MCMC burn-in iterations [20000]
    --ref )      ref="$2";      shift 2 ;; ## SNP panel used for interval partitioning, M/1/2/... [M]
    --dirPre )   dirPre="$2";   shift 2 ;; ## Prefix added to folders in multi-scenario, e.g., pre_multi_A_B [NULL]
    --method)    method="$2";   shift 2 ;; ## Breeding value estimation method: BLUP/GBLUP/ssGBLUP/mbBayesAB/bayesR [GBLUP]
    --nsnp_cor ) nsnp_cor="$2"; shift 2 ;; ## Number of QTLs in the interval when genetic correlation exists between breeds [10]
    --nbin_cor ) nbin_cor="$2"; shift 2 ;; ## Number of intervals when genetic correlation exists between breeds [10]
    --nsnp_sim ) nsnp_sim="$2"; shift 2 ;; ## Number of markers per interval when partitioning by fixed number in trait simulation [60]
    --bin_sim )  bin_sim="$2";  shift 2 ;; ## Interval partitioning method in trait simulation, binf path/win/chr [win]
    --nsnp_win ) nsnp_win="$2"; shift 2 ;; ## Number of markers per interval when partitioning by fixed number for multi-breed prediction [100]
    --nqtl )     nqtl="$2";     shift 2 ;; ## Number of QTLs in trait simulation [300]
    --nbin )     nbin="$2";     shift 2 ;; ## Approximate number of intervals when partitioning by fixed SNP number for multi-breed prediction [NULL]
    --min )      min="$2";      shift 2 ;; ## Minimum number of SNP markers in the interval for selecting QTLs [NULL]
    --tbv_col )  tbv_col="$2";  shift 2 ;; ## Column for true breeding values in phenotype file without correction [NULL]
    --bin )      bin="$2";      shift 2 ;; ## Interval partitioning method in multi-breed prediction, fix/frq/ld/lava/cubic [fix]
    --prior )    prior="$2";    shift 2 ;; ## Prior file for variance components in multi-breed prediction [NULL]
    --all_eff )  all_eff="$2";  shift 2 ;; ## Columns for fixed and random effects, e.g., "2 1" ["2 1"]
    --ran_eff )  ran_eff="$2";  shift 2 ;; ## Column for random effects, e.g., "1" [1]
    --binf )     binf="$2";     shift 2 ;; ## Interval partitioning file [NULL]
    --rep)       rep="$2";      shift 2 ;; ## Number of cross-validation repeats [1]
    --fold)      fold="$2";     shift 2 ;; ## Number of cross-validation folds [NULL]
    --gen)       gen="$2";      shift 2 ;; ## Generation of individuals with phenotypes used as validation population [ALL]
    --seed )     seed="$2";     shift 2 ;; ## Random seed for MCMC sampling and validation group partitioning [40296]
    --code )     code="$2";     shift 2 ;; ## Directory of script files, e.g., /BIGDATA2/cau_jfliu_2/liwn/code [NULL]
    --thread )   thread="$2";   shift 2 ;; ## Number of threads [NULL]
    --sim_dir )      sim_dir="$2";      shift 2 ;; ## Name of output folder for simulation results (folder must not already exist) [rep1]
    --nginds )       nginds="$2";       shift 2 ;; ## Number of individuals selected from each breed for output genotype population ["600 ..."]
    --seg_gens )     seg_gens="$2";     shift 2 ;; ## Number of generations each breed has been segregated from the historical population ["40 10"]
    --extentLDs )    extentLDs="$2";    shift 2 ;; ## Number of generations each breed has stabilized LD ["10 10"]
    --last_males )   last_males="$2";   shift 2 ;; ## Number of male individuals in the final stage of the population for each breed ["100 10"]
    --last_females ) last_females="$2"; shift 2 ;; ## Number of female individuals in the final stage of the population for each breed ["500 50"]
    --founder_sel )  founder_sel="$2";  shift 2 ;; ## Criteria for selecting from the historical population for each breed ["tbv /h,tbv /l"]
    --seg_sel )      seg_sel="$2";      shift 2 ;; ## Criteria for individual selection of each breed ["phen /h,phen /l"]
    --last_sel )     last_sel="$2";     shift 2 ;; ## Criteria for individual selection during the LD stabilization for each breed ["rnd,rnd"]
    --last_litters ) last_litters="$2"; shift 2 ;; ## Number of individuals per litter in the LD stabilization stage for each breed ["10 10"]
    --geno_gen )     geno_gen="$2";     shift 2 ;; ## Generation of output genotype individuals [8-10]
    --nchr )         nchr="$2";         shift 2 ;; ## Number of chromosomes [30]
    --nmloc )        nmloc="$2";        shift 2 ;; ## Number of markers per chromosome [300000]
    --nqloci )       nqloci="$2";       shift 2 ;; ## Number of QTLs per chromosome [100]
    --QMSim_h2 )     QMSim_h2="$2";     shift 2 ;; ## Broad-sense heritability of traits [0.3]
    --maf )          maf="$2";          shift 2 ;; ## Criteria for selecting simulated genotype markers within breed [0.01]
    --binDiv )       binDiv="$2";       shift 2 ;; ## Basis for interval partitioning when sampling SNPs, pos/frq [pos]
    --binThr )       binThr="$2";       shift 2 ;; ## Length of partition during SNP sampling, physical position (cM) or frequency step [10]
    --out )          out="$2";          shift 2 ;; ## Output filename [depends on analysis type]
    --nsnp )         nsnp="$2";         shift 2 ;; ## Number of markers to be selected in simulated genotype [50000]
    --debug )        debug=true;        shift   ;; ## Not running GEBV calculation program, usually used as a debugging test
    --dense )        dense=true;        shift   ;; ## Use dense matrix algorithm (multi-threaded) when using DMUAI module, i.e., method is 31
    --rg_local )     rg_local=true;     shift   ;; ## Different covariance for each genomic partition, for local ld, frq correlation coefficient
    --all_comb )     all_comb=true;     shift   ;; ## Evaluate all possible breed combinations in breeds
    --noCov )        noCov=true;        shift   ;; ## Constrain residual effect between traits to 0
    --overlap )      overlap=true;      shift   ;; ## QTL included in SNP markers
    --suffix )       suffix=true;       shift   ;; ## Add breed name suffix to union/blend/multi folder names, e.g., blend_YY_LL
    --evenly )       evenly=true;       shift   ;; ## Sample QTLs evenly across chromosomes
    --) shift; break ;;
    * ) echo "unknown option: $1"; exit 1 ;;
  esac
done


## Save the current working directory
workdir=$(pwd)

## Check if the required parameters are provided
if [[ ! -d ${proj} ]]; then
  echo "path ${proj} not found! "
  exit 1
fi

## Log folder
logp=${proj}/log
mkdir -p ${logp}

## Avoid warnings when running R scripts ("ignoring environment value of R_HOME")
unset R_HOME

## Folder containing the scripts
if [[ ${code} ]]; then
  [[ ! -d ${code} ]] && echo "${code} not exists! " && exit 5
else
  script_path=$(dirname "$(readlink -f "$0")")
  code=$(dirname "$script_path")
fi

## Scripts
pheno_sim=${code}/R/pheno_simulation.R
phe_adj=${code}/shell/dmu_get_pheno_adj.sh
GP_single=${code}/shell/GP_single_breed.sh
GP_multi=${code}/shell/GP_multi_breed.sh
run_QMSim=${code}/shell/run_QMSim.sh
genome_process=${code}/shell/QMSim_genome_process.sh
accuracy_summ=${code}/shell/accuracy_summary.sh
varComp_summ=${code}/shell/varcomp_summary.sh
time_summ=${code}/shell/time_summary.sh
block_define=${code}/shell/lava_cubic_bolck.sh
func=${code}/shell/function.sh

## Add program path to the environment variable
export PATH=${code}/bin:$PATH

## Load custom functions
[[ ! -s ${func} ]] && echo "Error: ${func} not found! " && exit 5
source ${func}

## Check if required programs are in the environment path and executable
check_command plink gmatrix mbBayesABLD LD_mean_r2 run_dmu4 run_dmuai

## Check if required script files exist and are executable
check_command $pheno_sim $GP_single $GP_multi $block_define $phe_adj

## Default parameters
method=${method:="GBLUP"}
all_eff=${all_eff:="2 1"}
ran_eff=${ran_eff:="1"}
seed=${seed:="8123"}
h2s=${h2s:="0.2"}
iter=${iter:="30000"}
burnin=${burnin:="20000"}
ref=${ref:="M"}
nqtl=${nqtl:="300"}
nsnp_cor=${nsnp_cor:="10"}
nbin_cor=${nbin_cor:="10"}
nsnp_sim=${nsnp_sim:="60"}
nsnp=${nsnp:="50000"}
nsnp_win=${nsnp_win:="100"}
bin_sim=${bin_sim:="win"}
type=${type:="sim"}
trait=${trait:="${traits}"}
bin=${bin:="fix"}
QMSim_rep=${QMSim_rep:="1"}
nchr=${nchr:="30"}
nmloc=${nmloc:="300000"}
nqloci=${nqloci:="100"}
QMSim_h2=${QMSim_h2:="0.3"}
QMSim_qtlh2=${QMSim_qtlh2:="${QMSim_h2}"}
bottleneck=${bottleneck:="250"}
prmpath=${prmpath:="${code}/prm"}
sim_dir=${sim_dir:="rep1"}
geno_gen=${geno_gen:="8-10"}
SLURM_JOB_ID=${SLURM_JOB_ID:="$RANDOM"}
binDiv=${binDiv:="pos"}
binThr=${binThr:="10"}
min=${min:="${nsnp_cor}"}
maf=${maf:=-0.01}

## Default parameters for different types of analysis
if [[ ${type} == "var" || ${type} == "accur" ]]; then
  ## Default parameters
  rg_sim=${rg_sim:=/}
  rg_dist=${rg_dist:=/}
  rep=${rep:=/}
  [[ ${out} ]] && out=" --out ${out} "
elif [[ ${type} == "psim" ]]; then
  rg_sim=${rg_sim:="0.2"}
  rg_dist=${rg_dist:=identical}
elif [[ ${type} == "within" ]]; then
  if [[ ! "ssGBLUP|bayesR|mbBayesAB|" =~ ${method} ]]; then
    echo "method can only be BLUP/GBLUP/ssGBLUP/mbBayesAB/bayesR."
    exit 2
  fi
  out=${out:="accur_${method}.txt"}
fi

## Extract information from parameters
read -ra breeds_array <<<"$breeds"
read -ra traits_array <<<"$traits"
read -ra trait_array <<<"$trait"
read -ra sim_dirs <<<"$sim_dir"
nbreed=${#breeds_array[@]}
ntrait=${#trait_array[@]}
[[ ${overlap} ]] && overlap=" --overlap "
[[ ${evenly} ]] && evenly=" --evenly "
[[ ${dense} ]] && dense=" --dense "
[[ ${debug} ]] && debug=" --debug "
[[ ${suffix} ]] && suffix=" --suffix "
[[ ${noCov} ]] && noCov=" --res_const "
[[ ${rg_local} ]] && rg_local=" --rg_local "
[[ ${gen} ]] && gen=" --gen ${gen} "
[[ ${binf} ]] && binf=" --binf ${binf} "
[[ ${pedf} ]] && pedf=" --pedf ${pedf} "
[[ ${tbv_col} ]] && tbv_col=" --tbv_col ${tbv_col} "
[[ ${nbin} ]] && nbin=" --nbin ${nbin} "
[[ ${rg} ]] && rg=" --rg ${rg} "
[[ ${dirPre} ]] && dirPre=" --prefix ${dirPre} "
[[ ${prior} ]] && prior=" --priorVar ${prior} "

## Default parameters determined by the number of breeds
nginds=${nginds:=$(printf "%${nbreed}s" | sed "s/ /600 /g" | sed 's/ *$//')}
last_litters=${last_litters:=$(printf "%${nbreed}s" | sed "s/ /10 /g" | sed 's/ *$//')}
last_females=${last_females:=$(printf "%${nbreed}s" | sed "s/ /200 /g" | sed 's/ *$//')}
founder_sel=${founder_sel:=$(printf "%${nbreed}s" | sed "s/ /rnd,/g" | sed 's/,$//')}
seg_sel=${seg_sel:=${founder_sel}}
last_sel=${last_sel:=${founder_sel}}
extentLDs=${extentLDs:=$(printf "%${nbreed}s" | sed "s/ /10 /g" | sed 's/ *$//')}
last_males=${last_males:=$(printf "%${nbreed}s" | sed "s/ /40 /g" | sed 's/ *$//')}
seg_gens=${seg_gens:=$(printf "%${nbreed}s" | sed "s/ /40 /g" | sed 's/ *$//')}

## When simulating phenotypes, set the number of traits to 1 if there are no trait names
if [[ ${ntrait} == "0" ]]; then
  ntrait=1
  trait_array=('/')
  traits=/
fi

## Modify the job name
if [[ ${SLURM_CPUS_ON_NODE} ]]; then
  thread=${thread:=${SLURM_CPUS_ON_NODE}}

  ## Modify the job name
  scontrol update jobid=${SLURM_JOB_ID} name="${type}"
fi

## Number of parallel jobs
[[ ! ${thread} && ${rep} && ${fold} ]] && thread=$((rep * fold))

## Change working directory
cd ${proj} || exit 5

## Extract breed genotype information from merged group (requires map/ped files to generate simulated phenotypes)
if [[ ${bfile} ]]; then
  check_plink "${bfile}" ${nchr}
  unset bfiles
  for b in "${breeds_array[@]}"; do
    echo ${b} >fid.txt
    [[ ! -s ${b}m.map ]] &&
      plink --bfile ${bfile} --keep-fam fid.txt --chr-set ${nchr} --freq --recode --out ${b}m &>/dev/null
    bfiles="${bfiles} ${proj}/${b}m"
  done
  [[ -s fid.txt ]] && rm fid.txt
fi

## Log file
logf=${logp}/${type}_${SLURM_JOB_ID}.log

if [[ ${type} == "bin" ]]; then
  ## Interval division
  $block_define \
    --bfile ${bfile} \
    --win ${nsnp_win} \
    --maf ${maf} \
    --minSize ${nsnp_sim} \
    --type ${bin} \
    --out ${out} >${logp}/${bin}_block_${SLURM_JOB_ID}.log
    [[ ! -s ${out} ]] && echo "error in bin definition! " && exit 1
elif [[ ${type} == "adj" ]]; then
  ## Calculate adjusted phenotypes (ebv+re)
  for ti in "${trait_array[@]}"; do
    ## Get the index of the trait to be evaluated among all traits
    pi=$(echo "${traits}" | tr ' ' '\n' | grep -n -w -m1 "${ti}" | cut -d':' -f1)
    phedir=${proj}/${ti}

    ## Modify job name
    [[ ${SLURM_CPUS_ON_NODE} ]] && scontrol update jobid=${SLURM_JOB_ID} name="${type}_${ti}"

    ## Folder for analyzing the specified phenotype
    mkdir -p ${phedir}
    cd "${phedir}" || exit

    [[ -s ${phedir}/phe_adj_BLUP.txt ]] && \
      echo "Warning: phe_adj_BLUP.txt exists! "

    for b in "${breeds_array[@]}"; do
      ## Folder for analyzing the specified breed's phenotype
      mkdir -p ${phedir}/${b}
      cd ${phedir}/${b} || exit

      if [[ ! -s ${b}_dmu.txt ]]; then
        ## Check if the genotype file exists and is in binary format
        check_plink "${bfile}" ${nchr}

        ## Check if the family id of the breed is present in the plink file
        if [[ $(grep -c "${b}" ${bfile}.fam) -eq 0 ]]; then
          echo "no family id ${b} in ${bfile}.fam file! "
          exit 1
        fi

        ## Extract the specified breed's genotype (might need to modify the genotype file, so it's also a copy of the genotype)
        echo ${b} >tmp_fid.txt
        plink --bfile ${bfile} --keep-fam tmp_fid.txt --chr-set ${nchr} --make-bed --out ${b} &>>${logf}
        rm tmp_fid.txt

        ## Breed's corresponding phenotype file
        awk 'FNR==NR{a[$2];next} $1 in a' ${b}.fam ${phef} >${b}_dmu.txt
        # grep ${b} ${phef} | awk '{$NF="";print}' > ${b}_dmu.txt # Last column of phenotype file is the breed

        ## Calculate adjusted phenotypes
        $phe_adj \
          --phereal ${pi} \
          --bfile ${b} \
          --DIR phe_adj_BLUP \
          --phef ${b}_dmu.txt \
          --all_eff "${all_eff}" \
          --ran_eff "${ran_eff}" \
          --out ${phedir}/phe_adj_BLUP.txt \
          --append &>>${logf}
      fi
    done
  done
elif [[ ${type} == "gsim" ]]; then
  for dir in "${sim_dirs[@]}"; do
    ## Modify job name
    [[ ${SLURM_CPUS_ON_NODE} ]] && scontrol update jobid=${SLURM_JOB_ID} name="${type}_${dir}"

    $run_QMSim \
      --proj "${proj}" \
      --breeds "${breeds}" \
      --thread "${thread}" \
      --sim_dir "${dir}" \
      --geno_gen "${geno_gen}" \
      --nmloc "${nmloc}" \
      --nchr "${nchr}" \
      --last_females "${last_females}" \
      --last_litters "${last_litters}" \
      --extentLDs "${extentLDs}" \
      --seg_gens "${seg_gens}" \
      --last_males "${last_males}" \
      --founder_sel "${founder_sel}" \
      --seg_sel "${seg_sel}" \
      --last_sel "${last_sel}" \
      --logf "${logf}"
  done
elif [[ ${type} == "geno" ]]; then
  $genome_process \
    --proj "${proj}" \
    --breeds "${breeds}" \
    --geno_gen "${geno_gen}" \
    --last_females "${last_females}" \
    --last_litters "${last_litters}" \
    --nginds "${nginds}" \
    --binDiv ${binDiv} \
    --binThr ${binThr} \
    --maf ${maf} \
    --nsnp ${nsnp}
elif [[ ${type} == "psim" ]]; then
  ## Simulate phenotypes
  if [[ ! -s pheno_sim.txt ]]; then
    $pheno_sim \
      --h2 "${h2s}" \
      --mean "${means}" \
      --rg "${rg_sim}" \
      --gt "${bfiles}" \
      --bin ${bin_sim} \
      --win ${nsnp_sim} \
      --nqtl ${nqtl} \
      --nbin_cor ${nbin_cor} \
      --nsnp_cor ${nsnp_cor} \
      --dist_cor ${rg_dist} \
      --seed ${seed} \
      --fid \
      --min ${min} \
      ${overlap} \
      ${evenly} \
      ${binf} \
      --qtlf qtl_info.txt \
      --out pheno_sim.txt &>>${logf}
    [[ ! -s pheno_sim.txt ]] && echo "phenotypes simulation error! " && exit 1

    ## Remove QTL
    if [[ ! ${overlap} ]]; then
      awk '{print $2}' qtl_info.txt >qtl_snpid.txt
    else 
      [[ -s qtl_snpid.txt ]] && rm qtl_snpid.txt
    fi
  fi
elif [[ ${type} == "within" ]]; then
  ## Within-population prediction
  for ti in "${trait_array[@]}"; do
    ## Get the index of the trait to be evaluated among all traits
    pi=$(echo "${traits}" | tr ' ' '\n' | grep -n -w -m1 "${ti}" | cut -d':' -f1)
    phedir=${proj}/${ti}

    ## File containing true breeding values
    if [[ -s ${phedir}/phe_adj_BLUP.txt ]]; then
      tbvf=${phedir}/phe_adj_BLUP.txt
    else
      tbvf=${phef}
    fi

    for b in "${breeds_array[@]}"; do
      ## Modify job name
      [[ ${SLURM_CPUS_ON_NODE} ]] && scontrol update jobid=${SLURM_JOB_ID} name="${type}_${ti}_${b}"

      ## Working folder
      mkdir -p ${phedir}/${b}
      cd ${phedir}/${b} || exit 5

      ## Genotype file
      if [[ -s ${phedir}/qtl_snpid.txt ]]; then
        ## Remove QTL from genotype file
        plink \
          --file ${phedir}/${b}m \
          --exclude ${phedir}/qtl_snpid.txt \
          --make-bed \
          --out ${phedir}/${b}mq &>>${logf}
        bfile=${phedir}/${b}mq
      else
        bfile=${proj}/${b}m
      fi

      ## Phenotype file
      if [[ -s ${phedir}/pheno_sim.txt ]]; then
        grep "${b}" ../pheno_sim.txt | awk '{print $2="", $0}' | awk '$2="1"' >${b}_dmu_pheno.txt
        phef=${b}_dmu_pheno.txt
      fi

      ## Accuracy of DMU calculation
      $GP_single \
        --label ${b} \
        --phef ${phef} \
        --bfile ${bfile} \
        --method ${method} \
        --seed ${seed} \
        --all_eff "${all_eff}" \
        --ran_eff "${ran_eff}" \
        --rep ${rep} \
        --fold ${fold} \
        --phereal ${pi} \
        --thread ${thread} \
        --tbvf ${tbvf} \
        --iter ${iter} \
        --burnin ${burnin} \
        ${pedf} \
        ${dense} \
        ${debug} \
        ${gen} \
        ${binf} \
        ${tbv_col} \
        --out ${out} &>>${logf}
    done
  done
elif [[ ${type} == "single" || ${type} == "multi" ]]; then
  for ti in "${trait_array[@]}"; do
    ## Get the index of the trait to be evaluated among all traits
    pi=$(echo "${traits}" | tr ' ' '\n' | grep -n -w -m1 "${ti}" | cut -d':' -f1)
    phedir=${proj}/${ti}

    ## Switch to working folder
    cd ${phedir} || exit 5

    ## Modify job name
    # mtype=$(basename "$(dirname "$phedir")")_$(basename "$phedir" | tr '/' '_')
    if [[ ${SLURM_CPUS_ON_NODE} ]]; then
      if [[ ${method} == "mbBayesAB" ]]; then
        name="${type}_${method}_${ti}_${bin}"
      else
        name="${type}_${method}_${ti}"
      fi
      scontrol update jobid=${SLURM_JOB_ID} name="${name}"
    fi

    ## File containing true breeding values
    if [[ -s ${phedir}/phe_adj_BLUP.txt ]]; then
      ## Adjusted phenotypes
      tbvf=${phedir}/phe_adj_BLUP.txt
    else
      ## Phenotype file contains "true" breeding values
      tbvf=${phef}
    fi

    ## Check if all breeds in the combination have completed within-breed prediction, otherwise skip
    for breed in $breeds; do
      [[ ! -d ${phedir}/${breed} ]] && echo "${phedir}/${breed} not found." && not_run=true && break
    done
    [[ ${not_run} ]] && unset not_run && continue

    $GP_multi \
      --pops "${breeds}" \
      --type ${type} \
      --method ${method} \
      --tbvf ${tbvf} \
      --phereal ${pi} \
      --nsnp_win ${nsnp_win} \
      --thread ${thread} \
      --iter ${iter} \
      --burnin ${burnin} \
      --ref ${ref} \
      --seed ${seed} \
      --bin ${bin} \
      ${prior} \
      ${dirPre} \
      ${rg_local} \
      ${binf} \
      ${rg} \
      ${suffix} \
      ${nbin} \
      ${tbv_col} \
      ${noCov} \
      ${dense} \
      ${debug} &>>${logf}
  done
elif [[ ${type} == "accur" ]]; then
  ## Accuracy statistics
  $accuracy_summ \
    --proj "${proj}" \
    --breeds "${breeds}" \
    --traits "${traits}" \
    --cor "${rg_sim}" \
    --rep "${rep}" \
    --dist "${rg_dist}" \
    ${out} &>>${logf}
elif [[ ${type} == "var" ]]; then
  ## Genetic parameter statistics
  $varComp_summ \
    --proj "${proj}" \
    --breeds "${breeds}" \
    --traits "${traits}" \
    --bin "${bin}" \
    --code "${code}" \
    --cor "${rg_sim}" \
    --rep "${rep}" \
    --dist "${rg_dist}" \
    ${out} &>>${logf}
elif [[ ${type} == "time" ]]; then
  ## Genetic parameter statistics
  $time_summ \
    --proj "${proj}" \
    --breeds "${breeds}" \
    --traits "${traits}" \
    --bin "${bin}" \
    --code "${code}" \
    --cor "${rg_sim}" \
    --rep "${rep}" \
    --dist "${rg_dist}" \
    ${out} &>>${logf}
fi

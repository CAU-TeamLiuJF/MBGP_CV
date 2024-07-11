#!/bin/usr/bash

########################################################################################################################
## Version:   1.2.0
## Author:    Liweining liwn@cau.edu.cn
## Orcid:     0000-0002-0578-3812
## Institute: College of Animal Science and Technology, China Agricul-tural University, Haidian, 100193, Beijing, China
## Date:      2024-07-07
##
## Function：
##  Used to reproduce the research results in the article "Multi-trait Bayesian models enhance the accuracy of genomic 
## prediction in multi-breed reference populations"
##
##
## Data sources ：
##  Xie, L., J. Qin, L. Rao, X. Tang and D. Cui et al., 2021 Accurate prediction and genome-wide 
##  association analysis of digital intramuscular fat content in longissimus muscle of pigs. Animal Genetics 
##  52: 633-644. https://doi.org/10.1111/age.13121
##
## [KELLER B](https://doi.org/10.3389/fpls.2022.830896), ARIZA-SUAREZ D, PORTILLA-BENAVIDES A E, et al. 
## Improving Association Studies and Genomic Predictions for Climbing Beans With Data From Bush Bean 
## Populations[J]. Frontiers in Plant Science, 2022,13
##
## License:
##  This script is licensed under the GPL-3.0 License.
##  See https://www.gnu.org/licenses/gpl-3.0.en.html for details.
########################################################################################################################



## Initialize the script, modify some paths in the script (path need to be modefied)
# /public/home/liujf/liwn/code/GitHub/mbBayesABLDLD/initialize.sh

## Path of main script
code=/public/home/liujf/liwn/code/GitHub/mbBayesABLD
GP_cross=${code}/shell/GP_cross_validation.sh
combination=${code}/R/array_combination.R

## 将程序路径加到环境变量中
export PATH=${code}/bin:$PATH


###########################################################################################
## 
##  Real date
## 
###########################################################################################

for source in Xie2021 Keller2022; do
  ## path of project
  pro=${code}/output/${source}
  mkdir -p ${pro}/log
  cd ${pro} || exit

  ## files
  data_path=${code}/data/Real/${source}
  seed=$(cat "${data_path}/random.seed")
  bfile=${data_path}/genotype
  phef=${data_path}/phenotype.txt
  IFS=" " read -r -a breeds <<< "$(cat ${data_path}/breeds.txt)"
  IFS=" " read -r -a traits_all <<< "$(cat ${data_path}/traits_all.txt)"
  IFS=" " read -r -a traits_cal <<< "$(cat ${data_path}/traits_analysis.txt)"

  ## Calculate corrected phenotype
  if [[ ${source} == "Xie2021" ]]; then
    ## effect setting
    all_eff="3 2 1"

    for trait in "${traits_cal[@]}"; do
      ## Parameters need to be adjusted according to the characteristics and requirements of Linux servers
      ## Note: If the Slurm workload manager is not installed in the system, please comment out the line where 'sbatch' is located in whole script
      sbatch -c1 --mem=4G --output=${pro}/log/phe_adj.log \
        $GP_cross \
        --proj ${pro} \
        --breeds "${breeds[*]}" \
        --traits "${traits_all[*]}" \
        --trait "${trait}" \
        --bfile "${bfile}" \
        --phef "${phef}" \
        --all_eff "${all_eff}" \
        --code "${code}" \
        --type adj
      sleep 5
    done
  else
    ## The phenotype from keller are BLUEs, and there is no need to calculate corrected phenotype
    tbv_col=" --tbv_col same "

    ## effect setting
    all_eff="2 1"
  fi

  ## Wait for the correction phenotype calculation to be completed
  while [[ $(wc -l 2>/dev/null <${pro}/${traits_all[-1]}/${breeds[-1]}/phe_adj_BLUP.SOL) -lt 10 ]]; do
    sleep 3
  done

  ## Calculate the accuracy of GBLUP model in within-breed prediction
  for trait in "${traits_cal[@]}"; do
    sbatch -c50 --mem=100G --output=${pro}/log/within.log \
      $GP_cross \
      --proj ${pro} \
      --breeds "${breeds[*]}" \
      --traits "${traits_all[*]}" \
      --trait "${trait}" \
      --bfile "${bfile}" \
      --phef "${phef}" \
      --seed "${seed}" \
      --all_eff "${all_eff}" \
      ${tbv_col} \
      --code "${code}" \
      --rep 10 \
      --fold 5 \
      --type within
    sleep 5
  done

  ## Calculate the accuracy of GBLUP model in multi-breed prediction
  for trait in "${traits_cal[@]}"; do
    for type in blend union; do
      while IFS= read -r breed_comb; do
        sbatch -c50 --mem=80G --output=${pro}/log/${type}.log \
          $GP_cross \
          --type "${type}" \
          --proj "${pro}" \
          --breeds "${breed_comb}" \
          --traits "${traits_all[*]}" \
          --trait "${trait}" \
          --code "${code}" \
          --phef "${phef}" \
          ${tbv_col} \
          --seed ${seed} \
          --thread 50 \
          --suffix
        sleep 5
      done <"${data_path}/breeds_combination.txt"
    done
  done

  ## Calculate the accuracy of Bayes multi-breed prediction
  for trait in "${traits_cal[@]}"; do
    while IFS= read -r breed_comb; do
      for bin in fix lava cubic; do
        sbatch -c50 --mem=80G --output=${pro}/log/${bin}.log \
          $GP_cross \
          --type multi \
          --proj ${pro} \
          --breeds "${breed_comb}" \
          --traits "${traits_all[*]}" \
          --trait "${trait}" \
          --code ${code} \
          ${tbv_col} \
          --seed ${seed} \
          --phef ${phef} \
          --thread 50 \
          --bin ${bin} \
          --suffix
        sleep 10
      done
    done <"${data_path}/breeds_combination.txt"
  done

  ## Calculate the accuracy of Bayesian model in within-breed prediction
  for trait in "${traits_cal[@]}"; do
    binf=$(find ${pro} -name "cubic_*.txt" -type f -print -quit)
    [[ ! -f ${binf} ]] && continue
    for b in "${breeds[@]}"; do
      sbatch -c50 --mem=100G --output=${pro}/log/within_bayes.log \
        $GP_cross \
        --proj ${pro} \
        --breeds "${b}" \
        --traits "${traits_all[*]}" \
        --trait "${trait}" \
        --bfile ${bfile} \
        --phef ${phef} \
        --seed ${seed} \
        --method BayesAS \
        --binf ${binf} \
        ${tbv_col} \
        --code ${code} \
        --rep 10 \
        --fold 5 \
        --type within \
        --out accur_BayesAS.txt
      sleep 5
    done
  done

  ## Statistical Accuracy and Variance Component Results
  for type in accur var; do
    $GP_cross \
      --type ${type} \
      --proj ${pro} \
      --breeds "${breeds[*]}" \
      --bin "fix lava cubic" \
      --traits "${traits_cal[*]}" \
      --code ${code}
  done
done


###########################################################################################
## 
##  Simulation of QMSim
## 
###########################################################################################

for scene in Two Three; do # scene=Three;r=1;dist=identical;cor=0.2;bin=lava
  for r in {1..20}; do
    data_path=${code}/data/Simulation/${scene}Breeds
    pro=${code}/output/${scene}Breeds
    proi=${code}/output/${scene}Breeds/rep${r}
    mkdir -p ${pro}
    cd ${pro} || exit

    ## Generating simulated populations using QMSim
    sed "s#%seed_file%#${data_path}/QMSim_rep${r}_seed.prv#g" ${data_path}/QMSim.prm >QMSim_rep${r}.prm
    sed -i "s#%nthread%#20#g" QMSim_rep${r}.prm
    srun -c20 --mem=100G \
      QMSim QMSim_rep${r}.prm

    ## Load parameters
    source ${data_path}/parameter.sh

    ## Screening individuals and SNP markers from simulated populations
    $GP_cross \
      --type geno \
      --code ${code} \
      --proj "${pro}/rep${r}" \
      --breeds "${breed_sim}" \
      --last_females "${last_females}" \
      --nginds "${nginds}" \
      --binDiv "pos" \
      --maf "0.01" \
      --binThr "10" \
      --geno_gen "8-10" \
      --nsnp "50000" \
      --out "merge"

    ## Random number seed
    seed=$(cat ${pro}/rep${r}/random.seed)

    ## genome partitioning
    if [[ ${scene} == "Two" ]]; then
      phe_geno="${pro}/rep${r}/merge"
      phe_bin="${pro}/rep${r}/cubic_M_50_psim.txt"
    else
      plink \
        --bfile ${pro}/rep${r}/Am \
        -bmerge ${pro}/rep${r}/Bm \
        --make-bed \
        --out ${pro}/rep${r}/ABm
      phe_geno="${pro}/rep${r}/ABm"
      phe_bin="${pro}/rep${r}/cubic_AB_50_psim.txt"
    fi
    $GP_cross \
      --type bin \
      --proj ${proi} \
      --bfile ${pro}/rep${r}/merge \
      --bin "cubic" \
      --nsnp_win "50" \
      --nsnp_sim "100" \
      --out "${phe_bin}"

    for dist in uniform; do
      for cor in ${cors}; do
        proi=${pro}/rep${r}/${dist}/cor${cor}
        mkdir -m 777 -p ${proi}

        ## Simulate phenotype
        $GP_cross \
          --type psim \
          --proj ${proi} \
          --bfile ${pro}/rep${r}/merge \
          --breeds "${breed_sim}" \
          --code ${code} \
          --means "${means}" \
          --h2s "${h2s}" \
          --rg_sim "${cor}" \
          --rg_dist ${dist} \
          --nqtl "400" \
          --nsnp_cor "10" \
          --nbin_cor "10" \
          --evenly \
          --min "30" \
          --binf ${pro}/rep${r}/cubic_AB_50_psim.txt \
          --seed ${seed}

        ## Calculate the accuracy of GBLUP model in within-breed prediction
        for b in ${breed_sim}; do
          sbatch -c50 --mem=100G --output=${pro}/rep${r}/log/within.log  \
            $GP_cross \
              --type within \
              --proj ${proi} \
              --breeds "${b}" \
              --bfile ${pro}/rep${r}/merge \
              --phef ${proi}/pheno_sim.txt \
              --code ${code} \
              --thread 50 \
              --tbv_col 6 \
              --seed ${seed} \
              --rep 5 \
              --fold 5
          sleep 5
        done

        ## wait for the GBLUP accuracy calculation to be completed
        while [[ $(find ${proi}/*/val5/rep5/pheno.txt 2>/dev/null | wc -l) -lt ${np} ]]; do
          sleep 3
        done

        ## Calculate the accuracy of GBLUP model in multi-breed prediction
        for type in blend union; do
          [[  -d ${proi}/${type}_${breed_sim// /_} ]] && continue
          sbatch -c50 --mem=100G $GP_cross \
            --type ${type} \
            --proj ${proi} \
            --breeds "${breed_sim}" \
            --bfile ${pro}/rep${r}/merge \
            --phef ${proi}/pheno_sim.txt \
            --code ${code} \
            --suffix \
            --thread 50 \
            --seed ${seed} \
            --tbv_col 6
          sleep 5
        done

        ## Calculate the accuracy of Bayes multi-breed prediction
        for bin in fix lava cubic; do
          [[ -s ${proi}/multi_${breed_sim// /_}/val1/rep1/EBV_${bin}_y1.txt ]] && continue
          sbatch -c50 --mem=100G $GP_cross \
            --type multi \
            --proj ${proi} \
            --breeds "${breed_sim}" \
            --bfile ${pro}/rep${r}/merge \
            --phef ${proi}/pheno_sim.txt \
            --thread 50 \
            --code ${code} \
            --seed ${seed} \
            --suffix \
            --tbv_col 6 \
            --bin ${bin}
          sleep 30
        done

        ## Calculate the accuracy of Bayes within-breed prediction
        binf=${proi}/multi_${breed_sim// /_}/cubic_M_50.txt

        ## 等待品种内验证群划分完毕
        while [[ ! -s ${binf} ]]; do
          sleep 3
        done

        for b in ${breed_sim}; do
          [[ -s "${proi}/${b}/accur_BayesAS.txt" ]] && continue
          sbatch -c50 --mem=100G $GP_cross \
            --type within \
            --proj ${proi} \
            --breeds "${b}" \
            --bfile ${pro}/rep${r}/merge \
            --phef ${proi}/pheno_sim.txt \
            --code ${code} \
            --binf ${binf} \
            --method BayesAS \
            --thread 50 \
            --tbv_col 6 \
            --seed ${seed} \
            --rep 5 \
            --fold 5
          sleep 5
        done
      done
    done
  done

  ## Statistical Accuracy and Variance Component Results
  for type in accur var; do # type=accur $(seq -s " " 1 10)
    $GP_cross \
      --type ${type} \
      --proj ${pro} \
      --rep "$(seq -s " " 1 20)" \
      --rg_dist "identical uniform" \
      --rg_sim "${cors}" \
      --bin "fix lava cubic" \
      --breeds "${breed_sim}" \
      --code ${code} \
      --out accuracy.txt
  done
done

#!/bin/bash
#SBATCH --job-name=accur_GBLUP

########################################################################################################################
## 版本: 1.1.1
## 作者: 李伟宁 liwn@cau.edu.cn
## 日期: 2023-07-05
## 
## 统计各种情形下的方差组分结果
## 
## 使用: ./varcomp_summary.sh --help
## 
## License:
##  This script is licensed under the GPL-3.0 License.
##  See https://www.gnu.org/licenses/gpl-3.0.en.html for details.
########################################################################################################################


###################  参数处理  #####################
####################################################
## NOTE: This requires GNU getopt.  On Mac OS X and FreeBSD, you have to install this
## 参数名
TEMP=$(getopt -o h --long code:,proj:,breeds:,rep:,dist:,cor:,traits:,bin:,dirPre:,out:,help \
              -n 'javawrap' -- "$@")
if [ $? != 0 ] ; then echo "Terminating..." >&2 ; exit 1 ; fi
eval set -- "$TEMP"
## 解析参数
while true; do
  case "$1" in
    --proj )     proj="$2";     shift 2 ;; ## 项目目录 [必要参数]
    --breeds )   breeds="$2";   shift 2 ;; ## 群体/品种标识符，如'YY DD' [必要参数]
    --traits )   traits="$2";   shift 2 ;; ## 性状名称，如"DF DPM" ["/"]
    --rep )      rep="$2";      shift 2 ;; ## 第几次重复 ["/"]
    --dist )     dist="$2";     shift 2 ;; ## 加性遗传相关服从的分布 ["/"]
    --cor )      cor="$2";      shift 2 ;; ## 加性遗传相关大小 ["/"]
    --dirPre )   dirPre="$2";   shift 2 ;; ## ebv输出文件夹增加的前缀 [""]
    --bin )      bins="$2";     shift 2 ;; ## 多品种评估时区间划分方法，fix/frq/ld/lava/cubic ["fix lava cubic"]
    --code )     code="$2";     shift 2 ;; ## 脚本文件所在目录，如/BIGDATA2/cau_jfliu_2/liwn/code [NULL]
    --out )      out="$2";      shift 2 ;; ## 准确性输出文件名 [accuracy_$date.txt]
  -h | --help)    grep ";; ##" $0 | grep -v help && exit 1 ;;
  -- ) shift; break ;;
  * ) shift; break ;;
  esac
done

## 检查必要参数是否提供
if [[ ! -d ${proj} ]]; then
  echo "${proj} not found! "
  exit 1
elif [[ ! ${breeds} ]]; then
  echo "para --breeds is reduired! "
  exit 1
fi

## 日期
today=$(date +%Y%m%d)

## 默认参数
out=${out:=${proj}/varcomp_${today}.txt}
bins=${bins:="fix lava cubic"}
dirPre=${dirPre:=""}
traits=${traits:="/"}
rep=${rep:="/"}
dist=${dist:="/"}
cor=${cor:="/"}

## 避免执行R脚本时的警告("ignoring environment value of R_HOME")
unset R_HOME

## DMU参数卡名称
DIR_full=phe_adj_BLUP
DIR_within=within
DIR_blend=blend
DIR_union=union

## 解析参数
read -ra breeds_array <<<"$breeds"
read -ra bins_array <<<"$bins"
read -ra traits_array <<<"$traits"
read -ra reps_array <<<"$rep"
read -ra dists_array <<<"$dist"
read -ra cors_array <<<"$cor"

##############  方差组分整理  ##############
###########################################
echo "模拟重复 相关分布 相关大小 模型 参考群 品种 性状 交叉重复 交叉折数 类型 值" >${out}
for p in "${traits_array[@]}"; do # p=${traits_array[0]};re=${reps_array[0]};d=${dists_array[0]};c=${cors_array[0]}
  for re in "${reps_array[@]}"; do # b=${breeds_array[0]};r=1;f=1;bin=fix
    for d in "${dists_array[@]}"; do
      for c in "${cors_array[@]}"; do
        path=${proj}/${p}

        ## 模拟情形下的路径设置
        [[ ${re} != "/" ]] && path=${path}/rep${re}
        [[ ${d} != "/" ]] && path=${path}/${d}
        [[ ${c} != "/" ]] && path=${path}/cor${c}

        ## 处理路径中存在多个斜杠的情况
        path=$(echo "$path" | sed 's#/\{2,\}#/#g; s#/$##')

        ## 判断文件夹是否存在
        [[ ! -d ${path} ]] && continue

        ## 品种内遗传力
        for b in "${breeds_array[@]}"; do
          ## 交叉验证参数
          rep=$(find ${path}/${b}/val1 -name "rep*" -type d | wc -l)
          fold=$(find ${path}/${b}/val* -name "rep1" -type d | wc -l)

          ## 每个子集方差组分
          for r in $(seq 1 ${rep}); do # r=1;f=1
            for f in $(seq 1 ${fold}); do
                ## GBLUP
                lst=${path}/${b}/val${f}/rep${r}/${DIR_within}.lst
                if [[ -s ${lst} ]]; then
                  h2_within=$(grep -A2 'Trait  correlation' ${lst} | tail -n 1 | awk '{print $2}')
                  # echo "${re} ${d} ${c} w-GBLUP ${b} ${b} ${p} ${r} ${f} h2 ${h2_within}" >>${out}
                fi

                ## Bayes
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

        ## 品种组合
        mapfile -t blends < <(find ${path} -name "blend*" -type d 2>/dev/null)

        for t in "${blends[@]}"; do # t=${blends[0]}
          paths=$(basename ${t})
          types=${paths/blend_/}
          IFS='_' read -r -a breeds_sub <<<"$types"
          n=${#breeds_sub[@]}

          ## types
          if [[ ${types} == "blend" ]]; then
            types=""
            comb=${breeds// /_}
          else
            types="_${types}"
            comb=$(IFS="_"; echo "${breeds_sub[*]}")
          fi

          ## 组合中每个品种
          for i in "${!breeds_sub[@]}"; do # i=0
            ## 子集
            for r in $(seq 1 ${rep}); do # r=5;f=5
              for f in $(seq 1 ${fold}); do
                ## 单性状GBLUP联合评估
                lst=${path}/blend${types}/val${f}/rep${r}/${DIR_blend}.lst
                if [[ -s ${lst} ]]; then
                  h2_blend=$(grep -A2 'Trait  correlation' ${lst} | tail -n 1 | awk '{print $2}')
                # else
                #   h2_blend=""
                #   echo "${lst} not found! "
                fi

                ## 双性状联合评估
                lst=${path}/union${types}/val${f}/rep${r}/${DIR_union}.lst
                if [[ -s ${lst} ]]; then
                  message="Correlation matrix for random"
                  h2_union=$(grep -A $((i + 2)) 'Trait  correlation' ${lst} | tail -n 1 | awk '{print $2}')

                  for j in $(seq $((i + 1)) $((n - 1))); do
                    rg_union=$(grep -A $((j + 2)) "${message}" ${lst} | tail -n 1 | awk -v col=$((i + 2)) '{print $col}')
                    echo "${re} ${d} ${c} MT-GBUP ${comb} ${breeds_sub[i]}_${breeds_sub[j]} ${p} ${r} ${f} rg ${rg_union}" >>${out}
                  done
                # else
                #   h2_union=""
                #   echo "${lst} not found! "
                fi

                for bin in "${bins_array[@]}"; do
                  # for dirPre in /; do
                    multi_path=${path}/multi${types}/val${f}/rep${r}

                    ebvfi=$(find ${multi_path} -name "EBV_${bin}_y$((i + 1)).txt" 2>/dev/null)
                    varfi=$(find ${multi_path} -name "var_${bin}*.txt" 2>/dev/null)
                    [[ ! -s ${ebvfi} || ! -s ${varfi} ]] && continue

                    ## 遗传力
                    varg=$(sed '1d' ${ebvfi} | awk '{sum+=$2; sumsq+=($2)^2} END {print (sumsq/NR-(sum/NR)^2)}')
                    vare=$(tail -n 1 ${varfi} | awk -v col=$((i * n + i + 1)) '{print $col}')
                    h2_multi=$(echo "scale=4; $varg / ($varg + $vare)" | bc)
                    echo "${re} ${d} ${c} ${bin} ${comb} ${breeds_sub[i]} ${p} ${r} ${f} h2 ${h2_multi}" >>${out}

                    ## 遗传相关
                    for j in $(seq $((i + 1)) $((n - 1))); do
                      ebvfj=$(find ${multi_path} -name "EBV_${bin}_y$((j + 1)).txt" 2>/dev/null)
                      [[ ! -s ${ebvfj} ]] && continue

                      # 计算A文件中第2列减去均值的值
                      meanA=$(awk 'NR>1{sum+=$2}END{print sum/(NR-1)}' ${ebvfi})
                      awk -v meanA="$meanA" 'NR>1{print $2-meanA}' ${ebvfi} >A.tmp
                      # 计算B文件中第2列减去均值的值
                      meanB=$(awk 'NR>1{sum+=$2}END{print sum/(NR-1)}' ${ebvfj})
                      awk -v meanB="$meanB" 'NR>1{print $2-meanB}' ${ebvfj} >B.tmp
                      # 计算协方差
                      cov=$(paste A.tmp B.tmp | awk '{sum+=($1*$2)}END{print sum/(NR - 1)}')
                      var2=$(sed '1d' ${ebvfj} | awk '{sum+=$2; sumsq+=($2)^2} END {print (sumsq/NR-(sum/NR)^2)}')

                      ## 计算遗传相关
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

                ## 写出到文件
                {
                  echo "${re} ${d} ${c} ST-GBUP ${comb} ${breeds_sub[i]} ${p} ${r} ${f} h2 ${h2_blend}"
                  echo "${re} ${d} ${c} MT-GBUP ${comb} ${breeds_sub[i]} ${p} ${r} ${f} h2 ${h2_union}"
                } >>${out}
              done
            done
          done
        done
      done
    done
  done
done

## 去除na值行
sed -i '/ $/d' ${out}
sed -i 's/ /\t/g' ${out}

[[ -s A.tmp ]] && rm A.tmp B.tmp

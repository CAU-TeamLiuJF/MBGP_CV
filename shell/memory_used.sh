#!/usr/bin/bash

cd /public/home/liujf/liwn/code/GitHub/MBGP_CV/output
out_file="simulation_efficience.csv"

# 表头
echo "scenario,sim_rep,dist,cor,breed,command,val,rep,type,CPU%,Elapsed,MaxRSS,AvgRSS" > "$out_file"

# 只匹配 single_*.log 或 multi_*.txt
find /public/home/liujf/liwn/code/GitHub/MBGP_CV/output/*/rep*/*/cor*/log \
  -type f \( -name "single_*.log" -o -name "multi_*.log" \) | while read -r log; do
  awk -v OFS="," -v out="$out_file" '
    function output() {
      if (scenario != "") {
        print scenario, sim_rep, dist, cor, breed, command, val, rep, type, cpu, elapsed, maxrss, avgrss >> out
      }
    }

    /Command being timed:/ {
      output()
      scenario=sim_rep=dist=cor=breed=command=val=rep=type=cpu=elapsed=maxrss=avgrss=""

      match($0, /"([^"]+)"/, m)
      cmdline = m[1]

      # 提取 scenario / sim_rep / dist / cor / breed
      if (match(cmdline, /output\/([^\/]+)\/(rep[0-9]+)\/([^\/]+)\/(cor[0-9.]+)\/([^\/]+)/, p)) {
        scenario=p[1]
        sim_rep=p[2]
        dist=p[3]
        cor=p[4]
        breed_dir=p[5]
        if (match(breed_dir, /single_([^\/]+)/, b)) breed=b[1]
        else if (match(breed_dir, /multi_([^\/]+)/, b)) breed=b[1]
        else breed=breed_dir
      }

      # 根据命令分类
      if (cmdline ~ /run_dmuai single/) {
        command="single_dmuai"
        if (match(cmdline, /val([0-9]+)\/rep([0-9]+)/, m2)) { val=m2[1]; rep=m2[2] }
        type="-"
      }
      else if (cmdline ~ /run_dmuai multi/) {
        command="multi_dmuai"
        if (match(cmdline, /val([0-9]+)\/rep([0-9]+)/, m2)) { val=m2[1]; rep=m2[2] }
        type="-"
      }
      else if (cmdline ~ /BayesR\.sh/) {
        command="single_bayesR"
        if (match(cmdline, /val([0-9]+)\/rep([0-9]+)/, m2)) { val=m2[1]; rep=m2[2] }
        type="bayesR"
      }
      else if (cmdline ~ /mbBayesABLD/) {
        command="mbBayesABLD"
        if (match(cmdline, /val([0-9]+)\/rep([0-9]+)/, m2)) { val=m2[1]; rep=m2[2] }
        if (match(cmdline, /--binf [^ ]*\/([^\/]+)\.txt/, bt)) type=bt[1]
      }
    }

    /Percent of CPU this job got:/ {
      if (match($0, /[0-9]+%/, m)) cpu=m[0]
    }
    /Elapsed \(wall clock\) time/ {
      if (match($0, /[0-9]+:[0-9]+(:[0-9]+)?(\.[0-9]+)?/, tm)) elapsed=tm[0]
    }
    /Maximum resident set size/ { maxrss=$NF }
    /Average resident set size/ { avgrss=$NF }

    END { output() }
  ' "$log"
done

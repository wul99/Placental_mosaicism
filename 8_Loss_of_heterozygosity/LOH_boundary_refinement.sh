#!/bin/bash

set -euo pipefail

snv_phase_AD_in_region(){

    local family=${1}
    local region=${2} 
    local chr=${region%:*}

    awk 'BEGIN{OFS="\t"}NR==FNR{a[$1,$2]=$0; next} ($1,$2) in a {print a[$1,$2],$4,$5,$6,$7,$8,$9,$10,$11}' \
    <(bcftools view -r ${region} ${shapeit_phased_bcf}|grep -v "#"|cut -f 1-2,10-17 ) \
    <(bcftools view -r ${region} ${origin_vcf} |grep -v "#"|cut -f 1-2,10-17 | awk -F '[:]' 'BEGIN{OFS="\t"}{print $1,$2,$7,$12,$17,$22,$27,$32,$37}') \
    |awk '$3==$6 && $6==$7 && $7==$8 && $8==$9 && $9==$10 {print $0}' \
    |awk '{ for (i=1; i<=NF; i++) if (length($i) - length(gensub(/,/, "", "g", $i)) > 1) next } 1' 

}

count_F_AF_in_region(){

    paste \
    <(cut -f 1,2,3 ${1} | awk 'BEGIN{OFS="\t"}{if ($3=="0|1" || $3=="1|0") print $1,$2}') \
    <(cut -f 3,11 ${1} | awk -F '[\t,]' 'BEGIN{OFS="\t"}{if ($1=="0|1") print $2/($2+$3);else if ($1=="1|0") print $3/($2+$3)}') \
    <(cut -f 6,14 ${1} | awk -F '[\t,]' 'BEGIN{OFS="\t"}{if ($1=="0|1") print $2/($2+$3);else if ($1=="1|0") print $3/($2+$3)}') \
    <(cut -f 7,15 ${1} | awk -F '[\t,]' 'BEGIN{OFS="\t"}{if ($1=="0|1") print $2/($2+$3);else if ($1=="1|0") print $3/($2+$3)}') \
    <(cut -f 8,16 ${1} | awk -F '[\t,]' 'BEGIN{OFS="\t"}{if ($1=="0|1") print $2/($2+$3);else if ($1=="1|0") print $3/($2+$3)}') \
    <(cut -f 9,17 ${1} | awk -F '[\t,]' 'BEGIN{OFS="\t"}{if ($1=="0|1") print $2/($2+$3);else if ($1=="1|0") print $3/($2+$3)}') \
    <(cut -f 10,18 ${1} | awk -F '[\t,]' 'BEGIN{OFS="\t"}{if ($1=="0|1") print $2/($2+$3);else if ($1=="1|0") print $3/($2+$3)}') \
    | sed -E 's/(\t+)(\t)/\1\1\2/g;s/^\t/\t\t/g' 

}

chr=${region%:*}
pos=${region##*:}
start=${pos%-*}
end=${pos##*-}
length=$((end - start))
start2=$((start - length))
end2=$((end + length))
region2=${chr}:${start2}-${end2}
echo ${start2} ${end2} ${region2}
snv_phase_AD_in_region ${family} ${region2} > ${combine_vcf_tmp}
count_F_AF_in_region ${combine_vcf_tmp} > ${phased_input}

python LOH_boundary_refinement.py -r ${region} -i ${phased_input}

rm ${combine_vcf_tmp} 


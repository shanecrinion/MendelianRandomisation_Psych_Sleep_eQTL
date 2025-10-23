#!/bin/bash
snp_list=()
for output in $(awk -F "\"*,\"*" '{if (NR!=1) {print $7}}' /home/shane/Desktop/files/data/eqtl/results_PSYCHENCODE_GWAS.Bonf.csv)
do
    #echo $output
    pos=$(echo $output | awk -F'-' '{print $3}')
    chr=$(echo $pos | awk -F':' '{print $1}')
    chr=$(echo chr$chr)
    pos=$(echo $pos | awk -F':' '{print $2}')
    echo $pos
    rsid=$(echo $output | awk -F'-' '{print $2}')
    echo $chr
    if [[ ! " ${snp_list[*]} " =~ " ${rsid} " ]] 
    then
        echo $rsid
        zcat '/home/shane/Desktop/files/data/eqtl/psychencode/psychencode_Full_hg19_cis-eQTL.txt.header.gz' | awk -v var_chr=$chr -v var_pos=$pos '{if($2==var_chr && $10 <= var_pos+500000 && $10 >= var_pos-500000){print}}' >> region_$rsid.txt
#       zcat '/home/shane/Desktop/files/data/eqtl/psychencode/psychencode_Full_hg19_cis-eQTL.txt.header.gz' | awk -v var_chr=$2 -v var_pos=${10} "{if( var_chr==$chr && var_pos < $pos+500000 && var_pos > $pos-500000){print}}" >> region_$rsid.txt
#       zcat /home/shane/Desktop/files/data/eqtl/psychencode/psychencode_Full_hg19_cis-eQTL.txt.header.gz | awk '{if($2==$chr && $10>$pos+500000 && $10<$pos-500000){print}}' >> region_$rsid.txt
        snp_list[${#snp_list[@]}]=$rsid
        echo $rsid >> psychencode_unique.txt
      # done
    fi
done
# for OUTPUT in $(Linux-Or-Unix-Command-Here)
# do
#     command1 on $OUTPUT
#     command2 on $OUTPUT
#     commandN
# done
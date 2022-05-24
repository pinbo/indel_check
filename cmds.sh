samtools tview -p chr1A_part1:345786964 170_S2_L002.proms.sorted.bam /group/dubcovskygrp6/references/iwgsc_refseqv1.0_all_chromosomes/161010_Chinese_Spring_v1.0_pseudomolecules_parts.ABonly.fasta

###### steps to check big indels and inversions
# step 1
cd /group/dubcovskygrp5/projects/sequence/promoter_cap/mapped/maps/base_group1_group2/junli_test
sbatch test_editcall_by_bam.sh

# step 2:find out uniq variants
time awk 'FNR>1{aa[$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6]++; bb[$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6]=FILENAME"\t"$7}END{for (i in aa){if (aa[i]==1) print bb[i]"\t"i}}' out_* | awk '$2>2' > potienal_uniq_indels.txt
#echo -e "File\tCoverage\tChrom\tref_start\tref_end\talt\tsize\ttype" > potienal_uniq_indels_sorted.txt
sort -k1,1 -k3,3 -k4,4n potienal_uniq_indels.txt >> potienal_uniq_indels_sorted.txt
sed -i 's/out_//g; s/.txt/.bam/g' potienal_uniq_indels_sorted.txt


# with alt, there could be point mutations makes GA and AA differnt, but the same insertion of GA
time awk 'FNR>1{aa[$1"\t"$2"\t"$3]++; bb[$1"\t"$2"\t"$3]=FILENAME"\t"$0}END{for (i in aa){if (aa[i]==1) print bb[i]}}' out_* | awk '$8>2' > potential_uniq_indels.txt
sort -k1,1 -k2,2 -k3,3n potential_uniq_indels.txt >> potential_uniq_indels_sorted.txt
sed -i 's/out_//g; s/.txt/.bam/g' potential_uniq_indels_sorted.txt


#awk 'NR>1{print "samtools depth ../"$1" -r "$3":"$4"-"$4" >> depth.txt"}'  potienal_uniq_indels_sorted.txt > for_depth.txt
# here it takes too long
#time cat for_depth.txt | while read -r line; do eval $line; done

# try bed files
# make bed file for each bam
awk '{print $3"\t"$4-1"\t"$4 >"test_"$1".bed"}' potienal_uniq_indels_sorted.txt
awk '{aa[$1]++}END{for (i in aa) print "samtools depth ../"i" -b test_"i".bed > depth_"i".txt"}'  potienal_uniq_indels_sorted.txt > for_depth_cmd.txt

#time cat for_depth_cmd.txt | while read -r line; do eval $line; done
# still too long, about 4 min for each bam, so I have to use sbatch
#time head -1 for_depth_cmd.txt | while read -r line; do eval $line; done

# use sbatch
sbatch get_depth.sh


## add depth result
echo -e "File\tmutCov\tChrom\tref_start\tref_end\talt\tsize\ttype\ttotalCov" > potienal_uniq_indels_sorted_with_total_coverage.txt
time awk 'FNR==NR{aa["depth_"$1".txt\t"$3"\t"$4]=$0;next}{aa[FILENAME"\t"$1"\t"$2]=aa[FILENAME"\t"$1"\t"$2]"\t"$3}END{for (i in aa) print aa[i]}' potienal_uniq_indels_sorted.txt depth_*.txt >> potienal_uniq_indels_sorted_with_total_coverage.txt

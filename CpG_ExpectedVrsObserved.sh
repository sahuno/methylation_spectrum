#date
#intersect cpg sites
path_ref=/juno/work/greenbaum/users/ahunos/methyl_SPECTRUM/data/ref/CpG_GRCh37.bed
path_CpG_HeaderLess=/juno/work/greenbaum/users/ahunos/methyl_SPECTRUM/data/ref/CpG_GRCh37_noHeader
tail -n +2 ${path_ref} > $path_CpG_HeaderLess
# cat -A $path_CpG_HeaderLess | head
# awk '!($3="")' $path_CpG_HeaderLess

awk '$2 !~ /_/' $path_CpG_HeaderLess | cut -f2- - > ${path_CpG_HeaderLess}.chr1_toY.bed


#remove `chr` prefix in naming 
sed -ie 's/^chr//' ${path_CpG_HeaderLess}.chr1_toY.bed
wc -l ${path_CpG_HeaderLess}.chr1_toY.bed

bedtools intersect -a ${path_CpG_HeaderLess}.chr1_toY.bed \
-b /juno/work/greenbaum/users/ahunos/methyl_SPECTRUM/data/processed/methylSeq_Agilent_Twist.interval_list0based_1bpWindow_annotations_BSSeq.bed \
-c > CpG_GRCh37_counts_expected_observed.bed #complete overlaps and counts

less CpG_GRCh37_counts_expected_observed.bed
wc -l CpG_GRCh37_counts_expected_observed.bed
#awk '{print NF}' CpG_GRCh37_counts_expected_observed.bed #how many columns
awk '{if ($12 > 0) { print } }' CpG_GRCh37_counts_expected_observed.bed | head #how cpgs 
#awk '$12 >= 1' CpG_GRCh37_counts_expected_observed.bed | head

#clean dir
rm CpG_GRCh37_noHeader.chr1_toX
rm CpG_GRCh37_noHeader
rm CpG_GRCh37_noHeader.chr1_toX

#plot cpg observed Expected
#CpG_GRCh37_counts_expected_observed.bed

Rscript /juno/work/greenbaum/users/ahunos/methyl_SPECTRUM/scripts/Plot_CpGs_Hg19_counts_overlaps.r \
-i /juno/work/greenbaum/users/ahunos/methyl_SPECTRUM/data/processed/CpG_GRCh37_counts_expected_observed.bed \
--output_dir /juno/work/greenbaum/users/ahunos/methyl_SPECTRUM/figures \
-s "BSeq_all3_commonCpG_SitesTest"
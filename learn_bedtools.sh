################
###### bedtools intersection 
################
bedtools intersect -a gr1.bed -b gr2.bed #correct intersections are made but each intersection loci has 3 entries of the same loci. guess due to width, strand, gc, percentage nethylation  


bedtools intersect -a gr1.bed -b gr2.bed gr3.bed -f 1 -u #Yes, desired unique intersections are produced with any metadata but metadata is carried on from 1st genomic interval option `-a` stated in the command
bedtools multiinter -header -names gr1 gr2 gr3 -i gr1.bed gr2.bed gr3.bed > collapsed_multiIntersect.bed
cat collapsed_multiIntersect.bed
awk -F'\t' 'BEGIN{FS=OFS="\t"}{print $0, $3-$2}' collapsed_multiIntersect.bed # append width of granges to file, but "wdith" header missing
awk -F'\t' 'NR==1{$(NF+1)="Width"}NR>1{$(NF+1)=$3-$2}1' collapsed_multiIntersect.bed # append with of granges to file, but feild sepearators are really bad
bedtools intersect -a gr1.bed -b gr2.bed gr3.bed -f 1 -c #complete overlaps and counts

#now use only start and end
cut -f 1-3 gr1.bed > gr1_chrStartEndOnly.bed
cut -f 1-3 gr2.bed > gr2_chrStartEndOnly.bed
cut -f 1-3 gr3.bed > gr3_chrStartEndOnly.bed

bedtools intersect -a gr1_chrStartEndOnly.bed -b gr2_chrStartEndOnly.bed -u #unique but not necessarily complete overlaps
bedtools intersect -a gr1_chrStartEndOnly.bed -b gr2_chrStartEndOnly.bed -f 1 #complete overlap with bedfiles specified as `b`, but not uniq. becomes problematic when `-b` is given multiple bed files
bedtools intersect -a gr1_chrStartEndOnly.bed -b gr2_chrStartEndOnly.bed -f 1 -u #complete overlap with `b` and unique this is what i want
bedtools intersect -a gr1_chrStartEndOnly.bed -b gr2_chrStartEndOnly.bed gr3_chrStartEndOnly.bed -f 1 -u #complete overlap with `b` and unique this is what i want
bedtools intersect -a gr1_chrStartEndOnly.bed -b gr2_chrStartEndOnly.bed -f 1 -c #complete overlaps and counts


awk -vFS="\t" -vOFS="\t" '{ $3 += 1; print $0; }' gr1.bed > gr1_0based.bed
awk -vFS="\t" -vOFS="\t" '{ $3 += 1; print $0; }' gr2.bed > gr2_0based.bed
awk -vFS="\t" -vOFS="\t" '{ $3 += 1; print $0; }' gr3.bed > gr3_0based.bed
bedops --merge gr1_0based.bed gr2_0based.bed gr3_0based.bed | bedops --chop 1 - | bedmap --echo --count --echo-map-id --delim '\t' - <(bedops --everything gr1_0based.bed gr2_0based.bed gr3_0based.bed) > answer.bed
cat answer.bed
bedtools intersect -a gr1_0based.bed -b gr2_0based.bed gr3_0based.bed
awk '{ if ($4 == 3) { print } }' answer.bed | head
awk '$4 ~ 3' answer.bed




######################## didn't work 
##bedops, try using bedops
# bedops -u /juno/work/greenbaum/users/ahunos/methyl_SPECTRUM/data/processed/nf_core_methylseq/nf_core_methylseq.interval_list \
# /juno/work/greenbaum/users/ahunos/methyl_SPECTRUM/data/processed/capture1/capture1.interval_list \
# /juno/work/greenbaum/users/ahunos/methyl_SPECTRUM/data/processed/capture2/capture2.interval_list \
#     | bedmap --echo --echo-map-id-uniq --fraction-both 0.5 - \
#     | awk -F"|" '(split($2, a, ";") > 1)' | cut -f1 -d'|' > almost_answer.bed


# bedops -u /juno/work/greenbaum/users/ahunos/methyl_SPECTRUM/data/processed/nf_core_methylseq/nf_core_methylseq.interval_list \
# /juno/work/greenbaum/users/ahunos/methyl_SPECTRUM/data/processed/capture1/capture1.interval_list \
# /juno/work/greenbaum/users/ahunos/methyl_SPECTRUM/data/processed/capture2/capture2.interval_list | bedmap --echo --echo-map-id-uniq --fraction-both 0.5 - | awk -F"|" '(split($2, a, ";") > 1)' \
#     | cut -f1 -d'|' | awk -F"|" '(split($2, a, ";") > 1)' | cut -f1 -d'|' > answer.bed




# bedops -u /juno/work/greenbaum/users/ahunos/methyl_SPECTRUM/data/processed/nf_core_methylseq/nf_core_methylseq.interval_list \
# /juno/work/greenbaum/users/ahunos/methyl_SPECTRUM/data/processed/capture1/capture1.interval_list \
# /juno/work/greenbaum/users/ahunos/methyl_SPECTRUM/data/processed/capture2/capture2.interval_list  > /juno/work/greenbaum/users/ahunos/methyl_SPECTRUM/bedops_intervals/methylSeq_Agilent_Twist_50perOv.interval_list.bed

# bedops -u /juno/work/greenbaum/users/ahunos/methyl_SPECTRUM/data/processed/nf_core_methylseq/nf_core_methylseq.interval_list \
# /juno/work/greenbaum/users/ahunos/methyl_SPECTRUM/data/processed/capture1/capture1.interval_list \
# /juno/work/greenbaum/users/ahunos/methyl_SPECTRUM/data/processed/capture2/capture2.interval_list  \
# | bedmap --echo --echo-map-id-uniq --fraction-both 1 - > /juno/work/greenbaum/users/ahunos/methyl_SPECTRUM/bedops_intervals/methylSeq_Agilent_Twist_50perOv2.interval_list.bed


# ##
# bedops --everything /juno/work/greenbaum/users/ahunos/methyl_SPECTRUM/data/processed/nf_core_methylseq/nf_core_methylseq.interval_list \
# /juno/work/greenbaum/users/ahunos/methyl_SPECTRUM/data/processed/capture1/capture1.interval_list \
# /juno/work/greenbaum/users/ahunos/methyl_SPECTRUM/data/processed/capture2/capture2.interval_list | bedmap --echo --echo-map-id-uniq --delim '\t' - > \
# /juno/work/greenbaum/users/ahunos/methyl_SPECTRUM/bedops_intervals/methylSeq_Agilent_Twist_50perOv2.interval_list.bed

# awk -vthreshold=3 '(split($5, ids, ";") >= threshold)' methylSeq_Agilent_Twist_50perOv2.interval_list.bed > thresholded.bed

# #$ bedmap --echo --exact union.bed thresholded.bed > union-thresholded.bed


# bedops -u /juno/work/greenbaum/users/ahunos/methyl_SPECTRUM/data/processed/nf_core_methylseq/nf_core_methylseq.interval_list \
# /juno/work/greenbaum/users/ahunos/methyl_SPECTRUM/data/processed/capture1/capture1.interval_list \
# /juno/work/greenbaum/users/ahunos/methyl_SPECTRUM/data/processed/capture2/capture2.interval_list  \
# | bedmap --echo --echo-map-id-uniq --fraction-both 1 - > almost.bed


# bedops -u master_gr.bed \
# gr_2.bed gr_3.bed both.bed | bedmap --echo --echo-map-id-uniq --fraction-both 0.5 -> almost.bed
# | awk -F"|" '(split($2, a, ";") > 1)' \ > almost.bed
# | cut -f1 -d'|' \
# > almost.bed

# bedtools intersect -a master_gr.bed -b gr_2.bed gr_3.bed -f 1 -r -c > both.bed


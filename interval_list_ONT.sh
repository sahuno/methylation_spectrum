#date - may 31st 2023
#goal interval list for methylation


#do for whole genome bisulphite, uncomment if that's what you are doings
dir_results=/work/greenbaum/projects/ont_pipeline/projects/SPECTRUM_MTHY/results/methylation/megalodon
mod_5mc=`find ${dir_results}/* -type f -iname "*.5mC.bed" | xargs ls`
mod_5hmc=`find ${dir_results}/* -type f -iname "*.5hmC.bed" | xargs ls`


#echo "${mod_5mc} ${mod_5hmc}"
filesPaths=`echo "${mod_5mc} ${mod_5hmc}"`
echo $filesPaths
cat ${mod_5mc} | head

#make parent dir to save
methyl_spect=/juno/work/greenbaum/users/ahunos/methyl_SPECTRUM/
mkdir -p ${methyl_spect}/data/processed/ONT


#recreate modified base to save space and time
#final file will have;
#chr, start, end, strand, #total_number_reads, percentage methylated, #Number_of_methylated_reads, $10-($10*($11/100))

for i in $filesPaths
do
#extract filenames
echo "processing file $i" #echo file being procesed
fileName=`echo "${i}" | grep -o 'Spectrum-OV.*' | sed -e 's/\//_/'`
sampleName=`echo "${i}" | grep -o 'Spectrum-OV[^_]*'`

echo "sampleName is ${sampleName}"
#make dir to save
#echo "mkdir -p ${methyl_spect}/data/processed/ONT/${sampleName}"
mkdir -p ${methyl_spect}/data/processed/ONT/${sampleName}

#reduce modified base file
awk -vFS="\t" -vOFS="\t" '{printf ("%s\t%s\t%s\t%s\t%.0f\t%.1f\t%.0f\t%.0f\n", $1,$2,$3,$6,$10,$11, $10*($11/100), $10-($10*($11/100)))}' $i > ${methyl_spect}/data/processed/ONT/${sampleName}/${fileName}
#chr, start, end, strand, #total_number_reads, percentage methylated, #Number_of_methylated_reads, $10-($10*($11/100))
done

#grep -Eom1 "[,| ]" file | head -1

#sanity checks to compare length of resultant file and orginal file 
# (base) bash-4.2$ wc -l /work/greenbaum/projects/ont_pipeline/projects/SPECTRUM_MTHY/results/methylation/megalodon/Spectrum-OV-009_N/modified_bases.5mC.bed 
# 60274036 /work/greenbaum/projects/ont_pipeline/projects/SPECTRUM_MTHY/results/methylation/megalodon/Spectrum-OV-009_N/modified_bases.5mC.bed
# (base) bash-4.2$ wc -l Spectrum-OV-009/Spectrum-OV-009_N_modified_bases.5mC.bed
# 60274036 Spectrum-OV-009/Spectrum-OV-009_N_modified_bases.5mC.bed

#sanity check , see there seperators for a file
#cat -t /juno/work/greenbaum/users/ahunos/methyl_SPECTRUM/data/processed/ONT/Spectrum-OV-044/Spectrum-OV-044_N_modified_bases.5mC.bed | head

#sanity check, see if each site is both measured for 5hmc and 5mc? can 5hmc and 5mc asssumed to occur on the same cpG site? 
bedtools intersect -a /juno/work/greenbaum/users/ahunos/methyl_SPECTRUM/data/processed/ONT/Spectrum-OV-044/Spectrum-OV-044_N_modified_bases.5mC.bed \
-b /juno/work/greenbaum/users/ahunos/methyl_SPECTRUM/data/processed/ONT/Spectrum-OV-044/Spectrum-OV-044_N_modified_bases.5hmC.bed -v | head

bedtools intersect -a /juno/work/greenbaum/users/ahunos/methyl_SPECTRUM/data/processed/ONT/Spectrum-OV-044/Spectrum-OV-044_N_modified_bases.5hmC.bed \
-b /juno/work/greenbaum/users/ahunos/methyl_SPECTRUM/data/processed/ONT/Spectrum-OV-044/Spectrum-OV-044_N_modified_bases.5mC.bed -v | head

#conclusion
#each cpg site - or surveyed site is measured for 5hmc or 5mc methylation. there is probabilities for each site being 5mc or 5hmc. there is number of reads supporting

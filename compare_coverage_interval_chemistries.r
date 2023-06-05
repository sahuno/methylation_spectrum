#may 15th 2023; samuel ahuno
library(magrittr)
library(data.table)
library(plyranges)
library(ggvenn)
library(purrr)
library(VennDiagram)
library(GenomicRanges)
library(data.table)
library(cowplot)
library(seqsetvis)
theme_set(cowplot::theme_cowplot())
library(RColorBrewer)
myCol <- brewer.pal(3, "Pastel2")

#set libraries to work with
path_methyl_spect <- "/juno/work/greenbaum/users/ahunos/methyl_SPECTRUM"
paths_intervals <- "/juno/work/greenbaum/users/ahunos/methyl_SPECTRUM/data/processed/" 
files_intervals <- list.files(paths_intervals, full.names = T, recursive = TRUE, pattern=".interval_list") #%>% as_tibble()



##define function to read files
readplus <- function(file_path_in, snippet_size = Inf, nThreads_spe=8){
data_df <- data.table::fread(file = file_path_in, sep = "\t", nrows = snippet_size, nThread=nThreads_spe) %>% mutate(Chemistry_file = basename(file_path_in))
}

#read interval lists 
methyl_intervals_from_coverage <- files_intervals %>% map_df(~readplus(.))
names(methyl_intervals_from_coverage) <- c("Chromosome", "Start", "End", "Chemistry_file")
head(methyl_intervals_from_coverage)
#unique(methyl_intervals_from_coverage$Chemistry_file)
#table(methyl_intervals_from_coverage$Chemistry_file)


##refactor and rename 
methyl_intervals_from_coverage %<>% mutate(Chemistry_file = case_when(Chemistry_file == "capture1.interval_list" ~ "Agilent",
                                                                    Chemistry_file == "capture2.interval_list" ~ "Twist",
                                                                    Chemistry_file == "nf_core_methylseq.interval_list" ~ "MethylSeqWGS",
                                                                    Chemistry_file == "methylSeq_Agilent_Twist.interval_list" ~ "Agilent_Twist_MethylSeqWGS")) %>% 
                                                                        mutate(Chemistry_file = fct_relevel(Chemistry_file,c("Agilent","Twist", "MethylSeqWGS","Agilent_Twist_MethylSeqWGS")))


#plot common cpg sites by  each chemistry
common_cpg_chemistries <- methyl_intervals_from_coverage %>% group_by(Chemistry_file) %>% summarise(counts = n())

plot_common_sites_CpG <- common_cpg_chemistries %>% 
                                ggplot(aes(Chemistry_file, counts, fill=Chemistry_file)) + 
                                    geom_col() + 
                                        labs(title = "common cpg sites in tumor/Normal for each methyl bisulphite chemistry") + 
                                            scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                                                labels = scales::trans_format("log10", scales::math_format(expr = 10^.x))) + 
                                                scale_fill_viridis_d() +
                                                    theme(legend.position = "none") +
geom_text(data=common_cpg_chemistries, aes(label = counts), vjust = -0.5)

ggsave(plot_common_sites_CpG, file = paste0(path_methyl_spect, "/figures/common_cpg_sites_tumorNormal_per_bisulphite_chemistries.pdf"))



methyl_intervals_from_coverage_toGR <- methyl_intervals_from_coverage
names(methyl_intervals_from_coverage_toGR) <- c("seqnames", "start", "end", "chemistry_file")
#methyl_intervals_from_coverage %<>% mutate(End = Start)
head(methyl_intervals_from_coverage_toGR) 

#convert file into list of data frames
methyl_intervals_from_coverage_list <- methyl_intervals_from_coverage_toGR %>% group_split(chemistry_file) %>% setNames(unique(methyl_intervals_from_coverage_toGR$chemistry_file))

methyl_intervals_from_coverage_list %>% map(~dim(.x))

#create function to add identifier to tables
setkey <- function(df){
    df <- df %>% mutate(key =paste0(seqnames,"_",start))
    return(df)
}
methyl_intervals_from_coverage_list_kys <- methyl_intervals_from_coverage_list %>% map(~setkey(.x)) #add keys 
methyl_intervals_from_coverage_list_kys_head <- methyl_intervals_from_coverage_list_kys %>% map(~head(.x,1000)) #add keys



res_join <- purrr::reduce(methyl_intervals_from_coverage_list_kys[c("Agilent","Twist","MethylSeqWGS")], dplyr::left_join, by = 'key')

print(res_join, n=100)
names(res_join)
# res_join %>% rename(chemistry_file.Agilent = chemistry_file.x,
#                     chemistry_file.Twist = chemistry_file.y, 
#                     chemistry_file.MethylSeqWGS = chemistry_file)

#install.packages("VennDiagram")

res_join_2plot <- res_join %>% select(c(key,chemistry_file.x,chemistry_file.y,chemistry_file))




methyl_intervals_from_coverage_gr_list <- methyl_intervals_from_coverage_list %>% map(~as_granges(.x))

methyl_intervals_from_coverage_gr_list 
methyl_intervals_from_coverage_gr_list[c("Agilent", "Twist")] %>% map(~find_overlaps(.x,.y))
#find_overlaps(x, y, maxgap, minoverlap, suffix = c(".x", ".y"))
methyl_overlaps_Agi_twist <- find_overlaps(methyl_intervals_from_coverage_gr_list[[c("Agilent")]],methyl_intervals_from_coverage_gr_list[[c("Twist")]])


methyl_overlaps_Agi_twist %>% group_by(seqnames) %>% summarize(n())







#create new granges list
methyl_int_list <- GRangesList(methyl_intervals_from_coverage_list[[c("Agilent")]],methyl_intervals_from_coverage_list[[c("Twist")]],methyl_intervals_from_coverage_list[[c("MethylSeqWGS")]])
head(methyl_int_list)
names(methyl_int_list) <- names(methyl_intervals_from_coverage_list)[!grepl("_",names(methyl_intervals_from_coverage_list))]

methyl_int_list <- GRangesList(methyl_intervals_from_coverage_list[[c("MethylSeqWGS")]], methyl_intervals_from_coverage_list[[c("Agilent")]],methyl_intervals_from_coverage_list[[c("Twist")]])
names(methyl_int_list) <- c("MethylSeqWGS","Agilent","Twist")

head(methyl_int_list) #sanity checks
#make overlaps of all subset
olaps = ssvOverlapIntervalSets(methyl_int_list, use_first = TRUE)
head(olaps)
#olaps %>% mutate(width = width(olaps))

#bar plots of overalaps
pdf(paste0(path_methyl_spect, "/figures/barplots_of_CpG_overlaps2.pdf"))
ssvFeatureBars(olaps)
dev.off()

#Venn Diagrams plots of overalaps
pdf(paste0(path_methyl_spect, "/figures/VennDiagram_of_CpG_overlaps2.pdf"))
ssvFeatureVenn(olaps)
dev.off()

#########################################################################################
#################### main script to plot venn diagrams of overlaps######################
########################################################################################
#merged interbal list
#merged_all_bed <- "/juno/work/greenbaum/users/ahunos/methyl_SPECTRUM/data/processed/methylSeq_Agilent_Twist.interval_list"
#"/juno/work/greenbaum/users/ahunos/methyl_SPECTRUM/data/processed/methylSeq_Agilent_Twist_commonSites_1bpWindow.interval_list"
#merged_all_bed <- "/juno/work/greenbaum/users/ahunos/methyl_SPECTRUM/data/processed/methylSeq_Agilent_Twist.interval_list"
merged_all_bed <- "/juno/work/greenbaum/users/ahunos/methyl_SPECTRUM/data/processed/methylSeq_Agilent_Twist.interval_list0Test_1bpWindow_annotations.bed"

df_bed <- data.table::fread(merged_all_bed,nThread=16)
names(df_bed) <- c("seqnames", "start", "end", "numb_files_overlaps","identity_files_overlapped", "methylseq", "Agilent", "Twist")

#intervals_gr <- df_bed %>% as_granges() #convert to granges
#intervals_gr_tile <- intervals_gr  %>% tile_ranges(width=1) #tile the genome

#qc, set key and use for overlaps
#intervals_gr_tile %<>% mutate(seqKey=paste0(seqnames,"_",start)) 
# intervals_gr_tile_SeqKey <- intervals_gr_tile %>% group_by(seqKey)
# intervals_gr_tile_SeqKey_filt <- filter(intervals_gr_tile_SeqKey, n() > 3)

#intervals_gr_tile_SeqKey  <- intervals_gr_tile %>% as.data.table()
#qc_dt_Counts2 <- intervals_gr_tile_SeqKey[,.N, by=seqKey][order(-N)]

#intervals_gr_tile %>% filter(seqKey == "1_10610")
##df_bed[(seqnames=="1" && start>=10468),][,head(.SD,1)]
#df_bed[,head(.SD,20)]

#intervals_gr_tile <- join_overlap_left(intervals_gr_tile, intervals_gr, suffix = c(".x", ".y")) #merge metadata
# intervals_gr_tile %>% select(-c(partition, numb_files_overlaps,identity_files_overlapped)) %>% head()

# pdf(file = paste0(path_methyl_spect, "/figures/common_cpg_sites_tumorNormal_all_bisulphite_chemistries.pdf"))
# ssvFeatureVenn(intervals_gr_tile %>% select(-c(partition, numb_files_overlaps,identity_files_overlapped)), group_names = c("methylseq", "Agilent", "Twist")) 
# dev.off()



##function to ven diagram
# get_vennData <- function(df, var){
#     var <- ensym(var)
#   df <- df %>% mutate(key =paste0(seqnames,"_",start)) %>% 
#                 as.data.frame() %>% dplyr::filter({{var}} == 1)
#     return(df)
# }

#mamke it data tabel
#intervals_dt_tile <- intervals_gr_tile %>% as.data.table()
# head(intervals_dt_tile)
# dim(intervals_dt_tile)
# head(intervals_dt_tile)[, .lapply(c("numb_files_overlaps","methylseq"), as.logical)]
#Compute on columns: 
# head(intervals_dt_tile)[, .(as.logical(methylseq),as.logical(Agilent),as.logical(Twist))]
# head(intervals_dt_tile)[, .(meth=as.logical(methylseq),ag=as.logical(Agilent),tw=as.logical(Twist))]


#set columns to pull for membership tables
intervals_dt_tile <- df_bed
cols_use <- c("methylseq","Agilent","Twist")
intervals_dt_tile_memberships <- intervals_dt_tile[, lapply(.SD, as.logical), .SDcols=cols_use]
#plot membership tabeles now
v_plot <- ggplot(intervals_dt_tile_memberships,                         # Apply geom_venn function
       aes(A = methylseq, B = Agilent, C = Twist)) + geom_venn() + 
       labs(title = "Common sites for methylation Bisulphite Seq") + scale_fill_brewer(palette = "Dark2")
ggsave(v_plot, filename = paste0(path_methyl_spect, "/figures/common_cpg_sites_tumorNormal_all_bisulphite_chemistries_percentages.pdf"))

saveRDS(intervals_dt_tile, file = paste0(path_methyl_spect, "/data/commonSitesBisulphites_methyl.rds"))


# List of items
# h_gr <- head(intervals_gr_tile %>% select(-c(partition, numb_files_overlaps,identity_files_overlapped)), 1000)
# library(ggvenn)
# as.data.frame(h_gr)

# #data_venn, create membership table
# v_plot <- ggplot(as.data.frame(h_gr) %>% select(seqnames,methylseq,Agilent,Twist),                         # Apply geom_venn function
#        aes(A = as.logical(methylseq), B = as.logical(Agilent), C = as.logical(Twist))) + geom_venn()
# ggsave(v_plot, filename = paste0(path_methyl_spect, "/figures/v_diadg_all_sites_test.pdf"))

# intervals_dt_tile[, Mean:=counts(identity_files_overlapped), by=list(identity_files_overlapped)]
# intervals_dt_tile[identity_files_overlapped == "methylseq,Twist", .N]

# ans <- intervals_dt_tile[, .(.N), by = .(identity_files_overlapped)]
# intervals_dt_tile[, counts:= .(.N), by = .(identity_files_overlapped)] #[, prop := N/N]





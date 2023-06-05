#samuel ahuno
#tumor normal plots
library(magrittr)
library(data.table)
library(plyranges)
library(patchwork)
library(RColorBrewer)
library(purrr)



#set paths for tumor and normal
tumPath = "/work/greenbaum/projects/methylation_projects/SPECTRUM/nf_core_methylseq_tum/OUTDIR/bismark/methylation_calls/methylation_coverage/SA123T_1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz"
normPath =  "/work/greenbaum/projects/methylation_projects/SPECTRUM/nf_core_methylseq/OUTDIR/bismark/methylation_calls/methylation_coverage/SA123_1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz"

# track type=bedGraph
# 1       10468   10469   75
# 1       10469   10470   100
# 1       10470   10471   50
# 1       10471   10472   50


#save tumor normal paths as list
tumNorm_ls <- list(tumor = tumPath,
                    normal = normPath)
#set libraries to work with
path_methyl_spect <- "/juno/work/greenbaum/users/ahunos/methyl_SPECTRUM"

#read and rename files
load_data <- function(file_path){
    df <- data.table::fread(file = file_path, sep = "\t", nrows = Inf, nThread=12)
    names(df) <- c("chromosome", "start", "end", "methylation_percentage", "count_methylated", "count_unmethylated")
    return(df)
}

#read files
tumNorm_data_ls <- tumNorm_ls %>% map(~load_data(.x))

# tum_df <- load_data(tumPath)
# Norm_df <- load_data(normPath)

#read interval lists 


#####
#function to plot
geom_text_size = 3
plotHisto <- function(df,x_var){
    x_var <- ensym(x_var) #quote variabe to plot

    stats_df_hist <- df %>% summarise(mean=round(mean({{x_var}}, na.rm = TRUE), 2), 
          median=round(median({{x_var}}, na.rm = TRUE),2),
          sd=round(sd({{x_var}}, na.rm = TRUE), 2),
          max=round(max({{x_var}}, na.rm = TRUE), 2),
          min=round(min({{x_var}}, na.rm = TRUE), 2),
          counts=n())


    plot_obj <- df %>% ggplot(aes({{x_var}})) + geom_histogram() +
    geom_text(x = 10, y = 5, aes(label = paste0("\nCounts: ", counts,
                                        "\nRange: ", min, "-", max,
                                        "\nmean: ", mean, 
                                      "\nmedian: ", median,
                                      "\nSD: ", sd)), color="blue",
                              data = stats_df_hist, size=geom_text_size) +
     scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                   labels = scales::trans_format("log10", scales::math_format(expr = 10^.x)))
    #plot histogram
#ggsave(plot_obj, filename = paste0("figures/",x_var,"_histogram_",fileTag,".pdf")) #save
return(plot_obj)
}

#save
#tumNorm_data_ls %>% map(~plotHisto(df=.x, x_var=methylation_percentage, fileTag=.y)) 
methy_percentage_ls <- tumNorm_data_ls %>% map(~plotHisto(df=.x, x_var=methylation_percentage)) 
methylation_percentage <- wrap_plots(methy_percentage_ls) + plot_annotation('methylation_percentage; Tumor vrs Normal')
ggsave(methylation_percentage, filename = paste0(path_methyl_spect,"/figures/","methylation_percentage_HistogramTumNorm_allSequencedSites_WGBS.pdf")) 

count_methylated_ls <- tumNorm_data_ls %>% map(~plotHisto(df=.x, x_var=count_methylated)) 
count_methylated_wrapPlot <- wrap_plots(count_methylated_ls) + plot_annotation('count_methylated; Tumor vrs Normal')
ggsave(count_methylated_wrapPlot, filename = paste0(path_methyl_spect,"/figures/","count_methylated_HistogramTumNorm_allSequencedSites_WGBS.pdf")) 

count_unmethylated_ls <- tumNorm_data_ls %>% map(~plotHisto(df=.x, x_var=count_unmethylated)) 
count_unmethylated_wrapPlot <- wrap_plots(count_unmethylated_ls) + plot_annotation('count_unmethylated; Tumor vrs Normal')
ggsave(count_unmethylated_wrapPlot, filename = paste0(path_methyl_spect,"/figures/","count_unmethylated_HistogramTumNorm_allSequencedSites_WGBS.pdf")) 










#####################################
### Step 2; plot common sites only
######################################
#previously needed to run this script; but no need anymore since all common sites are comming from bedtools
#/juno/work/greenbaum/users/ahunos/methyl_SPECTRUM/scripts/coverage_interval_list_compare.r

#intervals_dt_tile <- readRDS(file = "/juno/work/greenbaum/users/ahunos/methyl_SPECTRUM/data/commonSitesBisulphites_methyl.rds")
intervals_dt_tile <- data.table::fread(file = "/juno/work/greenbaum/users/ahunos/methyl_SPECTRUM/data/processed/methylSeq_Agilent_Twist.interval_list0based_1bpWindow_annotations_BSSeq.bed",nThread=16)
names(intervals_dt_tile) <- c("seqnames", "start", "end", "numb_files_overlaps","identity_files_overlapped", "methylseq", "Agilent", "Twist")

dim(intervals_dt_tile)

########### qc-begin, check if the entries are unique ###########
# qc_dt <- intervals_dt_tile[,.(seqKey=paste0(seqnames,"_",start),seqnames, start,   end)] #create key
# qc_dt_Counts <- qc_dt[,.N, by=seqKey][order(-N)] #sort by decensing 
# dim(qc_dt_Counts)
# qc_dt_Counts[N>1, .N] #are there any granges with multiple entries 
# head(intervals_dt_tile,50)
############# QC-end ###########


#convert commonsites interval list to granges 
class(intervals_dt_tile)
if(is.data.table(intervals_dt_tile)){
    commonSites_BSulphite_gr <- makeGRangesFromDataFrame(intervals_dt_tile, keep.extra.columns=TRUE,starts.in.df.are.0based=TRUE)
    #commonSites_BSulphite_gr <- makeGRangesFromDataFrame(intervals_dt_tile[numb_files_overlaps == 3,],  keep.extra.columns=TRUE)
}


#subset only  cpg common sites of tumor normal samples
tumNorm_gr_ls <- tumNorm_data_ls %>% map(~makeGRangesFromDataFrame(.x, keep.extra.columns=TRUE))
#function to subset using plyranges filter_by_overlaps
filter_ov <- function(gr, subject_gr=commonSites_BSulphite_gr){
    ov_gr <- filter_by_overlaps(x=gr, y=subject_gr)
    #ov_gr <- filter_by_overlaps(x=gr, y=subject_gr)
    return(ov_gr)
}

#commonSites_BSulphite_gr[!complete.cases(commonSites_BSulphite_gr),]



#do actual filtering for common sites across all 3 chemistries
tumNorm_ovlaps_comSites_ls <- tumNorm_gr_ls %>% map(~filter_ov(gr=.x, subject_gr=commonSites_BSulphite_gr))

#Function to combine and save plots 
# save_cowPlots <- function(data_TN_ls, fileTag){
# methy_percentage_ls <- data_TN_ls %>% map(~plotHisto(df=.x, x_var=methylation_percentage)) 
# methylation_percentage <- wrap_plots(methy_percentage_ls) + plot_annotation('methylation_percentage; Tumor vrs Normal')
# ggsave(methylation_percentage, filename = paste0("figures/","methylation_percentage_HistogramTumNorm_",fileTag,".pdf")) 

# count_methylated_ls <- data_TN_ls %>% map(~plotHisto(df=.x, x_var=count_methylated)) 
# count_methylated_wrapPlot <- wrap_plots(count_methylated_ls) + plot_annotation('count_methylated; Tumor vrs Normal')
# ggsave(count_methylated_wrapPlot, filename = paste0("figures/","count_methylated_HistogramTumNorm_",fileTag,".pdf")) 

# count_unmethylated_ls <- data_TN_ls %>% map(~plotHisto(df=.x, x_var=count_unmethylated)) 
# count_unmethylated_wrapPlot <- wrap_plots(count_unmethylated_ls) + plot_annotation('count_unmethylated; Tumor vrs Normal')
# ggsave(count_unmethylated_wrapPlot, filename = paste0("figures/","count_unmethylated_HistogramTumNorm_",fileTag,".pdf")) 
# }

#convert granges back to data.table for ploting
tumNorm_ovlaps_comSites_DT_ls <- tumNorm_ovlaps_comSites_ls %>% map(~as.data.table(.x))



###begin function for qc
# qc_ov <- function(dt){
# #dt <- as.data.table(gr) 
# qc_dt <- dt[,.(seqKey=paste0(seqnames,"_",start),seqnames, start,   end)] #create key
# qc_dt_Counts <- qc_dt[,.N, by=seqKey][order(-N)] #sort by decensing 
# return(qc_dt_Counts)
# }

# tumNorm_ovlaps_comSites_DT_ls_qc <- tumNorm_ovlaps_comSites_DT_ls %>% map(~qc_ov(as.data.table(.x)))

# qc_Tn <- function(dt){
# #dt <- as.data.table(gr) 
# qc_dt <- dt[,.(seqKey=paste0(chromosome,"_",start),seqnames, start,   end)] #create key
# qc_dt_Counts <- qc_dt[,.N, by=seqKey][order(-N)] #sort by decensing 
# return(qc_dt_Counts)
# }
# tumNorm_data_ls_qc <- tumNorm_data_ls %>% map(~qc_Tn(.x))
# ###end of qc

#tumNorm_data_ls_qc_join <- list(tumNorm_data_ls_qc[[1]],tumNorm_data_ls_qc[[2]]) %>% reduce(left_join, by = "seqKey")
# tumNorm_ovlaps_comSites_DT_ls_qc_merged <- Reduce(function(...) merge(..., by="seqKey", all=T), tumNorm_ovlaps_comSites_DT_ls_qc)
# tumNorm_data_ls_qc_merged <- Reduce(function(...) merge(..., by="seqKey", all=T), tumNorm_data_ls_qc)

# tumNorm_ovlaps_comSites_DT_ls_qc_merged %>% filter(seqKey == "1_68057")
# tumNorm_ovlaps_comSites_DT_ls_qc_merged %>% filter(is.na(N.x))
# tumNorm_ovlaps_comSites_DT_ls_qc_merged %>% filter(is.na(N.y))


#intersect tumor and normal
tumNorm_intersect_rng <- join_overlap_intersect(tumNorm_ovlaps_comSites_ls[["tumor"]], tumNorm_ovlaps_comSites_ls[["normal"]],suffix = c(".tumor", ".normal"))
tumNorm_intersect_df <- as.data.table(tumNorm_intersect_rng)

# head(tumNorm_ovlaps_comSites_DT_ls[[1]])
# str(head(tumNorm_ovlaps_comSites_DT_ls[[1]]))
# str(head(tumNorm_data_ls[[1]]))
#test run functions
# tumNorm_ovlaps_comSites_DT_ls %>% iwalk(~save_cowPlots(data_TN_ls=.x, fileTag="WGS-BS_comSites"))
# tumNorm_data_ls %>% walk(~save_cowPlots(data_TN_ls=.x, fileTag="WGS-BS_")) #actually didn't work

fileTag="WGS-BS_comSites"
methy_percentage_ls <- tumNorm_ovlaps_comSites_DT_ls %>% map(~plotHisto(df=.x, x_var=methylation_percentage)) 
methylation_percentage <- wrap_plots(methy_percentage_ls) + plot_annotation('methylation_percentage; Tumor vrs Normal')
ggsave(methylation_percentage, filename = paste0(path_methyl_spect,"/figures/","methylation_percentage_HistogramTumNorm_",fileTag,".pdf")) 

count_methylated_ls <- tumNorm_ovlaps_comSites_DT_ls %>% map(~plotHisto(df=.x, x_var=count_methylated)) 
count_methylated_wrapPlot <- wrap_plots(count_methylated_ls) + plot_annotation('count_methylated; Tumor vrs Normal')
ggsave(count_methylated_wrapPlot, filename = paste0(path_methyl_spect,"/figures/","count_methylated_HistogramTumNorm_",fileTag,".pdf")) 

count_unmethylated_ls <- tumNorm_ovlaps_comSites_DT_ls %>% map(~plotHisto(df=.x, x_var=count_unmethylated)) 
count_unmethylated_wrapPlot <- wrap_plots(count_unmethylated_ls) + plot_annotation('count_unmethylated; Tumor vrs Normal')
ggsave(count_unmethylated_wrapPlot, filename = paste0(path_methyl_spect,"/figures/","count_unmethylated_HistogramTumNorm_",fileTag,".pdf")) 


# pdf("test_wgs_methy_percent.pdf")
# tumNorm_data_ls %>% map(~plotHisto(df=.x, x_var=methylation_percentage))
# dev.off()
# save_cowPlot(data_TN_ls=, fileTag)

#########
##calculate Tumor Normal common sites
#list(x, y, z) %>% reduce(left_join, by = "i")



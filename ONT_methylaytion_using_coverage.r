library(data.table)
library(magrittr)
library(optparse)
library(tidyverse)

########################################################
######accept arguments from environment
########################################################
option_list_val <- list(make_option(c("-s","--sample"), type = "character", default=NULL, help = "unique identifer sample"),
                        make_option(c("-i","--input_file"), type="character", default=NULL, help = "per read modified base file"),
                        make_option(c("-h","--hmc_output_file"), type="character", default=NULL, help = "5hmC output file name .txt"),
                        make_option(c("-m","--mc_output_file"), type="character", default=NULL, help = "5mC output file name .txt"))
parseObject <- OptionParser(option_list = option_list_val)
opt <- parse_args(parseObject)

#print(opt)
########################################################
#sanity check
if (is.null(opt$input_file)){
  print_help(parseObject)
  stop("Please supply sample to work on", call.=FALSE)
}
opt$input_file <- "/work/greenbaum/projects/ont_pipeline/projects/SPECTRUM_MTHY/results/methylation/megalodon/Spectrum-OV-044_N/per_read_modified_base_calls.txt"
df_perReadModifiedBases <- data.table::fread(opt$input_file, nThread=16)
head(df_perReadModifiedBases)


head(df_perReadModifiedBases, 2) %>% mutate(across(c(mod_log_prob, can_log_prob), ~ exp(.x))) %>% 
            mutate(can_plus_mod_prob = (mod_log_prob + can_log_prob))
head(df_perReadModifiedBases, 2) %>% mutate(across(c(mod_log_prob, can_log_prob), ~ exp(.x))) %>% 
            mutate(can_plus_mod_prob = (mod_log_prob + can_log_prob))
head(df_perReadModifiedBases)[mod_base == "m", by=list(chrm,pos, strand, read_id)]
head(df_perReadModifiedBases)[mod_base == "m", .N,by=list(chrm,pos, strand, read_id)]
head(df_perReadModifiedBases, 1000)[mod_base == "m", .N,by=list(chrm,pos, strand, read_id)]
head(df_perReadModifiedBases, 1000)[mod_base == "m", .N,by=list(chrm,pos, strand, read_id)][order(N)]
head(df_perReadModifiedBases, 1000)[mod_base == "m", .N,by=list(chrm,pos, strand, read_id)][order(-N)]
head(df_perReadModifiedBases, 10000)[mod_base == "m", .N,by=list(chrm,pos, strand, read_id)][order(-N)]
head(df_perReadModifiedBases, 100000)[mod_base == "m", .N,by=list(chrm,pos, strand, read_id)][order(-N)]
head(df_perReadModifiedBases, 1000000)[mod_base == "m", .N,by=list(chrm,pos, strand, read_id)][order(-N)]
head(df_perReadModifiedBases, 1000000) %>% filter(mod_base == "m") %>% group_by(chrm, pos, strand,read_id) %>% summarise(n())
head(df_perReadModifiedBases, 1000000) %>% filter(mod_base == "m") %>% group_by(chrm, pos, strand,read_id) %>% summarise(numb = n()) %>% arrange(desc(numb))
head(df_perReadModifiedBases, 2)
head(df_perReadModifiedBases, 1000000)[mod_base == "m", list(mean(exp(mod_log_prob))), by=list(chrm, pos, strand,read_id)][order(-N)]
head(df_perReadModifiedBases, 1000000)[mod_base == "m", list(mean(exp(mod_log_prob))), by=list(chrm, pos, strand,read_id)]#[order(-N)]
head(df_perReadModifiedBases, 1000000)[mod_base == "m", list(.N,mean(exp(mod_log_prob))), by=list(chrm, pos, strand,read_id)]#[order(-N)]
head(df_perReadModifiedBases, 1000000)[mod_base == "m", list(exp_mod_log_prob = exp(mod_log_prob)), by=list(chrm, pos, strand,read_id)]#[order(-N)]
head(df_perReadModifiedBases, 1000000)[mod_base == "m", 
                                        list(number_reads=.N, mean_exp_mod_log_prob=mean(exp(mod_log_prob)), 
                                        median_exp_mod_log_prob=median(exp(mod_log_prob)), 
                                        min_exp_mod_log_prob=min(exp(mod_log_prob)), max_exp_mod_log_prob=max(exp(mod_log_prob))), 
                                        by=list(chrm, pos, strand,read_id)]


head(df_perReadModifiedBases, 1000000)[!chrm %like% "_"]
head(df_perReadModifiedBases, 1000000)[!chrm %like% "_"][mod_base == "m", 
                                        list(number_reads=.N, mean_exp_mod_log_prob=mean(exp(mod_log_prob)), 
                                        median_exp_mod_log_prob=median(exp(mod_log_prob)), 
                                        min_exp_mod_log_prob=min(exp(mod_log_prob)), max_exp_mod_log_prob=max(exp(mod_log_prob))), 
                                        by=list(chrm, pos, strand,read_id)]#[order(-N)]
df_perReadModifiedBases[!chrm %like% "_"][mod_base == "m", 
                                        list(number_reads=.N), 
                                        by=list(chrm, pos, strand,read_id)][order(-N)]
df_perReadModifiedBases[!chrm %like% "_"][mod_base == "m", 
                                        list(number_reads=.N), 
                                        by=list(chrm, pos, strand,read_id)]
df_perReadModifiedBases[!chrm %like% "_"][mod_base == "m", 
                                        list(number_reads=.N), 
                                        by=list(chrm, pos, strand,read_id)][number_reads>1]
df_perReadModifiedBases[!chrm %like% "_"][mod_base == "m", 
                                        list(number_reads=.N), 
                                        by=list(chrm, pos, strand)][number_reads>1]
head(df_perReadModifiedBases, 1000000)[!chrm %like% "_"][mod_base == "m", 
                                        list(exp_mod_log_prob=exp(mod_log_prob),
                                            number_reads=.N), 
                                        by=list(chrm, pos, strand,read_id)][list(mean_exp_mod_log_prob=mean(exp_mod_log_prob), 
                                        median_exp_mod_log_prob=median(exp_mod_log_prob), 
                                        min_exp_mod_log_prob=min(exp_mod_log_prob), max_exp_mod_log_prob=max(exp_mod_log_prob))]
head(df_perReadModifiedBases, 1000000)[!chrm %like% "_"][mod_base == "m", 
                                        list(exp_mod_log_prob=exp(mod_log_prob),
                                            number_reads=.N), 
                                        by=list(chrm, pos, strand,read_id)]
head(df_perReadModifiedBases, 1000000)[!chrm %like% "_"][mod_base == "m", 
                                        list(exp_mod_log_prob=exp(mod_log_prob),
                                            number_reads=.N), 
                                        by=list(chrm, pos, strand)]
                                        [list(mean_exp_mod_log_prob=mean(exp_mod_log_prob), 
                                        median_exp_mod_log_prob=median(exp_mod_log_prob), 
                                        min_exp_mod_log_prob=min(exp_mod_log_prob), 
                                        max_exp_mod_log_prob=max(exp_mod_log_prob))]


head(df_perReadModifiedBases, 1000000)[!chrm %like% "_"][mod_base == "m", 
                                        list(exp_mod_log_prob=exp(mod_log_prob),
                                            number_reads=.N), 
                                        by=list(chrm, pos, strand)][list(mean_exp_mod_log_prob=mean(exp_mod_log_prob), 
                                        median_exp_mod_log_prob=median(exp_mod_log_prob), 
                                        min_exp_mod_log_prob=min(exp_mod_log_prob), 
                                        max_exp_mod_log_prob=max(exp_mod_log_prob))]
head(df_perReadModifiedBases, 1000000)[!chrm %like% "_"][mod_base == "m", 
                                        list(exp_mod_log_prob=exp(mod_log_prob),
                                            number_reads=.N), 
                                        by=list(chrm, pos, strand)]
head(df_perReadModifiedBases, 1000000)[!chrm %like% "_"][mod_base == "m", 
                                        list(exp_mod_log_prob=exp(mod_log_prob),
                                            number_reads=.N), 
                                        by=list(chrm, pos, strand)][list(mean_exp_mod_log_prob=mean(exp_mod_log_prob), 
                                        median_exp_mod_log_prob=median(exp_mod_log_prob), 
                                        min_exp_mod_log_prob=min(exp_mod_log_prob), 
                                        max_exp_mod_log_prob=max(exp_mod_log_prob)),by=list(chrm, pos, strand)]
head(df_perReadModifiedBases, 1000000)[!chrm %like% "_"][mod_base == "m", 
                                        list(number_reads=.N, mean_exp_mod_log_prob=mean(exp(mod_log_prob)), 
                                        median_exp_mod_log_prob=median(exp(mod_log_prob)), 
                                        min_exp_mod_log_prob=min(exp(mod_log_prob)), max_exp_mod_log_prob=max(exp(mod_log_prob))), 
                                        by=list(chrm, pos, strand)]#[order(-N)]
dt_perReadModifiedBases_5hmC <- head(df_perReadModifiedBases, 1000000)[!chrm %like% "_"][mod_base == "h", 
                                        list(number_reads=.N, mean_exp_mod_log_prob=mean(exp(mod_log_prob)), 
                                        median_exp_mod_log_prob=median(exp(mod_log_prob)), 
                                        min_exp_mod_log_prob=min(exp(mod_log_prob)), max_exp_mod_log_prob=max(exp(mod_log_prob))), 
                                        by=list(chrm, pos, strand)]#[order(-N)]
dt_perReadModifiedBases_5hmC



head(df_perReadModifiedBases, 2) %>% mutate(across(c(mod_log_prob, can_log_prob), ~ exp(.x))) %>% 
            mutate(can_plus_mod_prob = (mod_log_prob + can_log_prob))
dt_perReadModifiedBases_5mC
#head(df_perReadModifiedBases, 1000000)[!chrm %like% "_"]
dt_perReadModifiedBases_5mC <- head(df_perReadModifiedBases, 1000000)[!chrm %like% "_"][mod_base == "m", 
                                        list(number_reads=.N, mean_exp_mod_log_prob=mean(exp(mod_log_prob)), 
                                        median_exp_mod_log_prob=median(exp(mod_log_prob)), 
                                        min_exp_mod_log_prob=min(exp(mod_log_prob)), max_exp_mod_log_prob=max(exp(mod_log_prob))), 
                                        by=list(chrm, pos, strand)]#[order(-N)]
dt_perReadModifiedBases_5mC
dt_perReadModifiedBases_5hmC
dt_perReadModifiedBases_5mC[dt_perReadModifiedBases_5hmC, on=.(chrm,pos, strand,number_reads), nomatch = NULL]
dt_perReadModifiedBases_5mC[dt_perReadModifiedBases_5hmC, on=.(chrm, pos, strand, number_reads), .(mean_exp_mod_log_prob, i.mean_exp_mod_log_prob), nomatch = NULL]
dt_perReadModifiedBases_5mC[dt_perReadModifiedBases_5hmC, on=.(chrm, pos, strand, number_reads), 
                .(chrm, pos, strand, number_reads, mean_exp_mod_log_prob_5mc=mean_exp_mod_log_prob,
                    median_exp_mod_log_prob_5mc=median_exp_mod_log_prob, min_exp_mod_log_prob_5mc=min_exp_mod_log_prob, max_exp_mod_log_prob_5mc=max_exp_mod_log_prob,
                    mean_exp_mod_log_prob_5hmc=i.mean_exp_mod_log_prob, median_exp_mod_log_prob_5hmc=i.median_exp_mod_log_prob, min_exp_mod_log_prob_5hmc=i.min_exp_mod_log_prob,
                    max_exp_mod_log_prob_5hmc=i.max_exp_mod_log_prob), nomatch = NULL]
dt_perReadModifiedBases_5mC[dt_perReadModifiedBases_5hmC, on=.(chrm, pos, strand, number_reads), 
                    .(chrm, pos, strand, number_reads, mean_exp_mod_log_prob_5mc=mean_exp_mod_log_prob,
                    median_exp_mod_log_prob_5mc=median_exp_mod_log_prob, min_exp_mod_log_prob_5mc=min_exp_mod_log_prob, max_exp_mod_log_prob_5mc=max_exp_mod_log_prob,
                    mean_exp_mod_log_prob_5hmc=i.mean_exp_mod_log_prob, median_exp_mod_log_prob_5hmc=i.median_exp_mod_log_prob, min_exp_mod_log_prob_5hmc=i.min_exp_mod_log_prob,
                    max_exp_mod_log_prob_5hmc=i.max_exp_mod_log_prob), nomatch = NULL][,methylationState:=ifelse(mean_exp_mod_log_prob_5mc>mean_exp_mod_log_prob_5hmc,"m","h")]
dt_perReadModifiedBases_5mC[dt_perReadModifiedBases_5hmC, on=.(chrm, pos, strand, number_reads), 
                    .(chrm, pos, strand, number_reads, mean_exp_mod_log_prob_5mc=mean_exp_mod_log_prob,
                    median_exp_mod_log_prob_5mc=median_exp_mod_log_prob, min_exp_mod_log_prob_5mc=min_exp_mod_log_prob, max_exp_mod_log_prob_5mc=max_exp_mod_log_prob,
                    mean_exp_mod_log_prob_5hmc=i.mean_exp_mod_log_prob, median_exp_mod_log_prob_5hmc=i.median_exp_mod_log_prob, min_exp_mod_log_prob_5hmc=i.min_exp_mod_log_prob,
                    max_exp_mod_log_prob_5hmc=i.max_exp_mod_log_prob), nomatch = NULL][methylationState:=ifelse(mean_exp_mod_log_prob_5mc>mean_exp_mod_log_prob_5hmc,"m","h")]
dt_perReadModifiedBases_merged5mC_5hmC <- dt_perReadModifiedBases_5mC[dt_perReadModifiedBases_5hmC, on=.(chrm, pos, strand, number_reads), 
                    .(chrm, pos, strand, number_reads, mean_exp_mod_log_prob_5mc=mean_exp_mod_log_prob,
                    median_exp_mod_log_prob_5mc=median_exp_mod_log_prob, min_exp_mod_log_prob_5mc=min_exp_mod_log_prob, max_exp_mod_log_prob_5mc=max_exp_mod_log_prob,
                    mean_exp_mod_log_prob_5hmc=i.mean_exp_mod_log_prob, median_exp_mod_log_prob_5hmc=i.median_exp_mod_log_prob, min_exp_mod_log_prob_5hmc=i.min_exp_mod_log_prob,
                    max_exp_mod_log_prob_5hmc=i.max_exp_mod_log_prob), nomatch = NULL]
dt_perReadModifiedBases_merged5mC_5hmC
dt_perReadModifiedBases_5mC[dt_perReadModifiedBases_5hmC, on=.(chrm, pos, strand, number_reads), 
                    .(chrm, pos, strand, number_reads, mean_exp_mod_log_prob_5mc=mean_exp_mod_log_prob,
                    median_exp_mod_log_prob_5mc=median_exp_mod_log_prob, min_exp_mod_log_prob_5mc=min_exp_mod_log_prob, max_exp_mod_log_prob_5mc=max_exp_mod_log_prob,
                    mean_exp_mod_log_prob_5hmc=i.mean_exp_mod_log_prob, median_exp_mod_log_prob_5hmc=i.median_exp_mod_log_prob, min_exp_mod_log_prob_5hmc=i.min_exp_mod_log_prob,
                    max_exp_mod_log_prob_5hmc=i.max_exp_mod_log_prob), nomatch = NULL][,methylationState=ifelse(mean_exp_mod_log_prob_5mc>mean_exp_mod_log_prob_5hmc,"m","h")]
dt_perReadModifiedBases_5mC[dt_perReadModifiedBases_5hmC, on=.(chrm, pos, strand, number_reads), 
                    .(chrm, pos, strand, number_reads, mean_exp_mod_log_prob_5mc=mean_exp_mod_log_prob,
                    median_exp_mod_log_prob_5mc=median_exp_mod_log_prob, min_exp_mod_log_prob_5mc=min_exp_mod_log_prob, max_exp_mod_log_prob_5mc=max_exp_mod_log_prob,
                    mean_exp_mod_log_prob_5hmc=i.mean_exp_mod_log_prob, median_exp_mod_log_prob_5hmc=i.median_exp_mod_log_prob, min_exp_mod_log_prob_5hmc=i.min_exp_mod_log_prob,
                    max_exp_mod_log_prob_5hmc=i.max_exp_mod_log_prob), nomatch = NULL][,list(methylationState=ifelse(mean_exp_mod_log_prob_5mc>mean_exp_mod_log_prob_5hmc,"m","h"))]
dt_perReadModifiedBases_merged5mC_5hmC[,list(methylationState=ifelse(mean_exp_mod_log_prob_5mc>mean_exp_mod_log_prob_5hmc,"m","h"))]
dt_perReadModifiedBases_merged5mC_5hmC
dt_perReadModifiedBases_merged5mC_5hmC[,list(methylationState=ifelse(mean_exp_mod_log_prob_5mc > mean_exp_mod_log_prob_5hmc & mean_exp_mod_log_prob_5mc >=0.6,"m","h"))]
dt_perReadModifiedBases_merged5mC_5hmC[,list(methylationState=ifelse(mean_exp_mod_log_prob_5mc > mean_exp_mod_log_prob_5hmc & mean_exp_mod_log_prob_5mc >=0.6,"m","h"))][]
dt_perReadModifiedBases_merged5mC_5hmC[,list(methylationState=ifelse(mean_exp_mod_log_prob_5mc > mean_exp_mod_log_prob_5hmc & mean_exp_mod_log_prob_5mc >=0.6,"m","h"))][]
dim(dt_perReadModifiedBases_merged5mC_5hmC)
dt_perReadModifiedBases_merged5mC_5hmC
dt_perReadModifiedBases_merged5mC_5hmC[,list(methylationState=ifelse(mean_exp_mod_log_prob_5mc > mean_exp_mod_log_prob_5hmc & mean_exp_mod_log_prob_5mc >=0.6,"m","h")),]
dt_perReadModifiedBases_merged5mC_5hmC[,list(methylationState=ifelse(mean_exp_mod_log_prob_5mc >= 0.50,"m",ifelse(mean_exp_mod_log_prob_5mc >= 0.5, "h")))]
dt_perReadModifiedBases_merged5mC_5hmC[,list(methylationState=ifelse(mean_exp_mod_log_prob_5mc >= 0.50,"m",ifelse(mean_exp_mod_log_prob_5mc >= 0.5, "h", "UnMeth")))]
dt_perReadModifiedBases_merged5mC_5hmC[,list(methylationState=ifelse(mean_exp_mod_log_prob_5mc >= 0.50,"m",ifelse(mean_exp_mod_log_prob_5hmc >= 0.5, "h", "UnMeth")))]
dt_perReadModifiedBases_merged5mC_5hmC[,list(methylationState;=ifelse(mean_exp_mod_log_prob_5mc >= 0.50,"m",ifelse(mean_exp_mod_log_prob_5hmc >= 0.5, "h", "UnMeth")))]
dt_perReadModifiedBases_merged5mC_5hmC[,list(methylationState=ifelse(mean_exp_mod_log_prob_5mc >= 0.50,"m",ifelse(mean_exp_mod_log_prob_5hmc >= 0.5, "h", "UnMeth")))]
dt_perReadModifiedBases_merged5mC_5hmC
dt_perReadModifiedBases_merged5mC_5hmC[,list(methylationState=(ifelse(mean_exp_mod_log_prob_5mc >= 0.50,"m",ifelse(mean_exp_mod_log_prob_5hmc >= 0.5, "h", "UnMeth"))))]
dt_perReadModifiedBases_merged5mC_5hmC[,methylationState:=(ifelse(mean_exp_mod_log_prob_5mc >= 0.50,"m",ifelse(mean_exp_mod_log_prob_5hmc >= 0.5, "h", "UnMeth")))]
dt_perReadModifiedBases_merged5mC_5hmC[,methylationState:=ifelse(mean_exp_mod_log_prob_5mc >= 0.50,"m",ifelse(mean_exp_mod_log_prob_5hmc >= 0.5, "h", "UnMeth"))]
dt_perReadModifiedBases_merged5mC_5hmC
dt_perReadModifiedBases_merged5mC_5hmC[,methylationState := ifelse(mean_exp_mod_log_prob_5mc >= 0.50,"m",ifelse(mean_exp_mod_log_prob_5hmc >= 0.5, "h", "UnMeth"))]
dt_perReadModifiedBases_merged5mC_5hmC[,methylationState:= ifelse(mean_exp_mod_log_prob_5mc >= 0.50,"m",ifelse(mean_exp_mod_log_prob_5hmc >= 0.5, "h", "UnMeth"))]
dt_perReadModifiedBases_merged5mC_5hmC[,methylationState:= ifelse(mean_exp_mod_log_prob_5mc >= 0.50,"m","h")]
dt_perReadModifiedBases_merged5mC_5hmC[, methylationState:= ifelse(mean_exp_mod_log_prob_5mc >= 0.50,"m","h")]
dt_perReadModifiedBases_merged5mC_5hmC[, methylationState := ifelse(mean_exp_mod_log_prob_5mc >= 0.50,"m","h")]
dt_perReadModifiedBases_merged5mC_5hmC
dt_perReadModifiedBases_merged5mC_5hmC <- dt_perReadModifiedBases_5mC[dt_perReadModifiedBases_5hmC, on=.(chrm, pos, strand, number_reads), 
                    .(chrm, pos, strand, number_reads, mean_exp_mod_log_prob_5mc=mean_exp_mod_log_prob,
                    median_exp_mod_log_prob_5mc=median_exp_mod_log_prob, min_exp_mod_log_prob_5mc=min_exp_mod_log_prob, max_exp_mod_log_prob_5mc=max_exp_mod_log_prob,
                    mean_exp_mod_log_prob_5hmc=i.mean_exp_mod_log_prob, median_exp_mod_log_prob_5hmc=i.median_exp_mod_log_prob, min_exp_mod_log_prob_5hmc=i.min_exp_mod_log_prob,
                    max_exp_mod_log_prob_5hmc=i.max_exp_mod_log_prob), nomatch = NULL]
dt_perReadModifiedBases_merged5mC_5hmC[, methylationState := ifelse(mean_exp_mod_log_prob_5mc >= 0.50,"m","h")]
dt_perReadModifiedBases_merged5mC_5hmC
dt_perReadModifiedBases_merged5mC_5hmC <- dt_perReadModifiedBases_5mC[dt_perReadModifiedBases_5hmC, on=.(chrm, pos, strand, number_reads), 
                    .(chrm, pos, strand, number_reads, mean_exp_mod_log_prob_5mc=mean_exp_mod_log_prob,
                    median_exp_mod_log_prob_5mc=median_exp_mod_log_prob, min_exp_mod_log_prob_5mc=min_exp_mod_log_prob, max_exp_mod_log_prob_5mc=max_exp_mod_log_prob,
                    mean_exp_mod_log_prob_5hmc=i.mean_exp_mod_log_prob, median_exp_mod_log_prob_5hmc=i.median_exp_mod_log_prob, min_exp_mod_log_prob_5hmc=i.min_exp_mod_log_prob,
                    max_exp_mod_log_prob_5hmc=i.max_exp_mod_log_prob), nomatch = NULL]
dt_perReadModifiedBases_merged5mC_5hmC[,list(methylationState=(ifelse(mean_exp_mod_log_prob_5mc >= 0.50,"m",ifelse(mean_exp_mod_log_prob_5hmc >= 0.5, "h", "UnMeth"))))]
dt_perReadModifiedBases_merged5mC_5hmC
dt_perReadModifiedBases_merged5mC_5hmC[,methylationState:= ifelse(mean_exp_mod_log_prob_5mc >= 0.50,"m",ifelse(mean_exp_mod_log_prob_5hmc >= 0.5, "h", "UnMeth"))]
dt_perReadModifiedBases_merged5mC_5hmC
dt_perReadModifiedBases_merged5mC_5hmC

# this script is to draw ROC curve using score based on relative rank about correlation coefficient between residual and RBP expression
# made 2022/12/13

# activate package for drawing ROC curve
library(Epi)

# import list of treiber's physical interaction
# this list is located at "https://github.com/Ryosuke-Hirota/20221122_ROC_curve_of_pvalue_after_cutoff"
setwd("C:/Rdata/20221129_ROC_curve_of_pvalue_after_cutoff_revised")
phy.list <-read.table("list_of_treiber_physical_interaction_with_not_considering_primary_transcript.txt",sep="\t",header = T,stringsAsFactors = F)

# edit list of physical interaction 
phy.list[,5] <-paste0(phy.list[,3],"_vs_",phy.list[,2])
colnames(phy.list)[5] <-"combination"
phy.list <-phy.list[,c(5,4,1)]

# make list of score by each cutoff 
# these lists are located at "\\fsw-q02\okamura-lab\Files_related_to_M1_Projects\Hirota\CCLE_plot\20221213_calculate_score_about_relative_rank_of_correlation_between_residual_and_RBP_expression"
score.list <-list.files(path="C:/Rdata/20221213_calculate_score_about_relative_rank_of_correlation_between_residual_and_RBP_expression",pattern = "table_of_score")

# change order of score list
score.list <-score.list[c(1,10,2:9,11:19)]

# activate package for editing strings
library(stringr)

# set cutoff value
cutoff <-seq(0,900,50)

setwd("C:/Rdata")
dir.create("20221213_ROC_curve_of_score_about_correlation_between_residual_and_RBP_exp")

for (i in 1:length(score.list)) {
  setwd("C:/Rdata/20221213_calculate_score_about_relative_rank_of_correlation_between_residual_and_RBP_expression")
  
  # import each table 
  score.df <-read.table(score.list[i],sep="\t",header = T,stringsAsFactors = F)
  score.df[,c(5,6)] <-NA
  colnames(score.df)[5:6] <-c("primary_miRNA","transcript")
  
  # split name of combination and make new name in order to merge
  sp.comb <-str_split(score.df[,1],pattern = "_",simplify = T)
  combinations <-paste0(sp.comb[,5],"_vs_",sp.comb[,2])
  score.df[,1] <-combinations
  score.df[,5] <-sp.comb[,1]
  score.df[,6] <-sp.comb[,3]
  
  # merge table of score and list of physical interaction
  score.data <-merge(score.df,phy.list,by="combination",all=T)
  score.data <-subset(score.data,!is.na(score.data[,2]))
  
  # minus score -> plus score
  minus <-score.data[,4]<0 
  score.data[minus,4] <-score.data[minus,4]*-1
  
  # annotate whether combinations match with physical interaction
  score.data[,c(9,10)] <-NA
  colnames(score.data)[9:10] <-c("match","physical_interaction")
  
  score.data[!is.na(score.data[,7]),9] <-"match"
  score.data[is.na(score.data[,7]),9] <-"no_match"
  
  # combination with physical interaction as 1, combination without physical interaction as 0
  score.data[,7] <-ifelse(is.na(score.data[,7]),0,score.data[,7])
  
  score.data[score.data[,7]>3,10] <-1
  score.data[score.data[,7]<=3,10] <-0
  
  # count number of combinations with/without physical interaction
  # caution : this number is duplicated. If you wanna know non-duplicated number, write additional script.
  count <-as.data.frame(table(score.data[,9],score.data[,10]))
  p <-count[3,3]
  np <-count[1,3]
  
  # draw ROC curve
  setwd("C:/Rdata/20221213_ROC_curve_of_score_about_correlation_between_residual_and_RBP_exp")
  pdf(paste0("ROC_curve_of_score_cutoff_",cutoff[i],".pdf"))
  ROC(test=score.data$score, stat=score.data$physical_interaction, plot="ROC")
  mtext(text = paste0("physical interaction : ",p," , no physical interaction : ",np))
  dev.off()
}

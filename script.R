rm(list = ls())  
library(stringr)
library(data.table)
library(TwoSampleMR)
library(MRPRESSO)
library(dplyr)
library(ggsci)
library(ggplot2)
library(TSMRhelper)
library(forestploter)
library(grid)
library(stringr)
library(grDevices)
library(ieugwasr)
options(ieugwasr_api = 'gwas-api.mrcieu.ac.uk/')
setwd("")
# 结局
#2.本地结局
id<-""
outcome_all<-fread("./GCST90435416.tsv",data.table = F)

# 获取免疫细胞暴露
exp_clumped<-fread("exp_clumped(immune).csv")
# 获取结局
if(id==""){
  snp<-unique(exp_clumped$SNP)
  outcome<-outcome_all[outcome_all$rs_id%in%snp,]
  outcome_dat<-format_data(outcome,
                           type = "outcome",
                           beta_col = "beta",
                           se_col = "standard_error",
                           effect_allele_col = "effect_allele",
                           other_allele_col = "other_allele",
                           snp_col = "rs_id",
                           chr_col = "chromosome",
                           pos_col = "base_pair_location",
                           pval_col = "p_value",
                           eaf_col = "effect_allele_frequency")
}else{
  # 这是在线结局要修改的的数据
  while(TRUE){
    message_to_next <<- TRUE
    error_to_next <<- FALSE
    try({withCallingHandlers(outcome_dat <- extract_outcome_data(snps = exp_clumped$SNP,
                                                                 outcomes = id,
                                                                 proxies = F,
                                                                 access_token = NULL
    ), 
    message = function(c) if (stringr::str_detect(as.character(c),"Failed to")) message_to_next <<- FALSE)
      error_to_next <<- TRUE})
    if(message_to_next == TRUE&error_to_next == TRUE) { break }
  }
}
write.csv(outcome_dat,file = "outcome.csv",row.names = F)

# 合并数据
dat<-harmonise_data(exposure_dat = exp_clumped,outcome_dat = outcome_dat)
dat<-dat[dat$mr_keep,]
write.csv(dat,file = "dat.csv",row.names = F)

# 跑mr
res<-generate_odds_ratios(mr(dat))
write.csv(res,file = "res.csv",row.names = F)

# 异质性分析
heterogeneity<-mr_heterogeneity(dat)
write.csv(heterogeneity,file = "heterogeneity.csv",row.names = F)

# 多效性分析
pleiotropy_test<-mr_pleiotropy_test(dat)
write.csv(pleiotropy_test,file = "pleiotropy_test.csv",row.names = F)



# 绘图
res_id<-read.csv("./res.csv")%>%get_sbeta_res()%>%dplyr::filter(method=="Inverse variance weighted",pval<0.05)%>%pull(id.exposure)
dat_total<-read.csv("./dat.csv")
for (idr in res_id) {
  dat<-subset(dat_total,id.exposure==idr)
  scatter_plot<-mr_scatter_plot(mr(dat),dat)[[1]]+
    scale_color_lancet()+
    scale_fill_lancet()+
    theme(axis.title.y = element_text(size = 20))+
    theme_bw(base_size = 16)+ theme(
      plot.margin = margin(0.5,0.5,0.5,0.5, unit = "cm")
    )
  ggsave(plot=scatter_plot,filename = paste0(idr,"_scatter_plot_plot.pdf"),device = "pdf",width = 10,height = 10)
  funnel_plot<-mr_funnel_plot(mr_singlesnp(dat,all_method=c("mr_egger_regression","mr_weighted_median","mr_ivw","mr_simple_mode","mr_weighted_mode")))[[1]]+
    scale_color_lancet()+
    scale_fill_lancet()+
    theme(axis.title.y = element_text(size = 20))+
    theme_bw(base_size = 16)+ theme(
      plot.margin = margin(0.5,0.5,0.5,0.5, unit = "cm")
    )
  ggsave(plot=funnel_plot,filename = paste0(idr,"_funnel_plot_plot.pdf"),device = "pdf",width = 10,height = 10)
  forest_plot<-mr_forest_plot(mr_singlesnp(dat))[[1]]+
    scale_color_lancet()+
    scale_fill_lancet()+
    theme_bw()+
    theme(legend.position = 'none',plot.margin = margin(0.5,0.5,0.5,0.5, unit = "cm"))
  ggsave(plot=forest_plot,filename = paste0(idr,"_forest_plot_plot.pdf"),device = "pdf",width = 10,height = 10)
  
  leaveoneout_plot<-mr_leaveoneout_plot(mr_leaveoneout(dat))[[1]]+
    scale_color_lancet()+
    scale_fill_lancet()+
    theme_bw()+
    theme(legend.position = 'none',plot.margin = margin(0.5,0.5,0.5,0.5, unit = "cm"),
          axis.title.x = element_text(size = 12))
  ggsave(plot=leaveoneout_plot,filename = paste0(idr,"_leaveoneout_plot_plot.pdf"),device = "pdf",width = 10,height = 10)  
}



# 森林图
res_sign<-read.csv("./res.csv")%>%dplyr::filter(id.exposure%in%res_id)
res_sign$exposure<- res_sign$exposure
res_sign<-res_sign[res_sign$method=="Inverse variance weighted",]
res_sign$`OR (95% CI)` <- sprintf("%.3f (%.3f - %.3f)",res_sign$or, res_sign$or_lci95, res_sign$or_uci95)


res_sign$outcome[duplicated(res_sign$outcome)]<-""

dt<-res_sign[,c(3:6,9,12:14,15)]

# 用空格调整列宽
dt$` ` <- paste(rep(" ", 30), collapse = " ")

dt$pval<-round(dt$pval,digits = 3)
dt$pval<-ifelse(dt$pval<0.001,"<0.001",dt$pval)

colnames(dt)[5]<-"italic(P)*-Value"

dt<-dt%>%dplyr::rename(Outcome=outcome,
                       `Trait`=exposure,
                       Method=method,
                       nSNP=nsnp,)

dt$Method<-ifelse(dt$Method=="Inverse variance weighted","IVW",dt$Method)

tm <- forest_theme(base_size = 8,colhead=list(fg_params = list(parse=TRUE)))
p=forest(dt[,c(2:5,10,9)],#按顺序选择数据1-4列、12-13列为绘图区域和OR列的位置（图放中间）、8-11列作为森林图元素（Q、Q_pval和多效性P值）
         est = dt$or,
         lower = dt$or_lci95,
         upper = dt$or_uci95,
         sizes = 0.6,
         ci_column = 5,
         ref_line = 1,
         xlim = c(0.75,1.25),
         ticks_at = c(0.75,1,1.25),
         theme = tm)

ggsave("免疫细胞森林图.pdf",plot=p,device = "pdf",width = 20,height = 20)
# 安装并加载必要的R包
install.packages("TwoSampleMR")
install.packages("MVMR")
library(TwoSampleMR)
library(MVMR)

# 读取暴露（X）的GWAS数据
exposure_dat <- read_exposure_data(
  filename = "exposure_gwas.txt",
  sep = "\t",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval"
)

# 筛选独立SNP（LD剪切）
exposure_dat <- clump_data(exposure_dat, clump_r2 = 0.001)

# 读取中介变量（M）的GWAS数据
mediator_dat <- read_outcome_data(
  filename = "mediator_gwas.txt",
  snps = exposure_dat$SNP,
  sep = "\t",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval"
)

# 读取结果（Y）的GWAS数据
outcome_dat <- read_outcome_data(
  filename = "outcome_gwas.txt",
  snps = exposure_dat$SNP,
  sep = "\t",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval"
)

# 两步MR分析

# 第一步：暴露（X）对中介变量（M）的因果效应
dat_step1 <- harmonise_data(exposure_dat, mediator_dat)
res_step1 <- mr(dat_step1, method_list = c("mr_ivw"))

# 为中介变量（M）提取独立工具变量
mediator_exposure <- read_exposure_data(
  filename = "mediator_gwas.txt",
  sep = "\t",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval"
)
mediator_exposure <- clump_data(mediator_exposure, clump_r2 = 0.001)

# 第二步：中介变量（M）对结果（Y）的因果效应
dat_step2 <- harmonise_data(mediator_exposure, outcome_dat)
res_step2 <- mr(dat_step2, method_list = c("mr_ivw"))

# 计算中介效应（间接效应）
indirect_effect <- res_step1$b * res_step2$b

# 多变量MR（MVMR）分析

# 准备MVMR数据
mvmr_dat <- mv_harmonise_data(exposure_dat, list(mediator_dat), outcome_dat)

# 运行MVMR分析
mvmr_res <- mv_multiple(mvmr_dat)

# 提取直接效应
direct_effect <- mvmr_res$coef[1, "Estimate"]

# 计算总效应和中介效应
total_effect <- mr(harmonise_data(exposure_dat, outcome_dat), method_list = "mr_ivw")$b
indirect_effect_mvmr <- total_effect - direct_effect

# 输出结果
cat("两步MR结果：\n")
print(res_step1)
print(res_step2)
cat("间接效应（Indirect Effect）:", indirect_effect, "\n\n")

cat("MVMR结果：\n")
print(mvmr_res)
cat("直接效应（Direct Effect）:", direct_effect, "\n")
cat("间接效应（MVMR）:", indirect_effect_mvmr, "\n")


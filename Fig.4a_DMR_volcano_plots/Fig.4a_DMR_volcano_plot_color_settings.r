#!/usr/bin/env Rscript
# gx@lawlab

################################# frequently used
#################################
siren.color="#EA0A8C"
hyperTE.color="#5CC9E2"
clsy34.color="#D0554F"
taohong="#8c232a"
xiangye="#dfa041"
qunqing="#103472"
chenxiang="#423735"
ningzhi="#eaebd9"

################################# Guanghui's color palette
#################################
my.color=list(
	control="grey", 
	Control="grey",
	WT="grey", 
	WT1="grey", 
	WT2="grey", 
	WT3="grey", 
	col0="grey", 
	Col0="grey",
	clsy3="#7475AE", 
	c3="#7475AE", 
	an_clsy3="#7475AE", 
	an_c3="#7475AE", 
	ov_clsy3="#7475AE", 
	ov_c3="#7475AE", 
	fl_clsy3="#7475AE", 
	fl_c3="#7475AE", 
	clsy4="#BB6599", 
	c4="#BB6599", 
	an_clsy4="#BB6599", 
	an_c4="#BB6599", 
	ov_clsy4="#BB6599", 
	ov_c4="#BB6599", 
	fl_clsy4="#BB6599", 
	fl_c4="#BB6599", 
	clsy34="#D0554F", 
	c34="#D0554F", 
	an_clsy34="#D0554F", 
	an_c34="#D0554F", 
	ov_clsy34="#D0554F", 
	ov_c34="#D0554F", 
	fl_clsy34="#D0554F", 
	fl_c34="#D0554F", 
	poliv="#1F4497", 
	nrpd1="#1F4497", 
	an_poliv="#1F4497", 
	an_nrpd1="#1F4497", 
	ov_poliv="#1F4497", 
	ov_nrpd1="#1F4497", 
	fl_poliv="#1F4497", 
	fl_nrpd1="#1F4497", 
	clsy12="#F5AA6A", 
	c12="#F5AA6A", 
	an_clsy12="#F5AA6A", 
	an_c12="#F5AA6A", 
	ov_clsy12="#F5AA6A", 
	ov_c12="#F5AA6A", 
	fl_clsy12="#F5AA6A", 
	fl_c12="#F5AA6A", 
	rem16="#54C2B4", 
	rim16="#54C2B4", 
	an_rem16="#54C2B4", 
	an_rim16="#54C2B4", 
	ov_rem16="#54C2B4", 
	ov_rim16="#54C2B4", 
	fl_rem16="#54C2B4", 
	fl_rim16="#54C2B4", 
	rem22="#018282", 
	rim22="#018282",
	an_rem22="#018282", 
	an_rim22="#018282",
	ov_rem22="#018282", 
	ov_rim22="#018282",
	fl_rem22="#018282", 
	fl_rim22="#018282",
	EMS1146="#018282",
	ems1146="#018282",
	an_EMS1146="#018282",
	an_ems1146="#018282",
	ov_EMS1146="#018282",
	ov_ems1146="#018282",
	fl_EMS1146="#018282",
	fl_ems1146="#018282",
	cmf1="#EA90BC", 
	rim_cr="#EA90BC", 
	rem12="#EA90BC", 
	an_cmf1="#EA90BC", 
	an_rim_cr="#EA90BC", 
	an_rem12="#EA90BC", 
	ov_cmf1="#EA90BC", 
	ov_rim_cr="#EA90BC", 
	ov_rem12="#EA90BC", 
	fl_cmf1="#EA90BC", 
	fl_rim_cr="#EA90BC", 
	fl_rem12="#EA90BC", 
	rem16cmf1="purple", 
	rem22cmf1="brown", 
	rem16rem22="#3B54A4", 
	rem16rem22cmf1="#EC2526", 
	rim_quint="#EC2526", 
	rim_quintuple="#EC2526", 
	rem_quint="#EC2526", 
	rem_quintuple="#EC2526", 
	an_rem16rem22cmf1="#EC2526", 
	an_rim_quint="#EC2526", 
	an_rim_quintuple="#EC2526", 
	an_rem_quint="#EC2526", 
	an_rem_quintuple="#EC2526", 
	ov_rem16rem22cmf1="#EC2526", 
	ov_rim_quint="#EC2526", 
	ov_rim_quintuple="#EC2526", 
	ov_rem_quint="#EC2526", 
	ov_rem_quintuple="#EC2526", 
	fl_rem16rem22cmf1="#EC2526", 
	fl_rim_quint="#EC2526", 
	fl_rim_quintuple="#EC2526", 
	fl_rem_quint="#EC2526", 
	fl_rem_quintuple="#EC2526", 
	hyperTE="#5CC9E2", 
	HyperTE="#5CC9E2", 
	siren="#EA0A8C", 
	Siren="#EA0A8C",
	hyperTE_w_siren="purple",
	siren_w_hyperTE="purple",	
	an_LR60_ins1=taohong,
	an_LR60_ins2=taohong,
	an_C3OExLR60=xiangye,
	input=qunqing,
	`3xFLAG_REM12`=taohong,
	CLSY3_3xFLAG=xiangye,
	an_LR34_ins2=qunqing,
	NRPD1=taohong,
	NRPD1_cmf1=xiangye,
	ClusterDown="blue",
	ClusterUp="red",
	ClusterNone="grey",
	other="gray30"
	)

################################# Assign color function 
#################################
#usage: col = as.character(assign_colors(newGroup_level, my.color))
fallbacks = c(taohong, ningzhi, xiangye, qunqing, chenxiang)
assign_colors <- function(genotypes, color_list) {
  used_fallbacks <- list()
  fallback_index <- 1
  
  sapply(genotypes, function(gt) {
    if (gt %in% names(color_list)) {
      return(color_list[[gt]])
    } else {
      if (!gt %in% names(used_fallbacks)) {
        used_fallbacks[[gt]] <<- fallbacks[fallback_index]
        fallback_index <<- (fallback_index %% length(fallbacks)) + 1
      }
      return(used_fallbacks[[gt]])
    }
  })
}



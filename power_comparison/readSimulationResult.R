# load("C:/Users/Huan/Dropbox/0My research/codes/power_comparison/PowerComp100s100iter100bt.rda")

prefix = c('test_')
postfix = c('cox_Tcox','pl_Tcox', 'cox_Tpl', 'pl_Tpl')

var_names = c()
var_names_nc = c()
i = 1
for(pre in prefix){
  for (post in postfix){
    var_names[i] = paste0(pre,post)
    # var_names_nc[i] = paste0(pre,post,'_nc_simple')
    i = i + 1
  }
}

for(vn in var_names){
  cmd = paste0(vn,'<- c()')
  eval(parse(text = cmd))
}

# for(vn in var_names_nc){
#   cmd = paste0(vn,'<- c()')
#   eval(parse(text = cmd))
# }



for(i in 1:100){
  for (vn in var_names){
    cmd = paste0(vn ,'[i] <- as.numeric(TSLC[[i]]$', vn, '$p_value < 0.05)')
    eval(parse(text = cmd))
  }
}

# for(i in 1:100){
#   for (vn in var_names){
#     cmd = paste0(vn ,'[i] <- as.numeric(TSLC[[i]]$', vn, '$T_obs < TSLC[[i]]$', vn, '$T_critical)')
#     eval(parse(text = cmd))
#   }
# }




# for(i in 1:500){
#   for (vn in var_names_nc){
#     cmd = paste0(vn ,'[i] <- as.numeric(BHT_CB[[i]]$', vn, '$T_obs < BHT_CB[[i]]$', vn, '$T_critical)')
#     eval(parse(text = cmd))
#   }
# }

#Type I error or power
for(vn in var_names){
  cmd = paste0('res_',vn,' <- mean(',vn,')')
  eval(parse(text = cmd))
}

# for(vn in var_names_nc){
#   cmd = paste0('res_',vn,' <- mean(',vn,')')
#   eval(parse(text = cmd))
# }



library(xtable)
mytable = data.frame(row.names = postfix)
for (i in 1:length(postfix)) {
  for (j in 1:length(prefix)) {
    cmd = paste0('mytable[', i , ',', j,'] = res_', prefix[j],postfix[i])
    #print(cmd)
     eval(parse(text = cmd))
  }
}
colnames(mytable) = 'censoring'
xtable(mytable)

# mytable_nc = data.frame(row.names = postfix)
# for (i in 1:length(postfix)) {
#   for (j in 1:length(prefix)) {
#     cmd = paste0('mytable_nc[', i , ',', j,'] = res_', prefix[j],postfix[i],'_nc_simple')
#     #print(cmd)
#     eval(parse(text = cmd))
#   }
# }
# colnames(mytable_nc) = 'no censoring'
# xtable(mytable_nc)

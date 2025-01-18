prefix = c('test0','test1','test5','testm1','testm2','testm3','testm4')
postfix = c('exp1','exp05','gp011','gp12','wb51','wb11')

var_names = c()
i = 1
for(pre in prefix){
  for (post in postfix){
    var_names[i] = paste0(pre,'_nc_',post)
    i = i + 1
  }
}

for(vn in var_names){
  cmd = paste0(vn,'<- c()')
  eval(parse(text = cmd))
}

# for(i in 1:500){
#   for (vn in var_names){
#     cmd = paste0(vn ,'[i] <- BHT_cnc[[i]]$', vn, '$p_value')
#     eval(parse(text = cmd))
#   }
# }

for(i in 1:500){
  for (vn in var_names){
    cmd = paste0(vn ,'[i] <- as.numeric(BHT_ncnc[[i]]$', vn, '$T_obs < BHT_ncnc[[i]]$', vn, '$T_critical)')
    eval(parse(text = cmd))
  }
}

# for(i in 1:500){
#   for (vn in var_names){
#     cmd = paste0(vn ,'[i] <- BHT_cnc[[i]]$', vn, '$reject')
#     eval(parse(text = cmd))
#   }
# }

#Type I error or power
for(vn in var_names){
  cmd = paste0('res_',vn,' <- mean(',vn,' )')
  eval(parse(text = cmd))
}


library(xtable)
mytable = data.frame(row.names = postfix)
for (i in 1:length(postfix)) {
  for (j in 1:length(prefix)) {
    cmd = paste0('mytable[', i , ',', j,'] = res_', prefix[j], '_nc_',postfix[i])
    #print(cmd)
    eval(parse(text = cmd))
  }
}
colnames(mytable) = prefix
xtable(mytable,digits=c(0,3,3,3,3, 3, 3, 3))

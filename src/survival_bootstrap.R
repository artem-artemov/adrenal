library(ggplot2)
library(survival)
library(cowplot)
library(survminer)
library(ggfortify)


"
dat=readRDS('bulk_survival/GSE49711.dat.rds')

exp=dat$dat
ano=dat$ano

ano$time=as.numeric(ano$OS)
ano$vital_status=ano$OS.event
ano$inss=ano$inss.stage

index=ano$inss!='4S'
ano=ano[index,]
exp=exp[,index]
"

surv=function(gs, exp, ano, stage=c(3,4), mycn=0){
  #gs='SOX10' 
  #stage=NULL  # filter data by disease stage  
  #mycn=NULL
  title=''
  

  print(gs)
    
      res=NULL
    
  if(startsWith(gs, 'fate_')){
      res = ano
      res$exp = ano[[gs]]
    #res=data.frame('exp'=as.numeric(ano[,gs]),ano)  
  }else{
    if(!(gs %in% rownames(exp)))
      return()
  
    res=data.frame('exp'=as.numeric(exp[gs,]),ano)
}
    
  res=res[!is.na(res$exp),]
  
  #res=res[res$mycn.status!=1,]
  # print(res[1,])
  
  if (!is.null(stage)){
    res=res[res$inss %in% stage,]
  }
  
  
  if (!is.null(mycn)){
    res=res[!is.na(res$mycn.status) & res$mycn.status==mycn,]
    print(table(res$mycn.status))
  }
  
  
  cut1=quantile(res$exp,0.75)
  cut2=quantile(res$exp,0.25)
  
  res1=res[res$exp>=cut1,]
  res2=res[res$exp<=cut2,]
  
  
  res2$x='Low-exp'
  res1$x='High-exp'
  res=rbind(res1,res2)
  
  print(dim(res))
    
    fit=NULL
      fit <- survfit(Surv(time, vital_status) ~ x, data = res)
  #pvalue=survdiff(Surv(time, vital_status)~x, data=res)
  #p2=round(1-pchisq(pvalue$chisq,1),4)
  
  pv = surv_pvalue(fit, data=res)
    print(pv)
  
  p0=NULL
  p0=ggsurvplot(fit, data=res, pval = T,linetype = c("solid", "dashed"), #
                palette = c("red","blue"),title=title, #paste(gs,p2),
                legend.title="",legend.labs=c("High-exp","Low-exp"),
                conf.int = F)
    #risk.table = F
  
  #,legend=c(0.7,0.2)
  #print(p0)
  # p0=ggsurvplot(fit, pval = TRUE)
  #  plot(p0)
    p0
}

surv_bootstrap = function(genes_by_fate, exp, ano, Nbootstrap=100, mycn=0){
    lres = list()
    for(fate in names(genes_by_fate)){
        
        ano[[paste0('fate_', fate)]] = as.vector(colSums(exp[rownames(exp) %in% genes_by_fate[[fate]],]))
        df = surv(paste0('fate_', fate), exp, ano, mycn=mycn)$data.survplot
        df$iteration=0
        df$set = 'all'
        #lres[[fate]][[0]] = 
        for(i in 1:Nbootstrap){
            selrows = which(rownames(exp) %in% genes_by_fate[[fate]])
            selrows_bootstrap = sample(selrows, replace = T)
            ano[[paste0('fate_bootstrap_', fate)]] = as.vector(colSums(exp[selrows_bootstrap,]))
            dfi = surv(paste0('fate_bootstrap_', fate), exp, ano, mycn=mycn)$data.survplot
            dfi$iteration=i
            dfi$set = 'bootstrap'
            df = rbind(df, dfi)
        }
        lres[[fate]] = df
    }
    lres
}

#lres = surv_bootstrap(genes_by_fate, exp, ano, Nbootstrap=30, mycn=mycn)
#p = ggplot(lres[["chromaffin"]], aes(time, surv, col=x))+
#        geom_step(aes(group = paste(x, iteration), 
#                      alpha=ifelse(set=='bootstrap', 0.8,1)))+
#        theme_bw()+theme(legend.position = "none")+ylab('Survival probability')+xlab('Time')
#        print(p)

tmp1 = function(){
genes_by_fate = list()
    for(fate in unique(env_medulla5$SR.markers.AUC.fate$cluster)){
        print(fate)
        genes_cluster = rownames(env_medulla5$SR.markers.AUC.fate[env_medulla5$SR.markers.AUC.fate$cluster==fate,])
        genes_by_fate[[fate]] = setdiff(genes_cluster, genes_cell_cycle)

        ano[[paste0('fate_', fate)]] = as.vector(colSums(exp[rownames(exp) %in% genes_by_fate[[fate]],]))#as.vector(colSums(exp[rownames(exp) %in% genes_cluster,]))
        #hist(ano[[paste0('fate_', fate)]], main=length(genes_by_fate[[fate]]))
    }

    fit=NULL
    p0 = NULL
    for(fate in unique(env_medulla5$SR.markers.AUC.fate$cluster)){
        p = surv(paste0('fate_', fate), mycn=0)
        #ggsave(paste0('figures/figS2_survival/survival.',gene,'.pdf'), print(p[[2]]))
        #ggsave(paste0('figures/figS2_survival/survival_fate.mycn0.',fate,'.pdf'), p$plot, width=5, height=5)
        ggsave(paste0('figures/figS2_survival/survival_fate.mycn0.noCC.',fate,'.pdf'), p$plot, width=5, height=5)
        #plot(p)
        p = surv(paste0('fate_', fate), mycn=1)
        ggsave(paste0('figures/figS2_survival/survival_fate.mycn1.noCC.',fate,'.pdf'), p$plot, width=5, height=5)
    }
}
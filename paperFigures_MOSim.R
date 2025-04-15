
library(ggplot2)
library(dplyr)
library(gridExtra)

## Results simulations

######## Results for MLR + EN ########

# Take the data from MLR

setwd("D:/resultados_definitivos/MLR/")
load('res_MLR.Rdata')

res = resultsf1acu
res$scale = as.character(res$scale)
res$scale<-factor(res$scale,levels = c("auto", "softBlock","hardBlock"))

df_summary <- res %>%
  group_by(n, scale,thres) %>%
  summarise(
    mean_f1score = mean(f1score),
    sd_f1score = sd(f1score),  
    n_obs = 12 
  ) %>%
  mutate(
    se_f1score = sd_f1score / sqrt(n_obs),  
    lower_ci = mean_f1score - 1.96 * se_f1score,  
    upper_ci = mean_f1score + 1.96 * se_f1score  
  )

ggplot(df_summary, aes(x = n, y = mean_f1score, color = scale, fill = scale, group=scale)) +
  geom_line(linewidth = 1) +  
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), alpha = 0.2) +  
  scale_y_continuous(breaks = c(0.65, 0.7, 0.75, 0.8)) +
  scale_fill_biostat(palette = "main") +
  scale_color_biostat(palette = "main") +
  labs(x = "Number of observations", y = "F1-score",title='Results MLR+EN-MF', color='scaleType',fill='scaleType') +
  theme_minimal() + facet_wrap(~thres)+
  theme(
    text = element_text(size = 14),
    legend.position = "right"  
  )

######## Results for MLR + ISGL - MF ########

setwd("D:/resultados_definitivos/isglresults_COR/resultados/")
load('res_tableISGLCOR.Rdata')

res = results
res$scale = as.character(res$scale)
res$scale<-factor(res$scale,levels = c("auto", "softBlock","hardBlock"))

df_summary <- res %>%
  group_by(n, scale,thres) %>%
  summarise(
    mean_f1score = mean(f1score),
    sd_f1score = sd(f1score), 
    n_obs = 12  
  ) %>%
  mutate(
    se_f1score = sd_f1score / sqrt(n_obs), 
    lower_ci = mean_f1score - 1.96 * se_f1score, 
    upper_ci = mean_f1score + 1.96 * se_f1score  
  )

ggplot(df_summary, aes(x = n, y = mean_f1score, color = scale, fill = scale, group=scale)) +
  geom_line(linewidth = 1) + 
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), alpha = 0.2) + 
  scale_fill_biostat(palette = "main") +
  scale_color_biostat(palette = "main") +
  labs(x = "Number of observations", y = "F1-score",title='Results MLR+ISGL-MF', color='scaleType',fill='scaleType') +
  theme_minimal() + facet_wrap(~thres)+
theme(
  text = element_text(size = 14),
  legend.position = "right" 
)

######## Results for MLR + ISGL - PCA ########

setwd("D:/resultados_definitivos/isglresults_PCA/resultados/")
load('res_tableISGLPCA.Rdata')

setwd("D:/resultados_definitivos/isglresults_PCA_lowthres/")
load('res_table_PCAlow.Rdata')

res = rbind(results,resultsf1acu)
res$scale = as.character(res$scale)
res$scale<-factor(res$scale,levels = c("auto", "softBlock","hardBlock"))

df_summary <- res %>%
  group_by(n, scale,thres) %>%
  summarise(
    mean_f1score = mean(f1score),
    sd_f1score = sd(f1score), 
    n_obs = 12  
  ) %>%
  mutate(
    se_f1score = sd_f1score / sqrt(n_obs),  
    lower_ci = mean_f1score - 1.96 * se_f1score, 
    upper_ci = mean_f1score + 1.96 * se_f1score  
  )

ggplot(df_summary, aes(x = n, y = mean_f1score, color = scale, fill = scale, group=scale)) +
  geom_line(linewidth = 1) +  
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), alpha = 0.2) +  
  scale_fill_biostat(palette = "main") +
  scale_color_biostat(palette = "main") +
  labs(x = "Number of observations", y = "F1-score",title='Results MLR+ISGL-PCA', color='scaleType',fill='scaleType') +
  theme_minimal() + facet_wrap(~thres)+
theme(
  text = element_text(size = 14),
  legend.position = "none" 
)

## Results for PLS1 ####

setwd("D:/resultados_definitivos/PLS/")
load('res_tablePLS.Rdata')

res = resultsf1acu
res$scale = as.character(res$scale)
res$scale<-factor(res$scale,levels = c("auto", "softBlock","hardBlock"))

res$pval = as.character(res$pval)
res$pval[which(res$pval=='perm')] = 'Permutation'
res$pval[which(res$pval=='jack')] = 'Jack-Knife'

df_summary <- res %>%
  group_by(n, scale,pval) %>%
  summarise(
    mean_f1score = mean(f1score),
    sd_f1score = sd(f1score),  
    n_obs = 12  
  ) %>%
  mutate(
    se_f1score = sd_f1score / sqrt(n_obs),  
    lower_ci = mean_f1score - 1.96 * se_f1score,  
    upper_ci = mean_f1score + 1.96 * se_f1score   
  )

ggplot(df_summary, aes(x = n, y = mean_f1score, color = scale, fill = scale, group=scale)) +
  geom_line(linewidth = 1) +  
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), alpha = 0.2) +  
  scale_fill_biostat(palette = "main") +
  scale_color_biostat(palette = "main") +
  labs(x = "Number of observations", y = "F1-score",title='Results PLS1', color='scaleType',fill='scaleType') +
  theme_minimal() + facet_wrap(~pval)+
theme(
  text = element_text(size = 14),
  legend.position = "right"  
)

## Results for PLS2 ####

setwd("D:/resultados_definitivos/PLS2_TF_mirna/resultados")
load('res_table_PLS2.Rdata')

res = resultsf1acu

setwd("D:/resultados_definitivos/PLS2_TF_mirna/PLS_jack/resultados")
load('res_table_PLS2Jack.Rdata')

res = rbind(res,resultsf1acu)

res$scale = as.character(res$scale)
res$scale<-factor(res$scale,levels = c("auto", "softBlock","hardBlock"))

res$f1score<-as.numeric(res$f1score)
res$R2<-as.numeric(res$R2)
res$time<-as.numeric(res$time)

res$ngenes<-res$n
res$ngenes <- as.character(res$ngenes)

res$ngenes[which(res$ngenes=='5')]<-'~2000'
res$ngenes[which(res$ngenes=='10')]<-'~4250'
res$ngenes[which(res$ngenes=='20')]<-'~8000'
res$ngenes[which(res$ngenes=='30')]<-'~9250'
res$ngenes[which(res$ngenes=='50')]<-'~10500'
res$ngenes<-factor(res$ngenes,levels = c('~2000', '~4250','~8000','~9250','~10500'))

#Relizar las gráficas comparativas mediante curvas

res$varsel = as.character(res$varsel)
res$varsel[which(res$varsel=='Perm')] = 'Permutation'
res$varsel[which(res$varsel=='Jack')] = 'Jack-Knife'

df_summary <- res %>%
  group_by(n, scale,varsel) %>%
  summarise(
    mean_f1score = mean(f1score),
    sd_f1score = sd(f1score), 
    n_obs = 12  
  ) %>%
  mutate(
    se_f1score = sd_f1score / sqrt(n_obs), 
    lower_ci = mean_f1score - 1.96 * se_f1score,  
    upper_ci = mean_f1score + 1.96 * se_f1score  
  )

ggplot(df_summary, aes(x = n, y = mean_f1score, color = scale, fill = scale, group=scale)) +
  geom_line(linewidth = 1) + 
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), alpha = 0.2) + 
  scale_fill_biostat(palette = "main") +
  scale_color_biostat(palette = "main") +
  labs(x = "Number of observations", y = "F1-score",title='Results PLS2', color='scaleType',fill='scaleType') +
  theme_minimal() + facet_wrap(~varsel)+
theme(
  text = element_text(size = 14),
  legend.position = "right"  
)

## Comparison MORE-PLS1 vs KiMONo-----------

setwd('D:/resultados_definitivos/kimono')
load('results_kim.Rdata')

setwd("D:/resultados_definitivos/PLS_comparar/")
load('results_PLScomp.Rdata')

res = rbind(results, resultsf1acu)

res$met[which(res$met=='MORE-PLS')]='MORE-PLS1'

res$n<-factor(res$n,levels = c('5', '10','20','30','50'))

res$ngenes<-res$n
res$ngenes <- as.character(res$ngenes)

res$ngenes[which(res$ngenes=='5')]<-'~2000'
res$ngenes[which(res$ngenes=='10')]<-'~4250'
res$ngenes[which(res$ngenes=='20')]<-'~8000'
res$ngenes[which(res$ngenes=='30')]<-'~9250'
res$ngenes[which(res$ngenes=='50')]<-'~10500'
res$ngenes<-factor(res$ngenes,levels = c('~2000', '~4250','~8000','~9250','~10500'))

res$f1score = as.numeric(res$f1score)
res$time = as.numeric(res$time)
res$R2 = as.numeric(res$R2)

df_summary <- res %>%
  group_by(n, met) %>%
  summarise(
    mean_f1score = mean(f1score),
    sd_f1score = sd(f1score),  
    n_obs = 12  
  ) %>%
  mutate(
    se_f1score = sd_f1score / sqrt(n_obs), 
    lower_ci = mean_f1score - 1.96 * se_f1score, 
    upper_ci = mean_f1score + 1.96 * se_f1score   
  )

p1<-ggplot(df_summary, aes(x = n, y = mean_f1score, color = met, fill = met, group = met)) +
  geom_line(linewidth = 1) + 
  geom_point(data = df_summary, aes(x = n, y = mean_f1score, color = met)) +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), alpha = 0.07, color = NA) +  
  scale_fill_biostat(palette = "main")+
  scale_color_biostat(palette = "main")+
  ggtitle("Comparison F1-score") +
  labs(x="Number of observations", y = "F1-score")+
  theme_minimal() +
  theme(axis.text = element_text(size = 14),  
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.position = "none"
  )

df_summary <- res %>%
  group_by(n, met) %>%
  summarise(
    mean_f1score = mean(R2),
    sd_f1score = sd(R2),  
    n_obs = 12 
  ) %>%
  mutate(
    se_f1score = sd_f1score / sqrt(n_obs),  
    lower_ci = mean_f1score - 1.96 * se_f1score,  
    upper_ci = mean_f1score + 1.96 * se_f1score   
  )

p2<-ggplot(df_summary, aes(x = n, y = mean_f1score, color = met, fill = met, group = met)) +
  geom_line(linewidth = 1) + 
  geom_point(data = df_summary, aes(x = n, y = mean_f1score, color = met)) +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), alpha = 0.07, color = NA) +  
  scale_fill_biostat(palette = "main")+
  scale_color_biostat(palette = "main")+
  ggtitle("Comparison R2") +
  labs(x="Number of observations", y = "R2")+
  theme_minimal() +
  theme(axis.text = element_text(size = 14),  
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.position = "none"
  )


df <- res %>%
  mutate(ngenes = as.numeric(gsub("~", "", ngenes)))

df_summary <- res %>%
  group_by(ngenes, met) %>%
  summarise(
    mean_time = mean(time),
    sd_time = sd(time), 
    n_obs = 12  
  ) %>%
  mutate(
    se_time = sd_time / sqrt(n_obs),  
    lower_ci = mean_time - 1.96 * se_time,  
    upper_ci = mean_time + 1.96 * se_time   
  )

p3<-ggplot(df_summary, aes(x = ngenes, y = log(mean_time), color = met, fill = met, group = met)) +
  geom_line(linewidth = 1) + 
  geom_point(data = df_summary, aes(x = ngenes, y = log(mean_time), color = met)) +
  geom_ribbon(aes(ymin = log(lower_ci), ymax = log(upper_ci)), alpha = 0.07, color=NA) + 
  scale_fill_biostat(palette = "main")+
  scale_color_biostat(palette = "main")+
  ggtitle("Comparison computational time") +
  labs(x="Number of genes", y = "Time in minutes", color ='Method', fill='Method')+ 
  theme_minimal() +
  theme(axis.text = element_text(size = 14),  
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.position = "none"
  ) + scale_y_continuous(labels = function(x) round(exp(x)))

grid.arrange(p1,p2, p3, nrow = 1, widths=c(1,1,1))



## Final comparison #####

#En la comparativa se añaden los mejores settings de cada approach: 

# - MLR+EN-MF con auto escalado y 0.7 de threshold
# - MLR+ISGL-MF con auto escalado y 0.7 de threshold
# - MLR+ISGL-PCA con auto escalado y 20% de variabilidad explicada
# - PLS1 con auto escalado y Jack-Knife como técnica para el cálculo de la significancia estadística


setwd("D:/resultados_definitivos/MLR_final")
load('res_MLR.Rdata')

res = resultsf1acu

setwd("D:/resultados_definitivos/isglresults_COR/resultados/")
load('res_tableISGLCORcomp.Rdata')

res = rbind(res, results)

setwd("D:/resultados_definitivos/isglresults_PCA_lowthres/")
load('res_table_PCAcomp.Rdata')

res = rbind(res, resultsf1acu)

setwd("D:/resultados_definitivos/PLS_comparar/")
load('results_PLScomp.Rdata')

resultsf1acu = resultsf1acu[,c(1,6,2,4,3,5,7)]
res = rbind(res, resultsf1acu)

res$ngenes<-res$n
res$ngenes <- as.character(res$ngenes)

res$ngenes[which(res$ngenes=='5')]<-'~2000'
res$ngenes[which(res$ngenes=='10')]<-'~4250'
res$ngenes[which(res$ngenes=='20')]<-'~8000'
res$ngenes[which(res$ngenes=='30')]<-'~9250'
res$ngenes[which(res$ngenes=='50')]<-'~10500'
res$ngenes<-factor(res$ngenes,levels = c('~2000', '~4250','~8000','~9250','~10500'))

res$f1score = as.numeric(res$f1score)
res$time = as.numeric(res$time)
res$R2 = as.numeric(res$R2)

###### Figure 2A ----

df_summary <- res %>%
  group_by(n, met) %>%
  summarise(
    mean_f1score = mean(f1score),
    sd_f1score = sd(f1score), 
    n_obs = 12 
  ) %>%
  mutate(
    se_f1score = sd_f1score / sqrt(n_obs), 
    lower_ci = mean_f1score - 1.96 * se_f1score,  
    upper_ci = mean_f1score + 1.96 * se_f1score  
  )



p1<-ggplot(df_summary, aes(x = n, y = mean_f1score, color = met, fill = met, group = met)) +
  geom_line(linewidth = 1) + 
  geom_point(data = df_summary, aes(x = n, y = mean_f1score, color = met)) +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), alpha = 0.07, color = NA) + 
  scale_fill_biostat(palette = "main")+
  scale_color_biostat(palette = "main")+
  ggtitle("Comparison F1-score") +
  labs(x="Number of observations", y = "F1-score")+
  theme_minimal() +
  theme(axis.text = element_text(size = 14),  
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
    legend.position = "none"
  )

df_summary <- res %>%
  group_by(n, met) %>%
  summarise(
    mean_R2 = mean(R2),
    sd_R2 = sd(R2),  
    n_obs = 12  
  ) %>%
  mutate(
    se_R2 = sd_R2 / sqrt(n_obs), 
    lower_ci = mean_R2 - 1.96 * se_R2, 
    upper_ci = mean_R2 + 1.96 * se_R2  
  )



p2<-ggplot(df_summary, aes(x = n, y = mean_R2, color = met, fill = met, group = met)) +
  geom_line(linewidth = 1) +  
  geom_point(data = df_summary, aes(x = n, y = mean_R2, color = met)) +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), alpha = 0.07, color = NA) + 
  scale_fill_biostat(palette = "main")+
  scale_color_biostat(palette = "main")+
  ggtitle("Comparison R2") +
  labs(x="Number of observations", y = "R2")+
  theme_minimal() +
  theme(axis.text = element_text(size = 14),  
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.position = "none"
  )

df <- res %>%
  mutate(ngenes = as.numeric(gsub("~", "", ngenes)))

df_summary <- res %>%
  group_by(ngenes, met) %>%
  summarise(
    mean_time = mean(time),
    sd_time = sd(time),  
    n_obs = 12  
  ) %>%
  mutate(
    se_time = sd_time / sqrt(n_obs), 
    lower_ci = mean_time - 1.96 * se_time,
    upper_ci = mean_time + 1.96 * se_time   
  )

p3<-ggplot(df_summary, aes(x = ngenes, y = log(mean_time), color = met, fill = met, group = met)) +
  geom_line(linewidth = 1) +  
  geom_point(data = df_summary, aes(x = ngenes, y = log(mean_time), color = met)) +
  geom_ribbon(aes(ymin = log(lower_ci), ymax = log(upper_ci)), alpha = 0.07, color=NA) +  
  scale_fill_biostat(palette = "main")+
  scale_color_biostat(palette = "main")+
  ggtitle("Comparison computational time") +
  labs(x="Number of genes", y = "Time in minutes", color ='Method', fill='Method')+ 
  theme_minimal() +
  theme(axis.text = element_text(size = 14),  
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
    legend.position = "none"
  ) + scale_y_continuous(labels = function(x) round(exp(x)))


grid.arrange(p1,p2, p3, nrow = 1, widths=c(1,1,1))

###### Figure 2B ----

df_summary <- res %>%
  group_by(p, met) %>%
  summarise(
    mean_f1score = mean(f1score),
    sd_f1score = sd(f1score),  
    n_obs = 30  
  ) %>%
  mutate(
    se_f1score = sd_f1score / sqrt(n_obs),  
    lower_ci = mean_f1score - 1.96 * se_f1score, 
    upper_ci = mean_f1score + 1.96 * se_f1score   
  )



p1<-ggplot(df_summary, aes(x = p, y = mean_f1score, color = met, fill = met, group = met)) +
  geom_line(linewidth = 1) +  
  geom_point(data = df_summary, aes(x = p, y = mean_f1score, color = met)) +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), alpha = 0.07, color = NA) +  
  scale_fill_biostat(palette = "main")+
  scale_color_biostat(palette = "main")+
  ggtitle("Comparison F1-score") +
  labs(x="% simulated significant regulations", y = "F1-score")+
  theme_minimal() +
  theme(axis.text = element_text(size = 14),  
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.position = "none"
  )

df_summary <- res %>%
  group_by(p, met) %>%
  summarise(
    mean_R2 = mean(R2),
    sd_R2 = sd(R2), 
    n_obs = 30
  ) %>%
  mutate(
    se_R2 = sd_R2 / sqrt(n_obs), 
    lower_ci = mean_R2 - 1.96 * se_R2,  
    upper_ci = mean_R2 + 1.96 * se_R2   
  )


p2<-ggplot(df_summary, aes(x = p, y = mean_R2, color = met, fill = met, group = met)) +
  geom_line(linewidth = 1) +  
  geom_point(data = df_summary, aes(x = p, y = mean_R2, color = met)) +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), alpha = 0.07, color = NA) +  
  scale_fill_biostat(palette = "main")+
  scale_color_biostat(palette = "main")+
  scale_y_continuous(limits = c(0.6, 1)) +
  ggtitle("Comparison R2") +
  labs(x="% simulated significant regulations", y = "R2")+
  theme_minimal() +
  theme(axis.text = element_text(size = 14),  
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.position = "none"
  )



df <- res %>%
  mutate(ngenes = as.numeric(gsub("~", "", ngenes)))

df_summary <- res %>%
  group_by(p, met) %>%
  summarise(
    mean_time = mean(time),
    sd_time = sd(time),
    n_obs = 30 
  ) %>%
  mutate(
    se_time = sd_time / sqrt(n_obs),
    lower_ci = mean_time - 1.96 * se_time, 
    upper_ci = mean_time + 1.96 * se_time  
  )

p3<-ggplot(df_summary, aes(x = p, y = log(mean_time), color = met, fill = met, group = met)) +
  geom_line(linewidth = 1) +  
  geom_point(data = df_summary, aes(x = p, y = log(mean_time), color = met)) +
  geom_ribbon(aes(ymin = log(lower_ci), ymax = log(upper_ci)), alpha = 0.07, color=NA) +  
  scale_fill_biostat(palette = "main")+
  scale_color_biostat(palette = "main")+
  ggtitle("Comparison computational time") +
  labs(x="% simulated significant regulations", y = "Time in minutes", color ='Method', fill='Method')+ 
  theme_minimal() +
  theme(axis.text = element_text(size = 14),  
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.position = "none"
  ) + scale_y_continuous(labels = function(x) round(exp(x)))


grid.arrange(p1,p2, p3, nrow = 1, widths=c(1,1,1))


###### Figure 2C ----

df_summary <- res %>%
  group_by(pDEG, met) %>%
  summarise(
    mean_f1score = mean(f1score),
    sd_f1score = sd(f1score),  
    n_obs = 30 
  ) %>%
  mutate(
    se_f1score = sd_f1score / sqrt(n_obs),  
    lower_ci = mean_f1score - 1.96 * se_f1score,
    upper_ci = mean_f1score + 1.96 * se_f1score   
  )


p1<-ggplot(df_summary, aes(x = pDEG, y = mean_f1score, color = met, fill = met, group = met)) +
  geom_line(linewidth = 1) +  
  geom_point(data = df_summary, aes(x = pDEG, y = mean_f1score, color = met)) +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), alpha = 0.07, color = NA) + 
  scale_fill_biostat(palette = "main")+
  scale_color_biostat(palette = "main")+
  ggtitle("Comparison F1-score") +
  labs(x="% simulated DEG", y = "F1-score")+
  theme_minimal() +
  theme(axis.text = element_text(size = 14),  
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.position = "none"
  )

df_summary <- res %>%
  group_by(pDEG, met) %>%
  summarise(
    mean_R2 = mean(R2),
    sd_R2 = sd(R2), 
    n_obs = 30  
  ) %>%
  mutate(
    se_R2 = sd_R2 / sqrt(n_obs),  
    lower_ci = mean_R2 - 1.96 * se_R2,  
    upper_ci = mean_R2 + 1.96 * se_R2   
  )

p2<-ggplot(df_summary, aes(x = pDEG, y = mean_R2, color = met, fill = met, group = met)) +
  geom_line(linewidth = 1) + 
  geom_point(data = df_summary, aes(x = pDEG, y = mean_R2, color = met)) +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), alpha = 0.07, color = NA) + 
  scale_fill_biostat(palette = "main")+
  scale_color_biostat(palette = "main")+
  ggtitle("Comparison R2") +
  labs(x="% simulated DEG", y = "R2")+
  theme_minimal() +
  theme(axis.text = element_text(size = 14),  
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.position = "none"
  )

df <- res %>%
  mutate(ngenes = as.numeric(gsub("~", "", ngenes)))

df_summary <- res %>%
  group_by(pDEG, met) %>%
  summarise(
    mean_time = mean(time),
    sd_time = sd(time), 
    n_obs = 30  
  ) %>%
  mutate(
    se_time = sd_time / sqrt(n_obs),  
    lower_ci = mean_time - 1.96 * se_time, 
    upper_ci = mean_time + 1.96 * se_time  
  )

p3<-ggplot(df_summary, aes(x = pDEG, y = log(mean_time), color = met, fill = met, group = met)) +
  geom_line(linewidth = 1) +  
  geom_point(data = df_summary, aes(x = pDEG, y = log(mean_time), color = met)) +
  geom_ribbon(aes(ymin = log(lower_ci), ymax = log(upper_ci)), alpha = 0.07, color=NA) + 
  scale_fill_biostat(palette = "main")+
  scale_color_biostat(palette = "main")+
  ggtitle("Comparison computational time") +
  labs(x="% simulated DEG", y = "Time in minutes", color ='Method', fill='Method')+ 
  theme_minimal() +
  theme(axis.text = element_text(size = 14),  
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.position = "none"
  ) + scale_y_continuous(labels = function(x) round(exp(x)))


grid.arrange(p1,p2, p3, nrow = 1, widths=c(1,1,1))


## Comparison of MORE-PLS1 + MORE-PLS2 with other tools -----------

setwd('D:/resultados_definitivos/PLS2_TF_mirna/resultados/')
load('res_table_PLS2.Rdata')
resultsf1acu = resultsf1acu[resultsf1acu$scale=='hardBlock',]

resultsf1acu$met = 'MORE-PLS2'
resultados = resultsf1acu

setwd('D:/resultados_definitivos')
load('res_kim_mix_more.Rdata')

res = res[res$met!='mixOmics',-8]

resultados = rbind(res, resultados[,c(1,3,5,4,6,2,7)])

setwd('D:/resultados_definitivos/mixOmics')
load('mix_table.RData')

results = results[,c(1,3,5,4,6,2,7)]

resultados = rbind(resultados, results)

res = resultados

res$pDEG[which(res$pDEG=='20')]='50'

res$ngenes<-res$n
res$ngenes <- as.character(res$ngenes)

res$ngenes[which(res$ngenes=='5')]<-'~2000'
res$ngenes[which(res$ngenes=='10')]<-'~4250'
res$ngenes[which(res$ngenes=='20')]<-'~8000'
res$ngenes[which(res$ngenes=='30')]<-'~9250'
res$ngenes[which(res$ngenes=='50')]<-'~10500'
res$ngenes<-factor(res$ngenes,levels = c('~2000', '~4250','~8000','~9250','~10500'))

res$f1score = as.numeric(res$f1score)
res$time = as.numeric(res$time)
res$R2 = as.numeric(res$R2)

res$met[which(res$met=='MORE-PLS')] = 'MORE-PLS1'

# Calculate the mean, min, and max of F1-score for each combination of ngenes and met

###### Figure 3 ------------

df_summary <- res %>%
  group_by(n, met) %>%
  summarise(
    mean_f1score = mean(f1score),
    sd_f1score = sd(f1score),  
    n_obs = 12  
  ) %>%
  mutate(
    se_f1score = sd_f1score / sqrt(n_obs),  
    lower_ci = mean_f1score - 1.96 * se_f1score,  
    upper_ci = mean_f1score + 1.96 * se_f1score   
  )

p1<-ggplot(df_summary, aes(x = n, y = mean_f1score, color = met, fill = met, group = met)) +
  geom_line(linewidth = 1) +  
  geom_point(data = df_summary, aes(x = n, y = mean_f1score, color = met)) +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), alpha = 0.07, color = NA) + 
  scale_fill_biostat(palette = "main")+
  scale_color_biostat(palette = "main")+
  ggtitle("Comparison F1-score") +
  labs(x="Number of observations", y = "F1-score")+
  theme_minimal() +
  theme(axis.text = element_text(size = 14),  
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.position = "none"
  )

# Calculate the mean, min, and max of F1-score for each combination of ngenes and met

df_summary <- res %>%
  group_by(n, met) %>%
  summarise(
    mean_R2 = mean(R2),
    sd_R2 = sd(R2),  
    n_obs = 12  
  ) %>%
  mutate(
    se_R2 = sd_R2 / sqrt(n_obs), 
    lower_ci = mean_R2 - 1.96 * se_R2,  
    upper_ci = mean_R2 + 1.96 * se_R2 
  )

p2<-ggplot(df_summary, aes(x = n, y = mean_R2, color = met, fill = met, group = met)) +
  geom_line(linewidth = 1) +  
  geom_point(data = df_summary, aes(x = n, y = mean_R2, color = met)) +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), alpha = 0.07, color = NA) + 
  scale_fill_biostat(palette = "main")+
  scale_color_biostat(palette = "main")+
  ggtitle("Comparison R2") +
  labs(x="Number of observations", y = "R2")+
  theme_minimal() +
  theme(axis.text = element_text(size = 14),  
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.position = "none"
  )

df <- res %>%
  mutate(ngenes = as.numeric(gsub("~", "", ngenes)))

df_summary <- res %>%
  group_by(ngenes, met) %>%
  summarise(
    mean_time = mean(time),
    sd_time = sd(time), 
    n_obs = 12 
  ) %>%
  mutate(
    se_time = sd_time / sqrt(n_obs), 
    lower_ci = mean_time - 1.96 * se_time,  
    upper_ci = mean_time + 1.96 * se_time   
  )

p3<-ggplot(df_summary, aes(x = ngenes, y = log(mean_time), color = met, fill = met, group = met)) +
  geom_line(linewidth = 1) +  
  geom_point(data = df_summary, aes(x = ngenes, y = log(mean_time), color = met)) +
  geom_ribbon(aes(ymin = log(lower_ci), ymax = log(upper_ci)), alpha = 0.07, color=NA) + 
  scale_fill_biostat(palette = "main")+
  scale_color_biostat(palette = "main")+
  ggtitle("Comparison computational time") +
  labs(x="Number of genes", y = "Time in minutes", color ='Method', fill='Method')+ 
  theme_minimal() +
  theme(axis.text = element_text(size = 14),  
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.position = "none"
  ) + scale_y_continuous(labels = function(x) round(exp(x)))

grid.arrange(p1, p2, p3, nrow = 1, widths=c(1,1,1))

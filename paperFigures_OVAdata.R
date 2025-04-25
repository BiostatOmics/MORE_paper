
# Analysis from the main paper -------
load('OVAdata.Rdata')

library(MORE)

## Run MORE ----

PLSresults = more(targetData = GeneExpresion,
                    associations=assoc,
                    regulatoryData=data.omics,
                    condition=edesign,
                    omicType = 0,
                    scaleType = 'auto',
                    interactions = TRUE,
                    minVariation=0, vip = 0.8, varSel = 'Jack',
                    method = 'PLS1')

saveRDS(PLSresults, file='results.rds') #Save the results in case they are necessary in the future

summary(PLSresults)

output_regpcond = RegulationPerCondition(PLSresults)

## Create the summary plots ---------
summaryPlot(PLSresults,output_regpcond, filterR2 = 0.5)
summaryPlot(PLSresults,output_regpcond, filterR2 = 0.5, byTargetF = FALSE)
output_regpcond = FilterRegulationPerCondition(PLSresults, output_regpcond, filterR2 = 0.5) #Filter to R2>0.5
differentialRegPlot(output = PLSresults, outputRegpcond = output_regpcond)

#Analize the regulations in each of the experimental groups under study.
output_regincond_dif = RegulationInCondition(output_regpcond, 'differentiated')
output_regincond_im = RegulationInCondition(output_regpcond, 'immunoreactive')
output_regincond_mes = RegulationInCondition(output_regpcond, 'mesenchymal')
output_regincond_pro = RegulationInCondition(output_regpcond, 'proliferative')

## Create the networks  ----

#Create the network for the differentiated subtype

networkMORE(outputRegpcond = output_regpcond, group1 = 'differentiated')

#Show the networks for BRCA2

output_regpcond_BRCA2 = output_regpcond[output_regpcond$targetF %in% c('BRCA2'),, drop=FALSE]

networkMORE(outputRegpcond = output_regpcond_BRCA2, group1 = 'differentiated') #the differentiated network
networkMORE(outputRegpcond = output_regpcond_BRCA2, group1 = 'proliferative') #the proliferative network
networkMORE(outputRegpcond = output_regpcond_BRCA2, group1 = 'differentiated', group2 = 'proliferative') # the differential network

## Perform GSEA ------

# Create the annotation file from AnnotationDBi

library(org.Hs.eg.db)
library(GO.db)
library(AnnotationDbi)
library(KEGGREST)

# Create the annotation from KEGG

head(keys(org.Hs.eg.db, keytype = 'SYMBOL'))
k = keys(org.Hs.eg.db, keytype = 'SYMBOL')
cl = c('PATH','SYMBOL')
Annotation = AnnotationDbi::select(org.Hs.eg.db,keys = k,columns = cl,keytype = 'SYMBOL')

pathw = keggList("pathway")
names(pathw) = gsub('map','',names(pathw))
pathw = as.data.frame(pathw)
pathw$PATH = rownames(pathw)

Annotation = merge(Annotation, pathw,by='PATH',all.x=TRUE)

Annotation = Annotation[,c(2,1,3)]

#Create the GSEA
y = gseaMORE(output_regincond_dif,output_regincond_pro, Annotation, alpha = 0.05, p.adjust.method= 'none')
dotplot(y,  split=".sign") +  facet_grid(.~.sign, labeller = as_labeller(c(activated = "proliferative", suppressed = "differentiated")))


# Analysis for the supplementary materials ------------

#Prior analysis of the hub genes and global regulators
#Hub genes

hub_genes = list(dif=output_regincond_dif$HubTargetF, im = output_regincond_im$HubTargetF,
                 mes = output_regincond_mes$HubTargetF, pro = output_regincond_pro$HubTargetF)

x <- list(Differentiated = output_regincond_dif$HubTargetF , Immunoreactive = output_regincond_im$HubTargetF,
          Mesenchymal = output_regincond_mes$HubTargetF, Proliferative = output_regincond_pro$HubTargetF)

library(ggVennDiagram)
ggVennDiagram(x , category.names = names(x),order.set.by = 'name',force_upset = TRUE)

#Upload data from top genes in Ovarian Cancer to see if we find any among the results found in MORE
library(readxl)
top_genes <- read_excel("top_genes.xlsx")

top_genes<-as.vector(top_genes[,1])[[1]]

intersect(top_genes, output_regincond_dif$HubTargetF)
intersect(top_genes, output_regincond_im$HubTargetF)
intersect(top_genes, output_regincond_mes$HubTargetF)
intersect(top_genes, output_regincond_pro$HubTargetF)

# Venn diagram representation of the global regulators per condition

#Global regulators
x <- list(Differentiated = output_regincond_dif$GlobalRegulators , Immunoreactive = output_regincond_im$GlobalRegulators,
          Mesenchymal = output_regincond_mes$GlobalRegulators, Proliferative = output_regincond_pro$GlobalRegulators)

ggVennDiagram(x , label_alpha=0,  set_color = color_palette[-1], set_size = 3.5)+scale_fill_gradient(low = "#FFFAFA", high = "#8B475D")
ggVennDiagram(x , force_upset = TRUE)

#Unique Global regulators

output_regincond_pro$GlobalRegulators[!(output_regincond_pro$GlobalRegulators %in% c(output_regincond_dif$GlobalRegulators,output_regincond_im$GlobalRegulators, output_regincond_mes$GlobalRegulators))]
output_regincond_im$GlobalRegulators[!(output_regincond_im$GlobalRegulators %in% c(output_regincond_dif$GlobalRegulators,output_regincond_pro$GlobalRegulators, output_regincond_mes$GlobalRegulators))]
output_regincond_mes$GlobalRegulators[!(output_regincond_mes$GlobalRegulators %in% c(output_regincond_dif$GlobalRegulators,output_regincond_pro$GlobalRegulators, output_regincond_im$GlobalRegulators))]

## Create the UpsetPlots -------

listRegincond = list(Differentiated = output_regincond_dif , Immunoreactive = output_regincond_im,
                     Mesenchymal = output_regincond_mes, Proliferative = output_regincond_pro)

upsetMORE(listRegincond) #UpsetPlot for the Hub genes
upsetMORE(listRegincond, byHubs = FALSE) #UpsetPlot for the Hub genes

## Create the plot of the profiles of the pointed regulators of BRCA2 -----

plotMORE(output = results_lowcv, targetF = 'BRCA2', regulator = 'TF-STAT4', xlab = '')
plotMORE(output = results_lowcv, targetF = 'BRCA2', regulator = 'TF-NR4A1', xlab = '')
plotMORE(output = results_lowcv, targetF = 'BRCA2', regulator = 'TF-FOXM1', xlab = '')
plotMORE(output = results_lowcv, targetF = 'BRCA2', regulator = 'TF-BRIP1', xlab = '')



library(tidyr)
library(dplyr)
library(ggplot2)
library(cowplot)

AllEffects <- coeffTabfinal %>% 
  mutate(VarEffect = rownames(selected.model$parameters$betaStandXWTable)) %>% 
  tidyr::separate(VarEffect, into = c("DepVar", "Treat"), sep = "_") %>%
  #tidyr::separate(DepVar, into = c("Phylum", "Class")) %>% 
  # mutate(DepVar = ifelse(DepVar == "HAl", "H_Al", DepVar)) %>% #View
  mutate(DepVar = factor(DepVar, levels = unique(DepVar)))

TaxInfo.T.All <- TaxInfo.T %>% 
  mutate(TaxAbbrGJAM = gsub(" ", ".", TaxAbbr)) %>% 
  mutate(TaxAbbrGJAM = gsub("_", "", TaxAbbrGJAM)) %>% 
  mutate(TaxAbbrGJAM = gsub("-", ".", TaxAbbrGJAM)) %>%
  mutate(TaxAbbrGJAM = gsub("\\(", ".", TaxAbbrGJAM))

labelsTabT <- data.frame(TaxAbbrGJAM = unique(AllEffects$DepVar)) %>% 
  left_join(TaxInfo.T.All, by = "TaxAbbrGJAM") %>% 
  mutate(finalLabel = ifelse(!is.na(Genus), as.character(Genus), 
                             ifelse(!is.na(Family), as.character(Family),
                                    ifelse(!is.na(Order), as.character(Order),
                                           ifelse(!is.na(Class), as.character(Class),as.character(Phylum))))))


SankeyTab <- AllEffects %>% 
  mutate(sig95 = factor(sig95), Treat = factor(Treat)) %>% 
  mutate(TaxAbbrGJAM = as.character(DepVar)) %>% 
  left_join(labelsTabT, by = "TaxAbbrGJAM") %>% 
  mutate(Groups = ifelse(DepVar %in% bact, "bact",
                         ifelse(DepVar %in% fungi, "fungi",
                                ifelse(DepVar %in% soil, "soil",
                                       ifelse(DepVar %in% biom, "biom", "Error"))))) %>% #View
  mutate(Groups = ifelse(Domain == "Archaea", "archaea", as.character(Groups))) %>% #View
  mutate(Groups = #ifelse(DepVar %in% bact, "bact",
           ifelse(DepVar %in% fungi, "fungi",
                  ifelse(DepVar %in% soil, "soil",
                         ifelse(DepVar %in% biom, "biom", as.character(Groups))))) %>% #View
  mutate(Groups = factor(Groups, levels = c("bact", "fungi","archaea", "biom", "soil"))) %>% 
  mutate(Treat = gsub("LandUse", "", Treat)) %>% 
  mutate(Treat = factor(Treat, levels = c("FLO", "CapA", "CapM", "CapB", "CapEn", "SAFC", "SAFQ"))) %>% 
  #filter(Treat == varSelected) %>% 
  mutate(rank = rank(Estimate, ties.method = "first")) %>%
  ungroup() %>% 
  mutate(rank2 = 1:nrow(.)) %>% 
  mutate(DepVarAbb = factor(finalLabel, levels = unique(finalLabel))) %>% 
  mutate(Effect = ifelse(Estimate > 0, "pos",
                         ifelse(Estimate < 0, "neg", "Error"))) %>% 
  mutate(EffectInt = paste0(sig95,Effect)) %>% #View
  select(Treat,Groups,Domain:Species,DepVar,EffectInt) %>% #View
  # unique() %>% 
  pivot_wider(names_from = Treat, values_from = EffectInt) %>% 
  mutate(Phylum = as.character(Phylum),Groups = as.character(Groups)) %>% 
  mutate(Phylum = ifelse(is.na(Phylum), Groups, Phylum))

View(SankeyTab)

testres1 <- SankeyDiagram(SankeyTab[,c(1,3,10:16)], 
              max.categories = length(unique(SankeyTab$Phylum)),#This is the number of Phylum I have, you may choose a different one according to your case
              label.show.varname = FALSE,
              label.show.counts = TRUE,
              link.color = "None", #variables.share.values = TRUE,
              hovertext.show.percentages = TRUE,
)

testres2 <- SankeyPlotMF(SankeyTab[,c(1,3,10:16)],
             max.categories = length(unique(SankeyTab$Phylum)),#This is the number of Phylum I have, you may choose a different one according to your case
             label.show.varname = FALSE,
             label.show.counts = TRUE,
             link.color = "First variable", #variables.share.values = TRUE,
             hovertext.show.percentages = TRUE,
)

testres3 <- SankeyPlotMF(SankeyTab[,c(1,3,10:16)],
                         max.categories = length(unique(SankeyTab$Phylum)),#This is the number of Phylum I have, you may choose a different one according to your case
                         label.show.varname = FALSE,
                         label.show.counts = TRUE,
                         link.color = "Source", #variables.share.values = TRUE,
                         hovertext.show.percentages = TRUE,
)

testres1
testres2
testres3
testres$x$options$colourScale
testres1$x$links
testres2$x$links
testres3$x$links
View(cbind(testres1$x$links,testres2$x$links,testres3$x$links))
testres2$x$nodes <- testres1$x$nodes
testres2$x$nodes$name <- testres2$x$nodes$group
testres2

testres3$x$links

testres2$x$options$colourScale

testres3$x$links$group[1:21]
length(testres2$x$links$group)

testres$x$options$colourScale
testres$sizingPolicy


testres2$x$links$group <- testres3$x$links$group

testres2

p <- plot_ly()
p <- add_trace(p, link=trace1$link, node=trace1$node, type=trace1$type, domain=trace1$domain, orientation=trace1$orientation, valueformat=trace1$valueformat)
p <- layout(p, font=layout$font, title=layout$title, height=layout$height)

  # mutate(DepVarAbb = factor(abbreviate(DepVar, 10), levels = unique(abbreviate(DepVar, 10)))) %>% 
  # ggplot(aes(x = reorder(DepVarAbb,Phylum), y = Estimate, color = Groups))+
  # geom_point()+
  # geom_errorbar(aes(ymin=CI_025, ymax=CI_975), width=.1)+
  # ylab(varSelected)+
  # geom_hline(yintercept=0, linetype="dashed")+
  # coord_flip()+
  # scale_color_manual(values = unique(colorVec))+
  # facet_grid(Groups+Phylum~Treat,scales = "free",labeller=LU.names, 
  #            space = "free")+
  # theme_cowplot(font_size = 8)

#Alluvial plot
install.packages("ggalluvial")
library(ggalluvial)
DataAlluv <- SankeyTab[,c(1,3,10:16)] %>% 
  group_by(Groups, Phylum,FLO,CapA,CapM,CapB,CapEn,SAFC,SAFQ
           ) %>% 
  tally %>% 
  mutate(Groups = factor(Groups, levels =c("biom", "soil","archaea", "bact" ,    "fungi")),
         # Phylum = factor(Phylum, levels = unique(Phylum)),
         FLO = factor(FLO, levels = c("*pos", "pos", "neg", "*neg")),
         CapA = factor(CapA, levels = c("*pos", "pos", "neg", "*neg")),
         CapM = factor(CapM, levels = c("*pos", "pos", "neg", "*neg")),
         CapB = factor(CapB, levels = c("*pos", "pos", "neg", "*neg")),
         CapEn = factor(CapEn, levels = c("*pos", "pos", "neg", "*neg")),
         SAFC = factor(SAFC, levels = c("*pos", "pos", "neg", "*neg")),
         SAFQ = factor(SAFQ, levels = c("*pos", "pos", "neg", "*neg"))
         )
levels(DataAlluv$Groups) <- c("biom", "soil","archaea", "bact" ,    "fungi")
levelsPhy <- c(unique(DataAlluv$Phylum[DataAlluv$Groups == "archaea"]),
  unique(DataAlluv$Phylum[DataAlluv$Groups == "bact"]),
  unique(DataAlluv$Phylum[DataAlluv$Groups == "fungi"]),
  unique(DataAlluv$Phylum[DataAlluv$Groups == "biom"]),
  unique(DataAlluv$Phylum[DataAlluv$Groups == "soil"])
  )

DataAlluv$Phylum <- factor(DataAlluv$Phylum, levels = levelsPhy)
library(ggrepel)
DataAlluv %>% 
  arrange(Groups, FLO,Phylum,n) %>% 
  ggplot(aes(y = n, axis1 = factor(Groups,levels = c("biom", "soil","archaea", "bact" ,    "fungi")), 
             axis2 = factor(Phylum, levels = levelsPhy),axis3 = FLO, axis4 = CapA,
             axis5 = CapM,axis6 =CapB,axis7 =CapEn,axis8 =SAFC,axis9 =SAFQ)) +
  geom_alluvium(aes(fill = Groups)# width = 1/12
                # width = 0, #knot.pos = 1/4, #reverse = FALSE
                ) +
  geom_stratum(width = 0.2, #alpha = 0.1, 
               fill = "gray", 
               color = "black") +
  # geom_label(stat = "stratum", aes(label = after_stat(stratum)), size = 3, hjust = -0.5) +
  geom_text_repel(stat = "stratum", aes(label = after_stat(stratum)), size = 3, hjust = -0.8)+
  scale_fill_manual(values = c("darkgreen", "brown", "red", "blue", "orange"))+
  scale_x_discrete(limits = colnames(DataAlluv)[1:9]) +
  # facet_grid(Groups~., scales = "free_y",space ="free_y")+
  theme(#axis.text.x = element_text(angle = 90, hjust = 1),
        panel.background = element_rect(fill = "white", colour = "grey50", linetype = "solid"))
  # scale_fill_brewer(type = "qual", palette = "Set1") +
  # ggtitle("UC Berkeley admissions and rejections, by sex and department")


#SankeyPlot for Coefficients
devtools::install_github("Displayr/flipPlots")
library(flipPlots)
library(tidyr)
library(dplyr)

tHF
#For H0F1
SankeyDFH0F1 <- tHF %>% 
  left_join(TaxInfo, by = "TaxVar") %>% #check if the TaxVar column in TaxInfo has the same name as in the tHF
  mutate(MeanDiff  = Mean - H0F1,
         SignH0F1 = ifelse(SignH0F1 == "", "ns", SignH0F1),
         ClassSign = ifelse(MeanDiff>0, "pos", ifelse(MeanDiff<0, "neg", "zero"))) %>% 
  mutate(ClassSign = paste0(SignH0F1,ClassSign)) %>% 
  dplyr::select(TaxVar, Fact1, Fact2,ClassSign,Kingdom, Phylum) %>% 
  mutate(Kingdom = ifelse(is.na(Kingdom), TaxVar, Kingdom), #If you have other types of variable you will need to add this group information
         Phylum = ifelse(is.na(Phylum), TaxVar, Phylum)) %>% 
  # dplyr::select(-SignH0F2) %>% 
  pivot_wider(names_from = Fact2, values_from = ClassSign) %>% 
  filter(Fact1 == lvF1[1]) #In case you want to filter the Factor 1 to show only one treatment at time
SankeyDiagram(SankeyDFH0F1[,-1],
              max.categories = length(unique(SankeyDF$Phylum)),#This is the number of Phylum I have, you may choose a different one according to your case
              label.show.varname = FALSE,
              label.show.counts = TRUE,
              hovertext.show.percentages = TRUE,
)

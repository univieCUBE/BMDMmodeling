#### Main R script for Computational modelling reveals and predicts flux differences in nucleotide metabolism in IL-4-treated Cyp27a1-KO bone marrow-derived macrophages ####

library(readxl)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(dplyr)
library(reshape2)
library(ggvenn)
library(VennDiagram)
library(gridExtra)
library(ggpubr)
library(ggallin)
library(stringr)
library(magrittr)
library(pheatmap)
require(grid)
library(tidyverse)


my.theme <- theme(axis.text = element_text(colour="black", size=15, face = "bold"),
                  text = element_text(size=16),
                  panel.background = element_rect(fill = 'gray99',
                                                  colour = "black",
                                                  linewidth=0.5),
                  axis.title.x=  element_text(vjust=-0.45, face = "bold"),
                  axis.title.y = element_text(vjust=1.2, face = "bold"),
                  axis.ticks = element_line(colour="black"),
                  axis.line = element_line(),
                  panel.grid.major = element_line(colour = "lightgray", linetype="dotted"),
                  panel.grid.minor = element_line(colour = "lightgray", linetype="dashed"),
                  #legend.title=element_text(),
                  legend.text = element_text(size = 14, face = "bold"),
                  legend.title = element_text(size = 14, face = "bold"))

#Load reaction objects:
ControlRXN = read.csv('/Users/rokosango/PhD/MetabModelling/MATLAB_scripts_workspaces/Immunometabolism_models/tINIT_models/OnlyReactionsFromModels/ControlModel.csv',
                      header = T)
LPSRXN = read.csv('/Users/rokosango/PhD/MetabModelling/MATLAB_scripts_workspaces/Immunometabolism_models/tINIT_models/OnlyReactionsFromModels/LPSModel.csv',
                  header = T)
IL4RXN = read.csv('/Users/rokosango/PhD/MetabModelling/MATLAB_scripts_workspaces/Immunometabolism_models/tINIT_models/OnlyReactionsFromModels/IL4Model.csv',
                  header = T)
reactions = read.csv('/Users/rokosango/PhD/MetabModelling/MATLAB_scripts_workspaces/reactions.csv')

#### Getting proper media compositions values for settings uptake/secretion rates ####
setwd("/Users/rokosango/PhD/MetabModelling/MATLAB_scripts_workspaces/GEM_Comparison_R")
options(scipen=999)
MediaCompDf = read.csv('MediaCompDf.csv', header = T, row.names = 1)

MediaCompDf = MediaCompDf %>%
  mutate(MetPercentage = (Concentration..mg.L. / sum(Concentration..mg.L.)) / 10 ) %>% #changing mg/L to g/L and then taking percentage in one step
  mutate(MetWeight = MetPercentage * 3) %>%
  mutate(FluxConstraint_mmol_per_mouse_day = (MetWeight / Molecular.Weight) * 1000)

write.csv(MediaCompDf, 'MediaCompDf.csv')


#Venn diagram of Constrained GEMs 
setwd("~/PhD/MetabModelling/MATLAB_scripts_workspaces/Immunometabolism_models/tINIT_models")


LPS = read.csv('/Users/rokosango/PhD/MetabModelling/MATLAB_scripts_workspaces/Immunometabolism_models/tINIT_models/Venn/LPS_GEM_rxns.csv',
               header = T)
IL4 = read.csv('/Users/rokosango/PhD/MetabModelling/MATLAB_scripts_workspaces/Immunometabolism_models/tINIT_models/Venn/IL4_GEM_rxns.csv',
               header = T)
Control = read.csv('/Users/rokosango/PhD/MetabModelling/MATLAB_scripts_workspaces/Immunometabolism_models/tINIT_models/Venn/Control_GEM_rxns.csv',
                   header = T)

myCol <- c("#56B4E9",
           "#E69F00",
           "#984EA3")


venn.plot = venn.diagram(
  x = list(as.character(LPS$NAME), as.character(IL4$NAME), as.character(Control$NAME)),
  category.names = c("LPS" , "IL4" , "Control"),
  filename = NULL,
  output=TRUE,
  # Circles
  lwd = 2,
  lty = 1,
  fill = myCol,
  # Numbers
  cex = 2,
  fontface = "bold",
  fontfamily = "sans",
  # Set names
  cat.cex = 2,
  cat.fontface = "bold",
  cat.col = myCol,
  cat.fontfamily = "serif",
  rotation = 1,
  scaled = TRUE
)

grid.draw(venn.plot)


x <- list(
  LPS = as.character(LPS$NAME), 
  IL4 = as.character(IL4$NAME), 
  Control = as.character(Control$NAME)
)

ggvenn(
  x, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 4) + ggtitle("Genes in GEMs")


#### ReporterMetabolites analysis ####

par(mfrow=c(3,2))

#### Curated FBA ####

file.choose()

options(scipen = 999)

path = '/Users/rokosango/PhD/MetabModelling/MATLAB_scripts_workspaces/Immunometabolism_models/tINIT_models/newFluxResults/'

LPSGapFilled = read_excel(paste0(path, "LPSFluxTable.xlsx"),
                          sheet = 1)
IL4GapFilled = read_excel(paste0(path, "IL4FluxTable.xlsx"),
                          sheet = 1)
CtrlGapFilled = read_excel(paste0(path, "ControlFluxTable.xlsx"),
                           sheet = 1)


prepFluxTables = function(FluxTable, allReactions, modelSpecificRxns) {
  
  # FluxTable = FluxTable %>%
  #   dplyr::select(ReactionID, Flux)
  
  
  FluxTable$ReactionName = allReactions$rxnRecon3DID[match(FluxTable$ReactionID, allReactions$rxns)]
  FluxTable$Flux = as.numeric(FluxTable$Flux)
  FluxTable$Subsystem = modelSpecificRxns$SUBSYSTEM[match(FluxTable$ReactionID, modelSpecificRxns$ID)]
  FluxTable$Subsystem = as.factor(FluxTable$Subsystem)
  FluxTable$Equation = modelSpecificRxns$EQUATION[match(FluxTable$ReactionID, modelSpecificRxns$ID)]
  
  return(FluxTable)
  
}

LPSGapFilled = prepFluxTables(LPSGapFilled, reactions, LPSRXN)
IL4GapFilled = prepFluxTables(IL4GapFilled, reactions, IL4RXN)
CtrlGapFilled = prepFluxTables(CtrlGapFilled, reactions, ControlRXN)


LPSGapFilled[6870, c("Subsystem", "Equation")] = list("Exchange/demand reactions", "NO[s] <=>")
LPSGapFilled[6869, c("Subsystem", "Equation")] = list("Transport reactions", "NO[c] <=> NO[s]")
LPSGapFilled[5153, "ReactionName"] = "NOS2"

#### RelativeContribution for Fig 2E ####

PathwaysToKeep = c('Glycolysis / Gluconeogenesis',
                   'Pentose phosphate pathway',
                   'Pyruvate metabolism',
                   'Nucleotide metabolism',
                   'Arginine and proline metabolism',
                   'Purine metabolism',
                   'Pyrimidine metabolism',
                   'Fatty acid biosynthesis',
                   'Pyrimidine metabolism',
                   'Fatty acid oxidation',
                   'Oxidative phosphorylation',
                   'Tricarboxylic acid cycle and glyoxylate/dicarboxylate metabolism')

FunRelativeContrib = function(LPSModel, IL4Model, ControlModel, title) {
  
LPSFluxSum = LPSModel  %>%
  dplyr::group_by(Subsystem) %>%
  dplyr::summarize(AbsSum = sum(abs(Flux), na.rm = T))

IL4FluxSum = IL4Model %>%
  dplyr::group_by(Subsystem) %>%
  dplyr::summarize(AbsSum = sum(abs(Flux), na.rm = T))

CtrlFluxSum = ControlModel %>%
  dplyr::group_by(Subsystem) %>%
  dplyr::summarize(AbsSum = sum(abs(Flux), na.rm = T))

CommonSubsystems = intersect(intersect(LPSFluxSum$Subsystem, IL4FluxSum$Subsystem),
                             CtrlFluxSum$Subsystem)

LPSFluxSum = LPSFluxSum  %>%
  dplyr::filter(Subsystem %in% CommonSubsystems) %>%
  dplyr::filter(!Subsystem == "NA")

IL4FluxSum = IL4FluxSum %>%
  dplyr::filter(Subsystem %in% CommonSubsystems) %>%
  dplyr::filter(!Subsystem == "NA")

CtrlFluxSum = CtrlFluxSum %>%
  dplyr::filter(Subsystem %in% CommonSubsystems) %>%
  dplyr::filter(!Subsystem == "NA")

TotalFlux = data.frame(Subsystem = LPSFluxSum$Subsystem,
                       LPSFlux = LPSFluxSum$AbsSum,
                       IL4Flux = IL4FluxSum$AbsSum,
                       CtrlFlux = CtrlFluxSum$AbsSum)


TotalFlux = TotalFlux %>%
  dplyr::filter(Subsystem %in% PathwaysToKeep)
 
TotalFlux$Subsystem = str_replace_all(TotalFlux$Subsystem,
                                         "Tricarboxylic acid cycle and glyoxylate/dicarboxylate metabolism",
                                         "TCA cycle")

TotalFlux$SumAllModels = apply(TotalFlux[c(2,3,4)], MARGIN = 1, sum)

TotalFlux = TotalFlux[!rowSums(TotalFlux[c(2,3,4)]) <= 1e-1,] #if the sum across different GEMs is < 1e-2
TotalFluxRelContribution = TotalFlux %>%
  mutate(LPSRelContribution = (LPSFlux * 1) / SumAllModels) %>%
  mutate(IL4RelContribution = (IL4Flux * 1) / SumAllModels) %>%
  mutate(CtrlRelContribution = (CtrlFlux * 1) / SumAllModels) %>%
  dplyr::select(Subsystem, LPSRelContribution, IL4RelContribution, CtrlRelContribution)

TotalFluxRelContribution = na.omit(TotalFluxRelContribution)

TotalFluxRelContribution = TotalFluxRelContribution %>%
  mutate(Order = case_when(
    LPSRelContribution > IL4RelContribution & LPSRelContribution > CtrlRelContribution ~ "1",
    IL4RelContribution > LPSRelContribution & IL4RelContribution > CtrlRelContribution ~ "2",
    CtrlRelContribution > LPSRelContribution & CtrlRelContribution > IL4RelContribution ~ "3"
  ))

TotalFluxRelContribution = na.omit(TotalFluxRelContribution)

meltedDf = melt(TotalFluxRelContribution)

meltedDf = meltedDf %>%
  arrange(Order)

meltedDf$Order = as.numeric(meltedDf$Order)

meltedDf$Subsystem = str_replace_all(meltedDf$Subsystem,
                                         pattern = "Tricarboxylic acid cycle and glyoxylate/dicarboxylate metabolism",
                                         replacement = "TCA Cycle")

relContribPlot = ggplot(meltedDf, aes(x=reorder(Subsystem, Order), y=value, fill=variable))+
  geom_bar(position='stack', stat="identity", color="black", linewidth=0.7) +
  coord_flip() +
  scale_fill_manual(values=c("LPSRelContribution" = "#56B4E9",
                             "IL4RelContribution" = "#E69F00",
                             "CtrlRelContribution" = "#984EA3"),
                    labels = c('LPS', 'IL4', 'Control')) +
  theme_minimal() +
  theme(axis.text = element_text(face = "bold", size = 11),
        axis.title = element_text(face = "bold", size = 10),
        legend.title = element_text(face = "bold", size = 10),
        legend.text = element_text(face = "bold", size = 10),
        title = element_text(face = "bold")) +
  labs(y = "Relative Flux Contribution (pFBA)", x = NULL) +
  guides(fill=guide_legend(title="Models:"))


FluxDistr = rbind(LPSModel, IL4Model, ControlModel)

FluxDistr$Model = rep(c("LPS", "IL4", "Control"), times = c(nrow(LPSModel),
                                                             nrow(IL4Model),
                                                            nrow(ControlModel)))
FluxSumarize = FluxDistr %>%
  group_by(Subsystem)
FluxSumarize$Order = TotalFluxRelContribution$Order[match(FluxSumarize$Subsystem,
                                                          TotalFluxRelContribution$Subsystem)]

FluxSumarize = FluxSumarize %>%
  dplyr::filter(!is.na(Order)) %>%
  dplyr::filter(!Flux == 0) %>%
  dplyr::filter(Subsystem %in% PathwaysToKeep)

FluxSumarize$Order = as.numeric(FluxSumarize$Order)

boxPlot = ggplot(FluxSumarize, aes(x=Flux, y=reorder(Subsystem, Order), fill=Model)) +
  geom_boxplot() +
  scale_fill_manual(values=c("LPS" = "#56B4E9",
                             "IL4" = "#E69F00",
                             "Control" = "#984EA3")) +
  theme_minimal() +
  theme(axis.text = element_text(face = "bold", size = 11),
        axis.title = element_text(face = "bold", size = 10),
        legend.title = element_text(face = "bold", size = 10),
        legend.text = element_text(face = "bold", size = 10),
        axis.text.y=element_blank()) +
  labs(x = "Absolute Flux (pFBA)", y = NULL)

fig = ggarrange(relContribPlot, boxPlot,
                ncol = 2, nrow = 1,
                common.legend = T,
                legend = "bottom")

fig = annotate_figure(fig, top = text_grob(title, face = "bold"))

 return(fig)

}

FunRelativeContrib(LPSGapFilled, IL4Model = IL4GapFilled, ControlModel = CtrlGapFilled, "WT")


#### Curated FVA ####
#load all 3 datasets containing fluxes for models (IL4RXN etc., plus "reactions" file)
options(scipen = 999)
path = '/Users/rokosango/PhD/MetabModelling/MATLAB_scripts_workspaces/Immunometabolism_models/tINIT_models/NewFVAResults/'

fvaprep = function(df, rxnlist, AllKnownRxnsFile) {
  
  df$Equation = rxnlist$EQUATION[match(df$ReactionID, rxnlist$ID)]
  df$ReactionName = AllKnownRxnsFile$rxnRecon3DID[match(df$ReactionID, AllKnownRxnsFile$rxns)]
  df = df %>%
    #dplyr::filter_if(., is.numeric, any_vars((.) != 0)) %>%
    dplyr::filter(MinFlux != 0 & MaxFlux !=0) %>%
    mutate(FluxDiff = abs(MinFlux - MaxFlux)) %>%
    relocate(FluxDiff, .after = MaxFlux) %>%
    relocate(ReactionName, .after = ReactionID) %>%
    #dplyr::filter(Compartment %in% InterestingSubsystems) %>%
    arrange(FluxDiff)
  
  
}


# LPS WT
LPS = read.csv(paste0(path, 'FVA_LPS_WT.csv'))
LPS = fvaprep(LPS, LPSRXN, reactions)
# LPS Cyp27a1 Knockout 
Cyp27a1_KO_FVA_LPS = read.csv(paste0(path, 'FVA_LPS_KO_CYP.csv'))
# LPS Gapdh Knockout
Gapdh_KO_LPS = read.csv(paste0(path, 'FVA_LPS_KO_GAPDH.csv'))
# IL4 WT
IL4Cpy_WT = read.csv(paste0(path, 'FVA_IL4_WT.csv'))
# IL4 Cyp27a1 KO
IL4Cpy_KO = read.csv(paste0(path, 'FVA_IL4_KO_CYP.csv'))
# Control
Control = read.csv(paste0(path, 'FVA_CONTROL.csv'))
Control = fvaprep(Control, ControlRXN, reactions)

IL4Cpy_WT$Compartment = IL4$Compartment[match(IL4Cpy_WT$ReactionID, IL4$ReactionID)]
IL4Cpy_KO$Compartment = IL4$Compartment[match(IL4Cpy_KO$ReactionID, IL4$ReactionID)]


#### knock-out Fluxes ####
setwd('/Users/rokosango/PhD/MetabModelling/MATLAB_scripts_workspaces/Immunometabolism_models/tINIT_models/GeneKnockoutFluxResults/')

LPSGeneDelFluxSolution = read_excel('LPSGeneDelFluxSolution.xlsx',
                                   sheet = 1)
IL4GeneDelFluxSolution = read_excel('IL4GeneDelFluxSolution.xlsx',
                                    sheet = 1)
CtrlGeneDelFluxSolution = read_excel('CtrlGeneDelFluxSolution.xlsx',
                                     sheet = 1)

#for each knockout flux, make a flux table with reaction information etc.

GeneDelFlux = function(df, gene) {
  
  df = as.data.frame(df)
  rownames(df) = df$Row
  
  df$Row = NULL
  
  df = df %>%
    dplyr::select(all_of(gene))
  
  names(df)[1] = "Flux"
  df = df %>%
    mutate(ReactionID = rownames(.))
  

  df$ReactionName = reactions$rxnRecon3DID[match(df$Reaction, reactions$rxns)]
  df$Flux = as.numeric(df$Flux)
  df$Subsystem = LPSRXN$SUBSYSTEM[match(df$ReactionID, LPSRXN$ID)]
  df$Subsystem = as.factor(df$Subsystem)
  df$Equation = LPSRXN$EQUATION[match(df$ReactionID, LPSRXN$ID)]
  

  return(df)
  
}

#GeneDelFlux(LPSGeneDelFluxSolution, "Cyp27a1")

#make a list for each model and store each respect Knockout and Flux Distribution

#-----------#
LPS_Knockout_List <- list()

for (i in names(LPSGeneDelFluxSolution)[2:length(names(LPSGeneDelFluxSolution))]) {
  
  name =  paste("LPS_KO", i, sep = "_")
  assign(name, 
         LPSGeneDelFluxSolution[i])
  
  LPS_Knockout_List[[i]] = GeneDelFlux(LPSGeneDelFluxSolution, i)
}

LPS_Knockout_List[["WT"]] = LPSGapFilled

#-----------#

IL4_Knockout_List <- list()

for (i in names(IL4GeneDelFluxSolution)[2:length(names(IL4GeneDelFluxSolution))]) {
  
  name =  paste("IL4_KO", i, sep = "_")
  assign(name, 
         IL4GeneDelFluxSolution[i])
  
  IL4_Knockout_List[[i]] = GeneDelFlux(IL4GeneDelFluxSolution, i)
}

IL4_Knockout_List[["WT"]] = IL4GapFilled

#----------#

Ctrl_Knockout_List <- list()

for (i in names(CtrlGeneDelFluxSolution)[2:length(names(CtrlGeneDelFluxSolution))]) {
  
  name =  paste("Ctrl_KO", i, sep = "_")
  assign(name, 
         CtrlGeneDelFluxSolution[i])
  
  Ctrl_Knockout_List[[i]] = GeneDelFlux(CtrlGeneDelFluxSolution, i)
}

Ctrl_Knockout_List[["WT"]] = CtrlGapFilled

for (i in names(LPS_Knockout_List)) {
  print(FunRelativeContrib(LPS_Knockout_List[[i]], IL4_Knockout_List[[i]], 
                     Ctrl_Knockout_List[[i]], paste(i, "knockout", sep = "_")))
  
}

FunRelativeContrib(LPS_Knockout_List[["Cyp27a1"]], IL4_Knockout_List[["Cyp27a1"]], 
                   Ctrl_Knockout_List[["Cyp27a1"]], "Cyp27a1 knockout")


#check contributions within each model, WT vs KO:

FunRelativeContribWithAllKnockouts = function(ListModelKnockouts, tissue) {
  
  CommonSubsystems = Reduce(intersect, list(ListModelKnockouts[["WT"]]$Subsystem,
                                            ListModelKnockouts[["Cyp27a1"]]$Subsystem))
  
  SummarizedListModelKnockouts = ListModelKnockouts
  #names(SummarizedListModelKnockouts) = names(ListModelKnockouts)
  
  for (i in names(ListModelKnockouts)) {
    
    SummarizedListModelKnockouts[[i]] = ListModelKnockouts[[i]] %>%
      dplyr::group_by(Subsystem) %>%
      dplyr::summarize(AbsSum = sum(abs(Flux), na.rm = T))
    
    SummarizedListModelKnockouts[[i]] = SummarizedListModelKnockouts[[i]] %>%
      dplyr::filter(!Subsystem == "NA") %>% 
      dplyr::filter(Subsystem %in% CommonSubsystems)
  }
  
  options(scipen = 999)
  
  TotalFlux = data.frame(Subsystem = SummarizedListModelKnockouts[["WT"]]$Subsystem,
                         WT = SummarizedListModelKnockouts[["WT"]][["AbsSum"]],
                         Cyp27a1 = SummarizedListModelKnockouts[["Cyp27a1"]][["AbsSum"]])
  
  
  TotalFlux$SumAllModels = rowSums(TotalFlux[2:ncol(TotalFlux)])
  TotalFlux$SD = apply(TotalFlux[3:ncol(TotalFlux)-1], 1, sd)
  TotalFlux = TotalFlux %>% arrange(desc(SD))
  
  TotalFlux = TotalFlux[!rowSums(TotalFlux[c(2:ncol(TotalFlux))]) <= 1e-2,] #if the sum across different GEMs is < 1e-2
  TotalFluxRelContribution = TotalFlux %>%
    mutate(WT_RelContribution = (WT * 1) / SumAllModels) %>%
    mutate(Cyp27a1_RelContribution = (Cyp27a1 * 1) / SumAllModels) %>%
    dplyr::select(Subsystem, WT_RelContribution, Cyp27a1_RelContribution)
  
  TotalFluxRelContribution = na.omit(TotalFluxRelContribution)
  
  meltedDf = melt(TotalFluxRelContribution)
  
  meltedDf$Subsystem = as.factor(meltedDf$Subsystem)
  
  # for clarity sake; otherwise too many pathways are plotted:
  PathwaysToRemove = c('Cysteine and methionine metabolism',
                       'Cholesterol metabolism',
                       'Carnitine shuttle (mitochondrial)',
                       'Carnitine shuttle (endoplasmic reticular)',
                       'Carnitine shuttle (cytosolic)',
                       'Butanoate metabolism',
                       'Biopterin metabolism',
                       'Bile acid recycling',
                       "Pool reactions",
                       "Transport reactions",
                       "Miscellaneous",
                       "Exchange/demand reactions",
                       "Inositol phosphate metabolism",
                       "Peptide metabolism",
                       "Prostaglandin biosynthesis",
                       "Valine, leucine, and isoleucine metabolism",
                       "Sulfur metabolism",
                       "Steroid metabolism",
                       "Starch and sucrose metabolism",
                       "Sphingolipid metabolism",
                       "ROS detoxification",
                       "Retinol metabolism",
                       "Protein modification",
                       "Porphyrin metabolism",
                       "Heparan sulfate degradation",
                       "Protein degradation",
                       "Arachidonic acid metabolism",
                       "Carnitine shuttle (peroxisomal)",
                       "Glycerolipid metabolism",
                       "Glycerophospholipid metabolism",
                       "Glutathione metabolism",
                       "Galactose metabolism",
                       "Aminoacyl-tRNA biosynthesis",
                       "Alanine, aspartate and glutamate metabolism",
                       "Amino sugar and nucleotide sugar metabolism",
                       "Artificial reactions",
                       "Isolated",
                       "Glycine, serine and threonine metabolism",
                       "Leukotriene metabolism",
                       "Eicosanoid metabolism",
                       "Pentose phosphate pathway", 
                       "Nicotinate and nicotinamide metabolism",
                       "Pentose and glucuronate interconversions",
                       "Omega-6 fatty acid metabolism",
                       "N-glycan metabolism",
                       "Phenylalanine metabolism",
                       "Lipoic acid metabolism",
                       "Acylglycerides metabolism",
                       "Drug metabolism",
                       "Lysine metabolism",
                       "Riboflavin metabolism",
                       "Terpenoid backbone biosynthesis",
                       "Cholesterol biosynthesis 1 (Bloch pathway)",
                       "Cholesterol biosynthesis 2",
                       "Ubiquinone synthesis",
                       "Heme synthesis",
                       "Phenylalanine, tyrosine and tryptophan biosynthesis",
                       "Nucleotide metabolism",
                       "Pantothenate and CoA biosynthesis",
                       "Arginine and proline metabolism",
                       "Folate metabolism",
                       "Protein assembly",
                       "Fructose and mannose metabolism",
                       "Heme synthesis")
                      
  
  meltedDf = meltedDf %>%
    dplyr::filter(!Subsystem %in% PathwaysToRemove)
  
  if (tissue == "LPS") {
    
    WT_col = "#56B4E9"
    col = "#084594"
    
  } else {
    
    WT_col = "#E69F00"
    col = "#FC4E2A"
    
    meltedDf = meltedDf %>%
      dplyr::filter(!Subsystem %in% 'Bile acid biosynthesis')
    
  }
  
  meltedDf = meltedDf[!grepl("Vitamin", meltedDf$Subsystem),] #for Poster plot
  meltedDf = meltedDf[!grepl("Fatty", meltedDf$Subsystem),] #for Poster plot
  meltedDf = meltedDf[!grepl("Glycos", meltedDf$Subsystem),] #for Poster plot
  
  
  meltedDf$Subsystem = str_replace_all(meltedDf$Subsystem, 
                                       "Tricarboxylic acid cycle and glyoxylate/dicarboxylate metabolism",
                                       "TCA cycle")
  
  meltedDf$Subsystem <- factor(meltedDf$Subsystem, levels=rev(c("Bile acid biosynthesis",
                                                            setdiff(meltedDf$Subsystem, "Bile acid biosynthesis"))))
                                                            
                                                            
  
  relContribPlot = ggplot(meltedDf, aes(x=Subsystem, y=value, fill=variable))+
    geom_bar(position='stack', stat="identity", color="black", linewidth=0.5, width=0.7) +
    coord_flip() +
    scale_fill_manual(values=c("WT_RelContribution" = WT_col,
                               "Cyp27a1_RelContribution" = col),
    labels = c('WT', 'Cyp27a1 KO')) +
    theme_minimal() +
    theme(axis.text.x = element_text(face = "bold", size = 14),
          axis.text.y = element_text(face = "bold", size = 14),
          axis.title = element_text(face = "bold", size = 14),
          legend.title = element_text(face = "bold", size = 10),
          legend.text = element_text(face = "bold", size = 10),
          legend.position = "right",
          legend.direction = "vertical", 
          title = element_text(face = "bold")) +
    labs(y = "Relative Flux Contribution (FBA)", x = NULL, title = NULL) +
    guides(fill=guide_legend(title="Models:"))
  
 FluxDistrBoxPlot = rbind(ListModelKnockouts[["WT"]], ListModelKnockouts[["Cyp27a1"]])
 
 FluxDistrBoxPlot$Model = rep(c("WT", "Cyp27a1 KO"), times = c(nrow(ListModelKnockouts[["WT"]]),
                                                             nrow(ListModelKnockouts[["Cyp27a1"]])))
 
FluxSumarize = FluxDistrBoxPlot %>%
  group_by(Subsystem)


FluxSumarize = FluxSumarize %>%
  dplyr::filter(!Flux == 0) %>%
  dplyr::filter(!Subsystem %in% PathwaysToRemove) %>%
  dplyr::filter(Subsystem %in% meltedDf$Subsystem | Subsystem %in% "Tricarboxylic acid cycle and glyoxylate/dicarboxylate metabolism")

  if (tissue == "IL4") {
  
    FluxSumarize = FluxSumarize %>%
    dplyr::filter(!Subsystem %in% 'Bile acid biosynthesis')
  
  }

FluxSumarize = na.omit(FluxSumarize)



FluxSumarize$Subsystem = str_replace_all(FluxSumarize$Subsystem,
                                     "Tricarboxylic acid cycle and glyoxylate/dicarboxylate metabolism",
                                     "TCA cycle")

FluxSumarize$Subsystem <- factor(FluxSumarize$Subsystem, levels=rev(c("Bile acid biosynthesis",
                                                            setdiff(meltedDf$Subsystem, "Bile acid biosynthesis"))))


boxPlot = ggplot(FluxSumarize, aes(x=Flux, y=Subsystem, fill=Model)) +
  geom_boxplot(width=0.7) +
  scale_fill_manual(values=c("WT" = WT_col,
                             "Cyp27a1 KO" = col),
                    labels = c('WT', 'Cyp27a1 KO')) +
  theme_minimal() +
  theme(axis.text = element_text(face = "bold", size = 14),
        axis.title = element_text(face = "bold", size = 14),
        legend.title = element_text(face = "bold", size = 10),
        legend.text = element_text(face = "bold", size = 10),
        legend.direction = "vertical",
        axis.text.y=element_blank()) +
  labs(x = "Absolute Flux (FBA)", y = NULL)


fig = ggarrange(relContribPlot, boxPlot,
                ncol = 2, nrow = 1,
                common.legend = T,
                legend = "bottom")



  res = list(TotalFlux, fig)
  names(res) = c("TableSummedFluxes", "CombinedPlot")

  return(res)

}


LPS_KO_Results = FunRelativeContribWithAllKnockouts(LPS_Knockout_List, 'LPS')
IL4_KO_Results = FunRelativeContribWithAllKnockouts(IL4_Knockout_List, 'IL4')
Ctrl_KO_Results = FunRelativeContribWithAllKnockouts(Ctrl_Knockout_List, 'Control')

#### LPS Knockout Fluxes: ####
# choose: 
# Bile acid biosynthesis, 
# Tricarboxylic acid cycle and glyoxylate/dicarboxylate metabolism
# Pyrimidine metabolism
# Heme synthesis

InterestingSubs = c("Bile acid biosynthesis", 
                    "Tricarboxylic acid cycle and glyoxylate/dicarboxylate metabolism",
                    "Pyrimidine metabolism", "Heme synthesis")

ProcessPlotFun = function(df1, df2, InterestingSubs, gene) {
  
   df1 = LPSGapFilled %>%
    dplyr::filter(Subsystem %in% InterestingSubs) %>%
    dplyr::group_by(Subsystem) %>%
    dplyr::summarize(AbsMean = mean(abs(Flux), na.rm = T),
              SDMean = sd(Flux, na.rm = T))
   
  df2 = LPS_Knockout_List[[gene]] %>%
    dplyr::filter(Subsystem %in% InterestingSubs) %>%
    dplyr::group_by(Subsystem) %>%
    dplyr::summarize(AbsMean = mean(abs(Flux), na.rm = T),
              SDMean = sd(Flux, na.rm = T))
  
  combined = rbind(df1, df2)
  
  combined$Subsystem = str_replace_all(combined$Subsystem, 
                                       "Tricarboxylic acid cycle and glyoxylate/dicarboxylate metabolism",
                                       "TCA cycle")
  
  combined$Model = rep(c("WT", "Cyp27a1"), each = 4)
  combined$Model = as.factor(combined$Model)
  combined$Subsystem = str_wrap(combined$Subsystem, width = 10)
  
  ggplot(combined, aes(x = Subsystem, y = AbsMean,
                 fill=Model)) +
    geom_bar(stat = "identity", color = "black", 
             position = position_dodge()) +
    scale_fill_manual(values=c("WT" = "#1B9E77", #D95F02
                               'Cyp27a1' = "#D95F02"), # 1B9E77
                      labels = c('Cyp27a1_KO', 'WT')) +
    geom_errorbar(aes(ymin=(AbsMean - SDMean), 
                      ymax=(AbsMean + SDMean)), 
                      width=.2,
                      position=position_dodge(.9)) +
    #coord_flip() +
    labs(x = NULL, y = "Flux") +
    theme_bw() +
    theme(axis.text = element_text(face = "bold", size = 15),
          axis.title = element_text(face = "bold", size = 15),
          legend.title = element_text(face = "bold", size = 10),
          legend.text = element_text(face = "bold", size = 10),
          legend.position = c(0.15, 0.9),
          legend.background = element_rect(fill = "white", color = "black"))
    
  
}

ProcessPlotFun(LPSGapFilled, LPS_Knockout_List[["Cyp27a1"]],
               InterestingSubs, 'Cyp27a1')



#### focus on specific subsystems and reactions within subsystems of interest:


grabFBAFVA = function(wtModel, koModel, fvaWtModel, 
                      fvaKoModel, reaction, WhichGene, Condition=NULL
                      ) {
  
  if (Condition == "LPS") {
    WT_col = "#56B4E9"
    col = "#084594"
  } else {
    WT_col = "#E69F00"
    col = "#FC4E2A"
    
  }
  
  WtFindRxns = intersect(wtModel$ReactionID, fvaWtModel$ReactionID)
  KoFindRxns = intersect(koModel$ReactionID, fvaKoModel$ReactionID)
   
   wtModel = wtModel %>%
       dplyr::filter(ReactionID %in% WtFindRxns)
   fvaWtModel = fvaWtModel %>%
     dplyr::filter(ReactionID %in% WtFindRxns)
 

  dfWt = cbind(wtModel, fvaWtModel)
  dfWt = dfWt[, !duplicated(colnames(dfWt))]
  dfWt = dfWt %>% dplyr::filter(ReactionID == reaction)
  

  koModel = koModel %>%
    dplyr::filter(ReactionID %in% KoFindRxns)
  
  fvaKoModel = fvaKoModel %>%
    dplyr::filter(ReactionID %in% KoFindRxns)

  # if (WhichGene == "Gapdh_KO") {
  # 
  #   fvaKoModel = fvaKoModel %>%
  #     dplyr::rename(minFlux = MinFlux) %>%
  #     dplyr::rename(maxFlux = MaxFlux)
  # }


  dfKo = cbind(koModel, fvaKoModel)
  dfKo = dfKo[, !duplicated(colnames(dfKo))]
  dfKo = dfKo %>% dplyr::filter(ReactionID == reaction)
  
  dfKo = dfKo[, names(dfWt)]
  
  
  dfCombined = rbind(dfWt, dfKo)
  dfCombined$Model = c("WT", WhichGene)
  dfCombined$Model = as.factor(dfCombined$Model)

  if (WhichGene == "Cyp27a1_KO") {

    ggplot(dfCombined, aes(x=Model, y=Flux, fill=Model, color = Model)) +
      geom_bar(stat="identity", position=position_dodge(), width = 0.2) +
      geom_errorbar(aes(ymin=minFlux, ymax=maxFlux), width=.5, linewidth = 2,
                    position=position_dodge(.9)) +
      scale_color_manual(values=c("WT" = WT_col,
                                  "Cyp27a1_KO" = col)) +
      scale_fill_manual(values=c("WT" = WT_col,
                                "Cyp27a1_KO" = col)) +
      labs(x = NULL, y = NULL) + my.theme + theme(axis.text.x = element_blank(),
                                                  axis.ticks.x = element_blank(),
                                                  legend.position = "none",
                                                  axis.text.y = element_text(size = 20)) +
      scale_y_continuous(position = "right")  +
      geom_hline(yintercept=0, linetype="dashed", color = "black")
  }
  else {

    ggplot(dfCombined, aes(x=Model, y=Flux, fill=Model, color = Model)) +
      geom_bar(stat="identity", position=position_dodge(), width = 0.2) +
      geom_errorbar(aes(ymin=MinFlux, ymax=MaxFlux), width=.5, linewidth = 2,
                    position=position_dodge(.9)) +
      scale_color_manual(values=c("WT" = WT_col,
                                 "Gapdh_KO" = col)) +
     scale_fill_manual(values=c("WT" = WT_col,
                                "Gapdh_KO" = col)) +
    labs(x = NULL, y = NULL) + my.theme + theme(axis.text.x = element_blank(),
                                                axis.ticks.x = element_blank(),
                                                legend.position = "none",
                                                axis.text.y = element_text(size = 20)) +
    scale_y_continuous(position = "right")  +
    geom_hline(yintercept=0, linetype="dashed", color = "black")
  }
}
  
grabFBAFVA(wtModel = IL4GapFilled, koModel = IL4_Knockout_List[["Cyp27a1"]],
           fvaWtModel = IL4Cpy_WT, fvaKoModel = IL4Cpy_KO, reaction = "MAR04002",
           WhichGene = "Cyp27a1_KO", Condition = "IL4")


grabFBAFVA(wtModel = LPSGapFilled, koModel = LPS_Knockout_List[["Gapdh"]],
           fvaWtModel = LPS, fvaKoModel = Gapdh_KO_LPS, reaction = "MAR04358", 
           WhichGene = "Gapdh_KO", Condition = "LPS")


#r0191 = PFK
#### Random sampling by ACHR ####

LPSSampledModelRxnID = read.table('~/PhD/MetabModelling/MATLAB_scripts_workspaces/SamplingResults/LPSModelSamplingRxns.txt', 
                                  col.names = "Reaction")

IL4SampledModelRxnID = read.table('~/PhD/MetabModelling/MATLAB_scripts_workspaces/SamplingResults/IL4ModelSamplingRxns.txt', 
                                  col.names = "Reaction")

KO_LPSSampledModelRxnID = read.table('~/PhD/MetabModelling/MATLAB_scripts_workspaces/SamplingResults/KO_LPSModelSamplingRxns.txt',
                                  col.names = "Reaction")

KO_IL4SampledModelRxnID = read.table('~/PhD/MetabModelling/MATLAB_scripts_workspaces/SamplingResults/KO_IL4ModelSamplingRxns.txt',
                                     col.names = "Reaction")



#CommonRxns = intersect(IL4SampledModelRxnID$Reaction,
#                       LPSSampledModelRxnID$Reaction)

CommonRxns = intersect(KO_IL4SampledModelRxnID$Reaction,
                       IL4SampledModelRxnID$Reaction)

filter = function(df) {
  df = df %>%
    dplyr::filter(rownames(df) %in% CommonRxns) %>%
    arrange(rownames(.))
  return(df)
}
                       


#LPS

path = '~/PhD/MetabModelling/MATLAB_scripts_workspaces/SamplingResults/LPS/'


# load up the data: reaction IDs & flux matrices
#find reactions that are shared across LPS & IL4
#arrange rows in both conditions (10 frames per model)
# for each output file (10), find mean of each sampled flux
# explore densities of those reactions (side by side densities per model)
#focus on reactions in glycolysis, tca, oxphos, ppp and arginine pathways
# use Kruskal–Wallis test to find significant reactions
LPSdata_files <- list.files(paste0(path))  # Identify file names
LPSdata_files


for(i in 1:length(LPSdata_files)) {                              # Head of for-loop
  assign(paste0("LPSSamplingFile_", i),                                   # Read and store data frames
         read.csv(paste0(path,
                          LPSdata_files[i]),
                   header = F, 
                   col.names = paste0(rep("sample", times = 1000), rep(1:1000)),
                   row.names =  LPSSampledModelRxnID$Reaction))
}


LPS_list = list(LPSSamplingFile_1, LPSSamplingFile_2, LPSSamplingFile_3, LPSSamplingFile_4,
                LPSSamplingFile_5, LPSSamplingFile_6, LPSSamplingFile_7, LPSSamplingFile_8,
                LPSSamplingFile_9, LPSSamplingFile_10)

SampleFileNames = paste0("LPSSamplingFile_10", 1:10)

for (i in 1:length(LPS_list)) {
  
  LPS_list[[i]] = filter(LPS_list[[i]])
  names(LPS_list)[i] = SampleFileNames[i]
  
}

LPSSamplingMeanDf = Reduce(`+`, LPS_list)/length(LPS_list)

#IL4


path = '~/PhD/MetabModelling/MATLAB_scripts_workspaces/SamplingResults/IL4/'

IL4data_files <- list.files(paste0(path))  # Identify file names
IL4data_files

for(i in 1:length(IL4data_files)) {                              # Head of for-loop
  assign(paste0("IL4SamplingFile_", i),                                   # Read and store data frames
         read.csv(paste0(path,
                         IL4data_files[i]),
                  header = F, 
                  col.names = paste0(rep("sample", times = 1000), rep(1:1000)),
                  row.names = IL4SampledModelRxnID$Reaction))
}


IL4_list = list(IL4SamplingFile_1, IL4SamplingFile_2, IL4SamplingFile_3, IL4SamplingFile_4,
                IL4SamplingFile_5, IL4SamplingFile_6, IL4SamplingFile_7, IL4SamplingFile_8,
                IL4SamplingFile_9, IL4SamplingFile_10)

SampleFileNames = paste0("IL4SamplingFile_", 1:10)

for (i in 1:length(IL4_list)) {
  
  IL4_list[[i]] = filter(IL4_list[[i]])
  names(IL4_list)[i] = SampleFileNames[i]
  
}

IL4SamplingMeanDf = Reduce(`+`, IL4_list)/length(IL4_list)

# LPS Knockout

path = '~/PhD/MetabModelling/MATLAB_scripts_workspaces/SamplingResults/KO_LPS/'

KO_LPSdata_files <- list.files(paste0(path))  # Identify file names
KO_LPSdata_files

for(i in 1:length(KO_LPSdata_files)) {                              # Head of for-loop
  assign(paste0("KO_LPSSamplingFile_", i),                                   # Read and store data frames
         read.csv(paste0(path,
                         KO_LPSdata_files[i]),
                  header = F, 
                  col.names = paste0(rep("sample", times = 1000), rep(1:1000)),
                  row.names = KO_LPSSampledModelRxnID$Reaction))
}

KO_LPS_list = list(KO_LPSSamplingFile_1, KO_LPSSamplingFile_2, KO_LPSSamplingFile_3, KO_LPSSamplingFile_4,
                   KO_LPSSamplingFile_5, KO_LPSSamplingFile_6, KO_LPSSamplingFile_7, KO_LPSSamplingFile_8,
                   KO_LPSSamplingFile_9, KO_LPSSamplingFile_10)

SampleFileNames = paste0("KO_LPSSamplingFile_", 1:10)

for (i in 1:length(KO_LPS_list)) {
  
  KO_LPS_list[[i]] = filter(KO_LPS_list[[i]])
  names(KO_LPS_list)[i] = SampleFileNames[i]
  
}

#do element-wise addition and divide by the number of the list elements. i.e., get the mean of
#each sampled flux for each cell across 10 list elements
KO_LPSSamplingMeanDf = Reduce(`+`, KO_LPS_list)/length(KO_LPS_list)

# IL4 Knockout

path = '~/PhD/MetabModelling/MATLAB_scripts_workspaces/SamplingResults/KO_IL4/'

KO_IL4data_files <- list.files(paste0(path))  # Identify file names
KO_IL4data_files

for(i in 1:length(KO_IL4data_files)) {                              # Head of for-loop
  assign(paste0("KO_IL4SamplingFile_", i),                                   # Read and store data frames
         read.csv(paste0(path,
                         KO_IL4data_files[i]),
                  header = F, 
                  col.names = paste0(rep("sample", times = 1000), rep(1:1000)),
                  row.names = KO_IL4SampledModelRxnID$Reaction))
}

KO_IL4_list = list(KO_IL4SamplingFile_1, KO_IL4SamplingFile_2, KO_IL4SamplingFile_3, KO_IL4SamplingFile_4,
                   KO_IL4SamplingFile_5, KO_IL4SamplingFile_6, KO_IL4SamplingFile_7, KO_IL4SamplingFile_8,
                   KO_IL4SamplingFile_9, KO_IL4SamplingFile_10)

SampleFileNames = paste0("KO_IL4SamplingFile_", 1:10)

for (i in 1:length(KO_IL4_list)) {
  
  KO_IL4_list[[i]] = filter(KO_IL4_list[[i]])
  names(KO_IL4_list)[i] = SampleFileNames[i]
  
}

KO_IL4SamplingMeanDf = Reduce(`+`, KO_IL4_list)/length(KO_IL4_list)


#LPS GAPDH knockout

KO_Gapdh_LPS_SamplingDf = read.csv('~/PhD/MetabModelling/MATLAB_scripts_workspaces/SamplingResults/Sampling_LPS_GAPDH_KO.csv')


#concat different sampling models to compare significance
concat_sampling_dfs = function(sampling_mean_df1, 
                               sampling_mean_df2,
                               sampling_df1_names,
                               sampling_df2_names,
                               model,
                               InterestingSubsystems) {
  
  options(scipen = 4)
  
  df = cbind(sampling_mean_df1,
             sampling_mean_df2)
  
  pvals <- apply(df, 1, function(x) {
    wilcox.test(x[1:1000], x[1001:2000])$p.value
  })
  

  adjpvals = p.adjust(pvals, method = "bonferroni")
  adjpvals = adjpvals[adjpvals < 1e-5]
  
  df = df[names(adjpvals),]
  names(df) = c(paste0(sampling_df1_names, 1:1000), paste0(sampling_df2_names, 1:1000))
  df$Subsystem = model$Subsystem[match(rownames(df),
                                              model$ReactionID)]
  df$ReactionName= model$ReactionName[match(rownames(df),
                                                   model$ReactionID)]
  
  df = df %>% 
    relocate(Subsystem, .before = IL4_Sample_1) %>%
    relocate(ReactionName, .after = Subsystem) %>%
    arrange(Subsystem) %>%
    mutate(AdjPvals = adjpvals[match(names(adjpvals),
                                     rownames(df))]) %>%
    relocate(AdjPvals, .after = ReactionName) %>%
    mutate(ReactionID = rownames(.)) %>%
    relocate(ReactionID, .after = AdjPvals) %>%
    dplyr::filter(Subsystem %in% InterestingSubsystems)
  
}

InterestingSubsystems = c("Bile acid biosynthesis",
                          "Heme synthesis",
                          "Pyrimidine metabolism",
                          "Tricarboxylic acid cycle and glyoxylate/dicarboxylate metabolism",
                          "Glycolysis / Gluconeogenesis")


combined_mean_dfs = concat_sampling_dfs(LPSSamplingMeanDf, KO_LPSSamplingMeanDf,
                         "LPS_Sample_", "KO_LPS_Sample_",
                         LPSGapFilled, InterestingSubsystems)


combined_mean_dfs_IL4 = concat_sampling_dfs(IL4SamplingMeanDf, KO_IL4SamplingMeanDf,
                                      "IL4_Sample_", "KO_IL4_Sample_",
                                      IL4GapFilled, c("Glycolysis / Gluconeogenesis",
                                                      "Tricarboxylic acid cycle and glyoxylate/dicarboxylate metabolism",
                                                      "Purine metabolism",
                                                      "Pyrimidine metabolism"))



#function for generating density plots from significant
#reactions sampled from 5 different subsystems

densityplotsfun = function(sampling_mean_df1, 
                           sampling_mean_df2,
                           reactionID, Model1Factor,
                           Model2Factor, Condition=NULL,
                           title) {
  
  if (Model2Factor == "Cyp27a1_KO") {
    col = "#D95F02"
  } else {
    col = "#A16864"
  }
  
  if (Condition == "LPS") {
    WT_col = "#56B4E9"
    col = "#084594"
  } else {
    WT_col = "#E69F00"
    col = "#FC4E2A"
    
  }
  
  combined = melt(rbind(sampling_mean_df1[reactionID, ],
                        sampling_mean_df2[reactionID, ]))

  combined$Model = rep(c(Model1Factor, Model2Factor))

  m = combined$value[seq(1, length(combined$value), 2)] #every step-wise first sampling model value
  n = combined$value[seq(2, length(combined$value), 2)] #every step-wise second sampling model value
  
  
  
  if (Model2Factor == "Cyp27a1_KO") {
    ggplot(combined, aes(x=value, fill=Model)) + 
      geom_density(color="black", alpha=0.9) + 
      scale_x_continuous(trans = pseudolog10_trans) +
      scale_fill_manual(values=c("WT" = WT_col,
                                 "Cyp27a1_KO" = col)) +
      theme_bw() +
      theme(axis.text = element_text(face = "bold", size = 12),
            axis.title = element_text(face = "bold", size = 16),
            legend.title = element_text(face = "bold", size = 10),
            legend.text = element_text(face = "bold", size = 10),
            legend.position = "none") +
      labs(x = NULL, y = NULL)
    
  } else {
    
    ggplot(combined, aes(x=value, fill=Model)) + 
      geom_density(color="black", alpha=0.9) + 
      scale_x_continuous(trans = pseudolog10_trans) +
      scale_fill_manual(values=c("WT" = WT_col,
                                 "Gapdh_KO" = col)) +
      theme_bw() +
      theme(axis.text = element_text(face = "bold", size = 12),
            axis.title = element_text(face = "bold", size = 16),
            legend.title = element_text(face = "bold", size = 10),
            legend.text = element_text(face = "bold", size = 10),
            legend.position = "none") +
      labs(x = NULL, y = NULL)
  }
  
}

test = IL4_Knockout_List[["Cyp27a1"]] %>% dplyr::filter(Subsystem %in% "Pyrimidine metabolism")


densityplotsfun(LPSSamplingMeanDf, KO_LPSSamplingMeanDf,
                "MAR04137", "WT", "Cyp27a1_KO", 
                Condition = "LPS", "Pyruvate dehydrogenase")


densityplotsfun(IL4SamplingMeanDf, KO_IL4SamplingMeanDf,
                "MAR04032", "WT", "Cyp27a1_KO", 
                Condition = "IL4", "MAR04574")


grabFBAFVA(wtModel = IL4GapFilled, koModel = IL4_Knockout_List[["Cyp27a1"]],
           fvaWtModel = IL4Cpy_WT, fvaKoModel = IL4Cpy_KO, 
           reaction = "MAR04032", WhichGene = "Cyp27a1_KO", Condition = "IL4")




# GAPDH KO

tKO_Gapdh_LPS_SamplingDf = t(KO_Gapdh_LPS_SamplingDf)
colnames(tKO_Gapdh_LPS_SamplingDf) = paste0(rep("sample", times = 1000), rep(1:1000))

checkGapdhDensities = function(wtModelSampling,
                               koModelSampling,
                               reaction) {

  test_lps = wtModelSampling[reaction,] %>%
    melt() %>%
    dplyr::select(value) %>%
    mutate(Model = rep("WT"))
  
  test_ko = koModelSampling[reaction,] %>%
    as.data.frame() %>%
    dplyr::rename(value = ".") %>%
    mutate(Model = rep("Gapdh_KO"))
  
  test_combined = rbind(test_lps, test_ko)
  
  wt_mean_value = wtModelSampling[reaction,] %>%
    melt() %>%
    dplyr::select(value) %>%
    mutate(Model = rep("WT")) %>%
    summarise(mean = mean(value)) %>%
    as.numeric()
  
  ko_mean_value = koModelSampling[reaction,] %>%
    as.data.frame() %>%
    dplyr::rename(value = ".") %>%
    summarise(mean = mean(value)) %>%
    as.numeric()
  
  #WT_col = "#56B4E9"
  #col = "#C6DBEF"
  
  ggplot(test_combined, aes(x = value, fill = Model)) +
    scale_fill_manual(values=c("WT" = "#56B4E9",
                               "Gapdh_KO" = "#084594")) + 
    geom_density(color="black", alpha=0.9) + 
    geom_vline(xintercept = wt_mean_value, linetype="dotted", 
               color = "#56B4E9", size=3) +
    geom_vline(xintercept = ko_mean_value, linetype="dotted", 
               color = "#084594", size=3) +
    labs(x = NULL, y = NULL, title = NULL) + my.theme +
    theme(#legend.position = "none",
          #legend.direction = "horizontal",
          axis.text = element_text(size = 20))
}

checkGapdhDensities(LPSSamplingMeanDf,
                    tKO_Gapdh_LPS_SamplingDf,
                    'MAR04394')

grabFBAFVA(wtModel = LPSGapFilled, koModel = LPS_Knockout_List[["Gapdh"]],
           fvaWtModel = LPS, fvaKoModel = Gapdh_KO_LPS, reaction = "MAR04358", WhichGene = "Gapdh_KO", Condition = "LPS")


#use this densityplotsfun to get sampling results from rxns that are changing
#between conditions in SpecFluxesFun

unique_rxns = factor(unique_rxns, levels = unique(rev(sort(unique_rxns))))

unique_rxns = TCAMet[[2]]$ReactionID


for (i in 1:length(unique_rxns)){
          
           plotlist[[i]] = print(densityplotsfun(LPSSamplingMeanDf,
                                        KO_LPSSamplingMeanDf,
                                        unique_rxns[i], 
                                        'WT', 'Cyp27a1_KO', 
                                        LPSGapFilled, LPS_Knockout_List[["Cyp27a1"]],
                                        unique_rxns[i]))
}


arrangedplots = function(sampling_mean_df1, 
          sampling_mean_df2,
          subsystem_name,
          subsystem_list, #subsystem = "HemeSynth", "TCAMet", "PyrMet" from SpecFluxesFun. No "BileAcid"
          model1Factor,
          model2Factor,
          FBAModel1, 
          FBAModel2,
          nrow) { 
  
          plotlist = list()
  
          #unique_rxns = unique(subsystem_list[[2]]$ReactionID)
          unique_rxns = unique(sort(subsystem_list[[2]]$ReactionID))
          
    for (i in 1:length(unique_rxns)){
    
    plotlist[[i]] = print(densityplotsfun(sampling_mean_df1,
                                          sampling_mean_df2,
                                          unique_rxns[i], 
                                          model1Factor, model2Factor, 
                                          FBAModel1, FBAModel2,
                                          unique_rxns[i]))
  }
  
  fig = ggarrange(plotlist=plotlist,
                  labels = letters[1:length(plotlist)],
                  ncol = 2, nrow = nrow,
                  common.legend = T)
  
  
  fig
  #require(grid)
  anotfig = annotate_figure(fig, 
                            left = textGrob("Density", rot = 90, vjust = 1, gp = gpar(cex = 0.8)),
                            bottom = textGrob("Sampled Flux", gp = gpar(cex = 0.8)),
                            top = textGrob(subsystem_name))
  
  return(anotfig)
  
}


arrangedplots(LPSSamplingMeanDf,KO_LPSSamplingMeanDf,"TCA Cycle",
              TCAMet, "WT", "Cyp27a1_KO", LPSGapFilled,
              LPS_Knockout_List[["Cyp27a1"]], 4)

arrangedplots(LPSSamplingMeanDf,KO_LPSSamplingMeanDf,"Pyrimidine Metabolism",
              PyrMet, "WT", "Cyp27a1_KO", LPSGapFilled,
              LPS_Knockout_List[["Cyp27a1"]], 7)

arrangedplots(LPSSamplingMeanDf,KO_LPSSamplingMeanDf, "Heme synthesis",
              HemeSynth, "LPS", "Cyp27a1_KO", LPSGapFilled,
              LPS_Knockout_List[["Cyp27a1"]], 1)



#from concatenated LPS + IL4 sampling dataset find reactions 
#that are also found in FVA dataset

FVAfiltSamplingDf = df %>%
  dplyr::filter(rownames(.) %in% FVACombined$ReactionID)

#PYK - MAR04358 - pyruvate kinase reaction:
#it is not in FVA dataset as it is only captured in LPS (not in Control or IL4) but with wide flux range {1, 1000}
# hence it is not in the FBA-FVA figure. but still it would be nice to show last step in
# glycolysis where ATP is being made:

FVAfiltSamplingDf['MAR04358', ] = df['MAR04358', ] #add reaction PYK
FVAfiltSamplingDf['MAR04379', ] = df['MAR04379', ] #add reaction PFK
FVAfiltSamplingDf['MAR04368', ] = df['MAR04368', ] #add reaction PGK
FVAfiltSamplingDf['MAR04394', ] = df['MAR04394', ] #add reaction HEX1

FVAfiltSamplingDf = FVAfiltSamplingDf[!FVAfiltSamplingDf$ReactionName %in% "r0173",]  #remove reaction r0173 (peroxisomal lactate dehydrogenase)


df[!(row.names(df) %in% c("1","2")),]

arrangedplots = function(x, subsystem, nrow) {
 
   plotlist = list()
   
   filtered = x %>%
     dplyr::filter(Subsystem == subsystem)
  
  for (i in 1:dim(filtered)[1]) {
    
    plotlist[[i]] = print(densityplotsfun(rownames(filtered)[i], filtered[i, "ReactionName"]))
    
  }
   
   fig = ggarrange(plotlist=plotlist,
                   labels = letters[1:length(plotlist)],
                   ncol = 2, nrow = nrow, common.legend = T)
   
   fig
   #require(grid)
   anotfig = annotate_figure(fig, 
                left = textGrob("Density", rot = 90, vjust = 1, gp = gpar(cex = 0.8)),
                bottom = textGrob("Sampled Flux", gp = gpar(cex = 0.8)),
                 top = textGrob(subsystem))
   
   return(anotfig)
}

FVAfiltSamplingDf %>% 
  group_by(Subsystem) %>%
  summarise(no_rows = length(Subsystem))


#Subsystem                                                        no_rows
#<fct>                                                                  <int>
#1 Arginine and proline metabolism                                        6
#2 Glycolysis / Gluconeogenesis                                           7
#3 Oxidative phosphorylation                                              2
#4 Pentose phosphate pathway                                              8
#5 Tricarboxylic acid cycle and glyoxylate/dicarboxylate metabolism       2

arrangedplots(FVAfiltSamplingDf, "Arginine and proline metabolism", 3)
arrangedplots(FVAfiltSamplingDf, "Glycolysis / Gluconeogenesis", 5)
arrangedplots(FVAfiltSamplingDf, "Pentose phosphate pathway", 4)
arrangedplots(FVAfiltSamplingDf, 
              "Tricarboxylic acid cycle and glyoxylate/dicarboxylate metabolism", 2)
arrangedplots(FVAfiltSamplingDf, "Oxidative phosphorylation", 2)



#this is just to compare means and variances LPS vs IL4
LPSSamplingMeanDf$MeanLPS = apply(LPSSamplingMeanDf, 1, mean)
IL4SamplingMeanDf$MeanIL4 = apply(IL4SamplingMeanDf, 1, mean)
LPSSamplingMeanDf$VarLPS = apply(LPSSamplingMeanDf, 1, var)
IL4SamplingMeanDf$VarIL4 = apply(LPSSamplingMeanDf, 1, var)

SamplingDfCombined = data.frame(cbind(
                           rownames(LPSSamplingMeanDf),
                           LPSSamplingMeanDf$MeanLPS,
                           IL4SamplingMeanDf$MeanIL4,
                           LPSSamplingMeanDf$VarLPS,
                           IL4SamplingMeanDf$VarIL4))

names(SamplingDfCombined) = c(
  "ReactionID", "MeanLPS", "MeanIL4", "VarLPS", "VarIL4"
)


rownames(SamplingDfCombined) = rownames(LPSSamplingMeanDf)


######### Cyp27a1 macs validation with metabolomics #########

options(ggrepel.max.overlaps = Inf)

path = '/Users/rokosango/PhD/Metabolomics/data'

RP = read_excel(paste0(path, "/Cyp27a1_KO_Metabolomics_first_experiment.xlsx"),
                sheet = 1)
HILIC = read_excel(paste0(path, "/Cyp27a1_KO_Metabolomics_first_experiment.xlsx"),
                sheet = 2)
RepeatHILIC = read_excel(paste0(path, "/Cyp27a1_KO_Metabolomics_repeat_experiment.xlsx"),
                         sheet = 1)
RepeatRP = read_excel(paste0(path, "/Cyp27a1_KO_Metabolomics_repeat_experiment.xlsx"),
                         sheet = 2)
GAPDH = read.csv(paste0(path, "/GAPDH_experiment.csv"))

AnalyzeFun = function(df, title) {
  
  dfNum = df %>%
    dplyr::mutate(Group = paste(df$Group, df$Phenotype, sep = "_")) %>%
    dplyr::arrange(Group) %>%
    dplyr::select(-Sample, -Phenotype, -Group)
  
  dfNum = dfNum[!grepl("nq", x = dfNum)]
  dfNum = dfNum[!grepl("n.d.", x = dfNum)]
  dfNum = dfNum[!grepl("n.q.", x = dfNum)]
  

  dfNum = as.matrix(dfNum)
  
  logm = log(dfNum, 10)
  
  logm = scale(logm, center = T, scale = T)
  row.names(logm) = paste(df$Group, df$Sample, sep = "_")
  row.names(dfNum) = paste(df$Group, df$Sample, sep = "_")
  

  pca <- prcomp(logm, center = F, scale = F)
  pca.var <- pca$sdev^2
  pca.var.per <- round(pca.var/sum(pca.var)*100, 1)

  barplot(pca.var.per, main="Scree Plot", xlab="Principal Component", ylab="Percent Variation")

  pcaresults <- summary(pca)

  scree.data <- as.data.frame(pcaresults$importance) # eigenvalues / standard deviations
  score.data <- as.data.frame(pcaresults$x) # coordinates of the samples (i.e., observations)
  loading.data <- as.data.frame(pcaresults$rotation) # loading scores for individual variables


  score.data <- score.data[, 1:2]
  score.data$Group <- df$Group
  score.data$Phenotype <- df$Phenotype

  perc.var = 100 * pcaresults[["importance"]][2, ]


  pca = ggplot(score.data, aes(x = PC1, y = PC2)) +
    geom_point(aes(shape = Group, colour = Phenotype)) +
    geom_label_repel(label = rownames(score.data)) +
    stat_ellipse(geom = "polygon", alpha = 1/2, aes(fill = Group)) +
    xlab(paste("PC1:", perc.var[1], "%")) +
    ylab(paste("PC2:", perc.var[2], "%")) +
    ggtitle(title) +
    my.theme


  order = sort(rownames(logm))
  logm = logm[order, ]
  t.logm = as.data.frame(t(logm))
  
  annot <- data.frame(Group = rep(c("Control-KO", "Control-WT",
                                    "IL-4-KO", "IL-4-WT"), each = 4))

  annot$Group <- factor(annot$Group, levels = c("Control-KO", "Control-WT",
                                                "IL-4-KO", "IL-4-WT"))

  rownames(annot) <- names(t.logm)

  p = pheatmap(t.logm,
           color = rev(colorRampPalette(c("#B2182B", "#FDDBC7", "#2166AC"))(n = 256)),
           border_color = "gray40",
           scale = "row",
           gaps_col = c(8,16),
           angle_col = 45,
           fontsize_row = 13,
           cluster_cols = F,
           cellwidth = 15,
           annotation_col =  annot,
           annotation_names_col = F,
           main = title
  )
  

# ANOVA #

  logm = t(t.logm)

  logdf = as.data.frame(logm) %>%
    dplyr::mutate(Group = annot$Group) %>%
    dplyr::relocate(Group, .before = names(logm)[1])

  metabolites.list = list()
  MetNames = c()
  
  

  for(i in 2:ncol(logdf)) { # if it doesn't work - comment out line 218 and 219 then run the for loop, then run it with removing comments - weird!

    column <- names(logdf[i]) # to print each ROI at the top of each ANOVA summary

    avz <- aov(logdf[,i] ~ Group, data = logdf) # Each ANOVA test iterating through each column of my dataframe

    result <- summary(avz) # summarize each ANOVA in a table
    hsd <- TukeyHSD(avz)

    print(column)
    print(result)
    print(hsd)

    metabolites.list[[i]] = hsd
    names(metabolites.list)[i] = column


    if (any(as.data.frame(metabolites.list[[i]]$Group)["IL-4-WT-IL-4-KO", 4] < 0.05,
            as.data.frame(metabolites.list[[i]]$Group)["Control-WT-Control-KO", 4] < 0.05)) {

      MetNames[i] = column
    }

  }


  MetNames = na.omit(MetNames)

  metabolites.list = metabolites.list[names(metabolites.list) %in% MetNames]

  logm = as.data.frame(logm)


  BoxPlotsDf = logm[names(logm) %in% names(metabolites.list)]

  BoxPlotsDf = logm
  BoxPlotsList = list()

  for (i in 1:ncol(logm)) {

    BoxPlotsDfSingle = matrix(logm[,i], ncol = 1, nrow = nrow(logm))
    colnames(BoxPlotsDfSingle)[1] = "ZScore"
    rownames(BoxPlotsDfSingle) = rownames(logm)
    BoxPlotsDfSingle = as.data.frame(BoxPlotsDfSingle)
    BoxPlotsDfSingle$Group = annot$Group

    BoxPlotsDfSingle = BoxPlotsDfSingle %>%
      dplyr::filter(Group %in% c("IL-4-WT", "IL-4-KO"))


    BoxPlotsList[[i]] = ggplot(BoxPlotsDfSingle, aes(x=Group, y=ZScore, fill=Group)) +
                        geom_boxplot(position=position_dodge(1)) +
                        geom_dotplot(binaxis='y', stackdir='center',
                        position=position_dodge(1), dotsize = 0.8) +
                        scale_fill_manual(
                                values=c("IL-4-KO" = "#FC4E2A",
                                 "IL-4-WT" = "#E69F00")) +
                        labs(x = NULL, y = NULL, title = names(logm)[i]) +
                        theme(legend.title = element_text(face = "bold"),
                              axis.text.x = element_blank(),
                              legend.position = "none") +
                        my.theme

    names(BoxPlotsList)[i] = names(logm)[i]

  }


  DatasetList = list(metabolites.list,
                     MetNames,
                     p,
                     pca,
                     BoxPlotsList,
                     )

  names(DatasetList) = c("PostHocStats", "PostHocMetNames", "Heatmap", "PCAPlot", "BoxPlotsList")

  return(DatasetList)


}


HILICList = AnalyzeFun(HILIC, "HILIC")
RPList = AnalyzeFun(RP, "RP")

RepeatHILICList = AnalyzeFun(RepeatHILIC, "RepeatHILIC")
RepeatRPList = AnalyzeFun(RepeatRP, "RepeatRP")


# plotting metabolites from CYP experiments with IL-4 focus (Figure 6):
PlotCombinedMets = function(Metabolite, whichDataset) {
  
  if (whichDataset == "HILIC") {
    
    HILICList[["BoxPlotsList"]][[Metabolite]][["data"]]$Experiment = rep("1st_experiment")
    RepeatHILICList[["BoxPlotsList"]][[Metabolite]][["data"]]$Experiment = rep("2nd_experiment")
    test_df = rbind(HILICList[["BoxPlotsList"]][[Metabolite]][["data"]], 
                    RepeatHILICList[["BoxPlotsList"]][[Metabolite]][["data"]])
    test_df$Metabolite = Metabolite
    
  } else {
    
    RPList[["BoxPlotsList"]][[Metabolite]][["data"]]$Experiment = rep("1st_experiment")
    RepeatRPList[["BoxPlotsList"]][[Metabolite]][["data"]]$Experiment = rep("2nd_experiment")
    test_df = rbind(RPList[["BoxPlotsList"]][[Metabolite]][["data"]], 
                    RepeatRPList[["BoxPlotsList"]][[Metabolite]][["data"]])
    test_df$Metabolite = Metabolite
  }
  
  ggplot(test_df, aes(x=Group, y=ZScore, fill=Group)) +
    geom_boxplot(position=position_dodge(1)) +
    geom_dotplot(binaxis='y', stackdir='center',
                 position=position_dodge(1), dotsize = 0.8) +
    scale_fill_manual(
      values=c("IL-4-KO" = "#FC4E2A",
               "IL-4-WT" = "#E69F00")) +
    labs(x = NULL, y = NULL, title = test_df$Metabolite[1]) +
    theme(legend.title = element_text(face = "bold"),
          axis.text.x = element_blank(),
          strip.text.x = element_text(
            size = 16, color = "black", face = "bold.italic"),
          legend.position = "none") +
    facet_wrap(vars(Experiment)) +
    my.theme
}
#check that it works:
PlotCombinedMets("Guanosine", "RP")
PlotCombinedMets("Serine", "HILIC")


############################
#### GAPDH metabolomics ####
############################



df = read.csv(paste0(path, "/GAPDH_experiment.csv"))

#df = df %>% dplyr::filter(Group == "IL-4")

dfNum = df %>%
  dplyr::select(-Sample, -Phenotype, -Group)

dfNum = dfNum[!grepl("nq", x = dfNum)]
dfNum = dfNum[!grepl("n.d.", x = dfNum)]
dfNum = dfNum[!grepl("n.q.", x = dfNum)]
dfNum = dfNum %>% mutate(across(where(is.character), as.numeric))


logm = log(dfNum, 10)

logm = scale(logm, center = T, scale = T)
row.names(logm) = paste(df$Group, df$Sample, sep = "_")
row.names(dfNum) = paste(df$Group, df$Sample, sep = "_")

pca <- prcomp(logm, center = F, scale = F)
pca.var <- pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100, 1)

barplot(pca.var.per, main="Scree Plot", xlab="Principal Component", ylab="Percent Variation")

pcaresults <- summary(pca)

scree.data <- as.data.frame(pcaresults$importance) # eigenvalues / standard deviations
score.data <- as.data.frame(pcaresults$x) # coordinates of the samples (i.e., observations)
loading.data <- as.data.frame(pcaresults$rotation) # loading scores for individual variables


score.data <- score.data[, 1:2]
score.data$Group <- df$Group
score.data$Phenotype <- df$Phenotype

perc.var = 100 * pcaresults[["importance"]][2, ]


pca = ggplot(score.data, aes(x = PC1, y = PC2)) +
  geom_point(aes(shape = Group, colour = Phenotype)) +
  geom_label_repel(label = rownames(score.data)) +
  stat_ellipse(geom = "polygon", alpha = 1/2, aes(fill = Phenotype)) +
  xlab(paste("PC1:", perc.var[1], "%")) +
  ylab(paste("PC2:", perc.var[2], "%")) +
  my.theme

pca

dfNum = df %>%
  dplyr::select(-Sample, -Phenotype, -Group)


t.logm = as.data.frame(t(logm))


pvals <- apply(t.logm, 1, function(x) {
  t.test(x[1:4], x[5:8])$p.value
})

pvals = p.adjust(pvals, method = "fdr")
pvals = sort(pvals[pvals < 0.05], decreasing = F)

annot <- data.frame(Group = rep(c("Control", "IL-4"), each = 8))
rownames(annot) <- names(t.logm)

ann_colors = list(
  Group = c(Control = "#56B4E9", `IL-4` = "#E69F00")
)

  pheatmap(t.logm,
             color = rev(colorRampPalette(c("#B2182B", "#FDDBC7", "#2166AC"))(n = 256)),
             border_color = "gray40",
             scale = "row",
             angle_col = 45,
             fontsize_row = 13,
             cluster_cols = F,
             cellwidth = 15,
             annotation_col =  annot,
             annotation_colors = ann_colors
)

# Supplementary Figure 5:
for (i in c("Glucose", "Glucose6Phosphate", "Fructose6Phosphate", "Phosphoglycerate", "Pyruvate", "PEP", "Lactate", "AcetylCoA", "Citrate")) {
  
  pheatmap(t.logm[i, ],
           color = rev(colorRampPalette(c("#B2182B", "#FDDBC7", "#2166AC"))(n = 256)),
           border_color = "gray40",
           scale = "row",
           angle_col = 45,
           show_rownames = F,
           show_colnames = F,
           cellwidth = 15,
           cellheight = 15,
           cluster_cols = F,
           cluster_rows = F,
           annotation_col =  annot,
           legend = T,
           annotation_colors = ann_colors,
           main = i
  )
}

p

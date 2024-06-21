file.choose()
library(xlsx)
library(dplyr)
library(readxl)

rp = read_excel("/Users/rokosango/PhD/Metabolomics/data/non-targ-metabolomics.xlsx",
                sheet = 6)
hilic = read_excel("/Users/rokosango/PhD/Metabolomics/data/non-targ-metabolomics.xlsx",
                   sheet = 5)


# cleaning the data for untargeted extracellular metabolomics.
# assumes Excel-level cleaning was already performed
# Purpose of this cleaning is to integrate this data with Genome scale metabolic models.
prep = function(x, outfile) {

setwd('/Users/rokosango/PhD/Metabolomics/data')
x = x %>%
  filter(Colour == "Green" | Colour == 'Yellow') %>%
  filter(Name != "Name")
x = x[!duplicated(x$Name),] 


x_num = x %>%
  dplyr::select(!c(Name, Colour))

x_num = x_num[,order(colnames(x_num))] # put all IL4 measurements together. same for LPS. makes it easy to do group calcs.

x_num = x_num %>%
  mutate_if(is.character, as.numeric)

x_num$meanIL4_24h = apply(x_num[c(1, 2, 3, 4)], 1, mean)
x_num$meanIL4_48h = apply(x_num[c(5, 6, 7, 8)], 1, mean)
x_num$meanLPS_24h = apply(x_num[c(9, 10, 11, 12)], 1, mean)
x_num$meanLPS_48h = apply(x_num[c(13, 14, 15, 16)], 1, mean)
x_num$meanCtrl_24h = apply(x_num[c(17, 18, 19, 20)], 1, mean)
x_num$meanCtrl_48h = apply(x_num[c(21, 22, 23, 24)], 1, mean)


x_num = x_num %>%
  select(meanIL4_24h, meanIL4_48h, meanLPS_24h,
         meanLPS_48h, meanCtrl_24h, meanCtrl_48h, 
         Padj_IL4_48_h_IL4_24_h, 
         Padj_LPS_48_h_LPS_24_h,
         Padj_Ctrl_48_h_Ctrl_24_h)
rownames(x_num) = x$Name
x_num = na.omit(x_num)

x_num$IL4_direction = ifelse(x_num$meanIL4_24h < x_num$meanIL4_48h,"uptake","secretion")
x_num$LPS_direction = ifelse(x_num$meanLPS_24h < x_num$meanLPS_48h,"uptake","secretion")

x_num = x_num %>%
  filter(Padj_IL4_48_h_IL4_24_h < 0.1 | Padj_LPS_48_h_LPS_24_h < 0.1 | Padj_Ctrl_48_h_Ctrl_24_h < 0.1) %>%
  dplyr::select(Padj_IL4_48_h_IL4_24_h, Padj_LPS_48_h_LPS_24_h, Padj_Ctrl_48_h_Ctrl_24_h,
                IL4_direction, LPS_direction, Ctrl_direction)

x_num$Padj_IL4_48_h_IL4_24_h = as.numeric(x_num$Padj_IL4_48_h_IL4_24_h)
x_num$Padj_LPS_48_h_LPS_24_h = as.numeric(x_num$Padj_LPS_48_h_LPS_24_h)
x_num$Padj_Ctrl_48_h_Ctrl_24_h = as.numeric(x_num$Padj_Ctrl_48_h_Ctrl_24_h)


x_num = x_num %>% mutate(WhichModel =
                             case_when(Padj_IL4_48_h_IL4_24_h < 0.1 & Padj_LPS_48_h_LPS_24_h < 0.1 ~ "Both", 
                                       Padj_IL4_48_h_IL4_24_h < 0.1 & Padj_LPS_48_h_LPS_24_h > 0.1 ~ "IL4",
                                       Padj_IL4_48_h_IL4_24_h > 0.1 & Padj_LPS_48_h_LPS_24_h < 0.1 ~ "LPS",
                                       Padj_IL4_48_h_IL4_24_h > 0.1 & Padj_LPS_48_h_LPS_24_h > 0.1 ~ "None")
                                
)

write.xlsx2(x_num, outfile)

}


prep(hilic, "hilic_num_full.xlsx")
prep(rp, "rp_num_full.xlsx")


# read and finalize the dataset for input to GEMs. Three datafiles, 1) LPS

RP = read.xlsx2("rp_num.xlsx", 
                sheetIndex = 1)
HILIC = read.xlsx2('hilic_num.xlsx', sheetIndex = 1)

Combined = rbind(HILIC, RP)

Combined = Combined %>%
  filter(Rxn_ID != "N/A") %>%
  arrange(WhichModel)

LPS = Combined %>%
  filter(WhichModel == "LPS" | WhichModel == "Both") %>%
  select(!c(Padj_IL4_48_h_IL4_24_h, IL4_direction)) %>%
  mutate(LB = ifelse(LPS_direction == "uptake", 0.5, 0)) %>%
  mutate(UB = ifelse(LPS_direction == "uptake", 0, 1000))
  

IL4 = Combined %>%
  filter(WhichModel == "IL4" | WhichModel == "Both") %>%
  select(!c(Padj_LPS_48_h_LPS_24_h, LPS_direction)) %>%
  mutate(LB = ifelse(IL4_direction == "uptake", 0.5, 0)) %>%
  mutate(UB = ifelse(IL4_direction == "uptake", 0, 1000))


write.xlsx2(Combined, "CombinedHilicRP.xlsx")
write.xlsx2(LPS, "ExMetForLPSGEM.xlsx")
write.xlsx2(IL4, "ExMetForIL4GEM.xlsx")


# filter LPS and IL4 datasets for MATLAB functions
file.choose()


LPS = read.xlsx2('/Users/rokosango/PhD/Metabolomics/data/ExMetForLPSGEM.xlsx',
                 sheetIndex = 1)


LPS = LPS %>%
  filter(Add_to_models == "Y")

write.xlsx2(LPS, "AddRxnsLPS.xlsx")


x = read_excel("/Users/rokosango/PhD/Metabolomics/data/non-targ-metabolomics.xlsx",
                sheet = 6)

prepControl = function(x, outfile) {
  
  setwd('/Users/rokosango/PhD/Metabolomics/data')
  
  x = x %>%
  filter(Colour == "Green" | Colour == 'Yellow') %>%
  filter(Name != "Name") %>%
  select(c(Name, contains("Ctrl"))) # selecting multiple columns (Name and all columns having "Ctrl"). Worth remembering!


x = x[!duplicated(x$Name),] 


x_num = x %>%
  dplyr::select(!Name)

x_num = x_num %>%
  mutate_if(is.character, as.numeric)

x_num$meanCtrl_24h = apply(x_num[c(2, 3, 4, 5)], 1, mean)
x_num$meanCtrl_48h = apply(x_num[c(6, 7, 8, 9)], 1, mean)

x_num = x_num %>%
  select(meanCtrl_24h, meanCtrl_48h,
         Padj_Ctrl_48_h_Ctrl_24_h)
x_num = as.data.frame(x_num)
rownames(x_num) = x$Name
x_num = na.omit(x_num)

x_num$Ctrl_direction = ifelse(x_num$meanCtrl_24h < x_num$meanCtrl_48h,"uptake","secretion")


x_num = x_num %>%
  filter(Padj_Ctrl_48_h_Ctrl_24_h < 0.1) %>%
  dplyr::select(Padj_Ctrl_48_h_Ctrl_24_h, Ctrl_direction)
x_num$Padj_Ctrl_48_h_Ctrl_24_h = as.numeric(x_num$Padj_Ctrl_48_h_Ctrl_24_h)

write.xlsx2(x_num, outfile)

}

prepControl(rp, "rpc.xlsx")
prepControl(hilic, "hilicc.xlsx")

# excel work done. now load them back up and combine
RP = read_excel("rpc.xlsx", 
                sheet = 1)
HILIC = read_excel('hilicc.xlsx', sheet = 1)

Combined = rbind(HILIC, RP)

Combined = Combined %>%
  filter(Rxn_ID != "N/A") %>%
  mutate(LB = ifelse(Ctrl_direction == "uptake", 0.5, 0)) %>%
  mutate(UB = ifelse(Ctrl_direction == "uptake", 0, 1000))


write.xlsx2(Combined, "Ctrl_combined.xlsx")


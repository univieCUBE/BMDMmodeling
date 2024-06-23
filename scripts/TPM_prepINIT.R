
options(java.parameters = "- Xmx1024m") # parameter needed for I/O big files (e.g., transcript files)

#########################################################################################################
################# Salmon output files (.sf) are needed for this step as they contain information ######## 
################# about TPM, and gene symbols which are needed for tINIT2.###############################
#########################################################################################################
library(stringr)
library(conflicted)
library(dplyr)
library(lattice)
library(ggplot2)
library(org.Mm.eg.db)
library(xlsx)
require(httr)
require(jsonlite)
setwd('/Users/rokosango/PhD/RNA-seq/quants') #where the salmon quantification files are

files = c('LPS_1.fastq_quant/quant.sf', #change to "LPS" and "Control"
          'LPS_2.fastq_quant/quant.sf',
          'LPS_3.fastq_quant/quant.sf',
          'LPS_4.fastq_quant/quant.sf')

for(i in 1:length(files)) {                           
  assign(paste0("LPS_", i),                                 
         read.delim(files[i]))

}

my_dplyr <- function(x) { 
  x %>%       
    dplyr::select(Name, TPM)
}

result = list(Control_1, Control_2, Control_3, Control_4) %>%
  lapply(my_dplyr)

for (i in 1:4) {
result[[i]][[1]] = substr(result[[i]][[1]], start = 22, stop = 39)
result[[i]][[1]] = str_remove(result[[i]][[1]], "[|]")

}

pick_dupl_max_val = function (df) {
  df = as.data.frame(df[with(df, ave(TPM, Name, FUN= mean)==TPM),]) #for each duplicate value in "Name", find max value in "TPM"
}
 
TPMS = data.frame( #keep changing this for each condition
      Control_1 = result[[1]][[2]],
      Control_2 = result[[2]][[2]],
      Control_3 = result[[3]][[2]],
      Control_4 = result[[4]][[2]])
TPMS = apply(TPMS, 1, median)
TPMS = cbind(result[[1]][[1]], TPMS)
TPMS = as.data.frame(TPMS)
names(TPMS)[1:2] = c("Name", "TPM")

TPMS = pick_dupl_max_val(TPMS)
TPMS = TPMS[!duplicated(TPMS), ]

###
# Multiple IDs to convert - use a POST request
###
url = "http://biotools.fr/mouse/ensembl_symbol_converter/"
ids = as.vector(TPMS$Name)
ids_json <- toJSON(ids)

body <- list(api=1, ids=ids_json)
r <- POST(url, body = body)

output = fromJSON( content(r, "text"), flatten=TRUE)

TPMS = TPMS %>% 
  mutate(GeneName = output) %>%
  filter(!GeneName == "NULL") %>%
  dplyr::select(Name, TPM, GeneName) %>%
  rename(genes = GeneName) %>%
  rename(IL4 = TPM)

setwd('/Users/rokosango/PhD/RNA-seq/quants/TPMSRefined')

options(max.print=999999)
capture.output(TPMS, file = "IL4MedianTPMS.csv")

table(TPMS$IL4 > 1) #how many genes with TPM > 1?

#Control
#FALSE  TRUE 
#39257 10243

#LPS
#FALSE  TRUE 
#39885  9615

#IL4
#FALSE  TRUE 
#39497 10003





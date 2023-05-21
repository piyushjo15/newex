library(reticulate)
Sys.setenv(RETICULATE_PYTHON = "/home/p541i/.conda/envs/pjexv2/bin/python3.7")
use_python("/home/p541i/.conda/envs/pjexv2/bin/python3.7")

## using module score of GCUBC to 
## obtain NMF, UMAP and slingshot trend components for GCUBC
## 1 load data -----
#rescale from negative to positive value
library(uwot)

to_pos <- function(x){
  a2 <- min(x)
  if(a2 <0){ 
    x <- x-a2
  }
  return(x)
  
}

#GCUBC mod score ---------
load("~/DATA/kons_analysis/indana/pdGCUBC_TFselClGEP_MSv2.RData")
st <-c("progenitor_RL_early", "progenitor_RL","GCP/UBCP","UBC_diff","UBC_Trpc3","UBC_Hcrtr2")
keep <- plot.data$subtype %in% st
plot.data <- plot.data[keep,]
# selx <- read.delim("MBnewana/TFselClGRNs_sel.txt")
# sel <- selx$GRNs
load("~/DATA/kons_analysis/indana/SCENICana/G34MB_nn/TFselcl_GEPv2.RData") #G34MB_nn GEP
sel<- names(GSEA)
##not used 
# cx <- c("ELF2","JUN","RLF","FOS",
#         "EPAS1","KDM5B","THRB","HIF3A","NR3C1")
# sel <- sel[!(sel %in% cx)]


head(sel)
pdUBC <- plot.data[,c(1:5)]
msnor <- plot.data[,sel]
msnor1 <- apply(msnor, 2, to_pos)
dim(msnor1)
## 2 using python code for NMF and UMAP stuff-----
repl_python()
import sklearn.decomposition as sk
import numpy as np
import umap
import random

model = sk.NMF(n_components=25, max_iter=10000, random_state=0, init="nndsvd")

random.seed(456)
model.fit(r.msnor1)
wn = model.transform(r.msnor1)

# modu = umap.UMAP(n_neighbors=10, random_state=4, min_dist=0.3)
# modu.fit(wn)
# random.seed(456)
# un = modu.transform(wn)

exit
###post python -----
nmf_n <- py$wn
#u_n <- py$un

set.seed(4)
u_n <- umap(nmf_n,n_neighbors=25, min_dist=0.3, ret_model = TRUE)

u_na <- u_n$embedding

pdUBC$UMAP1 <- u_na[,1]
pdUBC$UMAP2 <- u_na[,2]

save(nmf_n,pdUBC, file = "diffG34MBUBCref_r25_GEP.RData")
q()

##plot ------

library(ggplot2)
library(ggsci)
load("diffG34MBUBCref_r25_GEP.RData")
##color scale-----
st <-c("progenitor_RL_early", "progenitor_RL","GCP/UBCP","UBC_diff","UBC_Trpc3","UBC_Hcrtr2")
values1 <- read.delim("~/DATA/kons_analysis/indana/subtype_col.txt",row.names = 1)
valuesx1 <- values1[st,1]

dls <- c("progenitor","GCP/UBCP","UBC_diff","UBC_defined")
values2 <- read.csv("~/DATA/kons_analysis/indana/dev_state_colv3.csv",row.names = 1)
valuesx2 <- values2[dls,1]

##
pdUBC$subtype <- factor(pdUBC$subtype, levels = st)
pdUBC$dev_state <- factor(pdUBC$dev_state, levels = dls)
p <- ggplot(data=pdUBC, aes(x=UMAP1, y=UMAP2, color=dev_state)) +
  geom_point(size=.5 )+
  scale_color_manual(values=valuesx2)+
  theme_classic()+
  theme(legend.position='none',
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank())
tiff("UBC_MS_UMAP_DSr25GEP.tiff", width = 6, height = 4, units = "in", res = 300)
print(p)
dev.off()



##slingshot-----
library(slingshot)
cl <- as.character(pdUBC$subtype)
sl <- slingshot(u_na, cl, start.clus = "progenitor_RL_early",
                end.clus="UBC_Hcrtr2", reducedDim=NULL)

curve <- sl@curves$curve1
curveline <- curve$s
colnames(curveline) <- c("X","Y")
pdUBC <- cbind(pdUBC,curveline)

##combine plot 
p <- ggplot() +
  geom_point(data=pdUBC, aes(x=UMAP1, y=UMAP2, color=dev_state), size=.5 )+
  scale_color_manual(values=valuesx2)+
  theme_classic()+
  # add trend line plot
  geom_point(data=pdUBC, aes(x=X, y=Y), size=0.1)+
  theme(legend.position='none',
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank())
tiff(filename = "UBC_slingshot_DSr25GEP.tiff", width = 6, height = 4, units = "in", res = 300)
print(p)
dev.off()

library(ggridges)
pdUBC$Rank <- curve$lambda
p1 <- ggplot(pdUBC, aes(x=Rank, y=dev_state, fill=dev_state) ) +
  geom_density_ridges()+
  scale_fill_manual(values=valuesx2)+
  theme_classic()+
  theme(legend.position='none',
        axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank())
tiff(filename = "UBC_PST_DSr25GEP.tiff", width = 6, height = 4, units = "in", res = 300)
print(p1)
dev.off()
save(pdUBC,sl, file = "UBC_TFselClGEP_slingshotr25.RData")

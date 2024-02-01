library(ggplot2)
library(reshape2)
options(stringsAsFactors = F)
ggopt <-  theme_bw() +
  theme(strip.text = element_text(size=16, face = "bold"), 
        strip.background = element_rect(colour="white", fill="white"), 
        legend.position = 'bottom',  
        axis.title = element_text(size=16), 
        axis.text.y = element_text(size=13), 
        axis.text.x = element_text(size=13),
        plot.title = element_text(hjust = 0.5), 
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 16), 
        legend.title = element_blank()
    )

setwd("./DataAnalysis/Scerevisiae/")
colors = c("#FFCC00","#5599FF","#D40000")
complexity = read.table("antisense_prot_complexity.txt", header = T)
dkl = read.table("dKL_protein_antisenseprot.txt", header = T)
disorder = read.table("antisense_prot_iupred3.txt", header = T)[,c(1,2,4)]
colnames(disorder)[3] = "med.disorder"

fcombs = t(combn(unique(complexity$frame),2))
alldata = merge(disorder,merge(complexity,dkl, by=c("name","frame")), by=c("name","frame"))
metrics = c("dKL","ShannonH","lcomplexity","med.disorder")
wcxstats = do.call(rbind,lapply(c(1:3), function(z){
    do.call(rbind,lapply(metrics, function(m){
        xdata = alldata[alldata$frame==fcombs[z,1],m] 
        ydata = alldata[alldata$frame==fcombs[z,2],m]
        tail = if(median(xdata)>median(ydata)) 'g' else 'l'
        pval = wilcox.test(xdata,ydata,alternative = tail)$p.value
        data.frame(DataX = fcombs[z,1], DataY = fcombs[z,2], Metric = m, medX = median(xdata), medY = median(ydata), tail = tail, pvalue = pval)
    }))
}))

for(i in unique(wcxstats$Metric)){
    wcxstats[wcxstats$Variable==i,"pvalue"] = p.adjust(wcxstats[wcxstats$Variable==i,"pvalue"], method = 'bonf')
}

mdata = melt(alldata,id.vars = c("frame","name"), variable.name = "Metric")

plots = ggplot(mdata, aes(x=value, y = ..scaled.., color = frame)) + geom_density(size = 1.5) + ggopt + facet_wrap(~Metric, ncol = 2, scales = "free") + scale_color_manual(values = colors) + labs(y="Scaled Density", x = "Value")

ggsave(plots,filename = "sequence_composition_plots.pdf", width = 20, height = 20, units = "cm", device = "pdf")

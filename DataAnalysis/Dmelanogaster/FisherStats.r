options(stringsAsFactors = F)
setwd("../Dmelanogaster/")
obs = read.table("ORFcountdata.csv", header = T)

results = do.call(rbind, lapply(obs$Population, function(x){
    dx = obs[obs$Population==x,]
    ObsAS0 = c(dx$Obs_asORF0, dx$asRF_loci-dx$Obs_asORF0)
    ObsAS1 = c(dx$Obs_asORF1, dx$asRF_loci-dx$Obs_asORF1)
    ObsAS2 = c(dx$Obs_asORF2, dx$asRF_loci-dx$Obs_asORF2)
    ObsIG = c(dx$Obs_igORF, dx$igRF_loci-dx$Obs_igORF)

    ExpAS0 = c(dx$Exp_asORF0, dx$asRF_loci-dx$Exp_asORF0)
    ExpAS1 = c(dx$Exp_asORF0, dx$asRF_loci-dx$Exp_asORF1)
    ExpAS2 = c(dx$Exp_asORF0, dx$asRF_loci-dx$Exp_asORF2)
    ExpIG = c(dx$Exp_igORF, dx$igRF_loci-dx$Exp_igORF)

    f0I = fisher.test(rbind(ObsAS0,ObsIG), alternative = "l")
    f1I = fisher.test(rbind(ObsAS1,ObsIG), alternative = "l")
    f2I = fisher.test(rbind(ObsAS2,ObsIG), alternative = "l")

    f10 = fisher.test(rbind(ObsAS1,ObsAS0), alternative = "g")
    f12 = fisher.test(rbind(ObsAS1,ObsAS2), alternative = "g")

    f0x = fisher.test(rbind(ObsAS0,ExpAS0), alternative = "g")
    f1x = fisher.test(rbind(ObsAS1,ExpAS1), alternative = "g")
    f2x = fisher.test(rbind(ObsAS2,ExpAS2), alternative = "g")
    fIx = fisher.test(rbind(ObsIG,ExpIG), alternative = "l")
    alltests = list(f0I,f1I,f2I,f10,f12,f0x,f1x,f2x,fIx)
    do.call(rbind,lapply(c(1:9), function(z){
        dstr =  alltests[[z]]$data.name
        dpair = unlist(strsplit(substr(dstr,7,nchar(dstr)-1),", "))
        data.frame(Population = x, data1 = dpair[1], data2 = dpair[2],alternative = alltests[[z]]$alternative, odds = alltests[[z]]$estimate, pval = alltests[[z]]$p.value)
    }))
}))
rownames(results) <- NULL
results$pval = p.adjust(results$pval, method = 'bonf')

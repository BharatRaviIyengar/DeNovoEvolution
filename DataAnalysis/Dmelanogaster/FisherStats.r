options(stringsAsFactors = F)
setwd("../Dmelanogaster/")
obs = read.table("ORFcountdata.csv", header = T)

results = do.call(rbind, lapply(unique(obs$pop), function(x){
    xx = obs[obs$pop==x,]
	print(x)
    d1 = do.call(rbind,lapply(c(1:4), function(f){
	mat = rbind(c(xx[f,"ows"],xx[f,"loci"]-xx[f,"ows"]),c(xx[f,"expG"],xx[f,"loci"]-xx[f,"expG"]))
	if(xx[f,"ows"]>xx[f,"expG"]){
		tail = "g" 
		test = "Obs>Exp"
	} else{
		tail = "l"
		test = "Obs<Exp"
	}
	p = fisher.test(mat,alternative = tail)$p.value
	data.frame(Line = x, frame = xx[f,"frame"], test = test, pval = p)
	}))
    
	d2 = do.call(rbind,lapply(c(1:3), function(f){
	mat = rbind(c(xx[f,"ows"],xx[f,"loci"]-xx[f,"ows"]),c(xx[4,"ows"],xx[4,"loci"]-xx[4,"ows"]))
	if(xx[f,"ows"]/xx[f,"loci"] > xx[4,"ows"]/xx[4,"loci"]){
		tail = "g" 
		test = "as>ig"
	} else{
		tail = "l"
		test = "as<ig"
	}
	p = fisher.test(mat,alternative = tail)$p.value
	data.frame(Line = x, frame = xx[f,"frame"], test = test, pval = p)
	}))
	rbind(d1,d2)
}))
rownames(results) <- NULL
results$pval = p.adjust(results$pval, method = 'fdr')
   Line frame    test          pval
1    AK     0 Obs<Exp  1.000000e+00
2    AK     1 Obs<Exp  1.000000e+00
3    AK     2 Obs<Exp  1.426626e-03
4    AK    IG Obs>Exp 4.165805e-160
5    AK     0   as<ig 2.021329e-212
6    AK     1   as<ig 9.103113e-171
7    AK     2   as<ig 5.679093e-204
8    DK     0 Obs<Exp  1.000000e+00
9    DK     1 Obs<Exp  1.000000e+00
10   DK     2 Obs<Exp  2.015483e-12
11   DK    IG Obs>Exp 2.975263e-199
12   DK     0   as<ig 1.223371e-222
13   DK     1   as<ig 2.339881e-190
14   DK     2   as<ig 2.712468e-263
15   UM     0 Obs<Exp  2.214780e-02
16   UM     1 Obs>Exp  1.000000e+00
17   UM     2 Obs<Exp  1.099827e-10
18   UM    IG Obs>Exp 1.773450e-158
19   UM     0   as<ig 1.405953e-208
20   UM     1   as<ig 5.380123e-133
21   UM     2   as<ig 1.795443e-220
22   GI     0 Obs<Exp  4.016167e-02
23   GI     1 Obs<Exp  1.000000e+00
24   GI     2 Obs<Exp  1.086848e-02
25   GI    IG Obs>Exp 1.449433e-138
26   GI     0   as<ig 4.070331e-253
27   GI     1   as<ig 2.012301e-191
28   GI     2   as<ig 4.564539e-214
29   SW     0 Obs<Exp  3.887765e-01
30   SW     1 Obs<Exp  1.000000e+00
31   SW     2 Obs<Exp  8.139384e-09
32   SW    IG Obs>Exp 1.030224e-125
33   SW     0   as<ig 1.067335e-189
34   SW     1   as<ig 2.494980e-143
35   SW     2   as<ig 1.073527e-203
36   YE     0 Obs<Exp  4.784605e-01
37   YE     1 Obs>Exp  1.000000e+00
38   YE     2 Obs<Exp  1.004756e-14
39   YE    IG Obs>Exp 2.167557e-150
40   YE     0   as<ig 8.560770e-285
41   YE     1   as<ig 3.725396e-181
42   YE     2   as<ig  0.000000e+00
43   ZB     0 Obs<Exp  3.802174e-05
44   ZB     1 Obs<Exp  1.000000e+00
45   ZB     2 Obs<Exp  9.842928e-04
46   ZB    IG Obs>Exp 6.783207e-158
47   ZB     0   as<ig 7.788189e-222
48   ZB     1   as<ig 1.027642e-146
49   ZB     2   as<ig 1.190905e-171

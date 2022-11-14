cd(Base.source_dir())
include("nucleotidefuncts.jl")
stop = ["TAA","TAG","TGA"]

function nostop(codonset)
	codonset[codonset .∉ Ref(stop)]
end

gccontent = 0.53942 # all CDS #
nostopcodons = nostop(allcodons);

codonpairprob = (cpair) -> sum(nucprob.(cpair[:,1],gccontent) .* nucprob.(cpair[:,2],gccontent))

stopneighborsR1 = unique(vcat(frameXcpairs.(stop,-1)...), dims = 1);

stopNbrR1stopUP = stopneighborsR1[(stopneighborsR1[:,1] .∈ Ref(stop)),:];

stopNbrR1stopDN = stopneighborsR1[(stopneighborsR1[:,1].∉ Ref(stop)) .& (stopneighborsR1[:,2].∈ Ref(stop)),:];

stopNbrR1within = stopneighborsR1[(stopneighborsR1[:,1].∉ Ref(stop)) .& (stopneighborsR1[:,2].∉ Ref(stop)),:];

stopneighborsR2 = unique(vcat(frameXcpairs.(stop,-2)...), dims = 1);

stopNbrR2stopUP = stopneighborsR2[(stopneighborsR2[:,1] .∈ Ref(stop)),:];

stopNbrR2stopDN = stopneighborsR2[(stopneighborsR2[:,1].∉ Ref(stop)) .& (stopneighborsR2[:,2].∈ Ref(stop)),:];

stopNbrR2within = stopneighborsR2[(stopneighborsR2[:,1].∉ Ref(stop)) .& (stopneighborsR2[:,2].∉ Ref(stop)),:];

stopneighborsR0 = reverse.(comp.(stop));

p_stopR1within = codonpairprob(stopnbrR1within);
p_stopR1stopUP = codonpairprob(stopnbrR1stopUP);
p_stopR1stopDN = codonpairprob(stopnbrR1stopDN);

p_stopR2within = codonpairprob(stopnbrR2within);
p_stopR2stopUP = codonpairprob(stopnbrR2stopUP);
p_stopR2stopDN = codonpairprob(stopnbrR2stopDN);

probstopR0 = sum(nucprob.(stopneighborsR0,gccontent));

initneighborsR1 = unique(frameXcpairs("init",-1), dims = 1);
initNbrR1within = initneighborsR1[(initneighborsR1[:,1] .∉ Ref(stop)) .& (initneighborsR1[:,2].∉ Ref(stop)),:];

initneighborsR2 = unique(frameXcpairs("init",-2), dims = 1);
initNbrR2within = initneighborsR2[(initneighborsR2[:,1] .∉ Ref(stop)) .& (initneighborsR2[:,2].∉ Ref(stop)),:];

probinitR1 = sum(nucprob.(initneighborsR1[:,1],gccontent).*nucprob.(initneighborsR1[:,2],gccontent));
probinitR2 = sum(nucprob.(initneighborsR2[:,1],gccontent).*nucprob.(initneighborsR2[:,2],gccontent));


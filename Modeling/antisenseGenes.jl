cd(Base.source_dir())
include("nucleotidefuncts.jl")
stop = ["TAA","TAG","TGA"]

function nostop(codonset)
	codonset[codonset .∉ Ref(stop)]
end

gccontent = 0.53942 # all CDS #
nostopcodons = nostop(allcodons);

stopneighborsR1 = unique(vcat(frameXcpairs.(stop,-1)...), dims = 1);
stopneighborsR1 = stopneighborsR1[(stopneighborsR1[:,1].∉ Ref(stop)) .& (stopneighborsR1[:,2].∉ Ref(stop)),:];

stopneighborsR2 = unique(vcat(frameXcpairs.(stop,-2)...), dims = 1);
stopneighborsR2 = stopneighborsR2[(stopneighborsR2[:,1].∉ Ref(stop)) .& (stopneighborsR2[:,2].∉ Ref(stop)),:];

stopneighborsR0 = reverse.(comp.(stop));

probstopR1 = sum(nucprob.(stopneighborsR1[:,1],gccontent).*nucprob.(stopneighborsR1[:,2],gccontent));
probstopR2 = sum(nucprob.(stopneighborsR2[:,1],gccontent).*nucprob.(stopneighborsR2[:,2],gccontent));

probstopR0 = sum(nucprob.(stopneighborsR0,gccontent));

ATGneighborsR1 = unique(frameXcpairs("ATG",-1), dims = 1);
ATGneighborsR1 = ATGneighborsR1[(ATGneighborsR1[:,1].!= Ref(stop)) .& (ATGneighborsR1[:,2].∉ Ref(stop)),:];

ATGneighborsR2 = unique(frameXcpairs("ATG",-2), dims = 1);
ATGneighborsR2 = ATGneighborsR2[(ATGneighborsR2[:,1].∉ Ref(stop)) .& (ATGneighborsR2[:,2].∉ Ref(stop)),:];

probATGR1 = sum(nucprob.(ATGneighborsR1[:,1],gccontent).*nucprob.(ATGneighborsR1[:,2],gccontent));
probATGR2 = sum(nucprob.(ATGneighborsR2[:,1],gccontent).*nucprob.(ATGneighborsR2[:,2],gccontent));


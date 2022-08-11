cd(Base.source_dir())
include("nucleotidefuncts.jl")

blosum = readdlm("BLOSUM62N.csv", ',', Int, '\n'; header = true)
aasubmat = blosum[1][1:20,1:20]
aanames = Dict([(blosum[2][i],i) for i in 1:20])

allcodons = kmers(3)
nostopcodons = allcodons[allcodons .∉ Ref(["TAA","TAG","TGA"])]

f1cdn = permutedims(hcat(frameXcodons.(allcodons,1)...))
f2cdn = permutedims(hcat(frameXcodons.(allcodons,2)...))

r1cdn = reverse.(comp.(f1cdn))
r2cdn = reverse.(comp.(f2cdn))

gccontent = 0.42
function dfecdn(codon1, codon2)
        mp = nucprob(codon1,gccontent)*mprob(codon1,codon2)
        fe = aasubmat[aanames[transnucs(codon1)],aanames[transnucs(codon2)]]
        return [-fe mp]
end

dfedistv = vcat(reduce(vcat,[[dfecdn(y,x) for x in singlemutant(y)] for y in nostopcodons])...)

dfedistu = reduce(vcat,[[i sum(dfedistv[dfedistv[:,1] .== i,2])] for i in unique(dfedistv[:,1])])
dfedistu = dfedistu./sum(dfedistu, dims = 1)
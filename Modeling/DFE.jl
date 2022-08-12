cd(Base.source_dir())
include("nucleotidefuncts.jl")

aasubdata = readdlm("BLOSUM62N.csv", ',', Int, '\n'; header = true)
aasubmat = zeros(21,21)
aasubmat[1:20,1:20] = aasubdata[1][1:20,1:20]
aasubmat[21,1:20] .= -100
aasubmat[1:20,21] .= -90
aanames = Dict(vcat([(aasubdata[2][i],i) for i in 1:20],("O",21)))
allcodons = kmers(3)

function nostop(codonset)
        codonset[codonset .∉ Ref(["TAA","TAG","TGA"])]
end

nostopcodons = nostop(allcodons)

f1cdn = permutedims(hcat(frameXcodons.(nostopcodons,1)...))
f2cdn = permutedims(hcat(frameXcodons.(nostopcodons,2)...))

r1cdn = reverse.(comp.(f1cdn))
r2cdn = reverse.(comp.(f2cdn))

gccontent = 0.42

function dfecdn(codon1, codon2)
        mp = nucprob(codon1,gccontent)*mprob(codon1,codon2)
        fe = aasubmat[aanames[transnucs(codon1)],aanames[transnucs(codon2)]]
        return [-fe mp]
end

function dfedistgen(codonset1,codonset2)
        return vcat(reduce(vcat,[[dfecdn(y,x) for x in codonset2] for y in codonset1])...)
end

function combinedfe(dfeset)
        return reduce(vcat,[[i sum(dfeset[dfeset[:,1] .== i,2])] for i in unique(dfeset[:,1])])
end

function expdfe(dfeset)
        return sum(dfeset[:,1].*dfeset[:,2])/sum(dfeset[:,2])
end

dfedist0 = dfedistgen(nostopcodons,nostopcodons)

dfeXframeshift = zeros(61,4)

frames = [-2,-1,1,2]

sdiff2 = (s1,s2) -> s1[s1 .∉ Ref(s2)]

for i in 1:61
        for j in 1:4
               cdn = nostopcodons[i]
               frm = frames[j]
               fset0 = frameXcodons(cdn,frm)
               cset1 = singlemutant(cdn)
               fset1 = hcat([frameXcodons(x,frm) for x in cset1]...)
               # dfe1 = combinedfe(vcat) edit this
        end
end
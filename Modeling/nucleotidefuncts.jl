using Base.Iterators
using DelimitedFiles

gencode = readdlm(joinpath(Base.source_dir(),"../std_genetic_code.txt"), '\t', AbstractString, '\n')

nucs = "ATGC"
nucv = ['A','T','G','C']
nucnames = Dict('A'=>1,'T'=>2,'G'=>3,'C'=>4)
ncomp = Dict('A'=>'T', 'G'=>'C', 'T'=>'A', 'C'=> 'G')

aas = ["A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y"];

hydropathy = [1.8, 2.5, -3.5, -3.5, 2.8, -0.4, -3.2, 4.5, -3.9, 3.8, 1.9, -3.5, -1.6, -3.5, -4.5, -0.8, -0.7, 4.2, -0.9, -1.3]; # Kyte and Doolittle 

"""
# Generate all nucleotide sequences with length `k`
"""
function kmers(k)
    return vec(join.(product(repeated(nucs,k)...)))
end


"""
`frameXcodons(codon,n)`

Find all possibe codons in an alternate reading frame that overlap with a given codon
...

# Arguments
 - `codon::String`
 - `n::Int ∈ [-2,-1,1,2]`: reading frame (same sense)
# Output
 `x::Vector{String}`: list of 20 codons

# Examples
`julia> frameXcodons("ATG",1)`

`20-element Vector{String}:`

```"AAA"
"TAA"
"GAA"
"CAA"
"ATA"
"TTA"
"GTA"
"CTA"
"AGA"
"TGA"
"GGA"
"CGA"
"ACA"
"TCA"
"GCA"
"CCA"
"TGA"
"TGT"
"TGG"
"TGC"
 ```
"""
function frameXcodons(codon::String,frame::Int)
    if(frame ∉ [-2,-1,1,2])
        error("Frame should be ±1, or ±2")
    end
    frabs = abs(frame)
    z = 3-frabs
    down = codon[1+frabs:3].*vec(join.(product(repeated(nucs,frabs)...)))
    up = vec(join.(product(repeated(nucs,z)...))).*codon[1:frabs]
    both = vcat(up,down)
    if(frame>0)
        return both
    else
        return reverse.(comp.(both))
    end
end

"""
`comp(s)`
...
Find the complement of a nucleotide sequence, `s`
"""
function comp(s::String)
    join([ncomp[c] for c in s])
end



nsub = zeros(4,4)

nsub[1,2] = 0.056043956043956
nsub[1,3] = 0.242857142857143
nsub[1,4] = 0.074725274725275
nsub[3,1] = 0.438461538461539
nsub[3,2] = 0.075164835164835
nsub[3,4] = 1 - sum(nsub[1,:]) - sum(nsub[3,1:2])
nsub[2,1] = nsub[1,2]
nsub[2,3] = nsub[1,4]
nsub[2,4] = nsub[1,3]
nsub[4,3] = nsub[3,4]
nsub[4,2] = nsub[3,1]
nsub[4,1] = nsub[3,2]

mutrate = 7.8e-9

function µA(mutrate,nsub)
    return 2*mutrate*sum(nsub[1,:])
end

function µG(mutrate,nsub)
    return 2*mutrate*sum(nsub[3,:])
end


"""
`transnucs(codon)`
...
Find the amino acid encoded by a codon
"""
function transnucs(codon)
    return gencode[gencode[:,1].==codon,:][3]
end

function degen(consensus,position,subs)
    prefix = consensus[1:position-1]
    suffix = consensus[position+1:end]
    variants = prefix .* collect(subs) .* suffix
    return(variants)
end

function singlemutant(s)
    z = vcat([degen(s,x,nucs) for x in 1:3]...)
    return z[z.!=s]
end

function mprob(s1,s2)
    l = length(s1)
    if l != length(s2)
        error("length mismatch")
    end
    z = collect(s1) .!= collect(s2)
    if any(z)
        p = 1
        for i = (1:l)[z]
            p = p*nsub[nucnames[s1[i]], nucnames[s2[i]]]
        end
    else
        p = 0
    end
    return p*(2*mutrate)^sum(z)
end

function nucprob(x,gccontent)
    S = 0.5*gccontent
    W = 0.5 - S
    nw = count(r"[AT]", x)
    return W^nw * S^(length(x)-nw)
end

function numAT(x)
    return count(r"[AT]", x)
end

function featuregain(set1,set2,gccontent)
    v1 = sum([nucprob(s1,gccontent)*sum([mprob(s1,s2) for s2 in set2])
                for s1 in set1])
    return v1
end

function featurestay(set1,gccontent)
    lenseq =length(set1[1])
    v =  sum(
        [nucprob(s1,gccontent)*(
            sum([mprob(s1,s2) for s2 in set1]) +
            (1-µA(mutrate,nsub))^numAT(s1)*
            (1-µG(mutrate,nsub))^(lenseq-numAT(s1))
            ) for s1 in set1]
        )
    return v
end

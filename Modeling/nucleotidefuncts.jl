using Base.Iterators
using DelimitedFiles

gencode = readdlm(joinpath(Base.source_dir(),"../std_genetic_code.txt"), '\t', AbstractString, '\n')

nucs = "ATGC"
nucv = ['A','T','G','C']
nucnames = Dict('A'=>1,'T'=>2,'G'=>3,'C'=>4)
ncomp = Dict('A'=>'T', 'G'=>'C', 'T'=>'A', 'C'=> 'G')

aas = ["A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y"];

hydropathy = [1.8, 2.5, -3.5, -3.5, 2.8, -0.4, -3.2, 4.5, -3.9, 3.8, 1.9, -3.5, -1.6, -3.5, -4.5, -0.8, -0.7, 4.2, -0.9, -1.3]; # Kyte and Doolittle 

hydropathy2 = [0.17, -0.24, 1.23, 2.02, -1.13, 0.01, 0.96, -0.31, 0.99, -0.56, -0.23, 0.42, 0.45, 0.58, 0.81, 0.13, 0.14, 0.07, -1.85, -0.94
]; # Wimley and White

hydropathy3 = [0.5, -0.02, 3.64, 3.63, -1.71, 1.15, 2.33, -1.12, 2.8, -1.25, -0.67, 0.85, 0.14, 0.77, 1.81, 0.46, 0.25, -0.46, -2.09, -0.71
]; # Wimley and White
"""
# Generate all nucleotide sequences with length `k`
"""
function kmers(k)
    return vec(join.(product(repeated(nucs,k)...)))
end

normalize = x -> x/sum(x)

function randomize_nsub(nsm,sigma)
    nx = vcat(nsm[1,:],nsm[3,:]);
    nx = nx .+ randn(Float64,size(nx)).*sigma.*nx;
    nx = nx ./ sum(nx);
    nx[[1,7]] .= 0;
    smout = zeros(size(nsm))
    smout[1,:] = nx[1:4]
    smout[3,:] = nx[5:8]
    for i = [2,4]
        for j = 1:4
            if(j!=i)
                x = nucnames[ncomp[nucv[i]]];
                y = nucnames[ncomp[nucv[j]]];
                smout[i,j] = smout[x,y];
            end
        end
    end
    return smout
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

function frameXcpairs(codon::String,frame::Int)
    if(frame ∉ [-2,-1,1,2])
            error("Frame should be ±1, or ±2")
    end
    frabs = abs(frame)
    z = 3-frabs
    if(frame<0)
        codon = reverse(comp(codon))
    end
    down = codon[1+frabs:3].*vec(join.(product(repeated(nucs,frabs)...)))
    up = vec(join.(product(repeated(nucs,z)...))).*codon[1:frabs]
    return permutedims(hcat(collect.(product(up,down))...))
end


"""
`overlapXnucs(nucleotide.string,offset,is.complementary=true/false)`

For a given nucleotide sequence (length = x), find all overlapping nucleotide sequences of the same length (x), with an offset in the range (1-x,x-1)
...

# Examples
```
julia> frameXcodons("ATG",1)
4-element Vector{String}:
"TGA"
"TGT"
"TGG"
"TGC"

julia> overlapXnucs("ATG",-1)
4-element Vector{String}:
 "AAT"
 "TAT"
 "GAT"
 "CAT"
``` 
"""
function overlapXnucs(nstr::String,frame::Int, iscomp=false)
    slen = length(nstr)
    frabs = abs(frame)
    if(frabs>slen)
            error("offset ≥ nucleotide length")
    end
    if(frame == 0)
        error("offset is zero")
    end

    if(iscomp)
        nstr = reverse(comp(nstr))
    end

    if(frame>0)
        return nstr[1+frabs:slen].*vec(join.(product(repeated(nucs,frabs)...)))
    else
        z = slen-frabs
        return vec(join.(product(repeated(nucs,frabs)...))).*nstr[1:z]
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



nsub, nucsmbt = readdlm(joinpath(Base.source_dir(),"dmel_mutbias.txt"),'\t',header = true);

mutrate = 7.8e-9

mutratebac = 2e-8;

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
    if(isempty(s1) || isempty(s2))
        @warn "non-nucleotide string input, returning 0"
        return 0
    end
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

function mprob2(s1,s2)
    if(isempty(s1) || isempty(s2))
        @warn "non-nucleotide string input, returning 0"
        return 0
    end
    l = length(s1)
    if l != length(s2)
        error("length mismatch")
    end
    z = collect(s1) .!= collect(s2)
    p = 1
    if any(z)
        for i = (1:l)[z]
            p = p*nsub[nucnames[s1[i]], nucnames[s2[i]]]
        end
    end
    nmut = sum(z)
    unmut = join([s1[x] for x in findall(.!z)]);
    return (1-µA(mutrate,nsub))^numAT(unmut)*
    (1-µG(mutrate,nsub))^(length(unmut)-numAT(unmut))* p*(2*mutrate)^nmut
end

function nucprob(x,gccontent)
    if(isempty(x) || occursin(r"[^ATGC]",x))
        @warn "non-nucleotide string input, returning 0"
        return 0
    end
    S = 0.5*gccontent
    W = 0.5 - S
    nw = count(r"[AT]", x)
    return W^nw * S^(count(r"[GC]", x))
end

function numAT(x)
    return count(r"[AT]", x)
end

function featuregain(set1,set2,gccontent)
    if(isempty(set1) || isempty(set2))
        @warn "empty feature sets, returning 0"
        return 0
    end
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

nprob3 = (x,t3) -> t3[t3[:,1] .== x,2][1]
function nprob6(x,t6)
    if(isnothing(x))
        return 0
    else
        return  t6[t6[:,1] .== x,2][1]
    end
end

function featuregain3(set1,set2,t3)
    if(isempty(set1) || isempty(set2))
        @warn "empty feature sets, returning 0"
        return 0
    end
    v1 = sum([nprob3(s1,t3)*sum([mprob(s1,s2) for s2 in set2])
                for s1 in set1])
    return v1[1]
end

function featurestay3(set1,t3)
    v =  sum(
        [nprob3(s1,t3)*(
            sum([mprob(s1,s2) for s2 in set1]) +
            (1-µA(mutrate,nsub))^numAT(s1)*
            (1-µG(mutrate,nsub))^(3-numAT(s1))
            ) for s1 in set1]
        )
    return v[1]
end

function featuregain6(set1,set2,t6)
    if(isempty(set1) || isempty(set2))
        @warn "empty feature sets, returning 0"
        return 0
    end
    v1 = sum([nprob6(s1,t6)*sum([mprob(s1,s2) for s2 in set2])
                for s1 in set1])
    return v1[1]
end

function featurestay6(set1,t6)
    v =  sum(
        [nprob6(s1,t6)*(
            sum([mprob(s1,s2) for s2 in set1]) +
            (1-µA(mutrate,nsub))^numAT(s1)*
            (1-µG(mutrate,nsub))^(6-numAT(s1))
            ) for s1 in set1]
        )
    return v[1]
end
using DelimitedFiles

gencode = readdlm("/home/bravi/Documents/EBB/std_genetic_code.txt", '\t', AbstractString, '\n')
nucsv = ['A','T','G','C']

function transnucs(codon)
        return gencode[gencode[:,1].==codon,3]
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

muteff3 = permutedims(hcat([
                [string(x[1],x[2],i) for i in nucsv[nucsv.!=x[3]] ]
                for x in gencode[:,1]
        ]...))

muteff2 = permutedims(hcat([
                [string(x[1],i,x[3]) for i in nucsv[nucsv.!=x[2]] ]
                for x in gencode[:,1]
        ]...))

muteff1 = permutedims(hcat([
                [string(i,x[2],x[3]) for i in nucsv[nucsv.!=x[1]] ]
                for x in gencode[:,1]
        ]...))

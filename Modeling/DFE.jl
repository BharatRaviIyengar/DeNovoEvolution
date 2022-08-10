cd(Base.source_dir())
include("nucleotidefuncts.jl")

function frameXcodons(codon,frame)
        if(frame ∉ [1,2])
                error("Frame should be 1, or 2")
        end

        z = 3-frame
        u = codon[1:z].*vec(join.(product(repeated(nucs,frame)...)))
        d = vec(join.(product(repeated(nucs,z)...))).*codon[1:frame]
        return vcat(d,u)
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

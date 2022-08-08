include("nucleotidefuncts.jl")


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

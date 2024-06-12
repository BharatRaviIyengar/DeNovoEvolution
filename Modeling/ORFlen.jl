using Measures
using Plots
using Printf
using Distributed
using FileIO, JLD2

cm2pt = (cm) -> 28.3465*cm
figdir = joinpath(Base.source_dir(),"../Manuscripts/M3/Figures/");
colors = ["#FFCC00","#5599FF","#D40000","#754473","#000000"];
lstyles = [:solid,:dash,:dot]

default(linecolor = :black, linewidth = 2, tickfont = font(10,"Helvetica"), 
guidefont = font(13,"Helvetica"),framestyle = :box, legend = false);

withindels = false;

if withindels
    nprcs = 15
    fsuffix = "wtIndels"
else
    nprcs = 4;
    fsuffix = "noIndels"
end

addprocs(nprcs)

@everywhere cd(Base.source_dir())

@everywhere organism = "dmel"
# @everywhere organism = "scer"

@everywhere include("ORFlen_defs.jl")

gcrange = [0.3:0.1:0.6;];

@everywhere ncodons = [3:900;];

function markovcalc(gens,transmat,orfvalues,traj="s")
    xmt = transmat^gens;
    if traj=="e"
        xinit = orfvalues[xProb,:]'.*LowerTriangular(ones(size(transmat)));
        setindex!.(Ref(xinit), 0.0, eachindex(ncodons), eachindex(ncodons));
    elseif traj=="t"
        xinit = UpperTriangular(ones(size(transmat))).*orfvalues[xProb,:]';
        setindex!.(Ref(xinit), 0.0, eachindex(ncodons), eachindex(ncodons));
    else
        xinit = orfvalues[xProb,:].*Matrix(1.0*I,size(transmat));
    end
    return xinit*xmt
end

similarity90 = [poizero(mutrate,x) for x in Int.(ceil.(ncodons.*0.9))];

atgvals = hcat([ATGprobs(g) for g in gcrange]...);
stopvals = hcat([stopprobs(g) for g in gcrange]...);

orfvals = zeros(4,length(ncodons),length(gcrange));
tmat = zeros(length(ncodons),length(ncodons),length(gcrange));

# m3 = zeros(length(ncodons),length(ncodons),length(gcrange));
# m5 = zeros(length(ncodons),length(ncodons),length(gcrange));

atgvalsX = ATGprobsX(trimerfreq);
stopvalsX = stopprobsX(trimerfreq);
orfvalsX = hcat([orfprobs(atgvalsX,stopvalsX,k) for k in ncodons]...);

tmatX = pmap(j -> [tprobs(i,j,atgvalsX,stopvalsX,withindels) for i in ncodons], ncodons, on_error=identity);
tmatX = hcat(tmatX...);

m3x = pmap(j -> [T3(i,j,atgvalsX,stopvalsX,withindels) for i in ncodons], ncodons, on_error=identity);
m3x = hcat(m3x...);
m5x = pmap(j -> [T5(i,j,atgvalsX,stopvalsX) for i in ncodons], ncodons, on_error=identity);
m5x = hcat(m5x...);

m5g3x = log.(m5x./m3x);

for g in eachindex(gcrange)
    orfvals[:,:,g] = hcat([orfprobs(atgvals[:,g],stopvals[:,g],k) for k in ncodons]...);

    qp = pmap(j->[tprobs(i,j,atgvals[:,g],stopvals[:,g],withindels) for i in ncodons],ncodons; on_error=identity);

    tmat[:,:,g] = hcat(qp...);

    # qp = pmap(j->[T5(i,j,atgvals[:,g],stopvals[:,g]) for i in ncodons],ncodons; on_error=identity);

    # m5[:,:,g] = hcat(qp...);

    # qp = pmap(j->[T3(i,j,stopvals[:,g],withindels) for i in ncodons],ncodons; on_error=identity);

    # m3[:,:,g] = hcat(qp...);
    println(g)
end

# FileIO.save("DATA_ORFlen_"*fsuffix*".jld2","tmatX",tmatX,"tmat",tmat)

rmprocs(2:nprcs)

# mostlikelystart = zeros(Int32,length(ncodons),length(gcrange)+1);
# for g in eachindex(gcrange)
#     egt = eigen(tmat[:,:,g]);
#     vecs = orfvals[xProb,:,g].*egt.vectors; 
#     vecs2 = abs.(vecs)
#     mostlikelystart[:,g] = [ncodons[findlast(vecs2[:,x] .== maximum(vecs2[:,x]))] for x in eachindex(ncodons)];
# end

# egt = eigen(tmatX);
# vecs = orfvalsX[xProb,:].*egt.vectors; 
# vecs2 = abs.(vecs)
# mostlikelystart[:,end] = [ncodons[findlast(vecs2[:,x] .== maximum(vecs2[:,x]))] for x in eachindex(ncodons)];

gens = 2 .^[0:Int(floor(log2(100/mutrate)));]
(EgT, EgS, TgS) = [zeros(Int64,length(gens),length(gcrange)+1) for i=1:3];

(ev, tv, sv, Eop, Top, Sop) = [zeros(length(ncodons),length(gens),length(gcrange)+1) for i=1:6];

for r in eachindex(gens)
    for g in eachindex(gcrange)
        tmg = tmat[:,:,g]^gens[r]
        sameQ = markovcalc(1,tmg, orfvals[:,:,g],"s");
        extnQ = markovcalc(1,tmg, orfvals[:,:,g],"e");
        trncQ = markovcalc(1,tmg, orfvals[:,:,g],"t");
        (sv[:,r,g], ev[:,r,g], tv[:,r,g]) = [diag(x) for x in [sameQ, extnQ, trncQ]];
        # println([gcrange[g],gens[r]])
        if(any(ev[:,r,g] .> tv[:,r,g]))
            EgT[r,g] = ncodons[findfirst(ev[:,r,g] .> tv[:,r,g])]
        end    
        if(any(ev[:,r,g] .> sv[:,r,g]))
            EgS[r,g] = ncodons[findfirst(ev[:,r,g] .> sv[:,r,g])]
        end 
        if(any(tv[:,r,g] .> sv[:,r,g]))
            TgS[r,g] = ncodons[findfirst(tv[:,r,g] .> sv[:,r,g])]
        end 
        Eop[:,r,g] = [sum(tmg[x,x+1:end]) for x in eachindex(ncodons)]
        Top[:,r,g] = [sum(tmg[x,1:x-1]) for x in eachindex(ncodons)]
        Sop[:,r,g] = diag(tmg)
    end
    tmg = tmatX^gens[r]
    sameQ = markovcalc(1,tmg,orfvalsX,"s");
    extnQ = markovcalc(1,tmg,orfvalsX,"e");
    trncQ = markovcalc(1,tmg,orfvalsX,"t");
    (sv[:,r,end], ev[:,r,end], tv[:,r,end]) = [diag(x) for x in [sameQ, extnQ, trncQ]];
    if(any(ev[:,r,end] .> tv[:,r,end]))
        EgT[r,end] = ncodons[findfirst(ev[:,r,end] .> tv[:,r,end])]
    end    
    if(any(ev[:,r,end] .> sv[:,r,end]))
        EgS[r,end] = ncodons[findfirst(ev[:,r,end] .> sv[:,r,end])]
    end
    if(any(tv[:,r,end] .> sv[:,r,end]))
        TgS[r,end] = ncodons[findfirst(tv[:,r,end] .> sv[:,r,end])]
    end
    Eop[:,r,end] = [sum(tmg[x,x+1:end]) for x in eachindex(ncodons)]
    Top[:,r,end] = [sum(tmg[x,1:x-1]) for x in eachindex(ncodons)]
    Sop[:,r,end] = diag(tmg)
    println(r)
end


# fac = Int(5*10^(Int(floor(-log10(mutrate))) -2))
# eggen = findfirst(gens .>= fac)

# sameQ = markovcalc(fac,tmatX, orfvalsX,"s");
# extnQ = markovcalc(fac,tmatX, orfvalsX,"e");
# trncQ = markovcalc(fac,tmatX, orfvalsX,"t");
# (sv, ev, tv) = [diag(x) for x in [sameQ, extnQ, trncQ]];


# Plots #

# pEV = plot( 
#     xlabel = "Final ORF length (codons)",
#     ylabel = "Initial ORF length (codons)",
#     size = (width = cm2pt(11), height = cm2pt(10))
#     );

# for q in 1:5
#     plot!(pEV,ncodons,mostlikelystart[:,q],
#         linecolor = colors[q]
#     );
# end

mx = round(maximum(m5g3x[.!(isnan.(m5g3x))]),digits=2);
mn = round(minimum(m5g3x[.!(isnan.(m5g3x))]),digits=2);
rrx = mx - mn;
cr = cgrad(["#d40000",:black,"#66a3ff"], [-mn/rrx, 1])
setindex!.(Ref(m5g3x), NaN, eachindex(ncodons), eachindex(ncodons));

plt5g3x = heatmap(m5g3x, 
    c = cr,
    clim = (mn, mx),
    cb = :top,
    xrotation = 90,
    xticks = [100:100:900;],
    yticks = [100:100:900;],
    ylabel = "Initial ORF length",
    xlabel = "Final ORF length",
    size = (width = cm2pt(11.5), height = cm2pt(11)),
);

savefig(plt5g3x,figdir*"5g3_"*organism*".pdf")


TgE_gORF = hcat([ncodons[[findfirst(Top[:,x,g] .> Eop[:,x,g]) for x in eachindex(gens)]] for g =1:5]...);

pltTgE = plot( 
    xlabel = "Generations",
    ylabel = "Minimum ORF length\n(codons)",
    title = "P(T) > P(E)",
    xaxis = :log10,
    xticks = 10 .^[0:2:10;],
    yticks = [3:3:21;],
    size = (width = cm2pt(11.5), height = cm2pt(11))
    );

for q in 1:5
    plot!(pltTgE,gens,TgE_gORF[:,q],
        linecolor = colors[q]
    );
end

plot!(pltTgE,repeat([1/mutrate],2),vcat(ylims(pltTgE)...),
    linewidth = 0.8,
    linecolor = :grey,
    linestyle = :dash
)

savefig(pltTgE,figdir*"TgEvGen_fromORF_"*organism*"_"*fsuffix*".pdf")

pltEgT = plot( 
    xlabel = "Generations",
    ylabel = "Minimum ORF length\n(codons)",
    title = "P(E) > P(T)",
    xaxis = :log10,
    xticks = 10 .^[0:2:10;],
    size = (width = cm2pt(11.5), height = cm2pt(11))
    );

for q in 1:5
    plot!(pltEgT,gens,EgT[:,q],
        linecolor = colors[q]
    );
end


plot!(pltEgT,repeat([1/mutrate],2),vcat(ylims(pltEgT)...),
    linewidth = 0.8,
    linecolor = :grey,
    linestyle = :dash
)


savefig(pltEgT,figdir*"EgTvGen_"*organism*"_"*fsuffix*".pdf")

# SgE_gORF = hcat([ncodons[[findfirst(Sop[:,x,g] .> Eop[:,x,g]) for x in eachindex(gens)]] for g =1:5]...);

# pltSgE = plot( 
#     xlabel = "Generations",
#     ylabel = "Minimum ORF length\n(codons)",
#     title = "P(S) > P(E)",
#     xaxis = :log10,
#     xticks = 10 .^[0:2:10;],
#     size = (width = cm2pt(11.5), height = cm2pt(11))
#     );

# for q in 1:5
#     plot!(pltSgE,gens,SgE_gORF[:,q],
#         linecolor = colors[q]
#     );
# end

# savefig(pltTgE,figdir*"SgEvGen_fromORF_"*organism*".pdf")


mESGx = [findfirst(EgS[:,x] .>0) for x = 1:5];
mTSGx = [findfirst(TgS[:,x] .>0) for x = 1:5];

pltEgS = plot( 
    xlabel = "Generations",
    ylabel = "Minimum ORF length\n(codons)",
    title = "P(E) > P(C)",
    xaxis = :log10,
    xticks = 10 .^[floor(log10(gens[minimum(mESGx)])):floor(log10(maximum(gens)));],
    size = (width = cm2pt(11.5), height = cm2pt(11))
    );

for q in 1:5
    plot!(pltEgS,gens[mESGx[q]:end],EgS[mESGx[q]:end,q],
        linecolor = colors[q]
    );
end


plot!(pltEgS,repeat([1/mutrate],2),vcat(ylims(pltEgS)...),
    linewidth = 0.8,
    linecolor = :grey,
    linestyle = :dash
)


savefig(pltEgS,figdir*"EgSvGen_"*organism*"_"*fsuffix*".pdf")


pltTgS = plot( 
    xlabel = "Generations",
    ylabel = "Minimum ORF length\n(codons)",
    title = "P(T) > P(C)",
    xaxis = :log10,
    xticks = 10 .^[floor(log10(gens[minimum(mTSGx)])):floor(log10(maximum(gens)));],
    size = (width = cm2pt(11.5), height = cm2pt(11))
    );

for q in 1:5
    plot!(pltTgS,gens[mTSGx[q]:end],TgS[mTSGx[q]:end,q],
        linecolor = colors[q]
    );
end


plot!(pltTgS,repeat([1/mutrate],2),vcat(ylims(pltTgS)...),
    linewidth = 0.8,
    linecolor = :grey,
    linestyle = :dash
)


savefig(pltTgS,figdir*"TgSvGen_"*organism*"_"*fsuffix*".pdf")

# plt_eg = plot(ncodons[1:300],log10.(sv)[1:300], linecolor = :black, label = "Constant", xlabel = "ORF length (codons)", ylabel = "log10(Probability)")
# plot!(plt_eg, ncodons[1:300],log10.(ev)[1:300], linecolor = :red, label = "Extended");
# plot!(plt_eg,ncodons[1:300],log10.(tv)[1:300], linecolor = :blue, label = "Truncated");

# plot!(plt_eg, legend = :bottom, legendfont = font(11,"Helvetica"),
#     size = (width = cm2pt(13), height = cm2pt(11))
#     )

# savefig(plt_eg,figdir*organism*"_Prob_"* (@sprintf "%0.0e" fac)*"gens.pdf")

# allplots = []

# for i in Int.(floor.([1:length(gens)/4:length(gens);]))
#     plt_egOP = plot(ncodons[2:300],log10.(Eop[2:300,i,5]), linecolor = :red, title = string(gens[i])*" generations", label ="Extended");

#     plot!(plt_egOP,ncodons[2:300],log10.(Top[2:300,i,5]), linecolor = :blue, label ="Truncated")

#     if i==27
#         plot!(plt_egOP, legend = :bottom, legendfont = font(11,"Helvetica"),
#         size = (width = cm2pt(13), height = cm2pt(11))
#         )
#     end
#     push!(allplots,plt_egOP)
# end

# apx = plot(allplots..., layout = (2,2),size = (width = cm2pt(25), height = cm2pt(20)));

# savefig(apx,figdir*"Supp/TvE_fromORF"*organism*".pdf")


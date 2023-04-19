using Plots
using Measures
using DelimitedFiles
cd(Base.source_dir())

cm2pt = (cm) -> 28.3465*cm
figdir = joinpath(Base.source_dir(),"../../Manuscripts/Figures/M2_main/pdf/");
colors = ["#FFCC00","#5599FF","#D40000","#754473","#000000"];
lstyles = [:solid,:dash,:dot]

default(linecolor = :black, linewidth = 2, tickfont = font(10,"Helvetica"), 
guidefont = font(13,"Helvetica"),framestyle = :box, legend = false);

(framedata, header) =readdlm(joinpath(Base.source_dir(),"framecounts.txt"), '\t', header = true);

f100 = framedata[framedata[:,1] .== 100, [2,8,5,7]];

f100[:,3] = f100[:,3]/1000; 

pFcounts = plot(
    ylabel = "Value", size = (width = cm2pt(13), height = cm2pt(10))
);

bwf = 0.2
for j = 1:3
    bar!(pFcounts,[j+x*bwf for x in 1:3],f100[:,j+1],
        fill = colors[1:3],
        bar_width = bwf,
        linecolor = nothing
    )
end

xticks!(pFcounts,[1:3;] .+ 2*bwf, ["%ORFs","Total Length","Avg. Length"]);

savefig(pFcounts, figdir*"Framecounts_dmel.pdf")


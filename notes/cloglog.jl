using DataFrames, GLM, MixedModels
const df = readtable("artSim.csv", makefactors = true);
## compute proportions and totals
df[:Tot] = float.(df[:Dead]+df[:Alive]);
df[:prop] = df[:Dead] ./ df[:Tot];
## for CPU timing:
### Pkg.clone("https://github.com/schmrlng/CPUTime.jl.git")
## fit without timing, for compilation
## (Trt+0)/x formula notation doesn't work, need to expand it ourselves
g3 = fit!(glmm(@formula(prop ~ 0 + Trt + Trt&x + (1 + x|Rep)), df,
    Binomial(), CloglogLink(), wt=Array(df[:Tot])), fast = true);
e1 = @elapsed fit!(glmm(@formula(prop ~ 0 + Trt + Trt&x + (1 + x|Rep)), df,
    Binomial(), CloglogLink(), wt=Array(df[:Tot])), fast = true)
## non-fast (nAGQ>0)
e2 = @elapsed g4 = fit!(glmm(@formula(prop ~ 0 + Trt + Trt&x + (1 + x|Rep)), df,
    Binomial(), CloglogLink(), wt=Array(df[:Tot])))
## store timing info
open("julia_cloglog_timings.txt","w") do f
    write(f,string(e1)," ",string(e2),"\n")
end
## store fixed-effect results
ct = coeftable(g4)
outdf = DataFrame(variable = ct.rownms,
                 Estimate = ct.cols[1],
               StdError = ct.cols[2],
                z_val = ct.cols[3])
writetable("julia_cloglog_coefs.csv",outdf)
ct = coeftable(g3)
outdf = DataFrame(variable = ct.rownms,
                 Estimate = ct.cols[1],
               StdError = ct.cols[2],
                z_val = ct.cols[3])
writetable("julia_cloglog_coefs_nagq0.csv",outdf)

using DataFrames
using GLM
using MixedModels
df = readtable("artSim.csv");
## compute proportions and totals
df[:Tot] = df[:Dead]+df[:Alive];
df[:prop] = df[:Dead] ./ df[:Tot];
## convert Treatment to categorical
## https://stackoverflow.com/questions/43708635/indicator-matrix-for-categorical-data-in-glm-jl-with-dataframes-jl
pool!(df, [:Trt])
## basic (cloglog) GLM fit
e1 = @elapsed g1 = glm(@formula(prop ~ 0 + Trt), df, Binomial(), CloglogLink(), wts=float.(Array(df[:Tot])))
## GLMM with logit link
## note (Trt+0)/x formula notation doesn't work, need to expand it ourselves
e2 = @elapsed g2 = fit!(glmm(@formula(prop ~ 0+Trt+Trt&x + (x | Rep)), df, Binomial(),
           LogitLink(),
           wt = float.(Array(df[:Tot]))))
## GLMM with cloglog link
e3 = @elapsed g3 = fit!(glmm(@formula(prop ~ 0+Trt+Trt&x + (x | Rep)), df, Binomial(),
           CloglogLink(),
           wt = float.(Array(df[:Tot]))))
## store timing info
open("julia_cloglog_timings.txt","w") do f
    write(f,string(e1)," ",string(e2)," ",string(e3),"\n")
end
## store fixed-effect results
ct = coeftable(g3)
outdf = DataFrame(variable = ct.rownms,
                 Estimate = ct.cols[1],
               StdError = ct.cols[2],
                z_val = ct.cols[3])
writetable("julia_cloglog_coefs.csv",outdf)

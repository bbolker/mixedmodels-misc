using CSV, DataFrames, MixedModels;

"""
    outfun(m, outfn="output.csv")

output the coefficient table of a fitted model to a file
"""
outfun = function(m, outfn="output.csv")
    ct = coeftable(m)
    coef_df = DataFrame(ct.cols);
    rename!(coef_df, ct.colnms, makeunique = true)
    ## names!(coef_df, [ Symbol(nm) for nm in ct.colnms ]);
    coef_df[!, :term] = ct.rownms;
    CSV.write(outfn, coef_df);
end

data = CSV.read("gdata.csv", DataFrame);
data.Intervention = CategoricalArray(data.Intervention);
levels!(data.Intervention, ["Pre", "Post"]);

## DRY?
data_disagg = CSV.read("gdata_disagg.csv", DataFrame);
data_disagg.Intervention = CategoricalArray(data_disagg.Intervention);
levels!(data_disagg.Intervention, ["Pre", "Post"]);

## fit once for compilation
fit(MixedModels.GeneralizedLinearMixedModel,
               @formula(prop_succ~Intervention*Sex + (1|Subject) + (1|Image)),
    data, Binomial(), wts = data.tot);

@time(f1 = fit(MixedModels.GeneralizedLinearMixedModel,
               @formula(prop_succ~Intervention*Sex + (1|Subject) + (1|Image)),
               data, Binomial(), wts = data.tot));  ## 0.777
print(f1);
outfun(f1,"julia_agg.csv");

@time(f1_disagg = fit(MixedModels.GeneralizedLinearMixedModel, @formula(response~Intervention*Sex + (1|Subject) + (1|Image)),
                      data_disagg, Bernoulli()));
print(f1_disagg);
outfun(f1_disagg,"julia_disagg.csv");
## 73 seconds

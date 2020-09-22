using CSV, DataFrames, MixedModels;
import Pkg;
Pkg.status("MixedModels") ## now using devel version

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

cbpp = MixedModels.dataset(:cbpp);
describe(cbpp)

## fit once for compilation
fit(MixedModels.GeneralizedLinearMixedModel,
               @formula(incid~period + (1|herd)),
    cbpp, Binomial());

@time(f1 = fit(MixedModels.GeneralizedLinearMixedModel,
               @formula(prop~period 
               data, Binomial(), wts = data.tot));  ## 0.777
print(f1);
outfun(f1,"julia_agg.csv");

@time(f1_disagg = fit(MixedModels.GeneralizedLinearMixedModel, @formula(response~Intervention*Sex + (1|Subject) + (1|Image)),
                      data_disagg, Bernoulli()));
print(f1_disagg);
outfun(f1_disagg,"julia_disagg.csv");
## 73 seconds

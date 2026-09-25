"""statsmodels lane: dose-response of RBP count on genome length and host genus."""
import json
import pandas as pd
import statsmodels.formula.api as smf

def main():
    df = pd.read_csv("data/processed/phage_rbp_table.csv")
    df["host_genus"] = df.host_species.astype(str).str.split().str[0]
    df["log_len"] = (df.genome_length).map(lambda x: __import__("math").log10(x))
    m = smf.ols("n_rbps ~ log_len + C(host_genus)", data=df).fit()
    out = {"n": int(m.nobs), "r_squared": round(m.rsquared, 4),
           "log_len_coef": round(m.params.get("log_len", float("nan")), 3),
           "log_len_pvalue": float(m.pvalues.get("log_len", float("nan")))}
    json.dump(out, open("results/statsmodels_rbp_dose_response.json", "w"), indent=2)
    m.params.to_csv("results/statsmodels_rbp_dose_response_coefs.csv")
    print(out["n"], out["r_squared"], out["log_len_pvalue"])

if __name__ == "__main__":
    main()

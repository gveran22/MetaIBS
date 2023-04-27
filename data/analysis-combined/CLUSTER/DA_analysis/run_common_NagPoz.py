import pandas as pd
import anndata as ad
import numpy as np

import sccoda.util.comp_ana as mod

data_path = "/home/CLUSTER/DA_analysis/input"  # CHANGE THIS DIRECTORY ON YOUR COMPUTER
save_path = "/home/CLUSTER/DA_analysis/output"  # CHANGE THIS DIRECTORY ON YOUR COMPUTER

tax_levels = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
meta_col = ['host_disease', 'host_subtype', 'sequencing_tech', 'variable_region', 'sample_type', 'Collection',
            'host_sex', 'host_ID', 'host_age', 'host_bmi', 'collection_date',
            'Bristol', 'IBS_SSS', 'sex_flexibl', 'host_psy', 'ID_flexibl',
            'sample_storage_duration', 'sequencing_run', 'extraction_plate', 'author']


def read_shared_ASVs(a1, a2, data_path):
    raw = pd.read_csv(f"{data_path}/commonASV_{a1}-{a2}.csv", index_col=0)
    raw["Type"] = raw.loc[:, tax_levels].fillna("_").apply('*'.join, axis=1)

    raw_a1 = raw[raw["author"] == a1]
    counts_a1 = raw_a1.pivot("Sample", "OTU", "Abundance")
    meta_a1 = raw_a1.groupby("Sample").agg(dict([(x, "first") for x in meta_col]))
    tax_a1 = raw_a1.groupby("OTU").agg(dict([(x, "first") for x in tax_levels + ["Type"]]))
    data_a1 = ad.AnnData(X=counts_a1, obs=meta_a1, var=tax_a1)

    raw_a2 = raw[raw["author"] == a2]
    counts_a2 = raw_a2.pivot("Sample", "OTU", "Abundance")
    meta_a2 = raw_a2.groupby("Sample").agg(dict([(x, "first") for x in meta_col]))
    tax_a2 = raw_a2.groupby("OTU").agg(dict([(x, "first") for x in tax_levels + ["Type"]]))
    data_a2 = ad.AnnData(X=counts_a2, obs=meta_a2, var=tax_a2)

    if a2 == "Pozuelo":
        data_a2 = data_a2[data_a2.obs["Collection"] == "1st"]
        zero_sum = np.sum(data_a2.X, axis=0) > 0
        data_a2 = data_a2[:, zero_sum]

    counts_both = raw.pivot("Sample", "OTU", "Abundance")
    meta_both = raw.groupby("Sample").agg(dict([(x, "first") for x in meta_col]))
    tax_both = raw.groupby("OTU").agg(dict([(x, "first") for x in tax_levels + ["Type"]]))
    data_both = ad.AnnData(X=counts_both, obs=meta_both, var=tax_both)

    if a2 == "Pozuelo":
        data_both = data_both[data_both.obs["Collection"] == "1st"]

    return data_a1, data_a2, data_both


def run_model_shared(author1, author2, model, data_path, save_path, alpha=0.2, mode="all", total_scale=None):
    data_a1, data_a2, data_both = read_shared_ASVs(author1, author2, data_path)

    def model_run(d, m, alpha):

        if total_scale is not None:
            d.X = np.round(d.X/np.sum(d.X, axis=1, keepdims=True)*total_scale, 0).astype(int)

        if m == "sccoda":
            references = {
                "ASV": "Bacteria*Proteobacteria*Gammaproteobacteria*Burkholderiales*Sutterellaceae*Parasutterella*excrementihominis",
                "Genus": "Bacteria*Proteobacteria*Gammaproteobacteria*Burkholderiales*Sutterellaceae*Parasutterella",
                "Family": "Bacteria*Proteobacteria*Gammaproteobacteria*Burkholderiales*Sutterellaceae",
                "Order": "Bacteria*Proteobacteria*Gammaproteobacteria*Burkholderiales",
                "Class": "Bacteria*Proteobacteria*Gammaproteobacteria",
                "Phylum": "Bacteria*Proteobacteria",
            }

            ref = d.var[d.var["Type"] == references["ASV"]].index[0]
            print(ref)

            model = mod.CompositionalAnalysis(
                data=d,
                formula="C(host_disease, Treatment('Healthy'))",
                reference_cell_type=ref
            )
            result = model.sample_hmc(num_results=20000, num_burnin=5000)
            _, effect_df = result.summary_prepare(est_fdr=alpha)

            out = effect_df
            out["is_da"] = [True if x != 0 else False for x in out["Final Parameter"]]
            out.index = d.var.index

        else:
            return None

        res = d.var.merge(out, left_index=True, right_index=True)
        return res

    out = []

    if mode == "all" or mode == "a1":
        print(f"Running {author1}")
        if len(pd.unique(data_a1.obs["host_disease"])) < 2:
            print("Common ASVs were not found in both groups! Skipping...")
            da_a1 = pd.DataFrame()
        else:
            da_a1 = model_run(data_a1, model, alpha)
        da_a1.to_csv(save_path + f"/shared_{author1}{author2}_{model}_{total_scale}_a1.csv")
        out.append(da_a1)

    if mode == "all" or mode == "a2":
        print(f"Running {author2}")
        if len(pd.unique(data_a2.obs["host_disease"])) < 2:
            print("Common ASVs were not found in both groups! Skipping...")
            da_a2 = pd.DataFrame()
        else:
            da_a2 = model_run(data_a2, model, alpha)
        da_a2.to_csv(save_path + f"/shared_{author1}{author2}_{model}_{total_scale}_a2.csv")
        out.append(da_a2)

    if mode == "all" or mode == "combined":
        print(f"Running Combined data")
        da_both = model_run(data_both, model, alpha)
        da_both.to_csv(save_path + f"/shared_{author1}{author2}_{model}_{total_scale}_combined.csv")
        out.append(da_both)

    return out


author1 = "Nagel"
author2 = "Pozuelo"
sccoda_res = run_model_shared(author1, author2, "sccoda", data_path=data_path, save_path=save_path, total_scale=5161.0, mode="all")

import anndata as ad
import pandas as pd
import sccoda.util.comp_ana as mod

data_dir = "/home/CLUSTER/DA_analysis"  # CHANGE THIS DIRECTORY ON YOUR COMPUTER

tax_levels = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus"]


def agg_ibs_data(author, level, data_dir):

    # read data
    subdir_name = f"{data_dir}/input"
    file_name = subdir_name + f"/{author.lower()}_{level.lower()}-agg.csv"
    raw_data = pd.read_csv(file_name, index_col=0)

    # get taxonomic levels in the data
    tl = [x for x in tax_levels[:tax_levels.index(level) + 1]]

    # extract counts
    count_data = raw_data.pivot(index="Sample", columns=tl, values="Abundance")
    # get taxonomic tree (for data.var)
    tax_info = pd.DataFrame(index=count_data.columns).reset_index()
    tax_index = tax_info.apply('*'.join, axis=1)
    tax_index = [s.replace("(", "") .replace(")", "") for s in tax_index]
    tax_info.index = tax_index

    count_data.columns = tax_index

    # get metadata
    metadata_cols = raw_data.columns.drop(["Sample", "Abundance"] + tax_levels, errors="ignore")
    metadata = raw_data.groupby("Sample").agg(dict([(x, "first") for x in metadata_cols]))

    ret = ad.AnnData(X=count_data, obs=metadata, var=tax_info)
    return ret


def run_sccoda(author, level, data_dir, add=None, fdr_level=0.1):
    references = {
        "Genus": "Bacteria*Proteobacteria*Gammaproteobacteria*Burkholderiales*Sutterellaceae*Parasutterella",
        "Family": "Bacteria*Proteobacteria*Gammaproteobacteria*Burkholderiales*Sutterellaceae",
        "Order": "Bacteria*Proteobacteria*Gammaproteobacteria*Burkholderiales",
        "Class": "Bacteria*Proteobacteria*Gammaproteobacteria",
        "Phylum": "Bacteria*Proteobacteria",
    }

    data = agg_ibs_data(author, level, data_dir)
    if add is not None:
        data = data[data.obs[add[0]] == add[1]]

    model = mod.CompositionalAnalysis(
        data=data,
        formula="C(host_disease, Treatment('Healthy'))",
        reference_cell_type=references[level]
    )
    result = model.sample_hmc()
    _, effect_df = result.summary_prepare(est_fdr = fdr_level)

    return effect_df


def one_author_new(author, levels, adds, model, data_dir, alpha=0.1, run_no=None):

    for l in levels[1:]:
        print(l)

        for a in adds:
            print(a)

            if model == "sccoda":
                out = run_sccoda(author, l, data_dir, fdr_level=alpha, add=a)
            else:
                raise ValueError("Invalid model name!")

            filename = f"{author.lower()}_{l.lower()}_"
            if a is not None:
                filename += f"{a[1]}_"

            if run_no:
                out.to_csv(f"{data_dir}/output/{filename}{model}_alpha_{alpha}_{run_no}.csv")
            else:
                out.to_csv(f"{data_dir}/output/{filename}{model}_alpha_{alpha}.csv")


model = "sccoda"
alpha = 0.2
run_no = 2

author = "AGP"
one_author_new(author, tax_levels, [None], model, data_dir, alpha, run_no=run_no)

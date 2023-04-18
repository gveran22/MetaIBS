import os
import pandas as pd
import anndata as ad
import numpy as np
import toytree as tt
import toyplot
import toyplot.color

import sccoda.util.comp_ana as mod
import sccoda.model.other_models as om

from rpy2.robjects import numpy2ri, pandas2ri
import rpy2.robjects as rp
numpy2ri.activate()
pandas2ri.activate()
import rpy2.robjects.packages as rpackages

r_home = "/Library/Frameworks/R.framework/Resources" # CHANGE THIS DIRECTORY ON YOUR COMPUTER
r_path = r"/Library/Frameworks/R.framework/Resources/bin" # CHANGE THIS DIRECTORY ON YOUR COMPUTER

os.environ["R_HOME"] = r_home
os.environ["PATH"] = r_path + ";" + os.environ["PATH"]

tax_levels = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus"]


def agg_ibs_data(author, level, data_dir):

    # read data
    subdir_name = [x for x in os.listdir(data_dir) if x.startswith(author+"-")][0]
    file_name = data_dir + subdir_name + f"/{author.lower()}_{level.lower()}-agg.csv"
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


def run_sccoda(author, level, data_dir, add=None, fdr_level=0.1,):
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


def run_ancombc_model(author, level, data_dir, add=None, alpha=0.05):
    data = agg_ibs_data(author, level, data_dir)
    if add is not None:
        data = data[data.obs[add[0]] == add[1]]

    data.X[data.X == 0] = 0.5
    ac = om.ANCOMBCModel(data.copy(), covariate_column="host_disease")

    abc_full_dict = {
        "p_adj_method": "fdr",
        "zero_cut": 1,
        "lib_cut": 0,
        "struc_zero": "TRUE",
        "neg_lb": "TRUE",
        "tol": 1e-5
    }

    out_ = rp.r(f"""
                    library(ANCOMBC)
                    library(phyloseq)

                    #prepare phyloseq data format

                    counts = {pandas2ri.py2rpy_pandasdataframe(pd.DataFrame(ac.y, columns=ac.var.index)).r_repr()}

                    sample = {pandas2ri.py2rpy_pandasdataframe(pd.DataFrame(ac.x,
                                                                            columns=[ac.covariate_column])).r_repr()}

                    cell_types = colnames(counts)

                    OTU = otu_table(t(counts), taxa_are_rows = TRUE)

                    #create phyloseq data object
                    data = phyloseq(OTU, sample_data(sample))

                    ancombc_out = ancombc(phyloseq = data,
                                          formula = "{ac.covariate_column}",
                                          p_adj_method = "{abc_full_dict["p_adj_method"]}",
                                          zero_cut = {abc_full_dict["zero_cut"]},
                                          lib_cut = {abc_full_dict["lib_cut"]},
                                          group = "{ac.covariate_column}",
                                          struc_zero = {abc_full_dict["struc_zero"]},
                                          neg_lb = {abc_full_dict["neg_lb"]},
                                          tol = {abc_full_dict["tol"]},
                                          max_iter = 100,
                                          conserve = TRUE,
                                          alpha = {alpha},
                                          global = FALSE
                                          )

                    out = ancombc_out
                    out
                    """)
    x = pd.DataFrame(np.array(out_[1]))
    x.columns = ["struc_zero_0", "struc_zero_1"]
    x.index = list(out_[1].names[0])
    df = pd.DataFrame(dict(zip(out_[6].names, [[y[0] for y in x] for x in out_[6]])))
    df.index = list(out_[1].names[0])

    df = df.merge(x, right_index=True, left_index=True)

    # df["zero_ind"] = zero_ind
    out = df
    out.index = data.var.index
    out["is_da"] = [True if x == 1 else False for x in out["diff_abn"]]

    return out


def run_linda_model(author, level, data_dir, add=None, alpha=0.05, formula="host_disease"):
    data = agg_ibs_data(author, level, data_dir)
    if add is not None:
        data = data[data.obs[add[0]] == add[1]]

    meta = pandas2ri.py2rpy(data.obs)
    otus = pandas2ri.py2rpy(pd.DataFrame(data.X.T))

    linda = rpackages.importr("LinDA")

    lo = linda.linda(otus, meta, formula=f"~{formula}", alpha=alpha,
                     prev_cut=0, lib_cut=1)

    out = pd.DataFrame(lo[2][0])
    out.index = data.var.index

    return out


def read_authors_results(authors, data_dir, method, adds=None, alpha=None, run_no=None):
    out_dict = {}

    for l in tax_levels[1:]:
        ll = []
        for a in authors:
            subdir_name = [x for x in os.listdir(data_dir) if x.startswith(a + "-")][0]

            for add in adds[a]:
                name = f"{a.lower()}_{l.lower()}"
                if add:
                    name += f"_{add}"
                name += f"_{method}"
                if alpha:
                    name += f"_alpha_{alpha}"
                if run_no:
                    name += f"_{run_no}"
                name += ".csv"

                for f in os.listdir(data_dir + subdir_name):

                    if f == name:
                        res_ = pd.read_csv(data_dir + subdir_name + "/" + f, index_col=0)
                        if method == "sccoda":
                            res_ = res_.reset_index()
                            res_["Is credible"] = (res_["Final Parameter"] != 0)
                        elif method == "ANCOMBC":
                            res_ = res_.reset_index().rename(columns={
                                "index": "Cell Type",
                                "is_da": "Is credible"
                            })
                        elif method == "LinDA":
                            res_ = res_.reset_index().rename(columns={
                                "index": "Cell Type",
                                "padj": "p_adj"
                            })
                            res_["Is credible"] = [False if x == 0 else True for x in res_["reject"]]
                        res_["author"] = a
                        if add in ["stool", "sigmoid"]:
                            res_["source"] = add
                        elif a == "Mars":
                            res_["source"] = "sigmoid"
                        else:
                            res_["source"] = "stool"

                        tax = get_phylo_levels(res_, l)
                        res_ = pd.merge(res_, tax, left_index=True, right_index=True)
                        ll.append(res_)
        out_dict[l] = pd.concat(ll)

    return out_dict


def get_significances(out_dict, method):
    res_sigs = {}
    for l in tax_levels[1:]:

        res = out_dict[l]

        if method == "sccoda" or method == "scCODA":
            res_sig_gr = res.groupby("Cell Type").agg({
                "Kingdom": "first",
                "Phylum": "first",
                "Class": "first",
                "Order": "first",
                "Family": "first",
                "Genus": "first",
                "Is credible": "sum",
                "Final Parameter": [lambda x: (x > 0).sum(), lambda x: (x < 0).sum()],
            },
            )
            res_sig_gr.columns = ["Increase" if j == "<lambda_0>" else
                                  "Decrease" if j == "<lambda_1>" else
                                  "count" if i == "Is credible" else
                                  i for i, j in res_sig_gr.columns]
        elif method == "LinDA":
            res["eff"] = res.apply(lambda row: row["stat"] if row["reject"]==1 else 0, axis=1)
            res_sig_gr = res.groupby("Cell Type").agg({
                "Kingdom": "first",
                "Phylum": "first",
                "Class": "first",
                "Order": "first",
                "Family": "first",
                "Genus": "first",
                "Is credible": "sum",
                "eff": [lambda x: (x > 0).sum(), lambda x: (x < 0).sum()],
            },
            )
            res_sig_gr.columns = ["Increase" if j == "<lambda_0>" else
                                  "Decrease" if j == "<lambda_1>" else
                                  "count" if i == "Is credible" else
                                  i for i, j in res_sig_gr.columns]
        else:
            res_sig_gr = res.groupby("Cell Type").agg({
                "Kingdom": "first",
                "Phylum": "first",
                "Class": "first",
                "Order": "first",
                "Family": "first",
                "Genus": "first",
                "Is credible": "sum"
            }).rename(columns={"Is credible": "count"})

        res_sig_gr["model"] = method

        res_sigs[l] = res_sig_gr

    return res_sigs


def get_phylo_levels(results, level, col="Cell Type"):

    max_level_id = tax_levels.index(level)+1
    cols = tax_levels

    tax_table = pd.DataFrame(columns=cols, index=np.arange(len(results)))
    for i in range(len(results)):
        char = results.loc[i, col]
        split = char.split(sep="*")
        for j in range(max_level_id):
            try:
                tax_table.iloc[i, j] = split[j]
            except IndexError:
                tax_table.iloc[i, j] = None

    return tax_table


def traverse(df_, a, i, innerl):
    if i+1 < df_.shape[1]:
        a_inner = pd.unique(df_.loc[np.where(df_.iloc[:, i] == a)].iloc[:, i+1])

        desc = []
        for b in a_inner:
            desc.append(traverse(df_, b, i+1, innerl))
        if innerl:
            il = a
        else:
            il = ""
        out = f"({','.join(desc)}){il}"
    else:
        out = a

    return out


def df2newick(df, inner_label=True, tax_lev=tax_levels):

    df_tax = df.loc[:, [x for x in tax_lev if x in df.columns]]

    alevel = pd.unique(df_tax.iloc[:, 0])
    strs = []
    for a in alevel:
        strs.append(traverse(df_tax, a, 0, inner_label))

    newick = f"({','.join(strs)});"
    return newick


def build_fancy_tree(data, edge_color_dict=None, other_col="lightblue", tax_levels=["Kingdom", "Phylum", "Class", "Order", "Family", "Genus"], leaf_level="Genus", make_leaf_dict=False):

    data_all = pd.concat(data.values())

    if make_leaf_dict:
        data_gen = data_all
        for t in tax_levels:
            data_gen[t] = data_gen[t].fillna("_")
        data_gen["Cell Type"] = data_gen["Type"]
    else:
        data_gen = data[leaf_level]


    newick = df2newick(data_gen.reset_index(drop=True), tax_lev=tax_levels)
    tree = tt.tree(newick, tree_format=8)

    # Edge colors
    palette = toyplot.color.brewer.palette("Set1") + toyplot.color.brewer.palette("Set3")

    if edge_color_dict is None:
        edge_color_dict = {
            "Firmicutes": palette[0],
            "Proteobacteria": palette[1],
            "Bacteroidota": palette[2],
            "Actinobacteriota": palette[3],
            "Synergistota": palette[4],
            "Cyanobacteria": palette[5],

            "Fusobacteriota": palette[6],
            "Verrucomicrobiota": palette[7],

            "Desulfobacterota": palette[8],

            "Campilobacterota": palette[9],

            "Patescibacteria": palette[10],
            "Thermoplasmatota": palette[11],
            "Euryarchaeota": palette[12],
            "Spirochaetota": palette[13],

            "Elusimicrobiota": palette[14],

            "Deinococcota": palette[15],

            "Bdellovibrionota": palette[16],
            "Abditibacteriota": palette[17],
            "Acidobacteriota": palette[18],
        }

    markers = []

    max_height = np.max([n.height - 2 for n in tree.treenode.traverse()])
    c = 0
    for n in tree.treenode.traverse():
        if n.height == max_height:
            if n.name in edge_color_dict.keys():
                col = edge_color_dict[n.name]
                if type(col) == str:
                    col2 = col.lower()
                else:
                    col2 = '#%02x%02x%02x' % tuple([int(255 * x) for x in col.tolist()[:-1]])
                m = toyplot.marker.create(shape="o", size=8, mstyle={"fill": col2})
                markers.append((n.name, m))
            else:
                col = other_col
            n.add_feature("edge_color", col)
            for n_ in n.get_descendants():
                n_.add_feature("edge_color", col)



            c += 1
        elif n.height > max_height:
            n.add_feature("edge_color", "black")

    m = toyplot.marker.create(shape="o", size=8, mstyle={"fill": other_col})
    markers.append(("other", m))

    # assign taxonomic levels to nodes
    for n in tree.treenode.traverse():
        if n.height == "":
            l = tax_levels[-1]
        elif n.height >= len(tax_levels):
            l = ""
        else:
            l = tax_levels[-(int(n.height) + 1)]
        n.add_feature("tax_level", l)

    # add effects
    data_all["level"] = "Kingdom"
    for l in tax_levels[1:]:
        data_all.loc[pd.notna(data_all[l]), "level"] = l

    for n in tree.treenode.traverse():
        if n.tax_level == "":
            n.add_feature("effect", 0)
        else:
            l = data_all.loc[(data_all["level"] == n.tax_level) & (data_all[n.tax_level] == n.name), :]
            if len(l) == 1:
                n.add_feature("effect", l["Final Parameter"].values[0])
            elif len(l) == 0:
                n.add_feature("effect", 0)
            elif len(l) > 1:
                par_names = [n.name] + [m.name for m in n.get_ancestors()][:-1]
                par_names.reverse()
                full_name = '*'.join(par_names)
                ll = l[l["Cell Type"] == full_name]
                n.add_feature("effect", ll["Final Parameter"].values[0])


    # add node colors
    for n in tree.treenode.traverse():
        if np.sign(n.effect) == 1:
            n.add_feature("color", "black")
        elif np.sign(n.effect) == -1:
            n.add_feature("color", "white")
        else:
            n.add_feature("color", "cyan")

    return tree, markers

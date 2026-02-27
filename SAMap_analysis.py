import scanpy as sc
import pandas as pd
from samalg.gui import SAMGUI
from samap.mapping import SAMAP

from samap.analysis import (get_mapping_scores, GenePairFinder,
                           sankey_plot, chord_plot, CellTypeTriangles, 
                            ParalogSubstitutions, FunctionalEnrichment,
                            convert_eggnog_to_homologs, GeneTriangles)
from samalg import SAM
import pandas as pd
import numpy as np 
import os

# read in reciprocal blast tables
A=pd.read_csv('example_data/maps/hmsp/hm_to_sp.txt',sep='\t',index_col=0,header=None)
B=pd.read_csv('example_data/maps/hmsp/sp_to_hm.txt',sep='\t',index_col=0,header=None)


# define paths to single cell datasets
fn1 = 'ourspecies/h5ad/hm_germ_8.14.h5ad'
fn2 = 'ourspecies/h5ad/sp_germ_8.15.h5ad'
filenames = {'hm':fn1,'sp':fn2}

# run SAMap
sm = SAMAP(
        filenames,
        f_maps = 'example_data/maps/',
        save_processed=True, #if False, do not save the processed results to `*_pr.h5ad`
        resolutions = {"hm":0.8, "sp":1.2},
        keys = {'sp': 'cell_idents', 'hm':'cell_idents'}       
)

sm.run(pairwise=True, neigh_from_keys={'hm':True, 'sp': False})
samap = sm.samap

print("samap finished")
# save samap object
import dill

with open("samap_object_save.pkl", "wb") as f:
    dill.dump(sm, f)

print("dill saved")

folder = 'hm_germ_cell_experiment_8.14/SAMap_results_8.15/'
if not os.path.exists("hm_germ_cell_experiment_8.14/SAMap_results_8.15"):
    os.mkdir("hm_germ_cell_experiment_8.14/SAMap_results_8.15")

# Plot overlapping UMAP graphs
xy = sm.samap.adata.obsm["X_umap"]
# set masks for subsetting 
sp_mask = sm.samap.adata.obs["species"] == "sp"
hm_mask = sm.samap.adata.obs["species"] == "hm"
hm_07_mask = sm.samap.adata.obs["hm_orig.ident"] == "hmov_7"
hm_09_mask = sm.samap.adata.obs["hm_orig.ident"] == "hmov_9"
hm_10_mask = sm.samap.adata.obs["hm_orig.ident"] == "hmov_10"
hm_13_mask = sm.samap.adata.obs["hm_orig.ident"] == "hmov_13"
hm_16_mask = sm.samap.adata.obs["hm_orig.ident"] == "hmov_16"
hm_adult_mask = sm.samap.adata.obs["hm_stage"] == "hm_adult"

def plot_overlay(coords, species_mask_list, subset: str):
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(6, 6))

    color_map = {
        "sp": "#56B4E9",     # blue
        "hm": "#E69F00",     # orange
        "hm_wk7": "#E69F00",
        "hm_wk9": "#E69F00",
        "hm_wk10": "#E69F00",
        "hm_wk13": "#E69F00",
        "hm_wk16": "#E69F00",
        "hm_adult":"#E69F00"
    }

    for i, (label, mask) in enumerate(species_mask_list):
        alpha = 0.2 if label == "sp" else 0.5
        color = color_map.get(label, "#999999")  # fallback gray if label not in map

        ax.scatter(
            coords[mask, 0],
            coords[mask, 1],
            s=8,
            c=[color],
            alpha=alpha,
            linewidths=0,
            edgecolors="none",
            label=label,
            zorder=i
        )

    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")
    ax.set_title("Overlay_sp " + subset)
    ax.legend(title="Species", loc="upper right")

    plt.savefig(folder + "overlay_hmsp_" + subset + ".pdf", format="pdf", bbox_inches="tight")
    plt.close()
plot_overlay(xy, [("hm", hm_mask), ("sp", sp_mask)], "hm_whole")
plot_overlay(xy, [("hm_wk7", hm_07_mask), ("sp", sp_mask)], "hm_wk7")
plot_overlay(xy, [("hm_wk9", hm_09_mask), ("sp", sp_mask)], "hm_wk9")
plot_overlay(xy, [("hm_wk10", hm_10_mask), ("sp", sp_mask)], "hm_wk10")
plot_overlay(xy, [("hm_wk13", hm_13_mask), ("sp", sp_mask)], "hm_wk13")
plot_overlay(xy, [("hm_wk16", hm_16_mask), ("sp", sp_mask)], "hm_wk16")
plot_overlay(xy, [("hm_adult", hm_adult_mask), ("sp", sp_mask)], "adult")


# Plt cluster identities with UMAP
import numpy as np
from adjustText import adjust_text        # pip install adjustText

# existing imports
import seaborn as sns
import matplotlib.pyplot as plt
import scanpy as sc

import seaborn as sns
import matplotlib.pyplot as plt
import scanpy as sc


import matplotlib.colors as mcolors
adata = sm.samap.adata
# Ensure that the index is also a string 
adata.obs.index = adata.obs.index.astype(str)


# Save the modified AnnData object

hm_adata = adata[adata.obs['species'] == 'hm', :]
sp_adata = adata[adata.obs['species'] == 'sp', :]


hm_07 = adata[adata.obs["hm_orig.ident"] == "hmov_7", :]
hm_09 = adata[adata.obs["hm_orig.ident"] == "hmov_9", :]
hm_10 = adata[adata.obs["hm_orig.ident"] == "hmov_10", :]
hm_13 = adata[adata.obs["hm_orig.ident"] == "hmov_13", :]
hm_16 = adata[adata.obs["hm_orig.ident"] == "hmov_16", :]
hm_adult = adata[adata.obs["hm_stage"] == "hm_adult", :]

adata_list = [("hm_", hm_adata, "hm"), ("hm_",hm_07, "hm07"), ("hm_",hm_09, "hm09"), ("hm_", hm_10, "hm10"), ("hm_", hm_13, "hm13"), ("hm_",hm_16,"hm16"), ("hm_", hm_adult, "hm_adult"), ("sp_", sp_adata, "sp")]
print(adata.obs['hm_cell_idents'].unique())
def make_palette_57(seed=0):
    # Get colors from several categorical colormaps
    cmaps = [plt.get_cmap('tab20'),
             plt.get_cmap('tab20b'),
             plt.get_cmap('tab20c'),
             plt.get_cmap('tab10')]  
    colors = []
    for cmap in cmaps:
        colors.extend([mcolors.to_hex(cmap(i)) for i in range(cmap.N)])

    # Remove duplicates while preserving order
    seen = set()
    colors = [c for c in colors if not (c in seen or seen.add(c))]

    # If we still need more, fill from HSV space
    if len(colors) < 57:
        extra_needed = 57 - len(colors)
        rng = np.random.default_rng(seed)
        hsv_positions = np.linspace(0, 1, extra_needed, endpoint=False)
        rng.shuffle(hsv_positions)  # shuffle to avoid similar adjacents
        colors += [mcolors.to_hex(plt.cm.hsv(h)) for h in hsv_positions]

    return colors[:57]  # final length exactly 57

# 2) Collect unique category names
def collect_unique_categories(adata_list, category):
    names = set()
    for species, adata, _ in adata_list:
        key = species + category
        if key in adata.obs:
            names.update(adata.obs[key].astype('category').cat.categories)
    return sorted(names, key=str)  # deterministic

# 3) Build global mapping
def build_shared_palette(uniq_names, base_colors):
    if len(uniq_names) > len(base_colors):
        raise ValueError(f"You have {len(uniq_names)} categories but only {len(base_colors)} colors.")
    return {name: base_colors[i] for i, name in enumerate(uniq_names)}

# 4) Plotting function
def plot_umap(category, adata_list, folder):
    base_colors = make_palette_57()
    uniq_names = collect_unique_categories(adata_list, category)
    shared_palette = build_shared_palette(uniq_names, base_colors)
    print(shared_palette)

    for (species, adata, label) in adata_list:
        
        key = species + category
        if key not in adata.obs:
            continue

        adata.obs[key] = adata.obs[key].astype('category')
        cats = list(adata.obs[key].cat.categories)

        adata.uns[f"{key}_colors"] = [shared_palette[c] for c in cats]

        sc.pl.umap(adata, color=[key], legend_loc="right margin", show=False)

        ax = plt.gca()
        texts = []
        umap = adata.obsm['X_umap']

        for cat in cats:
            mask = (adata.obs[key] == cat).values
            if not np.any(mask):      # <- crucial: skip empty cats
                continue
            x, y = np.median(umap[mask, 0]), np.median(umap[mask, 1])
            texts.append(ax.text(x, y, cat, fontsize=11.5, weight='bold',
                                ha='center', va='center'))

        adjust_text(texts,
                    ax=ax,
                    arrowprops=dict(arrowstyle='-', lw=0.5, color='black'))

        plt.gcf().set_size_inches(14, 10)
        plt.tight_layout()
        plt.savefig(f"{folder}{label}_in_hmsp_umap_clusters.png", dpi=600, bbox_inches='tight')
        plt.close()


plot_umap('cell_idents', adata_list, folder)

# mappingtable analysis (cell-type conservation scores)
import pandas as pd
import os 
if not(os.path.isdir(folder + "mappingtables")):
    os.mkdir(folder + "mappingtables")


keys = {'hm': 'cell_idents', 'sp':'cell_idents'}
D,MappingTable = get_mapping_scores(sm,keys,n_top = 0)
MappingTable.to_csv(folder + "mappingtables/hmsp_whole_mappingtable_preannotated.csv",index=False)


adult_mask = sm.sams['hm'].adata.obs['stage'] == "hm_adult"

obs = sm.sams['hm'].adata.obs
if obs['orig.ident'].dtype.name == 'category':
   obs['orig.ident'] = obs['orig.ident'].cat.add_categories(['adult'])

# now assign
obs.loc[adult_mask, 'orig.ident'] = "adult"


# 2‑A‑1.  Build a composite key in the *human* AnnData
hm = sm.sams['hm'].adata            # convenience alias
hm.obs['cluster_tp'] = (
    hm.obs['cell_idents'].astype(str) + '_' +
    hm.obs['orig.ident'].astype(str)           # e.g. 'Granulosa_hmov_7'
).astype('category')
print(hm.obs)
# 2‑A‑2.  Add the same key to sm.samap.adata (needed by get_mapping_scores)
#        The combined object uses the original human indices
sm.samap.adata.obs['cluster_tp'] = hm.obs['cluster_tp']
keys = {'hm': 'cluster_tp', 'sp':'cell_idents'}


# 2‑A‑3.  Get mapping scores
keys = {'hm': 'cluster_tp', 'sp': 'cell_idents'}

D,MappingTable = get_mapping_scores(sm,keys,n_top = 0)


def clean_names(mtable):
    new_columns = []
    for col in mtable.columns:
        if "hm" in col:
            first_underscore = col.find("_")
            second_underscore = col.find("_",first_underscore + 1)
            new_col = "hm_" + col[first_underscore + 1 : second_underscore]
        else:
            new_col = col
        new_columns.append(new_col)
    mtable.columns = new_columns
        
hm7_columns = MappingTable.columns.str.contains('hmov_7') | MappingTable.columns.str.contains('sp')
MappingTable_hm7 = MappingTable.iloc[hm7_columns,hm7_columns]
hm9_columns = MappingTable.columns.str.contains('hmov_9') | MappingTable.columns.str.contains('sp')
MappingTable_hm9 = MappingTable.iloc[hm9_columns,hm9_columns]
hm10_columns = MappingTable.columns.str.contains('hmov_10') | MappingTable.columns.str.contains('sp')
MappingTable_hm10 = MappingTable.iloc[hm10_columns,hm10_columns]
hm13_columns = MappingTable.columns.str.contains('hmov_13') | MappingTable.columns.str.contains('sp')
MappingTable_hm13 = MappingTable.iloc[hm13_columns,hm13_columns]
hm16_columns = MappingTable.columns.str.contains('hmov_16') | MappingTable.columns.str.contains('sp')
MappingTable_hm16 = MappingTable.iloc[hm16_columns,hm16_columns]
hmadult_columns = MappingTable.columns.str.contains('adult') | MappingTable.columns.str.contains('sp')
MappingTable_adult = MappingTable.iloc[hmadult_columns,hmadult_columns]

clean_names(MappingTable_hm7)
clean_names(MappingTable_hm9)
clean_names(MappingTable_hm10)
clean_names(MappingTable_hm13)
clean_names(MappingTable_hm16)
clean_names(MappingTable_adult)

MappingTable_hm7.to_csv(folder + "mappingtables/hm_7_sp_mappingtable_preannotated.csv",index=False)
MappingTable_hm9.to_csv(folder + "mappingtables/hm_9_sp_mappingtable_preannotated.csv",index=False)
MappingTable_hm10.to_csv(folder + "mappingtables/hm_10_sp_mappingtable_preannotated.csv",index=False)
MappingTable_hm13.to_csv(folder + "mappingtables/hm_13_sp_mappingtable_preannotated.csv",index=False)
MappingTable_hm16.to_csv(folder + "mappingtables/hm_16_sp_mappingtable_preannotated.csv",index=False)
MappingTable_adult.to_csv(folder + "mappingtables/hm_adult_sp_mappingtable_preannotated.csv",index=False)


# gene pair analysis

import pandas as pd
import os 
keys = {'hm': 'cell_idents', 'sp': 'cell_idents'}

gpf = GenePairFinder(sm,keys=keys)


gene_pairs = gpf.find_all(align_thr=0.2)


gene_df = pd.DataFrame(gene_pairs)
if not os.path.exists(folder + "DGE"):
    os.mkdir(folder + "DGE")
gene_df.to_excel(folder + "DGE/gene_pair_table.xlsx", index=False)


import matplotlib.pyplot as plt
import scanpy as sc
import matplotlib.colors as mcolors
import os


# Define a colormap where 0 = light grey and expression values go from light blue to blue
# You can also use a custom colormap if you want more control
from matplotlib.colors import ListedColormap


def gene_pair_analysis(sm, gene1, gene2, folder):
    inner_folder = folder + gene1 + "+" + gene2
    os.makedirs(inner_folder, exist_ok=True)


# First color = light grey (for 0 expression), rest = blue gradient
    from matplotlib.colors import LinearSegmentedColormap

    # Create a smooth gradient from light grey (0) to light blue to dark blue (max)
    custom_cmap = LinearSegmentedColormap.from_list(
        "smooth_blue", ["lightgrey", "#add8e6", "#0000ff"], N=256
    )

    sc.pl.umap(
        sm.samap.adata,                       # or sm.samap.adata
        color=[gene1, gene2],   # replace with your gene name
        cmap=custom_cmap,
        na_color='lightgrey',       # fallback for missing data
        vmin=0,                 # ensures 0 stays mapped to lightgrey
        vmax=0.01,   
        ncols = 2,             # scale high values appropriately
        frameon=False,
        show=False 
    )
    plt.savefig(inner_folder + '/split.png', dpi=300, bbox_inches='tight')  # Save the figure
    plt.close()  # Close the figure to avoid overlap
    sm.plot_expression_overlap({'sp':gene1,'hm':gene2})
    plt.savefig(inner_folder + '/overlap.png', dpi=300, bbox_inches='tight')  # Save the figure
    plt.close()  # Close the figure to avoid overlap
    


sp_germline = pd.read_excel('sp_experiment/sp_seurat_markers_w_germ/sp_cluster_Germline .xlsx')
sp_Epithelial = pd.read_excel('sp_experiment/sp_seurat_markers_w_germ/sp_cluster_Epithelial.xlsx')
sp_Follicle = pd.read_excel('sp_experiment/sp_seurat_markers_w_germ/sp_cluster_Follicle.xlsx')
sp_Immune = pd.read_excel('sp_experiment/sp_seurat_markers_w_germ/sp_cluster_Immune.xlsx')
sp_Muscle = pd.read_excel('sp_experiment/sp_seurat_markers_w_germ/sp_cluster_Muscle.xlsx')
sp_Neuronal = pd.read_excel('sp_experiment/sp_seurat_markers_w_germ/sp_cluster_Neuronal.xlsx')
sp_Ribosomal = pd.read_excel('sp_experiment/sp_seurat_markers_w_germ/sp_cluster_Ribosomal.xlsx')

hm_Germline = pd.read_excel('hm_harmony_analysis/hm_seurat_markers_w_germ/hm_cluster_Germline.xlsx')
hm_PGC = pd.read_excel('hm_harmony_analysis/hm_seurat_markers_w_germ/hm_cluster_PGC.xlsx')
hm_Oocyte = pd.read_excel('hm_harmony_analysis/hm_seurat_markers_w_germ/hm_cluster_Oocyte.xlsx')
hm_Common_progenitor = pd.read_excel('hm_harmony_analysis/hm_seurat_markers_w_germ/hm_cluster_Common Progenitor.xlsx')
hm_Endothelial= pd.read_excel('hm_harmony_analysis/hm_seurat_markers_w_germ/hm_cluster_Endothelial.xlsx')
hm_Epithelial = pd.read_excel('hm_harmony_analysis/hm_seurat_markers_w_germ/hm_cluster_Epithelial.xlsx')
hm_Granulosa = pd.read_excel('hm_harmony_analysis/hm_seurat_markers_w_germ/hm_cluster_Granulosa.xlsx')
hm_Immune = pd.read_excel('hm_harmony_analysis/hm_seurat_markers_w_germ/hm_cluster_Immune.xlsx')
hm_Pre_granulosa = pd.read_excel('hm_harmony_analysis/hm_seurat_markers_w_germ/hm_cluster_Pre-granulosa.xlsx')
hm_Smooth_muscle = pd.read_excel('hm_harmony_analysis/hm_seurat_markers_w_germ/hm_cluster_Smooth muscle.xlsx')
hm_Stroma = pd.read_excel('hm_harmony_analysis/hm_seurat_markers_w_germ/hm_cluster_Stroma.xlsx')
hm_Perivascular = pd.read_excel('hm_harmony_analysis/hm_seurat_markers_w_germ/hm_cluster_Perivascular.xlsx')


sp_cell_dict = {
    "germline": sp_germline,
    "Germline": sp_germline,
    "Epithelial": sp_Epithelial,
    "Follicle": sp_Follicle,
    "Immune": sp_Immune,
    "Muscle": sp_Muscle,
    "Neuronal": sp_Neuronal,
    "Ribosomal": sp_Ribosomal
}

hm_cell_dict = {
    "PGC": hm_PGC,
    "Oocyte": hm_Oocyte,
    "Germline": hm_Germline,
    "germline": hm_Germline, 
    "Common progenitor" : hm_Common_progenitor,
    "Endothelial": hm_Endothelial,
    "Epithelial": hm_Epithelial,
    "Granulosa": hm_Granulosa,
    "Immune": hm_Immune,
    "Pre-granulosa": hm_Pre_granulosa,
    "Smooth muscle": hm_Smooth_muscle,
    "Stroma": hm_Stroma,
    "Perivascular": hm_Perivascular

}

def extract_cell_type(s):
    import re
    cell_ident = s.split('_')[1]
    match = re.search(r'\s(?=\d)', cell_ident)
    return cell_ident[:match.start()] if match else cell_ident


def check_in_marker_table(gene, marker_names):
    return gene in set(marker_names)

def find_specific(gp_table, sp_markers,hm_markers, DGE_folder_name):
    import os
    sp_mask = gp_table['sp'].apply(lambda gene: check_in_marker_table(gene, sp_markers['names']))
    gp_table = gp_table.loc[sp_mask,]
    hm_mask = gp_table['hm'].apply(lambda gene: check_in_marker_table(gene, hm_markers['names']))
    gp_table = gp_table.loc[hm_mask,]
    gp_table = gp_table.merge(sp_markers[['names', 'pct_ratio']], left_on='sp', right_on='names', how='left', suffixes=('', '_sp'))
    gp_table = gp_table.merge(hm_markers[['names', 'pct_ratio']], left_on='hm', right_on='names', how='left', suffixes=('', '_hm'))
    print(gp_table)

    # Compute a combined score, e.g., product or sum of both pct_ratios
    gp_table['combined_score'] = gp_table['pct_ratio'] * gp_table['pct_ratio_hm']

    # Sort by combined score (descending)
    gp_table = gp_table.sort_values(by='combined_score', ascending=False)

    # Loop over top 20 gene pairs
    count = 0
    for row in gp_table.itertuples(index=False):
        gene_pair_analysis(sm, "sp_" + row.sp, "hm_" + row.hm, DGE_folder_name)
        count += 1
        if count >= 20:
            break

    return(gp_table.iloc[:,0])



def trim_sp(col):
    col = col[col.notna()]
    return col.apply(lambda gp: gp.split(';')[1][3:])
    

def trim_hm(col):
    col = col[col.notna()]
    return col.apply(lambda gp: gp.split(';')[0][3:])



cols = []

def plot_top_pairs(gene_df, folder, hm_pct, sp_pct):
    filtered_gene_table = pd.DataFrame()
    for i in range(len(gene_df.columns) // 3): 
        col = gene_df.iloc[:,3*i]
        subfolder = os.path.join(folder, col.name)
        if not os.path.isdir(subfolder):
            os.makedirs(subfolder)
        

        # create dataframe of gene names 
        col = col[col.notna()]
        gp_df = pd.DataFrame(col)
        gp_df['sp'] = trim_sp(col)
        gp_df['hm'] = trim_hm(col)
        # map to corresponding marker tables
        print(col.name) 
        names = col.name.split(";")
        names = pd.Series(names)
        names = names.apply(extract_cell_type)
        hm_markers = hm_cell_dict[names[0]]
        sp_markers = sp_cell_dict[names[1]]
        hm_filter_mask = (hm_markers['pct_ratio'].isna()) | (hm_markers['pct_ratio'] > hm_pct)
        sp_filter_mask = (sp_markers['pct_ratio'].isna()) | (sp_markers['pct_ratio'] > sp_pct)
        hm_filtered = hm_markers[hm_filter_mask]
        sp_filtered = sp_markers[sp_filter_mask]
    
        interaction_column = find_specific(gp_df, sp_filtered, hm_filtered, subfolder + "/")
        interaction_column = interaction_column.reset_index(drop=True)
        cols.append(interaction_column)


       
        #filtered_gene_table[interaction_column.name] = interaction_column
    filtered_gene_table = pd.concat(cols, axis=1)
    filtered_gene_table.to_excel(folder + 'filtered_gene_table.xlsx')
    return(filtered_gene_table)



filtered_gene_table = plot_top_pairs(gene_df, folder + "DGE/",4,4)


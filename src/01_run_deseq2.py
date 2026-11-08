import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pydeseq2 as pyd 
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
import gseapy as gp
import os

os.makedirs("results", exist_ok=True)
os.makedirs("logs", exist_ok=True)

counts_path = "data/raw/counts.csv.gz"
metadata_path = "data/raw/metadata.csv.gz"
condition_column = "site"
plate_column = "plate"
condition_test = "Primary_Solid_Tumor"
condition_ref = "Solid_Tissue_Normal"
min_count = 10
min_samples = 20
alpha = 0.05
lfc_threshold = 1

print("[1/8] Leyendo datos de entrada (TCGA-BRCA)...")
counts_df = pd.read_csv(counts_path, index_col= 0)
metadata = pd.read_csv(metadata_path, index_col = 0)

print(counts_df.shape) #Comprobacion de las dimensiones de la matriz

print("[2/8] Alineando muestras entre conteos y metadatos...")
metadata = metadata.loc[counts_df.columns]
print(metadata.head()) 
print(counts_df.columns) #Comprobacion de orden entre metadatos y conteos

metadata[condition_column] = metadata[condition_column].astype('category')
metadata[plate_column] = metadata[plate_column].astype('category')
metadata[condition_column] = pd.Categorical(metadata[condition_column], categories= [condition_ref, condition_test], ordered= False)
print(counts_df.columns == metadata.index) #Comprobacion del alineamiento

print("[3/8] Filtrando genes de baja cuenta y creando objeto DESeq2...")
counts_df_filtered = (counts_df >= min_count).sum(axis=1) >= min_samples
counts_df_filtered = counts_df.loc[counts_df_filtered, :]
print(counts_df.info(), counts_df_filtered.info()) #Comprobacion del filtrado

counts_df_dds = counts_df_filtered.T
print(counts_df_dds.info()) #Comprobacion traspuesta

dds = DeseqDataSet(counts= counts_df_dds, metadata= metadata, design= f"~{plate_column} + {condition_column}")
print(dds)
print(dds.obsm["design_matrix"].head()) #Comprobacion dds

print("[4/8] Ejecutando DESeq...")
dds.deseq2()
plt.ion()
dds.plot_dispersions(save_path="results/plot_dispersions.png")
plt.title("Dispersion dds")
plt.close("all")

print("[5/8] Extrayendo resultados y guardando tablas")
result = DeseqStats(dds, contrast=[condition_column, condition_test, condition_ref], alpha= alpha)
result.summary()
result_df = result.results_df
print(result_df.head()) #Comprobar tabla resultados del contraste

result_df["gene_id"] = result_df.index
results_df= result_df.sort_values(by = "padj", ascending = True, na_position = "last")
print(results_df.head()) #Comprobacion de la reorganizacion

results_df.to_csv("results/deseq2_all_results.tsv", sep= "\t", index= True, header=True)

sig_df = results_df[results_df['padj'].notna() & (results_df['padj'] < alpha)]
de_df = results_df[results_df['padj'].notna() & (results_df['padj'] < alpha) & (results_df['log2FoldChange'].abs() > lfc_threshold)]
up_df = de_df[de_df["log2FoldChange"] > 0]
down_df = de_df[de_df["log2FoldChange"] < 0]

top10_up_tsv = up_df.head(10)
top10_down_tsv = down_df.head(10)

de_df["top10_label"] = "other"
de_df.loc[de_df['gene_id'].isin(top10_up_tsv['gene_id']), 'top10_label'] = 'top10_up'
de_df.loc[de_df['gene_id'].isin(top10_down_tsv['gene_id']), 'top10_label'] = 'top10_down'

sig_df.to_csv("results/deseq2_sig_padj_0.05.tsv", sep= "\t", index= True, header=True)
de_df.to_csv("results/deseq2_de_padj_0.05_lfc_1.tsv", sep= "\t", index= True, header=True)
top10_up_tsv.to_csv("results/top10_up_genes.tsv", sep= "\t", index= True, header=True)
top10_down_tsv.to_csv("results/top10_down_genes.tsv", sep= "\t", index= True, header=True)

print("[6/8] Representación gráfica de resultados...")

up_df_lfc = de_df[de_df["log2FoldChange"] > 0].sort_values("log2FoldChange", ascending=False)
down_df_lfc = de_df[de_df["log2FoldChange"] < 0].sort_values("log2FoldChange", ascending=True)
top10_up_plot = up_df_lfc.head(10)
top10_down_plot = down_df_lfc.head(10)

plt.ion()
fig, axs = plt.subplots(1, 2, figsize=(10, 4))
axs[0].barh( y = top10_down_plot["gene_id"], color= "#0A4281FF", width=  top10_down_plot["log2FoldChange"])
axs[0].set_title("Top 10 down")
axs[0].set_xlabel("log2FoldChange")
axs[1].barh( y = top10_up_plot["gene_id"], color= "#FB2A00FF", width=  top10_up_plot["log2FoldChange"])
axs[1].set_title("Top 10 up")
axs[1].set_xlabel("log2FoldChange")
plt.tight_layout()
plt.savefig("results/top_10_genes.png", dpi=300)
plt.close("all")

volcano_df = results_df.dropna(subset=['padj'])
volcano_df["status"] = "No significativo"
volcano_df.loc[volcano_df["gene_id"].isin(up_df["gene_id"]), "status"] = "Up"
volcano_df.loc[volcano_df["gene_id"].isin(down_df["gene_id"]), "status"] = "Down"

color_map = {'No significativo': "#676767FF", 'Up': '#FB2A00FF', 'Down': '#0A4281FF'}
volcano_df["color"] = volcano_df["status"].map(color_map)

plt.ion()
plt.scatter(x = volcano_df['log2FoldChange'], y = -np.log10(volcano_df['padj']), s=10, c=volcano_df["color"] )
plt.axvline(lfc_threshold, color='black', linestyle='--', linewidth=0.8)
plt.axvline(-lfc_threshold, color='black', linestyle='--', linewidth=0.8)
plt.axhline(-np.log10(alpha), color='black', linestyle='--', linewidth=0.8)
plt.xlabel('log2FoldChange')
plt.ylabel('-log10(padj)')
plt.title("Volcano plot", fontsize=14, fontweight='bold')
plt.tight_layout()

volcano_tpo10 = volcano_df.head(10)
for _, row in volcano_tpo10.iterrows():
  plt.text(
      row['log2FoldChange'],
      -np.log10(row['padj']),
      row['gene_id'],
      fontsize=8,
      ha='center',
      va='bottom', 
  )

plt.savefig('results/volcano_plot.png', dpi=300)
plt.close("all")

MA_df = results_df.dropna(subset=['padj', 'baseMean'])
MA_df["status"] = "No significativo"
MA_df.loc[MA_df["gene_id"].isin(up_df["gene_id"]), "status"] = "Up"
MA_df.loc[MA_df["gene_id"].isin(down_df["gene_id"]), "status"] = "Down"
MA_df["color"] = MA_df["status"].map(color_map)

plt.ion()
plt.scatter(x = np.log10(MA_df['baseMean']), y = MA_df["log2FoldChange"], s=10, c=MA_df["color"])
plt.axhline(0, color='black', linestyle='--', linewidth=0.8)
plt.xlabel('mean of normalized counts')
plt.ylabel('log2 Fold Change')
plt.title('MA Plot', fontsize=14, fontweight='bold')
plt.savefig('results/MA_plot.png', dpi=300)
plt.close("all")

Highlight = de_df.head(1)
Highlight = Highlight[['gene_id', 'log2FoldChange', 'padj']]
Highlight.to_csv("results/top_gene_context.txt", sep="\t", index=False)

print("[7/8] Enriquecimiento GO y Reactome...")
enrich_up_genes = up_df["gene_id"].tolist()
enrich_down_genes = down_df["gene_id"].tolist()
universe = de_df["gene_id"].tolist()

def enrichment (lista_gen, libreria, salidatsv, salidaplot, universo = universe, alfa = alpha):
  enrich = gp.enrichr(gene_list = lista_gen, gene_sets = libreria, background = universo, organism = "human", outdir=None)
  enrich_results = enrich.results[enrich.results["Adjusted P-value"] < alfa]
  enrich_results.to_csv(f"results/{salidatsv}.tsv", sep="\t", index=False)
  gp.dotplot(enrich_results, column="Adjusted P-value", top_term=10, title=f"{salidatsv}", ofname=f'results/{salidaplot}.png')

enrichment(enrich_up_genes, "GO_Biological_Process_2026","GO_Up","GO_Up_Dotplot")
enrichment(enrich_down_genes, "GO_Biological_Process_2026","GO_Down","GO_Down_Dotplot")
enrichment(enrich_up_genes, "Reactome_Pathways_2024","Rea_Up","Reactome_Up_Dotplot")
enrichment(enrich_down_genes, "Reactome_Pathways_2024","Rea_Down","Reactome_Down_Dotplot" )

print("[8/8] Guardando trazabilidad de sesion...")
import sys
save_session_info = [f"Python version: {sys.version}",f"numpy: {np.__version__}", f"pandas: {pd.__version__}",
                     f"matplotlib: {plt.matplotlib.__version__}", f"gseapy: {gp.__version__}", f"pydeseq2: {pyd.__version__}"]

open("logs/session_info.txt", "w").write("\n".join(save_session_info))

summary_lines = [f"n_samples: {counts_df.shape[1]}",f"n_genes_input: {counts_df.shape[0]}",
                 f"n_genes_after_filter: {dds.n_vars}",f"n_significant_padj_lt_{alpha}: {len(sig_df)}",
                 f"n_de_absLog2FC_gt_{lfc_threshold}_padj_lt{alpha}: {len(de_df)}",f"n_up: {len(up_df)}",
                 f"n_down: {len(down_df)}"]

open("logs/run_summary.txt", "w").write("\n".join(summary_lines))
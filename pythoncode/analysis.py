import scanpy as sc
import matplotlib
matplotlib.use('Agg')  # 服务器无GUI，必须加这行
import matplotlib.pyplot as plt
import os

sc.settings.verbosity = 3
os.makedirs('./figures', exist_ok=True)
sc.settings.figdir = './figures/'

# ========== 1. 读取6个样本 ==========
samples = {
    'A':   './matrix/A/',
    'B':   './matrix/B/',
    'C':   './matrix/C/',
    'A_2w': './matrix/A_2w/',
    'B_2w': './matrix/B_2w/',
    'C_2w': './matrix/C_2w/',
}

adatas = {}
for name, path in samples.items():
    # 自动找mtx所在的子目录
    for root, dirs, files in os.walk(path):
        if any(f.endswith('matrix.mtx.gz') or f == 'matrix.mtx' for f in files):
            adatas[name] = sc.read_10x_mtx(root, var_names='gene_symbols', cache=True)
            adatas[name].var_names_make_unique()
            adatas[name].obs['sample'] = name
            print(f"{name}: {adatas[name].n_obs} cells, {adatas[name].n_vars} genes")
            break

# ========== 2. 合并 ==========
import anndata
adata = anndata.concat(adatas, keys=list(adatas.keys()))
print(f"\n合并后：{adata.n_obs} cells, {adata.n_vars} genes")

# ========== 3. QC ==========
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], inplace=True)

sc.pl.violin(adata, ['n_genes_by_counts','total_counts','pct_counts_mt'],
             groupby='sample', jitter=0.4, rotation=45, save='_qc.png')

print("\nQC统计：")
print(adata.obs[['n_genes_by_counts','total_counts','pct_counts_mt']].describe())

# ========== 4. 过滤 ==========
n_before = adata.n_obs
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
adata = adata[adata.obs.n_genes_by_counts < 6000]
adata = adata[adata.obs.pct_counts_mt < 20]
print(f"\n过滤前：{n_before} cells → 过滤后：{adata.n_obs} cells")

# ========== 5. 归一化 ==========
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata

# ========== 6. 高变基因 + 降维 ==========
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5, batch_key='sample')
adata = adata[:, adata.var.highly_variable]
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata)
sc.pl.pca_variance_ratio(adata, n_pcs=50, save='_pca.png')

# ========== 7. 整合批次（Harmony） + UMAP + 聚类 ==========
try:
    sc.external.pp.harmony_integrate(adata, 'sample')
    use_rep = 'X_pca_harmony'
    print("使用Harmony整合批次")
except:
    use_rep = 'X_pca'
    print("Harmony不可用，使用普通PCA")

sc.pp.neighbors(adata, n_neighbors=10, n_pcs=30, use_rep=use_rep)
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution=0.5)

# ========== 8. 出图 ==========
sc.pl.umap(adata, color=['leiden'], save='_clusters.png')
sc.pl.umap(adata, color=['sample'], save='_samples.png')
sc.pl.umap(adata, color=['n_genes_by_counts','pct_counts_mt'], save='_qc_umap.png')

# ========== 9. 保存 ==========
adata.write('processed.h5ad')
print("\n完成！文件保存至 processed.h5ad")
print(f"聚类数：{adata.obs.leiden.nunique()}")

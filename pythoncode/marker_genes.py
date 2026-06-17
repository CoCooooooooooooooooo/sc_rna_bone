import scanpy as sc
import matplotlib
matplotlib.use('Agg')
import pandas as pd
import os

os.makedirs('./figures', exist_ok=True)
sc.settings.figdir = './figures/'

# 读取之前的结果
print("读取processed.h5ad...")
adata = sc.read_h5ad('processed.h5ad')
print(f"细胞数：{adata.n_obs}, 聚类数：{adata.obs.leiden.nunique()}")

# 恢复raw数据用于差异分析
adata_raw = adata.raw.to_adata()
adata_raw.obs = adata.obs  # 把聚类信息带过来

# 找每个聚类的marker基因
print("\n计算marker基因...")
sc.tl.rank_genes_groups(adata_raw, 'leiden', method='wilcoxon', n_genes=50)

# 出dot plot（最直观）
sc.pl.rank_genes_groups_dotplot(
    adata_raw, n_genes=5, save='_markers_dotplot.png')

# 出每个聚类top基因的热图
sc.pl.rank_genes_groups_heatmap(
    adata_raw, n_genes=5, save='_markers_heatmap.png')

# 保存marker基因表格
result = adata_raw.uns['rank_genes_groups']
markers_df = pd.DataFrame({
    cluster: result['names'][cluster]
    for cluster in adata.obs.leiden.unique()
})
markers_df.to_csv('marker_genes.csv', index=False)
print("marker基因保存至 marker_genes.csv")

# 同时保存带统计量的完整表格
markers_full = sc.get.rank_genes_groups_df(adata_raw, group=None)
markers_full.to_csv('marker_genes_full.csv', index=False)
print("完整marker基因统计保存至 marker_genes_full.csv")

print("\n完成！")

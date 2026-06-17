import scanpy as sc
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import os

os.makedirs('./figures/comparison', exist_ok=True)

print("读取数据...")
adata = sc.read_h5ad('processed_annotated.h5ad')

# ========== 1. 计算各组细胞比例 ==========
def get_proportion(adata, sample_list, label):
    sub = adata[adata.obs['sample'].isin(sample_list)]
    counts = sub.obs['cell_type'].value_counts()
    prop = (counts / counts.sum() * 100).round(2)
    prop.name = label
    return prop

groups = {
    'A_1w':  ['A'],
    'B_1w':  ['B'],
    'C_1w':  ['C'],
    'A_2w':  ['A_2w'],
    'B_2w':  ['B_2w'],
    'C_2w':  ['C_2w'],
}

prop_df = pd.DataFrame({k: get_proportion(adata, v, k) for k, v in groups.items()})
prop_df = prop_df.fillna(0)
prop_df.to_csv('cell_proportions.csv')
print("\n各组细胞比例（%）：")
print(prop_df.round(1).to_string())

# ========== 2. 画图函数 ==========
# 颜色方案
colors = plt.cm.tab20.colors

def plot_comparison(df, groups_to_plot, title, filename):
    """柱状图对比细胞比例"""
    sub = df[groups_to_plot]
    cell_types = sub.index.tolist()
    x = np.arange(len(cell_types))
    width = 0.8 / len(groups_to_plot)

    fig, ax = plt.subplots(figsize=(14, 6))
    for i, col in enumerate(groups_to_plot):
        ax.bar(x + i * width, sub[col], width, label=col,
               color=colors[i], alpha=0.85)

    ax.set_xticks(x + width * (len(groups_to_plot)-1) / 2)
    ax.set_xticklabels(cell_types, rotation=45, ha='right', fontsize=9)
    ax.set_ylabel('Cell proportion (%)')
    ax.set_title(title, fontsize=13)
    ax.legend(fontsize=10)
    plt.tight_layout()
    plt.savefig(f'figures/comparison/{filename}', dpi=150)
    plt.close()
    print(f"  已保存: {filename}")

def plot_heatmap(df, title, filename):
    """热图"""
    fig, ax = plt.subplots(figsize=(10, 8))
    im = ax.imshow(df.values, aspect='auto', cmap='YlOrRd')
    ax.set_xticks(range(len(df.columns)))
    ax.set_xticklabels(df.columns, rotation=45, ha='right')
    ax.set_yticks(range(len(df.index)))
    ax.set_yticklabels(df.index, fontsize=9)
    plt.colorbar(im, ax=ax, label='Cell proportion (%)')
    ax.set_title(title, fontsize=13)
    plt.tight_layout()
    plt.savefig(f'figures/comparison/{filename}', dpi=150)
    plt.close()
    print(f"  已保存: {filename}")

def plot_stacked_bar(df, groups_to_plot, title, filename):
    """堆叠柱状图"""
    sub = df[groups_to_plot].T
    fig, ax = plt.subplots(figsize=(10, 6))
    bottom = np.zeros(len(groups_to_plot))
    for i, ct in enumerate(sub.columns):
        ax.bar(groups_to_plot, sub[ct], bottom=bottom,
               label=ct, color=colors[i % len(colors)], alpha=0.85)
        bottom += sub[ct].values
    ax.set_ylabel('Cell proportion (%)')
    ax.set_title(title, fontsize=13)
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8)
    plt.tight_layout()
    plt.savefig(f'figures/comparison/{filename}', dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  已保存: {filename}")

# ========== 3. 总览热图 ==========
print("\n生成总览热图...")
plot_heatmap(prop_df, 'Cell type proportions across all groups', 'heatmap_all.png')
plot_stacked_bar(prop_df, list(prop_df.columns),
                 'Cell composition - all groups', 'stacked_all.png')

# ========== 4. C vs A 对比 ==========
print("\n生成 C vs A 对比图...")
plot_comparison(prop_df, ['A_1w','C_1w'],
                'C vs A at 1 week', 'CvsA_1w.png')
plot_comparison(prop_df, ['A_2w','C_2w'],
                'C vs A at 2 weeks', 'CvsA_2w.png')
plot_comparison(prop_df, ['A_1w','C_1w','A_2w','C_2w'],
                'C vs A: 1w and 2w comparison', 'CvsA_combined.png')
plot_stacked_bar(prop_df, ['A_1w','C_1w','A_2w','C_2w'],
                 'C vs A: cell composition', 'CvsA_stacked.png')

# ========== 5. C vs B 对比 ==========
print("\n生成 C vs B 对比图...")
plot_comparison(prop_df, ['B_1w','C_1w'],
                'C vs B at 1 week', 'CvsB_1w.png')
plot_comparison(prop_df, ['B_2w','C_2w'],
                'C vs B at 2 weeks', 'CvsB_2w.png')
plot_comparison(prop_df, ['B_1w','C_1w','B_2w','C_2w'],
                'C vs B: 1w and 2w comparison', 'CvsB_combined.png')
plot_stacked_bar(prop_df, ['B_1w','C_1w','B_2w','C_2w'],
                 'C vs B: cell composition', 'CvsB_stacked.png')

# ========== 6. C组时序变化 ==========
print("\n生成 C组时序变化图...")
plot_comparison(prop_df, ['C_1w','C_2w'],
                'C group: 1w vs 2w temporal change', 'C_temporal.png')
plot_stacked_bar(prop_df, ['C_1w','C_2w'],
                 'C group: temporal cell composition', 'C_temporal_stacked.png')

# ========== 7. 计算差值表 ==========
print("\n计算各组差异...")
diff_CvsA_1w = (prop_df['C_1w'] - prop_df['A_1w']).sort_values(ascending=False)
diff_CvsA_2w = (prop_df['C_2w'] - prop_df['A_2w']).sort_values(ascending=False)
diff_CvsB_1w = (prop_df['C_1w'] - prop_df['B_1w']).sort_values(ascending=False)
diff_CvsB_2w = (prop_df['C_2w'] - prop_df['B_2w']).sort_values(ascending=False)
diff_C_time  = (prop_df['C_2w'] - prop_df['C_1w']).sort_values(ascending=False)

diff_summary = pd.DataFrame({
    'C-A_1w': diff_CvsA_1w,
    'C-A_2w': diff_CvsA_2w,
    'C-B_1w': diff_CvsB_1w,
    'C-B_2w': diff_CvsB_2w,
    'C_2w-1w': diff_C_time,
})
diff_summary.to_csv('proportion_differences.csv')
print("\n比例差值（正值=C组更多，负值=C组更少）：")
print(diff_summary.round(2).to_string())

# ========== 8. 各组UMAP ==========
print("\n生成各组UMAP...")
sc.settings.figdir = './figures/comparison/'
for sample in ['A','B','C','A_2w','B_2w','C_2w']:
    sub = adata[adata.obs['sample'] == sample]
    sc.pl.umap(sub, color='cell_type', title=f'Sample {sample} ({sub.n_obs} cells)',
               legend_loc='right margin', legend_fontsize=7,
               save=f'_umap_{sample}.png')
    print(f"  {sample}: {sub.n_obs} cells")

print("\n全部完成！")

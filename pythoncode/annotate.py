import scanpy as sc
import matplotlib
matplotlib.use('Agg')
import os
os.makedirs('./figures', exist_ok=True)
sc.settings.figdir = './figures/'

adata = sc.read_h5ad('processed.h5ad')

cell_type_map = {
    '0':  'T cell',
    '1':  'Fibroblast',
    '2':  'Neutrophil',
    '3':  'Neutrophil/Monocyte',
    '4':  'Chondrocyte-related',
    '5':  'Fibroblast',
    '6':  'Proliferating cell',
    '7':  'Mature Neutrophil',
    '8':  'Erythroid precursor',
    '9':  'Macrophage/DC',
    '10': 'Antigen-presenting cell',
    '11': 'Neutrophil (defensin+)',
    '12': 'Erythrocyte',
    '13': 'B cell/DC',
    '14': 'Stromal/Perivascular cell',
    '15': 'B cell precursor',
    '16': 'Endothelial cell',
    '17': 'Dendritic cell',
    '18': 'Mast cell',
    '19': 'Chondrocyte',
}

adata.obs['cell_type'] = adata.obs['leiden'].map(cell_type_map)

sc.pl.umap(adata, color='cell_type', legend_loc='on data',
           legend_fontsize=8, save='_annotated.png')
sc.pl.umap(adata, color='cell_type', legend_loc='right margin',
           save='_annotated_legend.png')

adata.write('processed_annotated.h5ad', compression='gzip')
adata.obs[['sample','leiden','cell_type']].to_csv('cell_annotation.csv')
print("完成！")
print(adata.obs['cell_type'].value_counts())

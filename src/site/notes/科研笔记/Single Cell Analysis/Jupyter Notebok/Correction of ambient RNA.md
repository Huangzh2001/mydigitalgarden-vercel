---
{"dg-publish":true,"tags":["Notes","Single-Cell-Analysis","jupyter"],"permalink":"/科研笔记/Single Cell Analysis/Jupyter Notebok/Correction of ambient RNA/","dgPassFrontmatter":true}
---

# 1.导入包及数据
## (1)导入包及初始化参数
Seaborn是在Matplotlib的基础上开发的高级可视化库.有时候我们在Matplotlib上画一个好看的图需要输入很多参数,很费劲.但Seaborn相当于帮我们调好了参数,我们直接用就行了.


```python
import numpy as np
import scanpy as sc
import pandas as pd
import seaborn as sns
```

接着我们初始化scanpy的参数.


```python
## 设置日志等级
sc.settings.verbosity = 0
## 设置画图参数
sc.settings.set_figure_params(
    dpi = 80,
    facecolor = "white",
    frameon = False,
)
```

## (2)导入数据

数据我们可以在网站:'https://figshare.com/ndownloader/files/39546196'上下载.

文件名叫'filtered_feature_bc_matrix.h5'.


```python
## 数据集托管在Figshare上
adata = sc.read_10x_h5("D:/Document/DataSet/Single Cell Dataset/filtered_feature_bc_matrix.h5")
adata
```

    D:\Software\anaconda3\envs\scverse\lib\site-packages\anndata\_core\anndata.py:1899: UserWarning: Variable names are not unique. To make them unique, call `.var_names_make_unique`.
      utils.warn_names_duplicates("var")
    D:\Software\anaconda3\envs\scverse\lib\site-packages\anndata\_core\anndata.py:1899: UserWarning: Variable names are not unique. To make them unique, call `.var_names_make_unique`.
      utils.warn_names_duplicates("var")
    




    AnnData object with n_obs × n_vars = 16934 × 36601
        var: 'gene_ids', 'feature_types', 'genome'



读取数据后，scanpy 会显示一条警告，指出并非所有变量名称都是唯一的。
这表明某些变量（=基因）出现多次，这可能会导致下游分析任务出现错误或意外行为。
使用函数`var_names_make_unique()`来使变量名称唯一。


```python
adata.var_names_make_unique()
adata
```




    AnnData object with n_obs × n_vars = 16934 × 36601
        var: 'gene_ids', 'feature_types', 'genome'



从输出结果中我们可以获知:计数矩阵的形状为`n_obs x n_vars=16934 x 36601`.然后对基因还有三个注释:`gene_ids`,`feature_types`,`genome`.


```python
df = adata.to_df()
print(df)
```

                        MIR1302-2HG  FAM138A  OR4F5  AL627309.1  AL627309.3  \
    AAACAGCCAAGCTTAT-1          0.0      0.0    0.0         0.0         0.0   
    AAACAGCCATAGCTTG-1          0.0      0.0    0.0         0.0         0.0   
    AAACAGCCATGAAATG-1          0.0      0.0    0.0         0.0         0.0   
    AAACAGCCATGTTTGG-1          0.0      0.0    0.0         0.0         0.0   
    AAACATGCAACGTGCT-1          0.0      0.0    0.0         0.0         0.0   
    ...                         ...      ...    ...         ...         ...   
    TTTGTTGGTGGCTTCC-1          0.0      0.0    0.0         0.0         0.0   
    TTTGTTGGTTCTTTAG-1          0.0      0.0    0.0         0.0         0.0   
    TTTGTTGGTTGGCCGA-1          0.0      0.0    0.0         0.0         0.0   
    TTTGTTGGTTTACTTG-1          0.0      0.0    0.0         0.0         0.0   
    TTTGTTGGTTTGTGGA-1          0.0      0.0    0.0         0.0         0.0   
    
                        AL627309.2  AL627309.5  AL627309.4  AP006222.2  \
    AAACAGCCAAGCTTAT-1         0.0         0.0         0.0         0.0   
    AAACAGCCATAGCTTG-1         0.0         0.0         0.0         0.0   
    AAACAGCCATGAAATG-1         0.0         0.0         0.0         0.0   
    AAACAGCCATGTTTGG-1         0.0         0.0         0.0         0.0   
    AAACATGCAACGTGCT-1         0.0         0.0         0.0         0.0   
    ...                        ...         ...         ...         ...   
    TTTGTTGGTGGCTTCC-1         0.0         0.0         0.0         0.0   
    TTTGTTGGTTCTTTAG-1         0.0         0.0         0.0         0.0   
    TTTGTTGGTTGGCCGA-1         0.0         0.0         0.0         0.0   
    TTTGTTGGTTTACTTG-1         0.0         0.0         0.0         0.0   
    TTTGTTGGTTTGTGGA-1         0.0         0.0         0.0         0.0   
    
                        AL732372.1  ...  AC133551.1  AC136612.1  AC136616.1  \
    AAACAGCCAAGCTTAT-1         0.0  ...         0.0         0.0         0.0   
    AAACAGCCATAGCTTG-1         0.0  ...         0.0         0.0         0.0   
    AAACAGCCATGAAATG-1         0.0  ...         0.0         0.0         0.0   
    AAACAGCCATGTTTGG-1         0.0  ...         0.0         0.0         0.0   
    AAACATGCAACGTGCT-1         0.0  ...         0.0         0.0         0.0   
    ...                        ...  ...         ...         ...         ...   
    TTTGTTGGTGGCTTCC-1         0.0  ...         0.0         0.0         0.0   
    TTTGTTGGTTCTTTAG-1         0.0  ...         0.0         0.0         0.0   
    TTTGTTGGTTGGCCGA-1         0.0  ...         0.0         0.0         0.0   
    TTTGTTGGTTTACTTG-1         0.0  ...         0.0         0.0         0.0   
    TTTGTTGGTTTGTGGA-1         0.0  ...         0.0         0.0         0.0   
    
                        AC136616.3  AC136616.2  AC141272.1  AC023491.2  \
    AAACAGCCAAGCTTAT-1         0.0         0.0         0.0         0.0   
    AAACAGCCATAGCTTG-1         0.0         0.0         0.0         0.0   
    AAACAGCCATGAAATG-1         0.0         0.0         0.0         0.0   
    AAACAGCCATGTTTGG-1         0.0         0.0         0.0         0.0   
    AAACATGCAACGTGCT-1         0.0         0.0         0.0         0.0   
    ...                        ...         ...         ...         ...   
    TTTGTTGGTGGCTTCC-1         0.0         0.0         0.0         0.0   
    TTTGTTGGTTCTTTAG-1         0.0         0.0         0.0         0.0   
    TTTGTTGGTTGGCCGA-1         0.0         0.0         0.0         0.0   
    TTTGTTGGTTTACTTG-1         0.0         0.0         0.0         0.0   
    TTTGTTGGTTTGTGGA-1         0.0         0.0         0.0         0.0   
    
                        AC007325.1  AC007325.4  AC007325.2  
    AAACAGCCAAGCTTAT-1         0.0         0.0         0.0  
    AAACAGCCATAGCTTG-1         0.0         0.0         0.0  
    AAACAGCCATGAAATG-1         0.0         0.0         0.0  
    AAACAGCCATGTTTGG-1         0.0         0.0         0.0  
    AAACATGCAACGTGCT-1         0.0         0.0         0.0  
    ...                        ...         ...         ...  
    TTTGTTGGTGGCTTCC-1         0.0         0.0         0.0  
    TTTGTTGGTTCTTTAG-1         0.0         0.0         0.0  
    TTTGTTGGTTGGCCGA-1         0.0         0.0         0.0  
    TTTGTTGGTTTACTTG-1         0.0         0.0         0.0  
    TTTGTTGGTTTGTGGA-1         0.0         0.0         0.0  
    
    [16934 rows x 36601 columns]
    

我们从上面的输出结果可以看出行指标是一个个barcode;列指标是一个个基因编号,编号开头的字母表示不同来源的基因,比如来自线粒体的基因开头为"MT";来自核糖体的基因开头为"RPS", "RPL";来自血红蛋白的基因开头为"^HB[^(P)]".

# 2.校正环境RNA

## (1)导入包与初始化


```python
import os
os.environ['R_HOME'] = "D:/Software/R-4.3.1"
```


```python
import anndata2ri
import logging

import rpy2.rinterface_lib.callbacks as rcb
import rpy2.robjects as ro

rcb.logger.setLevel(logging.ERROR)
ro.pandas2ri.activate()
anndata2ri.activate()

%load_ext rpy2.ipython
```

    D:\Software\anaconda3\envs\scverse\lib\site-packages\rpy2\robjects\packages.py:367: UserWarning: The symbol 'quartz' is not in this R namespace/package.
      warnings.warn(
    

## (2)聚类并获取聚类结果

首先，我们生成 AnnData 对象的副本，对其进行标准化和 log1p 转换。


```python
adata_pp = adata.copy()
sc.pp.normalize_per_cell(adata_pp)
sc.pp.log1p(adata_pp)
```

接下来，我们计算数据的主成分以获得较低维的表示。然后使用该表示来生成数据的邻域图并在 KNN 图上运行leiden聚类。

我们将聚类结果储存在`soupx_groups`中


```python
sc.pp.pca(adata_pp)
sc.pp.neighbors(adata_pp)
sc.tl.leiden(adata_pp, key_added="soupx_groups")

soupx_groups = adata_pp.obs["soupx_groups"]
```

    D:\Software\anaconda3\envs\scverse\lib\site-packages\tqdm\auto.py:21: TqdmWarning: IProgress not found. Please update jupyter and ipywidgets. See https://ipywidgets.readthedocs.io/en/stable/user_install.html
      from .autonotebook import tqdm as notebook_tqdm
    

删除 AnnData 对象的副本以节省内存


```python
del adata_pp
```

## (3)利用R包SoupX进行环境mRNA校正

### 使用前的准备
我们保存细胞名称、基因名称和转置后的计数矩阵(SoupX 需要形状为 features x barcodes 的矩阵，因此我们必须转置).


```r
R -i data -i data_tod -i genes -i cells -i soupx_groups -o out`

上面代码中的`-i`表示input,我们把上面得到的data,data_tod,genes,cells,soupx_groups输入到R语言中;`-o`表示输出,即我们最后要把想要拿回python环境中的数据存储到out变量中.

SoupX 的输出是校正环境mRNA后的计数矩阵.


```r
%%R -i data -i data_tod -i genes -i cells -i soupx_groups -o out 

## 指定数据的行名和列名
rownames(data) = genes
colnames(data) = cells
## 确保计数表和液滴表的稀疏格式正确
data <- as(data, "sparseMatrix")
data_tod <- as(data_tod, "sparseMatrix")

## 为 SoupX 生成 SoupChannel 对象
sc = SoupChannel(data_tod, data, calcSoupProfile = FALSE)

## 向 SoupChannel 对象添加额外的元数据
soupProf = data.frame(row.names = rownames(data), est = rowSums(data)/sum(data), counts = rowSums(data))
sc = setSoupProfile(sc, soupProf)
# 在 SoupChannel 中设置集群信息
sc = setClusters(sc, soupx_groups)

# 估算污染分数
sc  = autoEstCont(sc, doPlot=FALSE)
# 将更正后的计数矩阵并将更正后矩阵的元素转变为整数
out = adjustCounts(sc, roundToInt = TRUE)
```


    581 genes passed tf-idf cut-off and 323 soup quantile filter.  Taking the top 100.
    Using 460 independent estimates of rho.
    Estimated global rho of 0.06
    Expanding counts from 16 clusters to 16934 cells.
    In addition: Warning message:
    In sparseMatrix(i = out@i[w] + 1, j = out@j[w] + 1, x = out@x[w],  :
      'giveCsparse' is deprecated; setting repr="T" for you
    



```python
## SoupX 输出校正后的计数矩阵，我们将其存储为附加层
adata.layers["counts"] = adata.X
adata.layers["soupX_counts"] = out.T
## 我们将.X用 soupX 校正后的计数矩阵覆盖
adata.X = adata.layers["soupX_counts"]
```


```python
print(f"Total number of genes: {adata.n_vars}")

# 我们另外过滤掉在细胞中未检测到的基因(至少 20 个)，因为这些基因没有提供信息
sc.pp.filter_genes(adata, min_cells=20)
print(f"Number of genes after cell filter: {adata.n_vars}")
```

    Total number of genes: 36601
    Number of genes after cell filter: 20755
    


```python

```

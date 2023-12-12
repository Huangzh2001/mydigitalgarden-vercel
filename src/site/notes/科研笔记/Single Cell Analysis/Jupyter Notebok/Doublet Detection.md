---
{"dg-publish":true,"tags":["Notes","Single-Cell-Analysis","jupyter"],"permalink":"/科研笔记/Single Cell Analysis/Jupyter Notebok/Doublet Detection/","dgPassFrontmatter":true}
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

# 2.检验双联体

## (1)导入包与初始化包


```python
import os
os.environ['R_HOME'] = "D:/Software/R-4.3.1" # 这个路径可以在RStudio中输入命令:Sys.getenv('R_HOME')得到
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

    The rpy2.ipython extension is already loaded. To reload it, use:
      %reload_ext rpy2.ipython
    


```r
R -i data_mat -o doublet_score -o doublet_class

set.seed(123)
sce = scDblFinder(
    SingleCellExperiment(
        list(counts=data_mat),
    ) 
)

doublet_score = sce$scDblFinder.score
doublet_class = sce$scDblFinder.class
```


    Creating ~13548 artificial doublets...
    Dimensional reduction
    Evaluating kNN...
    Training model...
    iter=0, 4464 cells excluded from training.
    iter=1, 4635 cells excluded from training.
    iter=2, 4687 cells excluded from training.
    Threshold found:0.46
    3711 (21.9%) doublets called
    In addition: Warning message:
    In .checkSCE(sce) :
      Some cells in `sce` have an extremely low read counts; note that these could trigger errors and might best be filtered out
    



```python
## 我们将得到的双联体得分和分类结果添加到adata.obs的注释中
adata.obs["scDblFinder_score"] = doublet_score
adata.obs["scDblFinder_class"] = doublet_class
# 我们可以查看有多少个单体和双联体
adata.obs.scDblFinder_class.value_counts()
```




    scDblFinder_class
    singlet    13223
    doublet     3711
    Name: count, dtype: int64




```python

```

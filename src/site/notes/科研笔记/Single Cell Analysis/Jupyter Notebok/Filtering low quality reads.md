---
{"dg-publish":true,"tags":["Notes","Single-Cell-Analysis","jupyter"],"permalink":"/科研笔记/Single Cell Analysis/Jupyter Notebok/Filtering low quality reads/","dgPassFrontmatter":true}
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

# 2.过滤低质量细胞

## (1)import包
scipy中的median_abs_deviation是用来计算中值绝对偏差的.


```python
from scipy.stats import median_abs_deviation
```

## (2)计算三个QC指标

这里的三个QC指标指的是:
1. 每个条形码(barcode)的计数数量（计数深度）
2. 每个条形码(barcode)的基因数量
3. 每个条形码(barcode)的线粒体基因计数分数

首先,我们根据基因编号开头的字母来给每个基因添加是否是线粒体基因,核糖体基因,血红蛋白基因的注释.
如果是则为Ture,否则为False.


```python
## 1.给基因添加是否是线粒体基因的注释
adata.var["mt"] = adata.var_names.str.startswith("MT-")
## 2.给基因添加是否是核糖体基因的注释
adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
## 3.给基因是否是血红蛋白基因进行注释
adata.var["hb"] = adata.var_names.str.contains(("^HB[^(P)]"))
```


```python
adata.var["mt"]
```




    MIR1302-2HG    False
    FAM138A        False
    OR4F5          False
    AL627309.1     False
    AL627309.3     False
                   ...  
    AC141272.1     False
    AC023491.2     False
    AC007325.1     False
    AC007325.4     False
    AC007325.2     False
    Name: mt, Length: 36601, dtype: bool



使用scanpy中的函数`sc.pp.calculate_qc_metrics`计算相应的QC指标.

计算得到的QC指标为一个dataframe矩阵,其中行指标为细胞编号;其中比较重要的列指标为:
* pct_counts_in_top_20_genes:每个细胞中前20个基因的表达量和占总的基因中表达量的百分比
* pct_counts_mt:每个细胞线粒体基因总表达量占所有基因表达量总和的百分比
* pct_counts_ribo:每个细胞核糖体基因总表达量占所有基因表达量总和的百分比
* pct_counts_hb:每个细胞血红蛋白基因总表达量占所有基因表达量总和的百分比
* total_counts:每个细胞的总的基因表达量
* n_genes_by_counts:每个细胞的基因表达量>0的基因数量

根据之前的源码讲解我们知道,当`inplace=True`的时候,计算得到的QC指标将被注释到adata中.


```python
sc.pp.calculate_qc_metrics(
    adata, qc_vars=["mt","ribo","hb"], inplace=True, percent_top=[20], log1p=True
)
adata
```




    AnnData object with n_obs × n_vars = 16934 × 36601
        obs: 'n_genes_by_counts', 'log1p_n_genes_by_counts', 'total_counts', 'log1p_total_counts', 'pct_counts_in_top_20_genes', 'total_counts_mt', 'log1p_total_counts_mt', 'pct_counts_mt', 'total_counts_ribo', 'log1p_total_counts_ribo', 'pct_counts_ribo', 'total_counts_hb', 'log1p_total_counts_hb', 'pct_counts_hb'
        var: 'gene_ids', 'feature_types', 'genome', 'mt', 'ribo', 'hb', 'n_cells_by_counts', 'mean_counts', 'log1p_mean_counts', 'pct_dropout_by_counts', 'total_counts', 'log1p_total_counts'



在上面的注释中,n_genes_by_counts, total_counts, pct_counts_mt为三个我们关心的 QC 指标.
它们分别代表了细胞检测到的基因数量较少、计数深度较低且线粒体计数较高.
下面我们做出三张与它们有关的图.

我们利用sns中的画图函数画出每个细胞的总的基因表达量的直方图.横轴代表基因表达数量,纵轴代表该基因表达数量的细胞个数;

使用scanpy中的函数sc.pl.violin画"每个细胞线粒体基因总表达量占所有基因表达量总和的百分比"的小提琴图;

画出以"每个细胞的总的基因表达量"为横坐标;以"每个细胞的基因表达量>0的基因数量"为纵坐标的散点图.


```python
p1 = sns.displot(adata.obs["total_counts"],bins=100,kde=False)
p2 = sc.pl.violin(adata, "pct_counts_mt")
p3 = sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt")
```

    D:\Software\anaconda3\envs\scverse\lib\site-packages\seaborn\_oldcore.py:1498: FutureWarning: is_categorical_dtype is deprecated and will be removed in a future version. Use isinstance(dtype, CategoricalDtype) instead
      if pd.api.types.is_categorical_dtype(vector):
    D:\Software\anaconda3\envs\scverse\lib\site-packages\seaborn\_oldcore.py:1119: FutureWarning: use_inf_as_na option is deprecated and will be removed in a future version. Convert inf values to NaN before operating instead.
      with pd.option_context('mode.use_inf_as_na', True):
    D:\Software\anaconda3\envs\scverse\lib\site-packages\seaborn\_oldcore.py:1498: FutureWarning: is_categorical_dtype is deprecated and will be removed in a future version. Use isinstance(dtype, CategoricalDtype) instead
      if pd.api.types.is_categorical_dtype(vector):
    D:\Software\anaconda3\envs\scverse\lib\site-packages\seaborn\_oldcore.py:1498: FutureWarning: is_categorical_dtype is deprecated and will be removed in a future version. Use isinstance(dtype, CategoricalDtype) instead
      if pd.api.types.is_categorical_dtype(vector):
    D:\Software\anaconda3\envs\scverse\lib\site-packages\seaborn\_oldcore.py:1498: FutureWarning: is_categorical_dtype is deprecated and will be removed in a future version. Use isinstance(dtype, CategoricalDtype) instead
      if pd.api.types.is_categorical_dtype(vector):
    D:\Software\anaconda3\envs\scverse\lib\site-packages\seaborn\_oldcore.py:1498: FutureWarning: is_categorical_dtype is deprecated and will be removed in a future version. Use isinstance(dtype, CategoricalDtype) instead
      if pd.api.types.is_categorical_dtype(vector):
    D:\Software\anaconda3\envs\scverse\lib\site-packages\seaborn\_oldcore.py:1498: FutureWarning: is_categorical_dtype is deprecated and will be removed in a future version. Use isinstance(dtype, CategoricalDtype) instead
      if pd.api.types.is_categorical_dtype(vector):
    D:\Software\anaconda3\envs\scverse\lib\site-packages\seaborn\_oldcore.py:1498: FutureWarning: is_categorical_dtype is deprecated and will be removed in a future version. Use isinstance(dtype, CategoricalDtype) instead
      if pd.api.types.is_categorical_dtype(vector):
    D:\Software\anaconda3\envs\scverse\lib\site-packages\seaborn\_oldcore.py:1498: FutureWarning: is_categorical_dtype is deprecated and will be removed in a future version. Use isinstance(dtype, CategoricalDtype) instead
      if pd.api.types.is_categorical_dtype(vector):
    D:\Software\anaconda3\envs\scverse\lib\site-packages\seaborn\_oldcore.py:1498: FutureWarning: is_categorical_dtype is deprecated and will be removed in a future version. Use isinstance(dtype, CategoricalDtype) instead
      if pd.api.types.is_categorical_dtype(vector):
    D:\Software\anaconda3\envs\scverse\lib\site-packages\seaborn\_oldcore.py:1498: FutureWarning: is_categorical_dtype is deprecated and will be removed in a future version. Use isinstance(dtype, CategoricalDtype) instead
      if pd.api.types.is_categorical_dtype(vector):
    D:\Software\anaconda3\envs\scverse\lib\site-packages\seaborn\_oldcore.py:1498: FutureWarning: is_categorical_dtype is deprecated and will be removed in a future version. Use isinstance(dtype, CategoricalDtype) instead
      if pd.api.types.is_categorical_dtype(vector):
    D:\Software\anaconda3\envs\scverse\lib\site-packages\seaborn\_oldcore.py:1119: FutureWarning: use_inf_as_na option is deprecated and will be removed in a future version. Convert inf values to NaN before operating instead.
      with pd.option_context('mode.use_inf_as_na', True):
    D:\Software\anaconda3\envs\scverse\lib\site-packages\seaborn\_oldcore.py:1119: FutureWarning: use_inf_as_na option is deprecated and will be removed in a future version. Convert inf values to NaN before operating instead.
      with pd.option_context('mode.use_inf_as_na', True):
    


    
![png](/img/user/科研笔记/Single Cell Analysis/Jupyter Notebok/Filtering low quality reads_files/Filtering low quality reads_20_1.png)
    



    
![png](/img/user/科研笔记/Single Cell Analysis/Jupyter Notebok/Filtering low quality reads_files/Filtering low quality reads_20_2.png)
    


    D:\Software\anaconda3\envs\scverse\lib\site-packages\scanpy\plotting\_utils.py:715: FutureWarning: Series.__getitem__ treating keys as positions is deprecated. In a future version, integer keys will always be treated as labels (consistent with DataFrame behavior). To access a value by position, use `ser.iloc[pos]`
      color = color[sort]
    


    
![png](/img/user/科研笔记/Single Cell Analysis/Jupyter Notebok/Filtering low quality reads_files/Filtering low quality reads_20_4.png)
    


## (3)过滤低质量细胞

下面我们定义一个函数`is_outlier()`来标记细胞的中值绝对偏差是否大于 MAD 的指定倍数

`is_outlier()`的输入参数为`adata`;`metric`为我们想要算中值绝对偏差的指标;`nmads`为指定MAD的倍数.

$$
|M-np.median(M)| > nmads * median\_abs\_deviation(M)
$$
也就是表示与中值的差大于 nmads 倍中值绝对偏差.

值得注意的是:
* (scipy中的median_abs_deviation是用来计算中值绝对偏差的.)
* `|`是位运算符.比如a = 1, b = 0, 则c = a | b = 1.
* 所以我们用`(M < np.median(M) - nmads * median_abs_deviation(M))`就得到以一个全由0,1组成的numpy数组,其中1代表$np.median(M)-M>nmads * median\_abs\_deviation(M)$
* 同理`(np.median(M) + nmads * median_abs_deviation(M) < M)`也得到以一个全由0,1组成的numpy数组,其中1代表$M-np.median(M)>nmads * median\_abs\_deviation(M)$
* 于是`(M < np.median(M) - nmads * median_abs_deviation(M)) | (np.median(M) + nmads * median_abs_deviation(M) < M)`得到的就是一个全由0,1组成的numpy数组,其中1代表$|M-np.median(M)| > nmads * median\_abs\_deviation(M)$,0代表$|M-np.median(M)| \leq nmads * median\_abs\_deviation(M)$


```python
def is_outlier(adata, metric: str, nmads: int):
    M = adata.obs[metric]
    outlier = (M < np.median(M) - nmads * median_abs_deviation(M)) | (
        np.median(M) + nmads * median_abs_deviation(M) < M
    )
    return outlier
```


```python
adata.obs["outlier"] = (
    is_outlier(adata, "log1p_total_counts", 5)
    | is_outlier(adata, "log1p_n_genes_by_counts", 5)
    | is_outlier(adata, "pct_counts_in_top_20_genes", 5)
)
adata.obs.outlier.value_counts()
```




    outlier
    False    16065
    True       869
    Name: count, dtype: int64




```python
adata.obs["mt_outlier"] = is_outlier(adata, "pct_counts_mt", 3) | (
    adata.obs["pct_counts_mt"] > 8
)
adata.obs.mt_outlier.value_counts()
```




    mt_outlier
    False    15240
    True      1694
    Name: count, dtype: int64




```python
print(f"Total number of cells: {adata.n_obs}")
adata = adata[(~adata.obs.outlier) & (~adata.obs.mt_outlier)].copy()

print(f"Number of cells after filtering of low quality cells: {adata.n_obs}")
```

    Total number of cells: 16934
    Number of cells after filtering of low quality cells: 14814
    


```python
p1 = sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt")
```

    D:\Software\anaconda3\envs\scverse\lib\site-packages\scanpy\plotting\_utils.py:715: FutureWarning: Series.__getitem__ treating keys as positions is deprecated. In a future version, integer keys will always be treated as labels (consistent with DataFrame behavior). To access a value by position, use `ser.iloc[pos]`
      color = color[sort]
    


    
![png](/img/user/科研笔记/Single Cell Analysis/Jupyter Notebok/Filtering low quality reads_files/Filtering low quality reads_26_1.png)
    



```python

```

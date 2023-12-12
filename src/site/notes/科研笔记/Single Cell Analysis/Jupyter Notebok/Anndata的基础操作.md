---
{"dg-publish":true,"tags":["Notes","Single-Cell-Analysis","jupyter"],"permalink":"/科研笔记/Single Cell Analysis/Jupyter Notebok/Anndata的基础操作/","dgPassFrontmatter":true}
---

## 1.首先我们导入包


```python
import numpy as np
import pandas as pd
import anndata as ad
from scipy.sparse import csr_matrix # 这个包是用来把矩阵存储为稀疏矩阵的
```

## 2.利用 AnnData 储存计数矩阵

我们使用随机泊松分布数据初始化一个 100*2000 的矩阵来模拟我们得到的计数矩阵.我们知道计数矩阵中最重要的数据有三个:细胞标识;基因标识;每个细胞的基因表达量.这个矩阵就表示有 100 个细胞, 2000 个基因,矩阵的(i,j)元就表示第i个细胞中第j个基因的表达量.

我们首先将其转化为稀疏矩阵储存(anndata支持稀疏矩阵导导入).

然后再导入到anndata中,生成一个anndata对象adata.


```python
counts = csr_matrix(np.random.poisson(1, size=(100,2000)), dtype=np.float32)
# csr_matrix的作用是将生成的矩阵转化为稀疏矩阵存储
adata = ad.AnnData(counts)
adata
```




    AnnData object with n_obs × n_vars = 100 × 2000



### (1)查看细胞标识

我们可以通过`adata.n_obs`查看细胞标识的个数(即有多少个细胞).

还可以通过`adata.obs_names`获取细胞标识,它是一个列表.


```python
print(len(adata.obs_names))
print(adata.n_obs)
adata.obs_names[:10]
```

    100
    100
    




    Index(['0', '1', '2', '3', '4', '5', '6', '7', '8', '9'], dtype='object')



为了好看,我们在每个细胞标识前面加上`Cell_`


```python
adata.obs_names = [f"Cell_{i:d}" for i in range(adata.n_obs)]
print(adata.obs_names[:10])
```

    Index(['Cell_0', 'Cell_1', 'Cell_2', 'Cell_3', 'Cell_4', 'Cell_5', 'Cell_6',
           'Cell_7', 'Cell_8', 'Cell_9'],
          dtype='object')
    

### (2)查看基因标识

我们可以通过`adata.n_vars`查看基因标识的个数(即有多少个基因).

还可以通过`adata.var_names`获取基因标识,它也是一个列表.


```python
print(len(adata.var_names))
print(adata.n_vars)
adata.var_names[:10]
```

    2000
    2000
    




    Index(['0', '1', '2', '3', '4', '5', '6', '7', '8', '9'], dtype='object')



同样地,我们给每个基因标识前面加上`Gene_`


```python
adata.var_names = [f"Gene_{i:d}" for i in range(adata.n_vars)]
print(adata.var_names[:10])
```

    Index(['Gene_0', 'Gene_1', 'Gene_2', 'Gene_3', 'Gene_4', 'Gene_5', 'Gene_6',
           'Gene_7', 'Gene_8', 'Gene_9'],
          dtype='object')
    

### (3)查看基因表达量

我们可以通过`adata.X`获取adata中的基因表达数据.

但由于它是用稀疏矩阵储存的,我们可以通过`adata.to_df()`命令将其转换为pandas的dataframe数据进行读取


```python
adata.X
```




    <100x2000 sparse matrix of type '<class 'numpy.float32'>'
    	with 126329 stored elements in Compressed Sparse Row format>




```python
df = adata.to_df()
print(df)
```

             Gene_0  Gene_1  Gene_2  Gene_3  Gene_4  Gene_5  Gene_6  Gene_7  \
    Cell_0      1.0     1.0     2.0     3.0     0.0     3.0     3.0     0.0   
    Cell_1      0.0     2.0     2.0     1.0     1.0     3.0     0.0     2.0   
    Cell_2      1.0     3.0     2.0     2.0     1.0     2.0     0.0     3.0   
    Cell_3      1.0     0.0     1.0     2.0     0.0     0.0     3.0     0.0   
    Cell_4      2.0     2.0     2.0     1.0     0.0     1.0     1.0     1.0   
    ...         ...     ...     ...     ...     ...     ...     ...     ...   
    Cell_95     4.0     1.0     1.0     1.0     1.0     0.0     1.0     1.0   
    Cell_96     2.0     1.0     4.0     1.0     0.0     0.0     0.0     0.0   
    Cell_97     1.0     1.0     0.0     4.0     0.0     2.0     2.0     2.0   
    Cell_98     2.0     0.0     2.0     0.0     0.0     1.0     0.0     1.0   
    Cell_99     0.0     0.0     2.0     1.0     0.0     2.0     1.0     0.0   
    
             Gene_8  Gene_9  ...  Gene_1990  Gene_1991  Gene_1992  Gene_1993  \
    Cell_0      1.0     1.0  ...        3.0        0.0        2.0        2.0   
    Cell_1      0.0     1.0  ...        1.0        1.0        1.0        1.0   
    Cell_2      0.0     2.0  ...        1.0        1.0        0.0        1.0   
    Cell_3      2.0     1.0  ...        1.0        4.0        0.0        3.0   
    Cell_4      0.0     0.0  ...        2.0        0.0        2.0        3.0   
    ...         ...     ...  ...        ...        ...        ...        ...   
    Cell_95     0.0     0.0  ...        0.0        0.0        2.0        4.0   
    Cell_96     0.0     1.0  ...        1.0        1.0        2.0        1.0   
    Cell_97     1.0     2.0  ...        2.0        1.0        0.0        3.0   
    Cell_98     0.0     1.0  ...        0.0        1.0        0.0        4.0   
    Cell_99     0.0     0.0  ...        1.0        1.0        1.0        0.0   
    
             Gene_1994  Gene_1995  Gene_1996  Gene_1997  Gene_1998  Gene_1999  
    Cell_0         1.0        1.0        1.0        0.0        2.0        2.0  
    Cell_1         0.0        0.0        2.0        1.0        0.0        0.0  
    Cell_2         0.0        1.0        3.0        0.0        3.0        0.0  
    Cell_3         1.0        1.0        0.0        1.0        2.0        0.0  
    Cell_4         0.0        1.0        2.0        0.0        1.0        0.0  
    ...            ...        ...        ...        ...        ...        ...  
    Cell_95        0.0        1.0        1.0        1.0        0.0        3.0  
    Cell_96        3.0        2.0        1.0        0.0        1.0        0.0  
    Cell_97        1.0        2.0        2.0        1.0        2.0        0.0  
    Cell_98        1.0        1.0        0.0        2.0        0.0        1.0  
    Cell_99        1.0        1.0        0.0        2.0        3.0        3.0  
    
    [100 rows x 2000 columns]
    

得到了dataframe数据之后,我们可以读取特定细胞的特定基因表达量


```python
## 1.读取单个细胞的单个基因表达量
print(df['Gene_100']['Cell_11'])

## 2.读取一个基因的所有数据
# print(df['Gene_1'])
## 或者用
# print(df.Gene_1)

## 3.读取一个细胞的所有数据
## 获取DataFrame中的一行数据时，不能直接用 data['行索引'] 或 data.行索引的方式
## 获取行数据也有两种方式，需要借助loc属性或iloc属性。
## loc属性基于行索引名获取数据，用法为 data.loc['行索引']
print(df.loc['Cell_1'])
## iloc属性基于数值索引获取数据，用法为 data.iloc[数值]
print(df.iloc[5])

## 4.切片
## 索引方法
print(df.loc[['Cell_1','Cell_6'],['Gene_5','Gene_23']])
## 数字方法
print(df.iloc[:4,2:4])
## 我们发现索引方法直观但无法索引一片,所以我们要用到如下格式转换
a = df.columns.get_indexer(['Gene_20'])
b = df.columns.get_indexer(['Gene_30'])
c = df.index.get_indexer(['Cell_5'])
d = df.index.get_indexer(['Cell_9'])
print(c)
print(df.iloc[c[0]:d[0], a[0]:b[0]])
```

    0.0
    Gene_0       0.0
    Gene_1       2.0
    Gene_2       2.0
    Gene_3       1.0
    Gene_4       1.0
                ... 
    Gene_1995    0.0
    Gene_1996    2.0
    Gene_1997    1.0
    Gene_1998    0.0
    Gene_1999    0.0
    Name: Cell_1, Length: 2000, dtype: float32
    Gene_0       0.0
    Gene_1       0.0
    Gene_2       1.0
    Gene_3       1.0
    Gene_4       1.0
                ... 
    Gene_1995    1.0
    Gene_1996    0.0
    Gene_1997    2.0
    Gene_1998    1.0
    Gene_1999    2.0
    Name: Cell_5, Length: 2000, dtype: float32
            Gene_5  Gene_23
    Cell_1     3.0      0.0
    Cell_6     0.0      2.0
            Gene_2  Gene_3
    Cell_0     2.0     3.0
    Cell_1     2.0     1.0
    Cell_2     2.0     2.0
    Cell_3     1.0     2.0
    [5]
            Gene_20  Gene_21  Gene_22  Gene_23  Gene_24  Gene_25  Gene_26  \
    Cell_5      1.0      1.0      1.0      0.0      0.0      3.0      1.0   
    Cell_6      1.0      0.0      0.0      2.0      1.0      1.0      1.0   
    Cell_7      2.0      0.0      1.0      0.0      2.0      0.0      1.0   
    Cell_8      1.0      1.0      0.0      1.0      1.0      1.0      1.0   
    
            Gene_27  Gene_28  Gene_29  
    Cell_5      2.0      0.0      1.0  
    Cell_6      1.0      1.0      0.0  
    Cell_7      1.0      0.0      1.0  
    Cell_8      0.0      2.0      1.0  
    

## 3.添加注释

### (1)添加一维注释

我们可以对计数矩阵中的细胞和基因添加注释,比如下面标记各个细胞的细胞类型.

当然,我们也可以给基因添加注释.


```python
ct = np.random.choice(["B", "T", "Monocyte"], size=(adata.n_obs,))
adata.obs["cell_type"] = pd.Categorical(ct)
# 当数据只有少数几种可能取值但有大量重复字符串字段，利用categorical类型可有效节省内存，提高数据分析的效率
adata.obs
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>cell_type</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>Cell_0</th>
      <td>Monocyte</td>
    </tr>
    <tr>
      <th>Cell_1</th>
      <td>Monocyte</td>
    </tr>
    <tr>
      <th>Cell_2</th>
      <td>B</td>
    </tr>
    <tr>
      <th>Cell_3</th>
      <td>B</td>
    </tr>
    <tr>
      <th>Cell_4</th>
      <td>B</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
    </tr>
    <tr>
      <th>Cell_95</th>
      <td>T</td>
    </tr>
    <tr>
      <th>Cell_96</th>
      <td>T</td>
    </tr>
    <tr>
      <th>Cell_97</th>
      <td>B</td>
    </tr>
    <tr>
      <th>Cell_98</th>
      <td>T</td>
    </tr>
    <tr>
      <th>Cell_99</th>
      <td>Monocyte</td>
    </tr>
  </tbody>
</table>
<p>100 rows × 1 columns</p>
</div>



我们可以提取出特定标记的细胞.


```python
bdata = adata[adata.obs.cell_type == "B"]
bdata.obs
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>cell_type</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>Cell_2</th>
      <td>B</td>
    </tr>
    <tr>
      <th>Cell_3</th>
      <td>B</td>
    </tr>
    <tr>
      <th>Cell_4</th>
      <td>B</td>
    </tr>
    <tr>
      <th>Cell_6</th>
      <td>B</td>
    </tr>
    <tr>
      <th>Cell_9</th>
      <td>B</td>
    </tr>
    <tr>
      <th>Cell_12</th>
      <td>B</td>
    </tr>
    <tr>
      <th>Cell_15</th>
      <td>B</td>
    </tr>
    <tr>
      <th>Cell_21</th>
      <td>B</td>
    </tr>
    <tr>
      <th>Cell_25</th>
      <td>B</td>
    </tr>
    <tr>
      <th>Cell_33</th>
      <td>B</td>
    </tr>
    <tr>
      <th>Cell_36</th>
      <td>B</td>
    </tr>
    <tr>
      <th>Cell_39</th>
      <td>B</td>
    </tr>
    <tr>
      <th>Cell_40</th>
      <td>B</td>
    </tr>
    <tr>
      <th>Cell_44</th>
      <td>B</td>
    </tr>
    <tr>
      <th>Cell_47</th>
      <td>B</td>
    </tr>
    <tr>
      <th>Cell_49</th>
      <td>B</td>
    </tr>
    <tr>
      <th>Cell_51</th>
      <td>B</td>
    </tr>
    <tr>
      <th>Cell_54</th>
      <td>B</td>
    </tr>
    <tr>
      <th>Cell_59</th>
      <td>B</td>
    </tr>
    <tr>
      <th>Cell_64</th>
      <td>B</td>
    </tr>
    <tr>
      <th>Cell_66</th>
      <td>B</td>
    </tr>
    <tr>
      <th>Cell_69</th>
      <td>B</td>
    </tr>
    <tr>
      <th>Cell_72</th>
      <td>B</td>
    </tr>
    <tr>
      <th>Cell_73</th>
      <td>B</td>
    </tr>
    <tr>
      <th>Cell_75</th>
      <td>B</td>
    </tr>
    <tr>
      <th>Cell_76</th>
      <td>B</td>
    </tr>
    <tr>
      <th>Cell_77</th>
      <td>B</td>
    </tr>
    <tr>
      <th>Cell_84</th>
      <td>B</td>
    </tr>
    <tr>
      <th>Cell_85</th>
      <td>B</td>
    </tr>
    <tr>
      <th>Cell_88</th>
      <td>B</td>
    </tr>
    <tr>
      <th>Cell_97</th>
      <td>B</td>
    </tr>
  </tbody>
</table>
</div>



在我们添加完注释之后,adata的结构就会发生相应的改变.


```python
adata
```




    AnnData object with n_obs × n_vars = 100 × 2000
        obs: 'cell_type'



### (2)添加二维以上注释

上面我们给细胞(或基因)添加的是一维的注释,我们当然也可以给细胞添加二维及以上的注释.

此时我们要用到 AnnData 中的`.obsm/.varm`属性.

但是要注意:
* `.obsm`矩阵的长度必须等于`.n_obs`(即细胞的数目)
* `.varm`矩阵的长度必须等于`.n_vars`(即基因的个数)

这是自然的,因为我们现在做的就是给每个细胞(或基因)进行注释.

比如我们要给每个细胞进行UMAP的注释.UMAP可以理解为给每个细胞在二维坐标系上标记一个坐标,从而我们可以很容易用图看出哪些细胞是相似的,哪些细胞是不同的.
更详细介绍可以看这里:

https://www.zhihu.com/question/455525769/answer/3138322711?utm_id=0


```python
adata.obsm["X_umap"] = np.random.normal(0, 1, size=(adata.n_obs, 2))
print(type(adata.obsm))
```

    <class 'anndata._core.aligned_mapping.AxisArrays'>
    

相应地,adata会进行更新


```python
adata
```




    AnnData object with n_obs × n_vars = 100 × 2000
        obs: 'cell_type'
        obsm: 'X_umap'



### (3)添加任意注释


```python

```

## 4.层数

最后，我们可以储存不同形式的原始数据，这些可以存储在 AnnData 的不同层中。比如，让我们对原始数据进行log转换并将其存储在新的一层中。


```python
adata.layers["log_transformed"] = np.log1p(adata.X)
adata
```




    AnnData object with n_obs × n_vars = 100 × 2000
        obs: 'cell_type'
        obsm: 'X_umap'
        layers: 'log_transformed'



我们可以像下面这样来读取不同层的数据为dataframe.


```python
adata.to_df(layer="log_transformed")
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>Gene_0</th>
      <th>Gene_1</th>
      <th>Gene_2</th>
      <th>Gene_3</th>
      <th>Gene_4</th>
      <th>Gene_5</th>
      <th>Gene_6</th>
      <th>Gene_7</th>
      <th>Gene_8</th>
      <th>Gene_9</th>
      <th>...</th>
      <th>Gene_1990</th>
      <th>Gene_1991</th>
      <th>Gene_1992</th>
      <th>Gene_1993</th>
      <th>Gene_1994</th>
      <th>Gene_1995</th>
      <th>Gene_1996</th>
      <th>Gene_1997</th>
      <th>Gene_1998</th>
      <th>Gene_1999</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>Cell_0</th>
      <td>0.693147</td>
      <td>0.693147</td>
      <td>1.098612</td>
      <td>1.386294</td>
      <td>0.000000</td>
      <td>1.386294</td>
      <td>1.386294</td>
      <td>0.000000</td>
      <td>0.693147</td>
      <td>0.693147</td>
      <td>...</td>
      <td>1.386294</td>
      <td>0.000000</td>
      <td>1.098612</td>
      <td>1.098612</td>
      <td>0.693147</td>
      <td>0.693147</td>
      <td>0.693147</td>
      <td>0.000000</td>
      <td>1.098612</td>
      <td>1.098612</td>
    </tr>
    <tr>
      <th>Cell_1</th>
      <td>0.000000</td>
      <td>1.098612</td>
      <td>1.098612</td>
      <td>0.693147</td>
      <td>0.693147</td>
      <td>1.386294</td>
      <td>0.000000</td>
      <td>1.098612</td>
      <td>0.000000</td>
      <td>0.693147</td>
      <td>...</td>
      <td>0.693147</td>
      <td>0.693147</td>
      <td>0.693147</td>
      <td>0.693147</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>1.098612</td>
      <td>0.693147</td>
      <td>0.000000</td>
      <td>0.000000</td>
    </tr>
    <tr>
      <th>Cell_2</th>
      <td>0.693147</td>
      <td>1.386294</td>
      <td>1.098612</td>
      <td>1.098612</td>
      <td>0.693147</td>
      <td>1.098612</td>
      <td>0.000000</td>
      <td>1.386294</td>
      <td>0.000000</td>
      <td>1.098612</td>
      <td>...</td>
      <td>0.693147</td>
      <td>0.693147</td>
      <td>0.000000</td>
      <td>0.693147</td>
      <td>0.000000</td>
      <td>0.693147</td>
      <td>1.386294</td>
      <td>0.000000</td>
      <td>1.386294</td>
      <td>0.000000</td>
    </tr>
    <tr>
      <th>Cell_3</th>
      <td>0.693147</td>
      <td>0.000000</td>
      <td>0.693147</td>
      <td>1.098612</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>1.386294</td>
      <td>0.000000</td>
      <td>1.098612</td>
      <td>0.693147</td>
      <td>...</td>
      <td>0.693147</td>
      <td>1.609438</td>
      <td>0.000000</td>
      <td>1.386294</td>
      <td>0.693147</td>
      <td>0.693147</td>
      <td>0.000000</td>
      <td>0.693147</td>
      <td>1.098612</td>
      <td>0.000000</td>
    </tr>
    <tr>
      <th>Cell_4</th>
      <td>1.098612</td>
      <td>1.098612</td>
      <td>1.098612</td>
      <td>0.693147</td>
      <td>0.000000</td>
      <td>0.693147</td>
      <td>0.693147</td>
      <td>0.693147</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>...</td>
      <td>1.098612</td>
      <td>0.000000</td>
      <td>1.098612</td>
      <td>1.386294</td>
      <td>0.000000</td>
      <td>0.693147</td>
      <td>1.098612</td>
      <td>0.000000</td>
      <td>0.693147</td>
      <td>0.000000</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>Cell_95</th>
      <td>1.609438</td>
      <td>0.693147</td>
      <td>0.693147</td>
      <td>0.693147</td>
      <td>0.693147</td>
      <td>0.000000</td>
      <td>0.693147</td>
      <td>0.693147</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>...</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>1.098612</td>
      <td>1.609438</td>
      <td>0.000000</td>
      <td>0.693147</td>
      <td>0.693147</td>
      <td>0.693147</td>
      <td>0.000000</td>
      <td>1.386294</td>
    </tr>
    <tr>
      <th>Cell_96</th>
      <td>1.098612</td>
      <td>0.693147</td>
      <td>1.609438</td>
      <td>0.693147</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.693147</td>
      <td>...</td>
      <td>0.693147</td>
      <td>0.693147</td>
      <td>1.098612</td>
      <td>0.693147</td>
      <td>1.386294</td>
      <td>1.098612</td>
      <td>0.693147</td>
      <td>0.000000</td>
      <td>0.693147</td>
      <td>0.000000</td>
    </tr>
    <tr>
      <th>Cell_97</th>
      <td>0.693147</td>
      <td>0.693147</td>
      <td>0.000000</td>
      <td>1.609438</td>
      <td>0.000000</td>
      <td>1.098612</td>
      <td>1.098612</td>
      <td>1.098612</td>
      <td>0.693147</td>
      <td>1.098612</td>
      <td>...</td>
      <td>1.098612</td>
      <td>0.693147</td>
      <td>0.000000</td>
      <td>1.386294</td>
      <td>0.693147</td>
      <td>1.098612</td>
      <td>1.098612</td>
      <td>0.693147</td>
      <td>1.098612</td>
      <td>0.000000</td>
    </tr>
    <tr>
      <th>Cell_98</th>
      <td>1.098612</td>
      <td>0.000000</td>
      <td>1.098612</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.693147</td>
      <td>0.000000</td>
      <td>0.693147</td>
      <td>0.000000</td>
      <td>0.693147</td>
      <td>...</td>
      <td>0.000000</td>
      <td>0.693147</td>
      <td>0.000000</td>
      <td>1.609438</td>
      <td>0.693147</td>
      <td>0.693147</td>
      <td>0.000000</td>
      <td>1.098612</td>
      <td>0.000000</td>
      <td>0.693147</td>
    </tr>
    <tr>
      <th>Cell_99</th>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>1.098612</td>
      <td>0.693147</td>
      <td>0.000000</td>
      <td>1.098612</td>
      <td>0.693147</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>...</td>
      <td>0.693147</td>
      <td>0.693147</td>
      <td>0.693147</td>
      <td>0.000000</td>
      <td>0.693147</td>
      <td>0.693147</td>
      <td>0.000000</td>
      <td>1.098612</td>
      <td>1.386294</td>
      <td>1.386294</td>
    </tr>
  </tbody>
</table>
<p>100 rows × 2000 columns</p>
</div>



## 5.AnnData 对象的读写

AnnData 的存储形式为：`h5ad`.


```python
adata.write("my_results.h5ad", compression="gzip")
```

然后读回来.


```python
adata_new = ad.read_h5ad("my_results.h5ad")
adata_new
```




    AnnData object with n_obs × n_vars = 100 × 2000
        obs: 'cell_type'
        obsm: 'X_umap'
        layers: 'log_transformed'



## 6.高效的数据访问

这里的知识点我没有整理,请看:https://www.sc-best-practices.org/introduction/analysis_tools.html#efficient-data-access

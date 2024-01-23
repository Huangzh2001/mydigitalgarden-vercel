---
{"dg-publish":true,"tags":["Notes","Single-Cell-Analysis"],"permalink":"/科研笔记/Single Cell Analysis/scverse 框架/","dgPassFrontmatter":true}
---

# 一.工具安装

首先我们在Ananconda中新创建一个虚拟环境:

```
conda create -n scverse python=3.10
```

然后我们激活这个环境:

```
conda activate scverse
```

> [!Tools]
> 我们总共要用到四个工具:
> * `anndata`:用于单峰数据存储
> * `scanpy`:用于单峰数据分析
> * `mudata`:存储多模态数据
> * `muon`:多模态数据分析

下面开始一一安装:

```
conda install -c conda-forge anndata
```

```
conda install -c conda-forge scanpy
```

```
conda install -c conda-forge mudata
```

```
pip install muon
```

最后我们输入

```
conda list
```

查看一下是否安装好了所有的文件.

但是,我在运行`import anndata`代码的时候,出现了错误.

怀疑从 conda 下载的 h5py 包与当前虚拟环境的python版本不匹配.

于是我们可以去网站[h5py · PyPI](https://pypi.org/project/h5py/#history)下载对应版本的whl包.比如我的虚拟环境python版本为python 3.10,于是我下载了"h5py-3.10.0-cp310-cp310-win_amd64.whl".

将Anaconda Prompt转移到下载的文件夹,并输入

```
pip install h5py-3.10.0-cp310-cp310-win_amd64.whl
```

安装完成后就成功解决了这个问题.
# 二.工具简介

## 1.anndata单峰数据存储

[[Anndata的基础操作\|Anndata的基础操作]]
## 2.scanpy单峰数据分析

在后续的项目学习中我们会慢慢学习
## 3.mudata多模态数据存储


## 4.muon多模态数据分析


# 三.其他工具

## 1.在python中应用R语言

我们知道,python和R语言在单细胞数据分析中有不同的优势,我们常常要两者一起来用.但是将数据在两个IDE里面传来传去非常麻烦,于是我们要装一些库,来使得在jupyter notebook中能同时使用python和R语言.

我们要安装:rpy2,anndata2ri

```ad-caution
注意:安装完了之后不能直接使用,你直接import的话会报错.比如说`_Kernel_ Restarting _The_ _kernel_ _for_ xxx._ipynb_ _appears_ _to_ _have_ _died_. _It_ _will_ _restart_ automatically.`

原因是python编辑器此时不知道你电脑里的R语言安装到了哪里.

这时有两个方法:
* 一是在import代码前加上
	```python
	import os
	os.environ['R_HOME'] = "D:/Software/R-4.3.1" # 这个路径可以在RStudio中输入命令:Sys.getenv('R_HOME')得到
	```
* 二是将R语言路径添加到系统变量,一劳永逸:[windows系统中安装、配置ryp2_哔哩哔哩_bilibili](https://www.bilibili.com/video/BV1Vb4y1C7T6/?spm_id_from=333.880.my_history.page.click&vd_source=06168f390bae49c4867767c52a20e87c)
```

使用方法:在使用R语言之前,请先用`python`换回python环境.

```ad-caution
注意:如果你调用电脑里的R语言没有的包,就会报错,所以要事先在RStudio中安装好.
```


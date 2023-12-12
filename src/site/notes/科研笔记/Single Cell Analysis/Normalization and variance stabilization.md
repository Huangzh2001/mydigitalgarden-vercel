---
{"tags":["Notes","Single-Cell-Analysis"],"dg-publish":true,"permalink":"/科研笔记/Single Cell Analysis/Normalization and variance stabilization/","dgPassFrontmatter":true}
---

我们将介绍三种Normalization and variance stabilization的方法,它们分别是:
* 移位对数变换
* scran归一化
* pearson残差的解析逼近

其中,前两种方法基于$\delta$方法,能够将细胞间基因表达量的方差差异控制在一定范围内;而第三种方法能有效消除采样技术导致的偏差.
# 1.方差稳定性

单细胞 RNA 测序 (RNA-seq) 得到的计数矩阵是异方差的.而我们后续所用到的方法大都假定是"同方差的".

因此,我们想要**将所有细胞基因表达量的方差变为是相同的,或是将其方差差异控制在一定范围内**.

我们将介绍两种方差稳定性方法,它们都是基于$\delta$方法的.下面我们先来介绍一下$\delta$方法.
## (1)$\delta$方法

我们有这样一个定理:

```ad-theorem
$\delta$ method

设 $k$维随机向量序列$T_n$ 满足如下的渐近正态性:
$$
\lim_{n\to\infty}\sqrt{n}(T_n-\theta)=\mathcal{N}(\boldsymbol{\mu},\boldsymbol{\Sigma}).
$$

(上面等号是指"依分布收敛")

若$k$ 变量多元向量值函数 $\phi(t)$ 在 $\theta$ 可微

$$
\left.\phi'_\theta=\left(∂ϕ1∂x1⋯∂ϕ1∂xk⋮⋱⋮∂ϕm∂x1⋯∂ϕm∂xk\right.\right)_{\boldsymbol{x}=\boldsymbol{\theta}}
$$

则成立

$$
\lim_{n\to\infty}\sqrt{n}(\phi(\boldsymbol{T}_n)-\phi(\boldsymbol{\theta}))=\mathcal{N}(\boldsymbol{\phi'}_\theta\boldsymbol{\mu},\boldsymbol{\phi'}_\theta\boldsymbol{\Sigma}(\boldsymbol{\phi'}_\theta)^T).
$$

(上面等号是指"依分布收敛")
```

```ad-proposition
一维情形的 $\delta$ 方法

当 $n\to\infty$ 时，若 $\sqrt{n}(T_n-\theta)\overset{w}{\operatorname*{\to}}\mathcal{N}(\mu,\sigma^2)$,又 $\phi(t)$ 为可微函数，则

$$
\lim_{n\to\infty}\sqrt{n}(\phi(T_n)-\phi(\theta))=\mathcal{N}(\phi'(\theta)\mu,[\phi'(\theta)]^2\sigma^2).
$$
```

从一维的$\delta$ 方法可以看出，如果我们将可微函数 $g$ 作用于具有均值 $\mu$ 的随机变量 $X$，则变换后的随机变量 $g(X)$ 的渐进方差可以近似为：

$$
Var[g(X)]\approx |g^{\prime}(\mu)|Var[X],
$$
现在考虑一组随机变量 $X_1, X_2, \cdots$,其方差和均值通过某个函数 $v$ 相互关联，即 $Var[X_{i}]=v(\mu_i)$.

我们找到一个函数$g$满足
$$
g^{\prime}(\mu)=\frac{\mathrm{const.}}{\sqrt{v(\mu)}},
$$
然后通过$g$这个变换我们就能将$Var[g(X_i)] = const.$.这不就达到了我们想要**将所有细胞基因表达量的方差变为是相同**的目的了吗.

```ad-info

从理论上和经验上建立的 UMI 数据模型都是服从 Gamma-Poisson 分布(也称为负二项分布),其方差满足
$$
Var[Y] = \mu + \alpha \mu^2
$$
即满足我们上面"方差和均值通过某个函数 $v$ 相互关联"的条件.
```

## (2)方差稳定性方法

### 移位对数变换

通过上面的讨论我们知道,可以适当地选取变换,使得通过该变换来将方差差异控制在一定的范围内.

经过研究发现,移位对数变换具有很好的效果.

假设$y$是原始计数矩阵上的值,则我们对$y$做如下变换:
$$
f(y) = \log(\frac{y}{s}+y_0)
$$
其中$s$称为尺寸因子;$y_0$称为伪计数(伪计数是为了防止log中的值为0).其中尺寸分子是跟细胞有关的,细胞$c$的尺寸分子$s_c$可以计算如下:
$$
s_c = \frac{\sum_g y_{gc}}{L}
$$
其中$g$代表不同基因;$y_{gc}$表示细胞$c$中$g$基因的表达量;$L$的值就很多样了,在 `scanpy` 默认值是数据集中原始计数深度的中位数,而许多分析模板使​​用固定得$L$值， 例如$L=10^5$， 或者$L=10^6$.

### scran归一化

Scran 归一化使用的也是移位对数变换,差异只是在于估计尺寸分子$s$的方法不同.

Scran 利用反卷积方法根据细胞池基因的线性回归来估计大小因子,能更好地解释数据集中所有细胞计数深度的差异.

Scran 需要聚类输入来提高估计尺寸因子的性能.

# 2.pearson残差的解析逼近


<div class="transclusion internal-embed is-loaded"><div class="markdown-embed">



## 🎈方法

我们观测到的细胞基因表达的计数矩阵实际上是受两个因素影响的:第一个因素是细胞本身的差异性表达;第二个因素是<span style="color: rgb(51, 51, 51)"><span style="background-color: rgb(255, 255, 255)">采样技术导致的.</span></span>

我们希望找到一个方法来挑选出受技术效应影响较小的数据.

### <span class="highlight" data-annotation="%7B%22attachmentURI%22%3A%22http%3A%2F%2Fzotero.org%2Fusers%2F13042945%2Fitems%2F5IJ4QJUI%22%2C%22pageLabel%22%3A%223%22%2C%22position%22%3A%7B%22pageIndex%22%3A2%2C%22rects%22%3A%5B%5B102.635%2C733.795%2C443.27%2C744.411%5D%5D%7D%2C%22citationItem%22%3A%7B%22uris%22%3A%5B%22http%3A%2F%2Fzotero.org%2Fusers%2F13042945%2Fitems%2FKIBCEQBX%22%5D%2C%22locator%22%3A%223%22%7D%7D" ztype="zhighlight"><a href="zotero://open-pdf/library/items/5IJ4QJUI?page=3">“Analytic Pearson residuals for normalization of UMI data”</a></span> <span class="citation" data-citation="%7B%22citationItems%22%3A%5B%7B%22uris%22%3A%5B%22http%3A%2F%2Fzotero.org%2Fusers%2F13042945%2Fitems%2FKIBCEQBX%22%5D%2C%22locator%22%3A%223%22%7D%5D%2C%22properties%22%3A%7B%7D%7D" ztype="zcitation">(<span class="citation-item"><a href="zotero://select/library/items/KIBCEQBX">Lause 等, 2023, p. 3</a></span>)</span>的原理

(1)在理想模型中,我们假设计数矩阵只受细胞自身的差异性表达所影响.

此时,假设细胞$c$ 中真实表达了总共$n_c$ 个mRNA,其中基因$g$ 的真实表达量占总表达量的比例为$p_g$ .

根据实验和理论,此时理论上观测到的细胞$c$ 中基因$g$ 的表达量(即UMI counts)$k_{cg}$ 服从如下负二项分布:

$k_{cg}\sim\mathrm{NB}(\mu_{cg},\theta)$

其中参数$\mu_{cg}=n_cp_g$ ,$\theta$ 是<span class="highlight" data-annotation="%7B%22attachmentURI%22%3A%22http%3A%2F%2Fzotero.org%2Fusers%2F13042945%2Fitems%2F5IJ4QJUI%22%2C%22pageLabel%22%3A%223%22%2C%22position%22%3A%7B%22pageIndex%22%3A2%2C%22rects%22%3A%5B%5B133.614%2C573.325%2C274.564%2C582.73%5D%5D%7D%2C%22citationItem%22%3A%7B%22uris%22%3A%5B%22http%3A%2F%2Fzotero.org%2Fusers%2F13042945%2Fitems%2FKIBCEQBX%22%5D%2C%22locator%22%3A%223%22%7D%7D" ztype="zhighlight"><a href="zotero://open-pdf/library/items/5IJ4QJUI?page=3">“inverse overdispersion parameter”</a></span> <span class="citation" data-citation="%7B%22citationItems%22%3A%5B%7B%22uris%22%3A%5B%22http%3A%2F%2Fzotero.org%2Fusers%2F13042945%2Fitems%2FKIBCEQBX%22%5D%2C%22locator%22%3A%223%22%7D%5D%2C%22properties%22%3A%7B%7D%7D" ztype="zcitation">(<span class="citation-item"><a href="zotero://select/library/items/KIBCEQBX">Lause 等, 2023, p. 3</a></span>)</span>.

<span class="highlight" data-annotation="%7B%22attachmentURI%22%3A%22http%3A%2F%2Fzotero.org%2Fusers%2F13042945%2Fitems%2F5IJ4QJUI%22%2C%22position%22%3A%7B%22rects%22%3A%5B%5D%7D%2C%22citationItem%22%3A%7B%22uris%22%3A%5B%22http%3A%2F%2Fzotero.org%2Fusers%2F13042945%2Fitems%2FKIBCEQBX%22%5D%7D%7D" ztype="zhighlight"><a href="zotero://open-pdf/library/items/5IJ4QJUI?page=NaN">“Note that θ in this formulation is shared between all genes; based on negative control UMI data,”</a></span> <span class="citation" data-citation="%7B%22citationItems%22%3A%5B%7B%22uris%22%3A%5B%22http%3A%2F%2Fzotero.org%2Fusers%2F13042945%2Fitems%2FKIBCEQBX%22%5D%7D%5D%2C%22properties%22%3A%7B%7D%7D" ztype="zcitation">(<span class="citation-item"><a href="zotero://select/library/items/KIBCEQBX">Lause 等, 2023</a></span>)</span> 请注意，此公式中的 θ 在所有基因之间共享；基于阴性对照 UMI 数据，

(2)在实际观测中,给定一个实际观测到的UMI计数矩阵,我们可以得出对于参数$\mu_{cg}$ 的一个极大似然估计:

$\hat{\mu}_{cg}=\frac{\sum_jk_{cj}\cdot\sum_ik_{ig}}{\sum_{ij}k_{ij}},$

我们计算该估计量的皮尔逊残差(观测值与模型预测之间的差异除以模型方差)为:

$R_{cg}=\frac{k_{cg}-\hat{\mu}_{cg}}{\sqrt{\hat{\mu}_{cg}+\hat{\mu}_{cg}^2/\theta}},$

皮尔逊残差衡量了估计出的参数与原模型的差异.差异越大,代表实际模型偏离原模型越远,这就代表了技术效应的影响越大.

此时,我们就可以根据残差大小进行选择,即选择那些残差较小的数据;去除残差较大的数据.

***


</div></div>






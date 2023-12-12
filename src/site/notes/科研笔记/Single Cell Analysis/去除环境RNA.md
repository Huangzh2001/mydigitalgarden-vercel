---
{"dg-publish":true,"tags":["Notes","Single-Cell-Analysis"],"permalink":"/科研笔记/Single Cell Analysis/去除环境RNA/","dgPassFrontmatter":true}
---


# 校正环境 RNA 的方法

## 方法 1 : SoupX

### (1)方法原理


<div class="transclusion internal-embed is-loaded"><div class="markdown-embed">



## 🎈方法

### 步骤简介

SoupX算法由下面3个步骤组成:

1.  从空液滴中估计环境mRNA的表达谱
2.  估计每个细胞的污染分数
3.  使用环境nRNA的表达谱和估计的污染来校正每个细胞的基因表达

### 方法原理

(1)首先我们得找出空液滴,我们明确假设所有总基因表达量$<N_{emp}$的UMIs都是空液滴.

(2)对于某个基因g,我们定义其背景表达比例(<span class="highlight" data-annotation="%7B%22attachmentURI%22%3A%22http%3A%2F%2Fzotero.org%2Fusers%2F13042945%2Fitems%2FBG32HPJZ%22%2C%22pageLabel%22%3A%222%22%2C%22position%22%3A%7B%22pageIndex%22%3A1%2C%22rects%22%3A%5B%5B130.974%2C283.291%2C270.47%2C291.159%5D%5D%7D%2C%22citationItem%22%3A%7B%22uris%22%3A%5B%22http%3A%2F%2Fzotero.org%2Fusers%2F13042945%2Fitems%2FAQ3VA93J%22%5D%2C%22locator%22%3A%222%22%7D%7D" ztype="zhighlight"><a href="zotero://open-pdf/library/items/BG32HPJZ?page=2">“The fraction of background expression”</a></span> <span class="citation" data-citation="%7B%22citationItems%22%3A%5B%7B%22uris%22%3A%5B%22http%3A%2F%2Fzotero.org%2Fusers%2F13042945%2Fitems%2FAQ3VA93J%22%5D%2C%22locator%22%3A%222%22%7D%5D%2C%22properties%22%3A%7B%7D%7D" ztype="zcitation">(<span class="citation-item"><a href="zotero://select/library/items/AQ3VA93J">Young 和 Behjati, 2020, p. 2</a></span>)</span>)为

$b_g=\frac{\sum_dn_{g,d}}{\sum_d\sum_gn_{g,d}},$

其中,$n_{g,d}$是基因$g$ 在液滴$d$ 中的计数数量;对$d$ 求和则取遍了所有UMIs$<N_{emp}$ 的液滴(即那些我们认为是空液滴的液滴).

(3)一个液滴 $c$ 中,某个基因$g$ 的来源有两个:<span style="color: rgb(64, 64, 64)"><span style="background-color: rgb(255, 255, 255)">细胞内的+背景污染的.于是我们就能得到下式:</span></span>

$n_{g,c}=m_{g,c}+o_{g,c},$

其中,对于某个基因$g$ ,$m_{g,c}$ 代表在液滴 $c$中来自细胞内的基因;而$o_{g,c}$ 代表着在液滴 $c$中来自背景污染的基因;$n_{g,c}$ 代表在液滴 $c$中$g$ 的总基因.我们想要求得的就是$m_{g,c}$ .

(4)我们假设来自背景污染的基因的相对丰度在细胞之间没有差异.于是我们就可以得出下式:

$o_{g,c}=N_{c}\rho_{c}b_{g},$

其中$N_c=\sum_gn_{g,c}$ 表示在液滴 $c$ 中的全部基因量;$\rho_c$ 是在液滴 $c$ 中的背景污染分数.

我们接着对某些基因$g$ 进行求和可得:

$\rho_cN_c\sum_gb_g=\sum_go_{g,c}$

$m_{g,c}=0$ 的意思是在细胞$c$ 中,基因$g$ 的表达量为0.此时我们称基因$g$ 在细胞类型$c$ 中是<span class="highlight" data-annotation="%7B%22attachmentURI%22%3A%22http%3A%2F%2Fzotero.org%2Fusers%2F13042945%2Fitems%2FBG32HPJZ%22%2C%22pageLabel%22%3A%222%22%2C%22position%22%3A%7B%22pageIndex%22%3A1%2C%22rects%22%3A%5B%5B452.049%2C608.524%2C540.882%2C616.392%5D%5D%7D%2C%22citationItem%22%3A%7B%22uris%22%3A%5B%22http%3A%2F%2Fzotero.org%2Fusers%2F13042945%2Fitems%2FAQ3VA93J%22%5D%2C%22locator%22%3A%222%22%7D%7D" ztype="zhighlight"><a href="zotero://open-pdf/library/items/BG32HPJZ?page=2">“strong negative markers”</a></span> .论文中举了一个例子:HBB基因只在血红细胞中表达,而不在其他任何细胞内表达,于是我们就说HBB基因在血红细胞中是<span class="highlight" data-annotation="%7B%22attachmentURI%22%3A%22http%3A%2F%2Fzotero.org%2Fusers%2F13042945%2Fitems%2FBG32HPJZ%22%2C%22pageLabel%22%3A%222%22%2C%22position%22%3A%7B%22pageIndex%22%3A1%2C%22rects%22%3A%5B%5B486.529%2C598.067%2C540.912%2C605.934%5D%2C%5B307.782%2C587.609%2C334.595%2C595.476%5D%5D%7D%2C%22citationItem%22%3A%7B%22uris%22%3A%5B%22http%3A%2F%2Fzotero.org%2Fusers%2F13042945%2Fitems%2FAQ3VA93J%22%5D%2C%22locator%22%3A%222%22%7D%7D" ztype="zhighlight"><a href="zotero://open-pdf/library/items/BG32HPJZ?page=2">“strong positive marker”</a></span> <span class="citation" data-citation="%7B%22citationItems%22%3A%5B%7B%22uris%22%3A%5B%22http%3A%2F%2Fzotero.org%2Fusers%2F13042945%2Fitems%2FAQ3VA93J%22%5D%2C%22locator%22%3A%222%22%7D%5D%2C%22properties%22%3A%7B%7D%7D" ztype="zcitation">(<span class="citation-item"><a href="zotero://select/library/items/AQ3VA93J">Young 和 Behjati, 2020, p. 2</a></span>)</span>;而在其他细胞中是strong negative markers.

于是我们在上面的求和中令$g$ 是在细胞$c$ 中strong negative markers的基因,此时$m_{g,c}=0$ ,于是我们就有

$\rho_cN_c\sum_gb_g=\sum_gn_{g,c}$

接着我们可以计算出液滴$c$ 的背景污染分数为:

$\rho_c=\frac{\sum_gn_{g,c}}{N_c\sum_gb_g},$

而上式中右边的$n_{g,c},N_c,b_g$ 我们都是已知的,于是我们就可以算出$\rho_c$ (在液滴 $c$ 中的背景污染分数).

然后我们就能很直观地算出我们心心念念的$m_{g,c}$ 如下:

$m_{g,c}=n_{g,c}-N_c\rho_cb_g,$

(5)最后只剩下一个问题,那就是如何找到那些$m_{g,c}=0$的基因.文章中给出了两种方法:

*   一种是人为指定.如果从已知的生物学知识(先验知识)知道哪些基因在细胞中是不会表达的,就可以令上式中的

    $g$

    为这些基因.(比如上面提到的HBB基因)

*   如果我们没有先验知识,文章也提供了一种自动化的方法.

<span class="highlight" data-annotation="%7B%22attachmentURI%22%3A%22http%3A%2F%2Fzotero.org%2Fusers%2F13042945%2Fitems%2FBG32HPJZ%22%2C%22pageLabel%22%3A%222%22%2C%22position%22%3A%7B%22pageIndex%22%3A1%2C%22rects%22%3A%5B%5B402.689%2C389.626%2C540.894%2C397.494%5D%2C%5B307.782%2C379.159%2C540.902%2C387.027%5D%2C%5B307.782%2C368.701%2C540.91%2C376.569%5D%2C%5B307.782%2C358.243%2C540.901%2C366.111%5D%2C%5B307.782%2C347.252%2C540.908%2C355.644%5D%2C%5B307.782%2C337.319%2C540.909%2C345.186%5D%2C%5B307.782%2C326.861%2C540.898%2C334.728%5D%2C%5B307.782%2C316.403%2C540.9%2C324.27%5D%2C%5B307.782%2C305.936%2C540.901%2C313.803%5D%2C%5B307.782%2C295.478%2C540.918%2C303.345%5D%2C%5B307.782%2C285.02%2C540.899%2C292.887%5D%2C%5B307.782%2C274.553%2C540.921%2C282.42%5D%2C%5B307.782%2C264.095%2C467.898%2C271.962%5D%5D%7D%2C%22citationItem%22%3A%7B%22uris%22%3A%5B%22http%3A%2F%2Fzotero.org%2Fusers%2F13042945%2Fitems%2FAQ3VA93J%22%5D%2C%22locator%22%3A%222%22%7D%7D" ztype="zhighlight"><a href="zotero://open-pdf/library/items/BG32HPJZ?page=2">“Where this is not known in advance, we provide an automated alternative to estimate the contamination fraction (see Supplementary Fig. S1). The automated approach first identifies markers of each cluster of cells in the data. For each strong marker, it is assumed that mg,c = 0 for all cells in clusters where the gene is not a marker and the contamination fraction is estimated (Supplementary Fig. S1). Performing this estimation across all strong marker genes provides a set of estimates of the contamination fraction. To obtain a final value, it is assumed that inaccurate estimates will have no preferred value while true estimates will cluster around the true value. The most common value is taken as the final estimate of the contamination fraction (see Fig. 1, Step 2.2).”</a></span> <span class="citation" data-citation="%7B%22citationItems%22%3A%5B%7B%22uris%22%3A%5B%22http%3A%2F%2Fzotero.org%2Fusers%2F13042945%2Fitems%2FAQ3VA93J%22%5D%2C%22locator%22%3A%222%22%7D%5D%2C%22properties%22%3A%7B%7D%7D" ztype="zcitation">(<span class="citation-item"><a href="zotero://select/library/items/AQ3VA93J">Young 和 Behjati, 2020, p. 2</a></span>)</span>

(这一段自动化的方法我没有看懂)

### 补充

<span style="color: rgb(51, 51, 51)"><span style="background-color: rgb(255, 255, 255)">SoupX 可以在没有聚类信息的情况下运行，但是，文章指出，如果提供聚类信息，结果会更好。</span></span>

<span class="highlight" data-annotation="%7B%22attachmentURI%22%3A%22http%3A%2F%2Fzotero.org%2Fusers%2F13042945%2Fitems%2FBG32HPJZ%22%2C%22pageLabel%22%3A%222%22%2C%22position%22%3A%7B%22pageIndex%22%3A1%2C%22rects%22%3A%5B%5B374.493%2C462.326%2C540.909%2C470.718%5D%2C%5B307.782%2C451.868%2C540.909%2C460.26%5D%2C%5B307.782%2C441.401%2C540.902%2C449.793%5D%2C%5B307.782%2C430.943%2C540.904%2C439.335%5D%2C%5B307.782%2C421.01%2C320.842%2C428.877%5D%5D%7D%2C%22citationItem%22%3A%7B%22uris%22%3A%5B%22http%3A%2F%2Fzotero.org%2Fusers%2F13042945%2Fitems%2FAQ3VA93J%22%5D%2C%22locator%22%3A%222%22%7D%7D" ztype="zhighlight"><a href="zotero://open-pdf/library/items/BG32HPJZ?page=2">“SoupX optionally uses clustering information to refine the set of cells for which it can be assumed that mg,c = 0. If it can be shown for any cell c in cluster P that mg,c > 0, then it is assumed that mg,c > 0 for all c ∈ P (see Supplementary Fig. S1).”</a></span> <span class="citation" data-citation="%7B%22citationItems%22%3A%5B%7B%22uris%22%3A%5B%22http%3A%2F%2Fzotero.org%2Fusers%2F13042945%2Fitems%2FAQ3VA93J%22%5D%2C%22locator%22%3A%222%22%7D%5D%2C%22properties%22%3A%7B%7D%7D" ztype="zcitation">(<span class="citation-item"><a href="zotero://select/library/items/AQ3VA93J">Young 和 Behjati, 2020, p. 2</a></span>)</span>

***


</div></div>

### (2)代码实现


### (3)项目实战

[[科研笔记/Single Cell Analysis/Jupyter Notebok/Correction of ambient RNA\|Correction of ambient RNA]]

## 方法 2 : DecontX

### (1)方法原理


### (2)代码实现


### (3)项目实战
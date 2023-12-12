---
{"tags":["Notes","Math","StochasticProcess"],"version":1,"dg-publish":true,"permalink":"/课堂笔记/随机过程/去相关与K-L展开/","dgPassFrontmatter":true}
---


# 去相关

设$X=\left(X_{1},\cdots,X_{n}\right)^{T}$是一组随机变量，它们之间的相关构成了一个矩阵$\left(E(X_iX_j)\right)_{ij}$，称为**相关矩阵$R_X$**。

我们希望找到一个线性变换$Y=AX,A\in \mathbb{R}^{n\times n}$，使得$E(Y_iY_j)=0,\forall i\neq j$。即消除$X_i$之间的相关性。

$Y$的相关矩阵为：
$$
R_{Y}=E(YY^{T})=E(AXX^{T}A^{T})=AE(XX^{T})A^{T}=AR_{X}A^{T}
$$

其中$R_Y$是一个对角阵，所以我们的问题就变成了：想要找到一个矩阵$A$，使得$AR_XA^T$是个对角阵。

由于$R_X$是个对称矩阵，所以可以正交对角化。即存在正交矩阵$U$，使得
$$
R_{X}=U\cdot\mathrm{diag}\left(\lambda_{1},\cdots,\lambda_{n}\right)U^{T}=\sum_{k=1}^{n}\lambda_{k}U_{k}U_{k}^{T}
$$
其中$U=(U_{1},\cdots,U_{n}),\lambda_i$为特征值，$U_i$为$\lambda_i$的特征向量。于是有$U^{T}R_X U=\mathrm{diag}(\lambda_{1},\cdots,\lambda_{n})$，所以我们取$A=U^T$即可。

于是有
$$
Y=U^{T}X\implies X=UY=(U_1,\cdots,U_{n})\left.\left(\begin{matrix}Y_{1}\\\vdots\\Y_{n}\end{matrix}\right.\right)=\sum_{k=1}^{n}U_{k}Y_{k}
$$

于是我们得到了：$X=\sum_{k=1}^{n}U_{k}Y_{k}$。这个分解被称为K-L分解(Karhunen Loeve Expansion)。我们观察这个式子，发现了以下几点：

```ad-info
(1) $U_k$是相关矩阵$R_X$的特征向量

(2) $U_i$和$U_j$是正交的（向量间的正交）

(3) $Y_i$和$Y_j$是正交的（随机变量间的正交，即$E(Y_iY_j)=0$）

因为既有$U_i$和$U_j$正交，又有$Y_i$和$Y_j$正交，所以这里被称为“**双正交**(Biorthogonal)”。

(4) $E(Y_iY_i)=\lambda_i$，其中$\lambda_i$是相关矩阵$R_X$的特征值
```

# K-L分解

我们想知道对于随机过程$X(t)$是否也存在类似的K-L分解？事实上，这是存在的。
$$
X(t)=\sum_{k=-\infty}^{+\infty}\alpha_{k}\phi_{k}(t)
$$

（具体的证明这里没有写，后面争取给出）

这个分解跟前面的K-L分解:$X=\sum_{k=1}^{n}U_{k}Y_{k}$是一一对应的，即：

```ad-info

(1) $\phi_{k}(t)$是相关函数$R_X(t,s)$的“特征向量”，即$\phi_{k}(t)$满足$\int_{T}R_{X}(t,s)\phi_{k}(s)ds=\lambda_{k}\phi_{k}(t)$

$\phi_{k}(t)$其实是$U_k$的推广。这怎么理解呢？由于$U_k$是相关矩阵$R_X$的特征向量$R_{X}U_{k}=\lambda_{k}U_{k}$，即
$$
\left.\mathrm{R_x}U_k=\lambda _kU_k\Longleftrightarrow\left(\begin{array}{ccc}\mathrm{Rx(1,1)}&\mathrm{Rx(1,2)}&\cdots\mathrm{Rx(1,n)}\\\vdots&\vdots&\vdots\\\vdots&\vdots&\vdots\\\mathrm{Rx(n,1)}&\mathrm{Rx(n,2)}&\cdots\mathrm{Rx(n,n)}\end{array}\right.\right)\left(\begin{array}{c}U_k(1)\\\vdots\\U_k(n)\end{array}\right)=\lambda_k\left(\begin{array}{c}U_k(1)\\\vdots\\U_k(n)\end{array}\right)
$$
写成分量形式就是
$$
\sum_{j}R_{X}\left(i,j\right)U_{k}\left(j\right)=\lambda_{k}U_{k}\left(i\right)
$$
如果我们将这里的连加号"$\sum_j$"换成"$\int_{L}$"就是$\int_{T}R_{X}\left(t,s\right)\phi_{k}\left(s\right)ds=\lambda_{k}\phi_{k}\left(t\right)$。

(2) $\phi_{k}(t),\phi_{l}(t)$相互正交，即$\int_{T}\phi_{k}\left(t\right)\phi_{l}\left(t\right)dt=0\left(k\neq l\right)$

这是$U_i,U_j$相互正交的推广

(3) $\alpha_m,\alpha_n$相互正交，即$E(\alpha_m\alpha_n)=\lambda_n\delta_{mn}$

我们可以把$\alpha_m$看成$Y_m$的推广。我们下面来证明这个命题：

首先我们先算出$\alpha_n$的表达式

由$X(t)=\sum_{k=-\infty}^{+\infty}\alpha_{k}\phi_{k}(t)$得：
$$
\begin{aligned}\int_{T}X(t)\phi_{n}(t)dt&=\int_{T}\sum_{k=-\infty}^{+\infty}\alpha_{k}\phi_{k}(t)\phi_{n}(t)dt\\&=\sum_{k=-\infty}^{+\infty}\alpha_{k}\int_{T}\phi_{k}(t)\phi_{n}(t)dt\\&=\alpha_n\end{aligned}
$$
所以$\alpha_{n}=\int_{T}X(t)\phi_{n}(t)dt$。

于是：
$$
\begin{array}{rcl}\mathrm{E}\left(\alpha_m\alpha_n\right)&=&\mathrm{E}\left(\int_TX(t)\phi_m(t)dt\right)\left(\int_TX(s)\phi_n(s)dt\right)\\\\&=&\int_T\int_T\mathrm{E}\left(X(t)X(s)\right)\phi_m(t)\phi_n(s)dtds\\\\&=&\int_T\int_TR_X(t,s)\phi_m(t)\phi_n(s)dtds\\\\&=&\int_T\lambda_n\phi_n(t)\phi_m(t)dt\\\\&=&\lambda_n\delta_{mn}\end{array}
$$
```

```ad-info
K-L展开的意义：它将随机和过程分开。

随机过程为什么难，是因为随机和过程搅合在了一起。K-L展开给我们提供了一个视角：我们可以把随机和过程分离开来，随机是随机（它只浓缩在系数$\alpha_k$里），过程是过程（它只取决于展开的基函数$\phi_k(t)$）。更何况K-L展开的性质实在是太好了，双正交。
```

于是我们就转到下一个问题了：如何求出这个K-L分解？事实上，我们只要求出了$\phi_k(t)$，就能求出整个分解式。那么，如何求出$\phi_k(t)$？

下面我们先给出求宽平稳(W.S.S)+周期(Periodic)过程的K-L分解的方法。

# 宽平稳周期过程的K-L分解与谱表示

## 宽平稳周期过程的K-L分解

根据W.S.S的性质$R_X(t,s)=R_X(t-s)$可得
$$
\int_{T}R_{X}(t-s)\phi_{k}(s)ds=\lambda_k\phi_{k}(t)
$$

在由于$X(t)$是周期的，所以$X(t)=X(t+T)$，于是我们有$E|X(t)-X(t+T)|^2=0$，这等价于$R_X(\tau)=R_X(\tau+T)$（因为我们之前证明一个命题的时候曾经搭了这样一座桥：$E|X(t)-X(t+T)|^2=0\Longleftrightarrow R_X(\tau)=R_X(\tau+T)$）

于是我们有$\phi_{k}(t)=\exp(j\frac{2k\pi}{T}t)$。我们不妨验证一下：

$$
\begin{aligned}\int_0^TR_X\left(t-s\right)\exp\left(j\frac{2k\pi}Ts\right)ds&=-\int_{t-T}^tR_X\left(s^{\prime}\right)\exp\left(-j\frac{2k\pi}Ts^{\prime}\right)ds^{\prime}\exp\left(j\frac{2k\pi}Tt\right)\\\\&=\int_0^TR_X\left(s^{\prime}\right)\exp\left(-j\frac{2k\pi}Ts^{\prime}\right)ds^{\prime}\exp\left(j\frac{2k\pi}Tt\right)\\\\
&=\lambda_{k}\exp(j\frac{2k\pi}{T}t)\end{aligned}
$$

（其中第一个等号做了变量替换：令$s'=t-s$；最后一个等号是因为$\int_{0}^{T}R_{X}(s^{\prime})\exp(-j\frac{2k\pi}{T}s^{\prime})ds^{\prime}$是个只依赖于$k$，不依赖于$t$的常数，于是我们把它设为$\lambda_k$）

于是我们就得到了宽平稳周期过程$X(t)$的K-L分解：
$$
\begin{aligned}
X(t)=\sum_{k=-\infty}^{+\infty}\alpha_{k}\phi_{k}(t)=\sum_{k=-\infty}^{+\infty}\alpha_{k}\exp{(j\frac{2k\pi}{T}t)},\\
~where~\alpha_{k}=\int_{T}X(t)\exp(j\frac{2k\pi}{T}t)dt
\end{aligned}
$$

## 随机过程的谱表示(Spectral Representation)

我们前面曾经对随机过程$X(t)$的相关函数$R_X(t)$进行傅里叶变换得到了功率谱密度(PSD)$S_X(\omega)$（PSD也被称为随机过程的相关函数的“谱表示”）。

那么我们自然会问：能否对随机过程$X(t)$自身做傅里叶变换，从而可以在频域上观察随机过程（即问：随机过程自身是否有谱表示？）。如果有，那么$X(t)$的谱表示和其相关函数的谱表示——PSD间存在什么关系？

```ad-info
对随机过程进行傅里叶变换即：
$$
\begin{aligned}\hat{X}(w)&=\int_{-\infty}^{+\infty}X(t)\exp{(-jwt)}dt\\X(t)&=\int_{-\infty}^{+\infty}\hat{X}(w)\exp{(jwt)}dw\end{aligned}
$$
（这里只是形式上的记号.上面的第一个式子是傅里叶变换；第二个式子是傅里叶逆变换）
```

很遗憾，往往是不能对$X(t)$这样做的。

为了能对随机过程$X(t)$进行傅里叶变换，我们先将积分推广：
$$
\begin{aligned}&(1)\int X(t)dt\quad\Longrightarrow\sum_{k}X(t_k)(t_k-t_{k-1})\\\\&(2)\int X(t)dg(t)\Longrightarrow\sum_{k}X(t_k)\left(g(t_k)-g(t_{k-1})\right)\end{aligned}
$$

数学家们证明了：对于宽平稳周期过程$X(t)$，存在谱表示(Spectral Representation)：
$$
X(t)=\int_{-\infty}^{+\infty}\exp(jwt)dF_{X}(w)
$$
```ad-info
如何理解呢？

我们前面求出了宽平稳周期过程$X(t)$的K-L分解：
$$
\begin{aligned}
X(t)=\sum_{k=-\infty}^{+\infty}\alpha_{k}\phi_{k}(t)=\sum_{k=-\infty}^{+\infty}\alpha_{k}\exp{(j\frac{2k\pi}{T}t)},
\end{aligned}
$$

我们发现它的形式中存在$\exp{(j\frac{2k\pi}{T}t)}$，于是我们定义$dF_X(\omega)=F_X(\omega+\Delta\omega)-F_X(\omega)=\alpha_k$，然后再令$T\to \infty$，就得到了：
$$
X(t)=\sum_{k}\alpha_{k}\exp(j\frac{2k\pi}{T}t)\xrightarrow{T\to\infty}X(t)=\int_{-\infty}^{+\infty}\exp(jwt)dF_{X(\omega)}
$$
如果这样理解，我们就能由$\alpha_k$间相互正交，推出$dF(\omega)$间相互正交，即：
$$
E(\alpha_n\overline{\alpha_m})=0\Longrightarrow E(dF_X(w_1)\overline{dF_X(w_2)})=0
$$
```

## 谱表示与功率谱密度(PSD)之间的关系

由于宽平稳周期过程的$\alpha_i$相互正交，可以推出$dF(\omega)$间相互正交，即：
$$
E(dF_{X}(w_{1})\overline{dF_{X}(w_{2})})=0,(w_{1}\neq w_{2})
$$
所以：
$$
\begin{aligned}R_X(t,s)&=\mathrm{E}\left(\left.X(t)\overline{X(s)}\right)\right.\\&=\mathrm{E}\left(\int_{-\infty}^{+\infty}\exp\left(jw_1t\right)\mathrm{d}F_X(w_1)\right)\left(\overline{\int_{-\infty}^{+\infty}\exp\left(jw_2s\right)\mathrm{d}F_X(w_2)}\right)\\\\&=\int_{-\infty}^{+\infty}\int_{-\infty}^{+\infty}\exp\left(jw_1t-jw_2s\right)\mathrm{E}\left(dF_X(w_1)\overline{dF_X(w_2)}\right)\\\\&=\int_{-\infty}^{+\infty}\exp(jw(t-s))E|dF_{X}(w)|^{2}\end{aligned}
$$

而根据功率谱密度(PSD)的表达式：
$$
R_{X}\left(t,s\right)=\frac{1}{2\pi}\int_{-\infty}^{+\infty}\exp\left(jw\left(t-s\right)\right)S_{X}\left(w\right)dw
$$

我们可得：
$$
\left.\left.\begin{array}{rrr}R_X(t,s)=\int_{-\infty}^{+\infty}\exp(jw(t-s))\mathrm{E}|\mathrm{d}F_X(w)|^2&\\\\R_X(t,s)=\frac{1}{2\pi}\int_{-\infty}^{+\infty}\exp(jw(t-s))\mathrm{S}_X(w)\mathrm{d}w&\end{array}\right.\right\}\Longrightarrow\mathrm{E}|\mathrm{d}F_X(w)|^2=\frac{1}{2\pi}S_X(w)\mathrm{d}w
$$

即随机过程的谱表示和其相关函数的谱表示（功率谱密度）之间的关系为：
$$
\mathrm{E}|\mathrm{d}F_X(w)|^2=\frac{1}{2\pi}S_X(w)\mathrm{d}w
$$
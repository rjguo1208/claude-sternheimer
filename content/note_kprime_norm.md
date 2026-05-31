# Implementation Note — the $k'$-sum normalization convention

*为什么 Layer-1 下折叠里的 rest 空间求和带一个 $1/N_k$ 因子,以及这如何把"逐通道相加过大"的问题一次性解决。EDT P2 → P3,详见 [Implementation Plan](plan.html)。*

> **结论(已定):** 对 $k'$ 的离散求和就是 BZ 积分的离散化,所以它带一个 **$1/N_k$**(BZ 测度)因子;这把朴素相加得到的 $\sim-70$ Ry 修正为物理的 $\sim-0.5$ Ry。此外,只要 $\Delta V$ 在超胞内收敛,矩阵元 $M$(进而 $\tilde V$)就**与超胞大小无关**,因此 $1/N_k$ 之外**没有**额外的 $N_{\rm sc}$ 因子(§4)。下面给出推导、与"逐通道闭合"的协调、超胞无关性,以及确认用的数值检验。

## 0. 待定的量

rest 自能(Feshbach 下折叠的二阶项):

$$
\Sigma_{mn}(\omega_0)=\big\langle m k_i\big|\,\Delta V\,G^R(\omega_0)\,\Delta V\,\big|n k_i\big\rangle
= \frac{1}{N_k}\sum_{k'}\sum_{r\in R(k')}\frac{M_{m,\,rk'}\;M_{rk',\,n}}{\omega_0-\varepsilon_{rk'}},
\qquad M_{m,rk'}\equiv\langle m k_i|\Delta V|r k'\rangle .
$$

那个 $1/N_k$ 就是下面要论证的测度因子($N_k$ = rest BZ 网格点数)。下折叠势 $\tilde V_{mn}=M_{mn}+\Sigma_{mn}$ 与 $M$ 同量纲、同归一化,Born 极限退回 $\tilde V\to M$。

## 1. 为什么是 $1/N_k$:$k'$ 求和 = BZ 积分

无穷晶体里的恒等分解是 BZ 积分:

$$
\mathbb 1=\sum_n\frac{\Omega_{\rm cell}}{(2\pi)^3}\int_{BZ}\!d^3k\;|\psi_{nk}\rangle\langle\psi_{nk}| .
$$

在 $N_k$ 点 Monkhorst–Pack 网格上离散化,$\dfrac{\Omega_{\rm cell}}{(2\pi)^3}\displaystyle\int_{BZ}d^3k\to\dfrac1{N_k}\sum_k$(因为 $V_{BZ}=(2\pi)^3/\Omega_{\rm cell}$),用 QE 的 **per–原胞归一** 态($\langle\psi_{nk}|\psi_{nk}\rangle_{\rm cell}=1$)即得

$$
\boxed{\;\mathbb 1=\frac{1}{N_k}\sum_{nk}|\psi_{nk}\rangle\langle\psi_{nk}|\;}
$$

因此预解算符与 rest Green 函数都带 $1/N_k$:

$$
G_0(\omega)=\frac{1}{N_k}\sum_{nk}\frac{|\psi_{nk}\rangle\langle\psi_{nk}|}{\omega-\varepsilon_{nk}},
\qquad
G^R(\omega_0)=\frac{1}{N_k}\sum_{k'}\sum_{r\in R}\frac{|rk'\rangle\langle rk'|}{\omega_0-\varepsilon_{rk'}} .
$$

**一致性**:Layer-2 的 $G^A(\omega)$ 同样带 $1/N_k$;rest 网格加密 $\frac1{N_k}\sum_{k'}\to\int_{BZ}$ 的收敛,正是早先"rest 必须覆盖全 BZ"那条修正的定量版;而 EDI 金规则里的 $\frac1{N_k}\sum_{k'}$(末态 BZ 平均)是**另一处**、同源的 $1/N_k$。——即:**每一处内部 $k$ 求和都带 $1/N_k$。**

## 2. 与"逐通道闭合"的协调(两件不同的事)

固定 $k'$、对**所有带**求和,是 $k'$ 扇区内的 Parseval,**不带因子**:

$$
\sum_{n'\,(\text{all})}\big|\langle n'k'|s\rangle\big|^2=\langle s(k')|s(k')\rangle_{\rm cell}=\sum_G|s(k',G)|^2,
\qquad |s\rangle=\Delta V|n k_i\rangle .
$$

把这个(无因子的)逐通道结果**带着 $1/N_k$** 对 $k'$ 求和,得到一个有限的、原胞尺度的量:

$$
\frac{1}{N_k}\sum_{k'}\sum_{n'\,(\text{all})}|M_{n,n'k'}|^2=\frac{1}{N_k}\sum_{k'}\langle s(k')|s(k')\rangle=\big\langle n k_i\big|\Delta V^{2}\big|n k_i\big\rangle .
$$

所以**两件事并不矛盾**:逐通道的带求和是扇区内 Parseval(无因子);$k'$ 之间的求和是 BZ 测度(带 $1/N_k$)。我此前把它们混为一谈,错误地下了"无 $1/N_k$"的结论。

## 3. 数量级:$-70$ Ry → $-0.5$ Ry

P2 的 Sternheimer 给出**逐通道** $\Delta\tilde V_{\rm chan}(k')\sim0.3$–$0.7$ Ry。朴素相加 144 个通道 $\sim-70$ Ry(非物理)。带上 $1/N_k$:

$$
\Sigma_{nn}=\frac{1}{N_k}\sum_{k'}\Delta\tilde V_{\rm chan}(k')\approx\frac{-70}{144}\approx-0.5\ \text{Ry},
$$

与 Born 自能 $M_{nn}\sim0.7$ Ry 同量级——是个**大但物理**的修正(空位是强散射,Born 本就该被显著修正,这正是 T-matrix 的用武之地)。

## 4. 超胞无关性 ⇒ 没有 $N_{\rm sc}$ 因子

**物理要求:** 只要 $\Delta V$ 在超胞内**收敛**(局域、在边界前衰减为零),电子–缺陷矩阵元 $M=\langle\psi_{nk_i}|\Delta V|\psi_{mk_f}\rangle$ 就**与超胞大小无关** —— 更大的超胞只是多出 $\Delta V=0$ 的体/真空区,不改变只在缺陷邻域有贡献的积分。

**实现一致:** EDI 的 $M=\frac1{N_{\rm nnr}}\sum_{r}u^*\,V^q_{\rm folded}\,u$ 用的是**原胞**测度 $\frac1{N_{\rm nnr}}$(不是超胞测度 $\frac1{N_{\rm nnr}^{\rm SC}}$),作用在把局域 $\Delta V$ 折回原胞的 $V^q_{\rm folded}$ 上;局域性 ⇒ 求和只在缺陷区有贡献(固定)⇒ $M$ **不随 $N_{\rm sc}$ 变**(T1 已证 $M$ 与 EDI 一致到 $2.3\times10^{-13}$)。单缺陷浓度由金规则里的 $n_d$ 处理,**不进 $M$**。

**推论:** $\Sigma=\frac1{N_k}\sum_{k'}\sum_r MM/(\omega_0-\varepsilon)$ 由超胞无关的 $M$ 构成 ⇒ $\Sigma$、$\tilde V$ 也**超胞无关**,前面估出的 $\Sigma_{nn}\approx-0.5$ Ry 是个 **per-defect 的物理量**。**所以 $1/N_k$ 之外没有额外的 $N_{\rm sc}$ 因子。**

**两条独立的收敛轴(别混):** (i) **超胞大小** —— 决定 $\Delta V$(进而 $M$)是否收敛、是否超胞无关;(ii) **rest BZ 网格 $N_k$** —— $\frac1{N_k}\sum_{k'}\to\int_{BZ}$ 的积分收敛($N_k$ 可比超胞折叠更密,用任意 $q$ 的实空间折叠取 $\Delta V(q)$)。

## 5. 确认用的检验(P3 起点)

$1/N_k$ 由 BZ 测度定下、$M$ 超胞无关(§4)⇒ 归一化已完全确定。以下检验**确认实现无误**:

- **逐通道闭合**:单个 $k'$,$\sum_{n'=1}^{150}|M_{n,n'k'}|^2\overset?=\sum_G|s(k',G)|^2$(差额=高带不完备,应与 Sternheimer 高带尾一致)。
- **全和 $=\langle\Delta V^2\rangle$**:$\dfrac1{N_k}\sum_{k'}\sum_{n'}|M|^2\overset?=\dfrac1\Omega\int_{\rm cell}|\psi_{nk_i}|^2\,\Delta V^2$;两边都应超胞无关。
- **超胞无关性**:若有两个超胞尺寸,$M$(及 $\Sigma$)应不变;EDI 输运已隐含验证(收敛后迁移率与超胞无关)。
- **Born 极限锚定(金标准)**:关掉 rest/active 重求和 ⇒ $\tilde V\to M$、$T_{PP}\to M$,跑 transport,**必须逐位复现 EDI 迁移率**(含 $n_d$、$\frac1{N_k}$)。

> **小结**:归一化已完全定下 —— $k'$ 求和带 BZ 测度 $1/N_k$;$M$ 因 $\Delta V$ 局域而**超胞无关**(无 $N_{\rm sc}$ 因子);逐通道带求和无因子,$\frac1{N_k}\sum_{k'}\sum_{n'}|M|^2=\langle\Delta V^2\rangle$。两条收敛轴:超胞大小($\Delta V$)与 rest BZ 网格($N_k$)。最终由 Born 极限复现 EDI 迁移率作总锚定。

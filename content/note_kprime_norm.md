# Implementation Note — the $k'$-sum normalization convention

*为什么 Layer-1 下折叠里的 rest 空间求和必须先定下一个归一化,才能报出物理的 $\tilde V$;以及如何把它钉死。EDT P2 → P3,详见 [Implementation Plan](plan.html)。*

> 背景:Sternheimer 求解器(P2)已验证给出**逐通道**的 rest dressing(单通道 $\sim0.3$–$0.7$ Ry),但把全 BZ 的通道朴素相加会过大。这条 note 解释症结在哪、以及用什么**确定性判据**把 $k'$ 求和的归一化因子定下来。

## 0. 待定的量

rest 自能(Feshbach 下折叠的二阶项)是

$$
\Sigma_{mn}(\omega_0)=\big\langle m k_i\big|\,\Delta V\,G^R(\omega_0)\,\Delta V\,\big|n k_i\big\rangle
= c\sum_{k'}\sum_{r\in R(k')}\frac{M_{m,\,rk'}\;M_{rk',\,n}}{\omega_0-\varepsilon_{rk'}},
\qquad M_{m,rk'}\equiv\langle m k_i|\Delta V|r k'\rangle .
$$

那个常数 $c$ 就是 $k'$ 求和前**未定的测度因子**。下折叠势 $\tilde V_{mn}=M_{mn}+\Sigma_{mn}$ 必须与 $M$ **同量纲、同归一化**,才能塞进金规则、并在 Born 极限退回 $\tilde V\to M$。

## 1. 三个相互独立的归一化因子

**(a) Bloch 态归一化。** QE / EDI 用 **per–原胞归一**:$\langle\psi_{nk}|\psi_{nk}\rangle_{\rm cell}=\sum_G|c_{nk}(G)|^2=1$。在整个 Born–von Kármán 晶体($N_k$ 个原胞)里,同一个态的范数是 $N_k$,不是 $1$ —— 所以"全晶体的恒等分解"里每个投影子要除以它的范数。

**(b) $k'$ 完备性求和的测度。** 固定 $k'$、对**所有带**求和,是 $k'$ 扇区里的 Parseval 恒等式,**不带任何因子**:

$$
\sum_{n'\,(\text{all bands})}|n'k'\rangle\langle n'k'| = \mathbb 1_{k'} .
$$

不同 $k'$ 的扇区相互正交,直和即全空间:$\sum_{k'}\mathbb 1_{k'}=\mathbb 1$。所以**单从预解算符** $G_0=(\omega-H_0)^{-1}$ 看,

$$
G_0(\omega)=\sum_{n'k'}\frac{|n'k'\rangle\langle n'k'|}{\omega-\varepsilon_{n'k'}}\qquad(\text{形式上不带 }1/N_k).
$$

**关键提醒:** 这与金规则里的 $\tfrac1{N_k}\sum_{k'}$ 来源**完全不同** —— 金规则的 $1/N_k$ 是**末态 BZ 平均(相空间)**,而这里是**完备性**。两者不可混为一谈。

**(c) 超胞 $\Delta V$ ↔ 单缺陷浓度($N_{\rm sc}$)。** 这是最隐蔽、也最可能造成"数过大"的因子。$\Delta V$ 是**一个**缺陷放在 $N_{\rm sc}$ 个原胞的超胞里(缺陷浓度 $=1/N_{\rm sc}$ 每原胞),而 $M$ 用**原胞**波函数配合 EDI 的 `V_folded` 折叠算出。折叠 $V^{q}_{\rm folded}(r)=\sum_R\Delta V(r+R)e^{iqR}$ 在"取出缺陷的傅里叶分量"时,会带进 $N_{\rm sc}$(或 $\sqrt{N_{\rm sc}}$)的因子。

## 2. 决定性判据 —— 闭合(closure)求和规则

由 (b) 的 Parseval,**对所有带求和**(不仅是 rest)在单个通道上应严格满足

$$
\boxed{\;\sum_{n'\,(\text{all})}\big|\langle n'k'|s\rangle\big|^2 = \langle s(k')|s(k')\rangle_{\rm cell} = \sum_G|s(k',G)|^2\;}
$$

其中 $|s\rangle=\Delta V|n k_i\rangle$ 是已验证(T1 $=2.3\times10^{-13}$)的源 ket。再对通道求和:

$$
\sum_{k'}\sum_{n'\,(\text{all})}|M_{n,n'k'}|^2 = \sum_{k'}\langle s(k')|s(k')\rangle = \big\langle n k_i\big|\Delta V^{2}\big|n k_i\big\rangle,
$$

是一个**有限的、原胞尺度的数**(不正比于 $N_k$)。由此得到两个**可直接测量、无歧义**的检验:

- **逐通道闭合:** 单个 $k'$ 上,$\sum_{n'=1}^{150}|M_{n,n'k'}|^2$ 应等于 $\sum_G|s(k',G)|^2$(差额 = 高带不完备,应与 Sternheimer 的高带尾一致)。→ 判定"逐通道**不带**因子"。
- **全和 $=\langle\Delta V^2\rangle$:** $\sum_{k'}\sum_{n'}|M|^2$ 应等于直接实空间算的 $\frac1\Omega\int_{\rm cell}|\psi_{nk_i}|^2\,\Delta V^2$(用 `V_folded^{q=0}` 即可)。→ 判定 **$k'$ 求和的测度**(及里面是否藏着 $N_{\rm sc}$)。

**为什么这能解释"过大":** 当前单通道 $\Delta\tilde V_{\rm chan}\sim0.3$–$0.7$ Ry,144 个通道朴素相加 $\sim70$ Ry。闭合检验会立刻区分两种情况:

1. 若逐通道闭合**成立**(无因子)且 $\langle\Delta V^2\rangle$ 实测真的 $\sim$ 几十 Ry²,那 $\sim70$ Ry 就是**真实**的 —— 意味着空位是**极强散射、Born 彻底失效**(这恰是需要 T-matrix 的理由,重求和后的 $T_{PP}$ 才是物理量);
2. 若闭合给出的逐通道范数比 $0.5$ Ry **小约 $N_{\rm sc}$ 倍**,说明 $M$ 多带了 $\sqrt{N_{\rm sc}}$ 的**超胞折叠因子**,$k'$ 求和(或 $M$ 本身)需除以 $N_{\rm sc}$。

## 3. 金标准 —— Born 极限复现 EDI 迁移率

无论 (a)(b)(c) 如何组合,**唯一能一次性锁死所有因子**的是 observable 层面的一致性:

- 关掉 rest dressing 与 active 重求和 $\Rightarrow$ $\tilde V\to M$、$T_{PP}\to M$;
- 跑 transport,**必须逐位复现当前 EDI 的迁移率**(含 $n_d$、$\tfrac1{N_k}$ 这些 EDI 已验证的因子)。

T1 已在**矩阵元层面**证明 EDT 的 $M$ 与 EDI 完全一致($2.3\times10^{-13}$),所以只要 $k'$ 求和测度定对,Born 极限就会自动复现 EDI —— 这一步把 (a)(b)(c) 全部钉死。

## 4. 判断与建议

- **(b) 纯完备性:** $k'$ 求和**不带** $1/N_k$(正交扇区直和)—— 这点较确定。
- **(c) 超胞因子:** 最可能是"过大"的来源;$1/N_{\rm sc}$(或 $1/\sqrt{N_{\rm sc}}$)需**闭合检验实测**,不能手算拍板。
- **金规则的 $1/N_k$** 属于末态平均,**不是** rest 求和该用的 —— 切勿混用。

**下一步(P3 起点):** 先跑两个便宜的检验 —— (i) 逐通道闭合 $\sum_{n'}|M|^2 \overset?= \lVert s\rVert^2$;(ii) 全和 $\overset?=\langle\Delta V^2\rangle$。它们只复用已验证的 source / overlap 机器,几分钟出结果,直接把测度因子 $c$ 钉死;随后再用 Born 极限复现 EDI 迁移率作总锚定。

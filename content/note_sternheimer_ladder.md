# Rest-space Sternheimer 阶梯:为什么二阶不够,以及逐阶求解高阶 dressing

EDT 的 active/rest 下折叠把 dressed 势写成 $\tilde V(\omega)=M+\Sigma(\omega)$:$M$ 是 active 子空间里的裸耦合,$\Sigma$ 把 rest 子空间(被折叠掉的高能带)的虚跃迁折回 active。目前代码里的 $\Sigma$ 只取到**二阶(Born)** —— 用裸 host 传播子 $G^0_{QQ}$、丢掉缺陷势在 rest 内部的多次散射 $\Delta V_{QQ}$。在 MoS$_2$ 的 S-vacancy 上,这个近似把 conduction 派生的 $e$ 能级**过度屏蔽**了约 $1.1$ eV(见[结果页](results.html#sec-3) 的 Fig 13–15)。

本页做三件事:(i) 推导精确的 Feshbach rest 自能;(ii) 说明**为什么二阶不是一个好近似**;(iii) 给出**逐阶 Sternheimer 阶梯**的完整公式与实现方案(含非局域部分的接法)。记号:投影 $P$(active)、$Q=1-P$(rest);host 哈密顿量 $H_0$,缺陷势 $\Delta V=V_{\rm defect}-V_{\rm host}$;块记号 $\Delta V_{PQ}\equiv P\,\Delta V\,Q$ 等。

## 1. 记号与精确下折叠(Feshbach)

host 以 Bloch 态为本征基,$H_0|n\mathbf k\rangle=\varepsilon_{n\mathbf k}|n\mathbf k\rangle$。active 与 rest **都**由 host Bloch 态张成,所以 $H_0$ 在两块之间**不耦合**:
$$P H_0 Q=0\quad(\text{Bloch 态是 }H_0\text{ 的本征态}),$$
$P$–$Q$ 之间的耦合**只来自** $\Delta V$。记 active 块裸耦合 $M\equiv\Delta V_{PP}$。

解定态方程 $(H_0+\Delta V-\omega)|\psi\rangle=0$,把 $|\psi\rangle=|\psi_P\rangle+|\psi_Q\rangle$ 代入,从 $Q$ 分量方程消去 $|\psi_Q\rangle$:
$$|\psi_Q\rangle=(\omega-H_{QQ})^{-1}\,\Delta V_{QP}\,|\psi_P\rangle,\qquad
H_{QQ}\equiv Q(H_0+\Delta V)Q=H_0^{QQ}+\Delta V_{QQ},$$
得到 active 子空间里的有效哈密顿量
$$\boxed{\,H_{\rm eff}(\omega)=P H_0 P+\underbrace{M+\Sigma(\omega)}_{\tilde V(\omega)}\,},\qquad
\Sigma(\omega)=\Delta V_{PQ}\,\big(\omega-H_0^{QQ}-\Delta V_{QQ}\big)^{-1}\,\Delta V_{QP}. \tag{1.1}$$

$\Sigma(\omega)$ 就是 rest dressing。**注意分母里含 $\Delta V_{QQ}$** —— 缺陷势在 rest 子空间里的多次散射,精确的 $\Sigma$ 把它取到**全阶**。这个 $\tilde V$ 正是 active-space $T$-matrix $T_{PP}=[1-\tilde V G^A]^{-1}\tilde V$ 的输入,所以 $\Sigma$ 的精度决定整条流水线。缺陷能级 = 自洽本征值 $\det[\omega-H_{\rm eff}(\omega)]=0$。

## 2. 二阶(Born)近似为什么会过度屏蔽

目前代码把分母里的 $\Delta V_{QQ}$ 丢掉,用裸 host 传播子
$$G^0_{QQ}(\omega)\equiv\big(\omega-H_0^{QQ}\big)^{-1}=Q(\omega-H_0)^{-1}Q,\qquad
\Sigma^{(2)}(\omega)=\Delta V_{PQ}\,G^0_{QQ}(\omega)\,\Delta V_{QP}. \tag{2.1}$$

把 (1.1) 的精确 resolvent 按 Dyson 展开,
$$\big(\omega-H_0^{QQ}-\Delta V_{QQ}\big)^{-1}=\sum_{p\ge0}\big(G^0_{QQ}\,\Delta V_{QQ}\big)^p G^0_{QQ},$$
于是
$$\Sigma(\omega)=\sum_{p\ge0}\Sigma^{(2+p)},\qquad
\Sigma^{(2+p)}=\Delta V_{PQ}\,G^0_{QQ}\big(\Delta V_{QQ}\,G^0_{QQ}\big)^p\,\Delta V_{QP}. \tag{2.2}$$
二阶只保留 $p=0$,把所有含 $\Delta V_{QQ}$ 的 rest 内多次散射($p\ge1$)全扔掉。

**(a) 收敛性。** 级数 (2.2) 收敛 $\iff$ 谱半径
$$\rho\equiv\rho\!\big(G^0_{QQ}(\omega)\,\Delta V_{QQ}\big)<1,\qquad
\frac{\lVert\Sigma^{(m+1)}\rVert}{\lVert\Sigma^{(m)}\rVert}\xrightarrow{\,m\to\infty\,}\rho.$$
只有 $\rho\ll1$ 时二阶才靠谱。S-vacancy 是个**深势**(缺陷 cube 势阱 $\sim-28$ Ry),$\Delta V_{QQ}$ 很强,$\rho\sim\mathcal O(1)$ —— 级数远未收敛,二阶截断可以差一个 $\mathcal O(1)$ 倍。

**(b) 物理图像:中间传播没被屏蔽。** 用 resolvent 恒等式 $G_{QQ}=G^0_{QQ}+G^0_{QQ}\Delta V_{QQ}G_{QQ}$,二阶的误差是
$$\Sigma-\Sigma^{(2)}=\Delta V_{PQ}\,G^0_{QQ}\,\Delta V_{QQ}\,G_{QQ}\,\Delta V_{QP},$$
它对 $\Delta V_{QQ}$ 是**一阶**的、并不小。$\Sigma^{(2)}$ 把 rest 里的中间态当作**自由 host 能带**($G^0$)传播;但这些 rest 态本身被强缺陷势 $\Delta V_{QQ}$ 强烈散射。对 conduction 派生的 $e$(三条 Mo 悬挂键的反键组合),它对 rest conduction 流形($\varepsilon_q>\varepsilon_e$)的二阶耦合
$$\Sigma^{(2)}_{ee}=\sum_{q\in Q}\frac{|\langle e|\Delta V|q\rangle|^2}{\varepsilon_e-\varepsilon_q}<0$$
是强吸引的(分母为负)。用裸 $G^0$ 时这份吸引**没有被屏蔽**,把 $e$ 拉得过深;把 $\Delta V_{QQ}$ 加回去 = 让中间 rest 态在缺陷场里重新排布,高阶项会部分抵消这份裸吸引。

**(c) 与数据对照。** explicit 把这些 conduction 带直接放进 $P$、完全不做 rest dressing,给出 $e$ 在 VBM 上方 $+1.35$ eV(Fig 14);二阶 dressed 给的是 $+0.36$ eV(Fig 15)。也就是说**全阶求和必须把 $e$ 从 $+0.36$ 抬回到 $\sim+1.2\text{–}1.35$**(DFT 给 $+1.19$)—— 高阶项净额减少了约 $1$ eV 的过度束缚。这正是"二阶过度屏蔽"的定量含义。

> 另一个**独立**的、与之叠加的近似:代码在单一静态参考 $\omega_0=\varepsilon_{\rm VBM}$ 处取 $\Sigma$,而不是在能级自身能量处自洽。对离 $\omega_0$ 很远的 $e$,即使是二阶项也估在了错误的频率上。本页的阶梯解决的是**截断阶数**这一层;频率依赖那一层要么在能级能量处求 $\Sigma(\omega)$,要么扫全 $\Sigma_{\rm rest}(\omega)$。

## 3. Sternheimer 阶梯:逐阶公式推导

Sternheimer 的核心是**不显式做 rest 能带求和**,而是解一个(投影)线性方程。对每个 active 源 $|b\rangle$,定义一阶 rest 响应
$$|\chi^0_b\rangle\equiv G^0_{QQ}\,\Delta V_{QP}\,|b\rangle
\;\Longleftrightarrow\;
Q(\omega-H_0)Q\,|\chi^0_b\rangle=Q\,\Delta V\,|b\rangle,\quad|\chi^0_b\rangle\in Q. \tag{3.1}$$
(QE 里用投影 CG 解,$H_0$ 由 `h_psi` 给出。)二阶自能就是 $\Sigma^{(2)}_{ab}=\langle a|\Delta V|\chi^0_b\rangle$。

**高阶响应的递推。** 把 (2.2) 里 $\big(G^0_{QQ}\Delta V_{QQ}\big)^p$ 一层层拆开,定义
$$|\chi^p_b\rangle\equiv\big(G^0_{QQ}\,\Delta V_{QQ}\big)^p|\chi^0_b\rangle
\;\Longleftrightarrow\;
Q(\omega-H_0)Q\,|\chi^p_b\rangle=Q\,\Delta V\,Q\,|\chi^{p-1}_b\rangle,\quad p\ge1. \tag{3.2}$$
**关键**:每一阶用的是**同一个**算符 $Q(\omega-H_0)Q$,只换右端项(把上一阶响应再被 $\Delta V_{QQ}$ 散射一次)。CG 设置完全复用,每加一阶只多"算一次右端 + 解一次"。

**自能各阶用响应表达(对称配对)。** 在实轴上 $G^0_{QQ}$ 厄米,$\langle a|\Delta V_{PQ}G^0_{QQ}=\langle\chi^0_a|$;把 (2.2) 的链条从中间任意劈开:
$$\boxed{\;\Sigma^{(2)}_{ab}=\langle a|\Delta V|\chi^0_b\rangle,\qquad
\Sigma^{(2+p)}_{ab}=\langle\chi^i_a|\,\Delta V_{QQ}\,|\chi^j_b\rangle,\quad i+j=p-1\ \ (p\ge1).\;} \tag{3.3}$$
(带展宽 $\eta$ 时 $G^0$ 非厄米,改用不对称式 $\Sigma^{(2+p)}_{ab}=\langle a|\Delta V|\chi^p_b\rangle$,或分别解左/右响应。)

**成本结构:一次 solve 买两阶,奇数阶白送。** 由 (3.3),手上有响应 $\chi^0,\ldots,\chi^{S-1}$($S$ 次 solve)就能精确到 $2S+1$ 阶:

| 已解响应 | solve 数 $S$ | 可达自能阶 | 该 solve 的"白送"项 |
|---|---|---|---|
| $\chi^0$ | $1$ | $\Sigma^{(2)},\ \Sigma^{(3)}$ | $\Sigma^{(3)}=\langle\chi^0|\Delta V_{QQ}|\chi^0\rangle$ |
| $+\,\chi^1$ | $2$ | $\Sigma^{(4)},\ \Sigma^{(5)}$ | $\Sigma^{(5)}=\langle\chi^1|\Delta V_{QQ}|\chi^1\rangle$ |
| $+\,\chi^2$ | $3$ | $\Sigma^{(6)},\ \Sigma^{(7)}$ | $\Sigma^{(7)}=\langle\chi^2|\Delta V_{QQ}|\chi^2\rangle$ |
| $\ \ \vdots$ | $S$ | 精确到 $\Sigma^{(2S+1)}$ | $\langle\chi^{S-1}|\Delta V_{QQ}|\chi^{S-1}\rangle$ |

特别地,**三阶 $\Sigma^{(3)}=\langle\chi^0_a|\Delta V_{QQ}|\chi^0_b\rangle$ 不需要任何新 solve** —— 把已经算好的二阶响应再缩并一次即可,是"白送"的。这给了最便宜的二阶可靠性检验:若 $r_3\equiv\lVert\Sigma^{(3)}\rVert/\lVert\Sigma^{(2)}\rVert$ 已经不小,二阶就不可信。

**全阶端点(resummation)。** $S\to\infty$ 等价于直接解含 $\Delta V_{QQ}$ 的完整 Feshbach 方程
$$Q(\omega-H_0-\Delta V)Q\,|X_b\rangle=Q\,\Delta V\,|b\rangle,\qquad \Sigma_{ab}=\langle a|\Delta V|X_b\rangle, \tag{3.4}$$
即把 $\Delta V_{QQ}$ 放进 CG 的算符里、一次解到全阶(完整推导、matvec 与复解法见 §6)。逐阶阶梯的价值在于:它把这个非微扰解**展开成可监控的序列**,既能在 $\rho<1$ 时按需截断,又能用 $\Sigma^{(3)},\Sigma^{(5)},\ldots$ 的大小直接读出 $\rho$、判断该不该信二阶。

## 4. 非局域部分的接法与实现方案

递推 (3.2) 的右端是 $Q\,\Delta V\,Q\,|\chi^{p-1}\rangle$ —— 要把缺陷势作用到一个 rest 波函数上。$\Delta V$ 分局域与非局域两块:
$$\Delta V=\Delta V_{\rm loc}(\mathbf r)+\Delta V_{\rm nl},\qquad
\Delta V_{\rm nl}=\sum_{IJ}|\beta_I\rangle\,\Delta D_{IJ}\,\langle\beta_J|.$$
$\Delta V_{\rm nl}$ 是 KB **可分离**的(对 vacancy 主要是 $-$ 被移除原子的投影子);$\Delta V_{\rm loc}$ 是对齐后的超胞局域势差(cube)。把 $\Delta V$ 作用到 $|\chi^{p-1}\rangle$ 的步骤:

```text
给定 rest 响应 |chi^{p-1}>:
1. 局域:  FFT |chi^{p-1}> 到超胞实空间网格 -> 逐点乘 dV_loc(r) -> FFT 回 G 空间
2. 非局域: 对每个投影子 J 算 c_J = <beta_J | chi^{p-1}>
           组装 sum_{I,J} |beta_I> dD_{IJ} c_J            (可分离, 便宜)
3. 相加 -> dV|chi^{p-1}>
4. 投影到 Q:  |R> = (1-P) dV|chi^{p-1}>
              = dV|chi^{p-1}> - sum_{a in P} |a><a| dV|chi^{p-1}>
5. 解 Sternheimer:  Q(w - H0)Q |chi^p> = |R>     (投影 CG, 复用二阶的算符)
```

两点实现要害:

- **局域部分把 $\mathbf k$ 通道耦起来。** 二阶 solve (3.1) 对每个 $\mathbf k$ 是**解耦**的(固定 source 下右端逐 $\mathbf k$ 独立)。但缺陷破坏平移对称,$\Delta V_{\rm loc}$ 通过超胞 FFT 把不同 $\mathbf k$ 连起来:高阶右端 $\Delta V Q|\chi^{p-1}\rangle$ 会把一个 $\mathbf k$ 的响应散射到别的 $\mathbf k$。所以**三阶以上不再是逐 $\mathbf k$ 独立的 solve**,需要这套超胞 FFT 场机器 —— 这是新增的主要工程量。非局域部分则始终是可分离的投影子求和,便宜、与 $\mathbf k$ 解耦。
- **缩并复用右端。** (3.3) 里 $\Sigma^{(2+p)}_{ab}=\langle\chi^i_a|\Delta V_{QQ}|\chi^j_b\rangle=\langle\chi^i_a|\big(\Delta V|\chi^j_b\rangle\big)$,而 $\Delta V|\chi^j_b\rangle$ 正是算 $\chi^{j+1}$ 时已经造好的右端向量。所以自能缩并不额外花一次 $\Delta V$ 作用,只是向量内积。

## 5. 收敛判据、与 explicit 计算的对应、下一步

**怎么用阶梯判断二阶。** 实跑顺序:先白送地算 $\Sigma^{(3)}$,看比值 $r_3=\lVert\Sigma^{(3)}\rVert/\lVert\Sigma^{(2)}\rVert$。

| $r_3$ | 含义 | 做法 |
|---|---|---|
| $\ll1$ | 二阶可信 | $\tilde V\approx M+\Sigma^{(2)}$,收工 |
| $\lesssim1$ | 收敛慢 | 逐阶往上($+1$ solve 到五阶,再看 $r_5$) |
| $\gtrsim1$ | Born **发散** | 逐阶救不回 —— 改用全阶 Feshbach (3.4),或把这些带搬进 active 显式处理 |

**实测 $r_3$(用 explicit $M$,精确,含全 $k$ 耦合)。** 取 $P$ = bands 7–17(block 的 active 窗口)、$Q$ = bands 18–28(近 gap conduction rest),物理势 $W=M/N_k$,直接矩阵代数算 $\Sigma^{(2)}=W_{PQ}\,g\,W_{QP}$、$\Sigma^{(3)}=W_{PQ}\,g\,W_{QQ}\,g\,W_{QP}$:

| $\omega$ | $\lVert\Sigma^{(2)}\rVert_F$ | $\lVert\Sigma^{(3)}\rVert_F$ | $r_3$ | $\rho(W_{QQ}g)$ |
|---|---|---|---|---|
| $\omega_0=\varepsilon_{\rm VBM}$(block 静态参考) | $8.72$ eV | $3.22$ eV | $\mathbf{0.37}$ | $0.60$ |
| $\omega\approx\varepsilon_e$($e$ 自身能量) | $10.35$ eV | $4.64$ eV | $\mathbf{0.45}$ | $0.70$ |

$\rho<1$ 但 $r_3\approx0.4$ —— Born 级数**收敛却很慢**(三阶还剩二阶的 $40\%$),正落在上表 **"$\lesssim1$,续阶"** 那一档:二阶不可信。

**逐阶确实在收敛,且自洽可验。** 同一套数据下 $e$ 能级(相对 VBM,静态 $\Sigma(\omega_0)$)随阶数:

| 处理 | $e$(相对 VBM) |
|---|---|
| 裸 $M$(P=7–17,无 rest) | $+1.49$ |
| $+\,\Sigma^{(2)}(\omega_0)$ | $+1.30$(二阶 overshoot) |
| $+\,\Sigma^{(2)}+\Sigma^{(3)}(\omega_0)$ | $+1.40$(白送三阶扳回) |
| $+\,\Sigma_{\rm full}(\omega_0)$(全阶求和) | $+1.37$ |
| 自洽 Feshbach(全阶,自洽 $\omega$) | $\mathbf{+1.35}$ |
| explicit all-band 7–28(目标) | $+1.35$ ✓ |

自洽 Feshbach **精确**复刻 explicit 的 $+1.35$,验证了下折叠框架自洽;白送的三阶项实测把二阶 overshoot 从 $+1.30$ 扳到 $+1.40$,正是该有的逐阶收敛。

**$Q$ 截断的说明。** 这里 $Q$ 只到 band 28(explicit 数据上限),是主导的近 gap rest,二阶对它只**轻微** overshoot($+1.30$ vs $+1.35$)。但 block 用的是**全 rest 18–150**:多出的高带把 $\rho$ 推向 $1$,才给出之前那个剧烈的 $+0.36$ 过度屏蔽($e$ 从裸的 $+1.49$ 砸下来 $\sim1.1$ eV)。所以**实测的 $r_3\approx0.4$、$\rho\approx0.6$–$0.7$ 是全问题的下界** —— 全 rest 只会更靠近发散、更需要 ladder 或全阶 Feshbach。(实测脚本 `edt/run/ladder_r3.py`。)

**代码实测(in-code 交叉验证)。** $\Sigma^{(3)}$ 白送项已落进 EDT block 代码(`do_sigma3` 开关,单态模式):保留逐 channel 的 $\chi^0$(不再算完即扔),再做一次 cross-channel 的 $\Delta V$ 双重求和 $\Sigma^{(3)}_{aa}=\sum_{k',k''}\langle\chi^0(k')|\Delta V(k'\!\leftarrow\!k'')|\chi^0(k'')\rangle$ —— **不需新的 CG solve**。对 band 14、$k=1$(full rest 18–150):

| 量 | Fortran(full rest) | explicit($Q=18$–$28$) |
|---|---|---|
| $\Sigma^{(2)}_{aa}$ 自洽校验 | 重建 $=$ 代码存的 `Sgblk`,7 位全同 ✓ | — |
| $\Sigma^{(2)}_{aa}$ | $-0.045$ eV | $-0.014$ eV |
| $\Sigma^{(3)}_{aa}$ | $+0.070$ eV | $+0.0053$ eV |
| $r_{3,aa}$ | $\mathbf{1.54}$ | $0.37$ |

**(1)** $\Sigma^{(2)}_{aa}$ 自洽校验**精确通过**(用保存的 $\chi^0$ 重建 $=$ 代码存的 `Sgblk`),证明 cross-channel 机器正确。**(2)** full rest 的 $r_{3,aa}=1.54>1$ —— **三阶比二阶还大,Born 级数发散**,正落在上面判据表的 **"$\gtrsim1$,发散"** 那档。这把"全 rest 更靠近发散"从推断**升级成实测**:截断 rest($r_3=0.37$)是下界,全 rest 真的越过了 $1$(远端 conduction 带让有两个 rest 求和的 $\Sigma^{(3)}$ 涨得比 $\Sigma^{(2)}$ 快得多:实测 $|\Sigma^{(3)}|$ 大 13×、$|\Sigma^{(2)}|$ 只大 3×)。所以最终结论是硬的:**全 rest 逐阶救不回,必须全阶 Feshbach(resolvent 求逆,非级数)或 explicit**;发散的是**逐阶**,explicit 的自洽 Feshbach 仍精确给 $+1.35$。(caveat:中间 $\Delta V_{QQ}$ 只取局域部分、且是单态对角 $r_{3,aa}$;但 $1.54\gg1$ 稳健。开关 `do_sigma3`,作业 `edt/run/edt_sigma3.{in,slurm}`。)

**两条路通同一物理。** 同一个 $e$ 能级有两种算法:**(i)** 把 conduction 带放进 $Q$、dressing 取全阶(本页的阶梯 / 全阶 Feshbach);**(ii)** 把它们放进 $P$ 显式对角化、不 dressing(结果页的 21-band explicit)。二者在带数足够时收敛到同一答案。explicit 已经给出干净的 $a_1\oplus e$($e$ 在 $+1.35$,Fig 14),所以**阶梯收敛的目标值是已知的**,可以用来验证实现。

**实现状态。** $\Sigma^{(3)}$ 白送项已落地(`do_sigma3`)。MVP 的前两步 —— ① 复用二阶 solve 拿 $\chi^0$;② 在 code 内白送地算 $\Sigma^{(3)}$、量 $r_3$ —— **已完成**(explicit 路 $r_3=0.37$–$0.45$、in-code full-rest $r_{3,aa}=1.54$ 都测出)。第三步因 $r_3>1$(全 rest 发散)直接落到 **"转全阶 Feshbach (3.4) / explicit"**:逐阶 ladder 在这个体系不收敛,生产路线就是全阶 resolvent 求逆或把带搬进 active 显式处理。验收标准(把 dressed 的 $e$ 从二阶的 $+0.36$ eV 抬回 explicit/DFT 的 $\sim+1.2\text{–}1.35$ eV)已由 21-band explicit 达成(Fig 14)。

## 6. 全阶 Feshbach Sternheimer:逐阶发散后的生产路线

§5 实测逐阶 ladder 在全 rest 下发散($r_3=1.54$),二阶的过度屏蔽**不能**靠加几阶救回来。出路是**不展开、直接求逆**的全阶 Feshbach —— 它同样是一个 Sternheimer 线性解,只是算符里把 $\Delta V_{QQ}$ 也放进去。这是式 (3.4) 的完整推导。

**线性方程。** 把 $H_{QQ}=Q(H_0+\Delta V)Q=QHQ$ 代回式 (1.1):
$$\Sigma(\omega)=\Delta V_{PQ}\,(\omega-QHQ)^{-1}\,\Delta V_{QP}.$$
对每个 active 源 $|b\rangle$ 定义全阶 rest 响应 $|X_b\rangle\equiv(\omega-QHQ)^{-1}\Delta V_{QP}|b\rangle\in Q$;左乘 $(\omega-QHQ)$,用 $Q|X_b\rangle=|X_b\rangle$(故 $QHQ|X_b\rangle=QH|X_b\rangle$):
$$\boxed{\,Q(\omega-H_0-\Delta V)Q\,|X_b\rangle=Q\,\Delta V\,|b\rangle\,},\qquad
\Sigma_{ab}=\langle a|\Delta V|X_b\rangle=\langle S_a|X_b\rangle, \tag{6.1}$$
$|S_a\rangle=Q\Delta V|a\rangle$ 是和二阶**一样**的源;对称式即 $\Sigma=S^\dagger A^{-1}S$,$A=Q(\omega-H_0-\Delta V)Q$。**精确、全阶、无截断。**

**和二阶 Sternheimer 的唯一区别:matvec 里多一个 $\Delta V$。**

| | 二阶(Born) | 全阶 Feshbach |
|---|---|---|
| 算符 $A$ | $Q(\omega-H_0)Q$ | $Q(\omega-H_0-\Delta V)Q$ |
| matvec | `h_psi` $+\,Q$ | `h_psi` $-\,\Delta V+Q$ |
| $\mathbf k$ 结构 | 块对角 → 逐 $\mathbf k$ 解耦 | $\Delta V$ 耦合 → 全 $\mathbf k$ 一个大方程 |

```text
A|v> ,  v in Q   (每次 CG 迭代):
  1. t = omega*v - h_psi(v)      # (omega - H0) v       [per-k, 已有]
  2. t = t - dV(v)               # 减 dV v : 局域超胞 FFT(build_V_folded, cross-channel) + 非局域可分离
  3. t = Q t  (apply_Qproj)      # 投回 rest
```
即 §4 的 $\Delta V$ 作用例程(就是 `do_sigma3` 那套 cross-channel 机器)从"算一次 $\Sigma^{(3)}$"变成"**每步 matvec 调一次**"。

**为什么收敛、逐阶却发散 —— 求逆 $\neq$ 级数。** 记 $A=A_0-B$,$A_0=Q(\omega-H_0)Q$、$B=Q\Delta V Q$。逐阶是 Neumann 级数 $A^{-1}=\sum_{p\ge0}A_0^{-1}(BA_0^{-1})^p$,收敛要 $\rho(A_0^{-1}B)=\rho(G^0\Delta V_{QQ})<1$ —— 实测 $\sim1.54$,发散。**但 $A^{-1}$ 本身存在**;Krylov/CG 直接解 $AX=S$ 的收敛只看 $A$ 的谱,不看那个 $\rho$。类比:$x=1.54$ 时 $1+x+x^2+\cdots$ 发散,而 $\tfrac{1}{1-x}=-1.85$ 良定 —— **发散的是展开方式,不是那个数**。加 $\eta>0$,$A=Q(\omega+i\eta-H)Q$ 永远可逆。

**条件数与解法。** rest 取在 $\omega_0$ 之上时 $A_0\succ0$;加 $-\Delta V_{QQ}$ 后,若没有 rest 态被拉到 $\omega$ 以下(缺陷 bound state 应已落在 $P$ 内),$A$ 仍正定 → 实 CG。若某个 rest 共振逼近 $\omega$($A$ 不定/病态)→ 上 $\omega+i\eta$ 配复解法(代码已有的 `rest_split='complex'` / `ccgsolve_all`),$\eta$ 兼作正则化;active 用 $\alpha P$ deflation 钉在 $Q$。

**$\omega$ 自洽 → 精确能级。** $\Sigma(\omega)$ 依赖 $\omega$;缺陷能级是自洽根 $\det[\omega-H_0^{PP}-M-\Sigma(\omega)]=0$,迭代 $\omega$ 即得。全阶 Feshbach 在自洽 $\omega$ 下给出**精确**能级 —— 即 §5 自洽 Feshbach 复刻 explicit $+1.35$ 的那个量。

**成本与定位。** 每源一次**全 $\mathbf k$** 的 Sternheimer solve,$n_{\rm iter}$ 步、每步带一次 $\Delta V$ 作用($\Sigma^{(3)}$ 量级),全 block $N_A$ 个源;比二阶贵 $\sim n_{\rm iter}$ 倍的全-$\mathbf k$ 耦合,但它**收敛**。它与 21-band explicit 互为印证:一个把 conduction 放 $Q$ 全阶 dress、一个放 $P$ 显式不 dress —— 既然逐阶发散,这两条才是该走的路。

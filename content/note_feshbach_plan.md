# 全阶 Feshbach Sternheimer:rest-space 求解的实现计划(v2)

[Ladder 页](sternheimer-ladder.html) §5 实测逐阶展开对全 rest **发散**($r_3=1.54$),§6 推导了不展开、直接求逆的全阶 Feshbach 路线。本页是它的**代码实现计划(v2)**:对 rest contribution 取**静态参考**($\omega\to\omega_0$,取在 gap 内),全阶 $\Sigma$ 化成一个全 $\mathbf k$ 耦合的 Sternheimer 线性解

$$A\,|X_b\rangle=|S_b\rangle,\qquad A=Q\,(H_0+\Delta V-\omega_0)\,Q+\alpha P,\qquad
|S_b\rangle=Q\,\Delta V\,|b\rangle,\qquad \Sigma_{ab}(\omega_0)=-\langle S_a|X_b\rangle .$$

唯一的新物理机器是 matvec 里的 **cross-channel $\Delta V$ 作用子**(`apply_dV`,超胞实现)。第一版(v1)已编码并跑过**隔离验收 —— 抓到一个 $1/N_k^2$ 量级的约定 bug**。本页先做 post-mortem,再给修正后的设计:以"对已验证的 $M$ 矩阵**逐元素**回归"为基石的验收阶梯、零通信的 source-parallel 布局、对算符不定性的诚实分析,以及实测锚定的成本。

## 1. v1 验收的 post-mortem:三个结构性根因

v1(预存 $\Delta V$ 场 + 超胞 FFT 作用子)编译、运行都干净,隔离验收用 `do_sigma3` 的 folded $\Sigma^{(3)}$ 做参考:

| 量 | 值 |
|---|---|
| $\Sigma^{(3)}$,folded 参考(raw) | $1.526\times10^{4}$ |
| $\Sigma^{(3)}$,v1 超胞 FFT(raw) | $6.70\times10^{-1}$ |
| **ratio** | $4.4\times10^{-5}\approx0.91/N_k^2$(对照 $1/N_k^2=4.82\times10^{-5}$,$N_k{=}144$) |

干净的 $\sim1/N_k^2$ 加 ~9% 残差 ⇒ 不是少乘一个因子,是**约定层面的结构错误**。三个根因:

1. **丢了胞内 Bloch 相位 $e^{i\mathbf k\cdot\mathbf r}$。** v1 的展开/折回用 $\sum_k u_k(\mathbf r)\,e^{i\mathbf k\cdot\mathbf R}$,但真实波函数是 $\psi(\mathbf r{+}\mathbf R)=\sum_k e^{i\mathbf k\cdot(\mathbf r+\mathbf R)}u_k(\mathbf r)$ —— 胞内与胞间相位**都要**。`build_V_folded` 的相位约定(`arg=irx*d1`,位置含胞内分数)正是完整的 $\mathbf q\cdot(\mathbf r{+}\mathbf R)$;v1 只配了胞间一半,相干性被破坏。
2. **$\Delta V$ 超胞场的组装方式结构错。** v1 从 144 个 folded 势 $V_f(\mathbf q)$ 做 $\mathbf q$-逆变换 —— 但 $V_f$ 的相位含连续 $\mathbf r$,$\mathbf q$-和不是干净的 $\delta$,产生 aliasing(那 ~9%)。**正确做法根本不经过 $V_f(\mathbf q)$**:$\Delta V$ 的超胞场就是已在内存里的 `V_colin`(对齐后的 cube 数据),只需一次网格重排。
3. **支撑域弄错。** $\Delta V$ 的支撑是 **6×6 cube 区域(36 个原胞)**,不是 12×12 BvK 的 144 个胞(BvK 晶胞 = 2×2 个缺陷超胞,势只在其中一个超胞上非零)。v1 把场摊到 144 胞 —— 既错又贵 4 倍。

**方法论教训:** 标量级验收(一个 $\Sigma^{(3)}$ 数)只能判错、不能定位。v2 的基石改成**逐元素**回归(§4 的 V0)。

## 2. 归一化契约(一次定死,只在文件边界出现)

BvK 归一基 $|kG\rangle$ 下,态 = evc 系数向量按通道堆叠,$\langle k'G'|\Delta V|kG\rangle$ **自带 $1/N_k$**。物理作用子:

$$\big[\Delta V_{\rm phys}X\big]_{k'}=\frac1{N_k}\sum_{\mathbf k}V_f^{\,\mathbf k-\mathbf k'}\!\circ u_{\mathbf k},
\qquad V_f^{\mathbf q}(\mathbf r)=\sum_{\mathbf R\in\text{cube}}\Delta V(\mathbf r{+}\mathbf R)\,e^{i\mathbf q\cdot(\mathbf r+\mathbf R)} .$$

一致性链(全部可验):$M_{\rm phys}=M_{\rm raw}/N_k$(项目既有契约,$H_{\rm eff}$ 已对 DFT 验证);$\Sigma^{(2)}_{\rm phys}={\rm Sgblk}/N_k$(已验);folded 路的 $Y$ **不带** $1/N_k$ ⇒ §4 V1 验收的 ratio 目标是 $1/N_k$ 而非 $1$。**求解器内部一律用 physical 量**($S^{\rm phys}$、$\Delta V_{\rm phys}$、$\Sigma_{\rm phys}=-\langle S^{\rm phys}_a|X_b\rangle$,无游离因子);只在写 `.dat` 时乘回 $N_k$,兼容既有后处理(`H_eff = diag(ε) + (M+Σ)/N_k`)。

## 3. 件 B v2:cube-anchored `apply_dV`

**预备(一次性,~30 行):** 把 `V_colin`(240×240×300)重排成 `dV_blk(nnr_p, 36)`(原胞网格 40×40×300 × 36 胞;启动 assert 公度 $240=6\times40$、$300=300$、k-grid $12\times12\times1$);相位表 $E(k,m)=e^{i\mathbf k\cdot\mathbf R_m}$($144\times36$)。锚点由构造对齐:cube 原点 = 原胞原点、$m=0..5$,胞内+胞间相位合成后逐点等于 `build_V_folded` 的 `irx*d1` 约定。

**每次作用(单源向量,~6–8 s):**

```text
1. 每通道:  ũ_k(r) = e^{i k·r} ∘ invfft(X(:,k))            # 144 个原胞 FFT + 胞内相位
2. 展开:    Psi_blk(nnr,36) = ũ(nnr,144) · E(144,36)        # ZGEMM -> cube 区域的真实波函数
3. 乘势:    W_blk = dV_blk ∘ Psi_blk                         # 逐点 (17M)
4. 折回:    ṽ(nnr,144) = (1/N_k) · W_blk · E^H               # ZGEMM
5. 每通道:  Y(:,k') = gather_G[ fwfft( e^{-i k'·r} ∘ ṽ_k' ) ] # 144 个 FFT, igk_all 映射
6. 非局域:  Y(k') += Σ_I |β_I(k')> D_I · (1/N_k) Σ_k <β_I(k)|X(k)>   # KB 可分离; get_betavkb/make_coeff 复用
7. 逐通道 apply_Qproj
```

**非局域必须在**(它是 $\Delta V$ 的一部分,(7a) 的源就含 nl;v1 的 folded 对照只比局域是验收口径,不是物理口径)。

## 4. 验收阶梯:V0(对 $M$ 块逐元素)是基石

| 级 | 测试 | 判据 | 成本 |
|---|---|---|---|
| **V0(基石)** | `apply_dV`(单位源 $=|b\rangle$ 在其 home-$k$ 通道)与全部 $\langle a|$ 缩并 | $=M_{\rm blk}(a,b)/N_k$ **逐元素**(若干 $b$ 扫全列;局域+非局域、相位、锚点、归一化一次全验) | 秒级,in-run 用内存里的 $M_{\rm blk}$ |
| V1 | 对 folded $\Sigma^{(3)}$(关 nl,局域对局域) | ratio $=1/N_k$ **精确** | 复用 `do_sigma3` harness |
| V2 | 求解器回归:matvec 关掉 $\Delta V$ 项 | $\Sigma$ 逐列 $={\rm Sgblk}/N_k$(对已存 2nd-order block 直接回归) | 小子集 |
| V3 | 物理验收:全解后 $H_{\rm eff}$ | **$e$ 回到 $+1.2$–$1.4$**、$a_1$ 贴 VBM;与 explicit / Koster–Slater 对照 | 全跑 |
| V4 | $\max\lvert\Sigma-\Sigma^\dagger\rvert$;两个 $\omega_0$ 抽查 | Hermitian;$\omega_0$ 行为合理 | 免费 |

V0 不过,后面全不动 —— 相位/锚点类 bug 在 V0 是逐元素的明牌,不再像 v1 那样只看到一个错的标量。

## 5. 求解器与并行布局(两处实质修订)

**(1) source-parallel,matvec 零通信。** 2nd-order 是 pool-per-channel;全阶 matvec 耦合所有通道,沿用通道并行就要每步 `mp_sum` 一个 ~280 MB 的 cube 场 —— 不可接受。改为:**每 rank 持有 $N_A/36\approx44$ 个源的完整全通道向量**($\approx15$ MB/源),`h_psi`(逐通道、按源批量)与 `apply_dV` 全部 rank 内完成,只在最后 `mp_sum` 拼 $\Sigma$。源天然并行 ⇒ 多节点线性加速。新 plumbing:`hpsi_setup_globalk(kg)` —— 用已 gather 的 `igk_all/xkc` 覆写本地 `igk_k` 槽 + 手算 `g2kin` + `get_betavkb` 填 `vkb`(带 save/restore,~40 行),让任意 rank 对任意全局 $k$ 跑 `h_psi`。启动自检:每 rank 对**全部 144 通道**各做一次 $\langle\psi|H_0|\psi\rangle=\varepsilon$ gate(现有 gate 的全局版),立刻暴露 igk/vkb/g2kin 错配。

**(2) 算符不定性:不再赌正定,先测谱再定案。** 之前"rest 在 $\omega_0$ 之上 ⇒ $A\succ0$ ⇒ 实 CG"的断言**不严谨**,两处风险:

- **深带 bands 1–6**(Mo 4s4p 半芯 $\sim-30$ 到 $-60$ eV、S 3s $\sim-12.5$ eV)在 $\omega_0$ 之下,$(H_0-\omega_0)$ 在其上为**负** —— 其实现有 2nd-order 算符就已轻度不定;它一直收敛且全部验证通过,是**实践幸存而非定理**。
- 加上 $\Delta V$(vacancy 核心区 $+11.8$ Ry)后,$Q(H_0{+}\Delta V)Q$ 可能有本征值被推近 $\omega_0$(rest-space 共振)⇒ 近奇异。

策略:沿用验证过的 **deflated block-CG 骨架 + breakdown 监控**($p^\dagger Ap$ 变号即报);smoke run 顺带用 Lanczos 系数**实测** $A|_Q$ 的极端 Ritz 值;若有负方向作怪 ⇒ 换 **MINRES**(对称不定的正解,~100 行、同存储);若近共振病态 ⇒ $\omega_0+i\eta$ 复解法。预条件子沿用每通道 Jacobi $1/\max(g_{\rm kin}^2,1)$。$\alpha=2(\omega_0-{\rm win_{min}})$ 不变。

## 6. 成本与内存(实测锚定)

| 项 | 锚 | 值 |
|---|---|---|
| `apply_dV` / 源·迭代 | 288 个原胞 FFT $\approx2$ s + 2 个 ZGEMM($4.8\times10^5\times144\times36$)$\approx4$ s + nl $\approx1$ s | **~6–8 s** |
| `h_psi` / 源·迭代(批量摊薄) | 2nd-order block(1 h 59 m)的逐通道成本 | ~1 s |
| $n_{\rm iter}$ | 待 smoke 实测($\Delta V$ 核心 $+11.8$ Ry,预条件子只含动能,可能高于 2nd-order) | ~100(假设) |
| **全块 1584 源** | $1584\times100\times7\,{\rm s}/36\ {\rm ranks}$ | **$\approx$8.6 h(1 node)**;零通信 ⇒ $\div N$(4 nodes $\approx$2 h) |

内存/rank:$X$ 批(8 源 × 5 个 CG 向量 × 15 MB)$\approx0.6$ GB + cube 缓冲 $2\times0.28$ GB + `evc_act_all` 165 MB $\approx$ **1.5 GB** ✓。ZGEMM 按小批源合并、`vkb` 每通道每迭代缓存,可再省 1.5–2×。

## 7. 阶段计划、风险与定位

| 阶段 | 内容 | 验收 | 工作量 |
|---|---|---|---|
| **P-I** | `dV_blk` 重排 + `apply_dV` v2(含 nl)+ V0/V1 harness | V0 逐元素过;V1 ratio $=1/N_k$ | ~1–2 天 |
| P-II | `hpsi_setup_globalk` + matvec + 全通道 gate | V2 回归 $={\rm Sgblk}/N_k$ | ~1 天 |
| P-III | block 求解(deflation+监控)→ smoke(288 列,band 13–14)→ 全 1584 列 @ $\omega_0{=}$VBM,2–4 nodes | Hermiticity;$n_{\rm iter}$/谱实测 | ~1 天 + 2–4 h 机时 |
| P-IV | $H_{\rm eff}$ 能级 + 上页 | $e$ 回 $+1.2$–$1.4$,与 explicit/KS 对照 | 半天 |

**风险表:** ① 相位/锚点(v1 的杀手)→ V0 逐元素钉死;② 不定性/迭代数 → 先测谱,MINRES / $+i\eta$ 备选;③ nl 漏项或重复计 → V0 含 nl;④ `h_psi` 全局 $k$ 的模块状态副作用 → save/restore + 全通道 gate;⑤ cube 缓冲内存 → 源批 $\le$4 时上限受控。

**定位。** 若只要缺陷能级,[Koster–Slater](koster-slater.html)(分钟级)与 explicit(46 min)已给出同样的 $a_1{+}e$ —— 本计划的价值在于把 rest dressing **治本**(全阶 $\Delta V_{QQ}$,消除二阶的过度屏蔽/发散),并产出谱函数级的生产对象:静态 $\Sigma(\omega_0)$ block 一次 ~半天/node;将来要整条 $\Sigma(\omega)$ 时,逐 $\omega$ 重复 P-III 即可(多节点线性摊)。

## 8. 执行清单:任务分解与验收判据

§7 的阶段表展开成带验收判据的任务分解(编号 #11–#26 对应工作跟踪器;勾掉一项的条件是其"验收判据"达成,**而不是**代码写完)。依赖链:$11\to12\to13\to14$;$15$ 可与 P-I 并行;$\{14,15\}\to16\to17\to18\to19\to20\to21\to\{22,23\}$;$22\to24$;#25/#26 独立副线。**关键路径是 $11\to12\to13$(V0)** —— V0 过了,其余是接线与跑机。

### P-I 件 B(`apply_dV`)+ 验收 —— 唯一的新物理机器(~1–2 天)

| # | 任务 | 验收判据 |
|---|---|---|
| 11 | `dV_blk` 预备:`V_colin` 重排成 $(40{\times}40{\times}300,\,36\ \text{胞})$ 分块 + 相位表 $E(k,m)=e^{i\mathbf k\cdot\mathbf R_m}$ + 公度 assert($240=6\times40$ 等) | 任取一胞与 `V_colin` 原序逐点一致 |
| 12 | `apply_dV` v2:胞内相位 $\to$ 展开 ZGEMM $\to$ 乘 `dV_blk` $\to$ 折回(含 $1/N_k$)$\to$ 逐通道 FFT/gather $\to$ 非局域(带 nl-off 开关)$\to$ Qproj | 编译过,走 #13/#14 验收 |
| 13 | **V0 基石**:单位源($=|b\rangle$ 在 home-$k$ 通道)经 `apply_dV` 与全部 $\langle a|$ 缩并,~8 个跨 band/$k$ 的 $b$ | **逐元素** $=M_{\rm blk}(a,b)/N_k$,相对误差 ~$10^{-10}$;不过则回 #12,后续冻结 |
| 14 | V1:nl-off 作用于保留的 $\chi^0$,对 folded $\Sigma^{(3)}$ 缩并 | ratio $=1/N_k=6.944\times10^{-3}$ **精确**(v1 失败值 $4.4\times10^{-5}$) |

### P-II 全局 matvec(~1 天)

| # | 任务 | 验收判据 |
|---|---|---|
| 15 | `hpsi_setup_globalk`:igk 槽覆写 + 手算 $g^2_{\rm kin}$ + `get_betavkb` 填 vkb,save/restore(**可与 P-I 并行开工**) | 每 rank 全 144 通道 $\langle\psi_a|H_0|\psi_a\rangle=\varepsilon_a$ gate ~$10^{-9}$ eV |
| 16 | 全阶 matvec 组装:$A=Q(H_0{+}\Delta V{-}\omega_0)Q+\alpha P$,源批 $\le8$/rank,cube 缓冲复用 | 随机 Q-向量对 $\langle Y|AX\rangle=\langle AY|X\rangle^*$ ~$10^{-12}$ |
| 17 | V2 回归:matvec 关掉 $\Delta V$ 项,小子集源解 $\Sigma$ | 逐列 $={\rm Sgblk}/N_k$ 到 CG 阈 —— 新布局(source-parallel + 全局 $k$)的端到端回归 |

### P-III 求解器 + 跑(~1 天 + 机时)

| # | 任务 | 验收判据 |
|---|---|---|
| 18 | `solve_feshbach_block`:block-CG + deflation + **breakdown 监控**($p^\dagger Ap$ 变号即报)+ Ritz 谱诊断;MINRES / $\omega_0{+}i\eta$ 仅在监控触发时上 | 收敛 + 谱报告(定不定性) |
| 19 | `do_feshbach` driver:源构建、source-parallel 分发($\approx44$ 源/rank)、$\Sigma_{ab}=-\langle S^{\rm phys}_a|X_b\rangle$ 组装、写 `.dat`(边界乘回 $N_k$,同格式) | 内存 ~1.5 GB/rank 核实 |
| 20 | smoke run:288 列(band 13–14 全 $k$),1 node | 实测 $n_{\rm iter}$(假设 ~100;$>300$ 则改预条件/上 $+i\eta$)、Ritz 谱定解法、墙钟校准 |
| 21 | 生产跑:全 1584 列 @ $\omega_0{=}$VBM,2–4 个 **exclusive** node(36 ranks/node,不共享) | Hermitize 前 $\max\lvert\Sigma-\Sigma^\dagger\rvert$ 与 2nd-order 的 $5.9\times10^{-14}$ 同量级;~2–4 h |

### P-IV 物理验收 + 上页(~1 天)

| # | 任务 | 验收判据 |
|---|---|---|
| 22 | $H_{\rm eff}$ 能级 + ΔDOS 第四 panel(与 Fig 15 同网格同 $\eta$) | **$e$ 从二阶的 $+0.36$ 回到 $+1.2$–$1.4$**(explicit $+1.35$ / DFT $+1.19$),$a_1$ 贴 VBM、简并度 $1{+}2$ —— 整条线的最终验收 |
| 23 | (可选)第二参考 $\omega_0'\approx-4.6$ 生产跑 + static-$\omega_0$ 误差评估 | 两个 $\omega_0$ 的 $e$ 差 $<0.1$ eV $\Rightarrow$ static 近似定案 |
| 24 | 上页收尾:结果图、Fig 13/15 叙事补全阶结论、本页 Drafted$\to$实测数字、catalog 翻 prod、记忆更新、push | 构建自检 + 发布集干净 |

### 副线:Koster–Slater(独立,可随时并行)

| # | 任务 | 验收判据 |
|---|---|---|
| 25 | host GF 加带收敛:扩大 Wannier 带空间,跑 $e(N_{\rm bands})$ 收敛曲线 | $e$ 从 $+1.50\to+1.35$/$+1.19$;Krein–Friedel 计数 $2.06\to3$ |
| 26 | (可选)$C_{3v}$ 特征标:$\sigma_v$ 字符严格定 $a_1$ vs $a_2$ + irrep 分块 | 撤掉 [KS 页](koster-slater.html) "by analogy" 的 caveat;顺带 $6$–$12\times$ 对称性加速 |

**状态(2026-06-09):** 全部 pending;当前未被阻塞、可立即开工的是 **#11、#15**(Feshbach 主线的并行起点)与 **#25、#26**(副线)。

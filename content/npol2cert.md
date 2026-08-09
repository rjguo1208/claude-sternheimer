# npol=2 认证:翻倍测试的裁决与谱系

## 0. 摘要

MODE C 下折叠链的 npol=2(旋量)移植在翻倍测试(noncolin=T、lspinorb=F、标量赝势)中
未能复现"能级精确 ×2",触发了一次完整的司法式裁决。两幕结论:第一幕,**EDT 的
npol=2 实现无罪**——三个独立规范不变判据全部机器精度通过;第二幕(按对仪器化转储后
的大反转),**EDI 的计算引擎同样无罪(逐对 $10^{-12}$ 精确),真凶是文件接口:
EDI-direct 的 .bin 在每个 $(k_i,k_f)$ 块内写 $(\mathrm{bra},\mathrm{ket})$,
而 EDT 读取约定是 $(\mathrm{ket},\mathrm{bra})$——纯带轴转置**,数值判决
$\max|E_\mathrm{EDI}-T_\mathrm{EDT}^{T(\mathrm{bands})}|=8.7\times10^{-14}$。
修复:EDI 写入端转置(格式版本 20260702→20260809);EDT 侧则把 born 块 $M_{AA}$
移植进本体(npol=2),SOC 生产线自给自足。终判:旋量对简并 **0.0000 meV**,
跨宿主 VBM 门控 ≤4.9 meV(宿主孪生底噪),链-vs-born 闭环 $4.6\times10^{-11}$,
修复后 EDI-vs-EDT 同波函数全矩阵 $\sim10^{-13}$。

## 1. 翻倍测试与最初的失败

翻倍测试构造:同一 MoS$_2$ V$_\mathrm{S}$ 6×6 体系,宿主用 noncolin=T、lspinorb=F 与
**标量赝势**重算(d66nc)。此时 $H$ 自旋对角、无磁化,物理与标量计算完全相同:
下折叠有效哈密顿量的每个能级都应精确出现两次。首次运行(EDI-direct 提供 born 块
$M_{AA}$)结果 FAIL:能级既不成对也不对应,链的算符恒等测试
$\max|P^\dagger(\Delta V\,\Psi)-M|$ 高达 $8.7\times10^{-1}$(标量情形为 $10^{-14}$ 量级)。

## 2. 三方分解:局域臂与 KB 臂分开审

用两个探针钩子把 apply_dV 拆开:`EDT_NOKB=1` 关闭 KB(非局域)通道得纯局域投影,
再以独立编写的 numpy 平面波积分器(裁判,零依赖重实现:补零 ifftn、$e^{iq\cdot r}$ 相位、
双旋量分量求和)对 160 个复矩阵元逐一对账:

| 审计对象 | 判据 | 结果 |
|---|---|---|
| EDT 旋量局域 fold | $\|P_\mathrm{nokb}\|/\|裁判\|$ | 52/52 元素比值 **1.000**(全部 3% 内) |
| EDT 旋量 KB | $\|P_\mathrm{full}-P_\mathrm{nokb}\|/\|M_\mathrm{EDI}-裁判\|$ | 中位 0.36,p10=0,散布 0.13–1.3 |

局域臂当场无罪(且裁判与 EDT 互验)。KB 臂的"罪证"随后被证明是冤案——分母
$M_\mathrm{EDI}-裁判$ 本身被污染。

## 3. 关键证据:混合旋量与简并对平行性

d66nc 的本征态是**重度混合旋量**(纯度 0.51–0.99:Davidson 在简并对内任意 SU(2) 混合,
纯旋量假设不成立)。对无 SO 的自旋对角算符,真实矩阵元满足

$$\langle m|\mathcal{O}|n\rangle = \langle\chi_m|\chi_n\rangle\,\mathcal{O}_{\mathrm{orb}}(m,n),$$

于是对简并对 $(m_1,m_2)$(同轨道、正交旋量),**KB 部分的 2-矢量
$[K(m_1,n),K(m_2,n)]$ 必须与局域 2-矢量 $[L(m_1,n),L(m_2,n)]$ 复比例平行**
(都 $\propto[\langle\chi_{m_1}|\chi_n\rangle,\langle\chi_{m_2}|\chi_n\rangle]$)。
这是规范不变的终极判据——对 per-k 幺正与简并混合完全免疫:

| 被审 KB | 平行度 $|\cos|$ | 通过率 (>0.99) |
|---|---|---|
| EDT($P_\mathrm{full}-P_\mathrm{nokb}$) | 中位 **1.0000**,p10 1.000 | **25/25** |
| EDI($M_\mathrm{EDI}-裁判$) | 中位 0.977,p10 0.759 | 44% |

同时,"EDI 总残差"$M_\mathrm{EDI}^{\rm file}-(裁判+K_\mathrm{EDT})$ 高达 0.50
(而 $\max|M|=0.42$)且与局域幅度相关(0.69)——此刻的证据链指向"EDI 矩阵元
局域+KB 皆不可信",据此把 born 块移植进了 EDT(§4)。**但这一定罪的机制随后被
推翻**(§7):按对仪器化转储证明 EDI 引擎逐对精确,坏的是文件接口的带轴约定。
(审读中顺带发现的 `ed_coarse_full_q` noncolin 分支 $\sigma_1$-only 缩并是
真实的代码错误并已修复,但它属于插值路径,不在本案传播链上。)

**方法论教训(记录在案)**:(1) 此前对 EDI 的 SV"平反"(中位数 6.5e-3)之所以
成立又不可信,根本原因是**奇异值对转置不变**——SV 类判据对这类接口错误天然全盲;
(2) 全矩阵 Hermiticity 使 $|M^T|$ 与 $|M|$ 幅度谱相同,"中位比值 1.0000 而
p10/p90 大散布"正是转置的指纹;(3) **永远看全分布,永远补逐元素复数级对账。**

## 4. 修复:born 块进 EDT,noncolin 全流程去 EDI 化

把 `vtilde_block_mpi`(born 路径,`born_only=.true.`)移植 npol=2:源矢量 $S$、
bec/coeff 数组加自旋维,fold 加分量循环($\Delta V_\mathrm{loc}$ 自旋对角),
KB 走与 apply_dV 相同的 `make_coeff_nc`(dvan 自旋对角 / dvan_so 2×2),
bra 侧逐分量收缩;npol=1 逐位不变。验证梯(kbval):

| 关卡 | 结果 |
|---|---|
| 标量 born 回归(r9 vs 存档) | **逐位相同**(8,296,152 B) |
| nc born(40 带 × 36 k,$N_A=1440$) | Hermiticity $1.1\times10^{-14}$ |
| 链-vs-born 闭环 unit test | $\mathbf{4.6\times10^{-11}}$(对 EDI edmat 时 0.87) |
| SV-翻倍审计,**全部** 1260 离对角块 | 中位 3.5e-3,p90 7.8e-3,max 1.1e-2(双 SCF 宿主差水平,无坏尾) |

## 5. 终判:双门框架

原"0.2 meV 内精确 ×2"隐含**同一哈密顿量**假设。实际上标量宿主(dout66f)与
noncolin 宿主(d66nc)是两次独立 SCF,QE 的 noncolin-XC($m\to0$ 路径)使原始本征值
差 3–52 meV、VBM 差 10.8 meV——异宿主上 0.2 meV 原则上不可达。正确框架:

| 门 | 性质 | 结果 |
|---|---|---|
| 旋量对简并 $\max|e_{2j}-e_{2j-1}|$ | 同宿主,**精确** | **0.0000 meV — PASS** |
| 各自 VBM 门控后 scalar-vs-nc | 跨宿主,受孪生底噪限 | ≤4.86 meV — PASS |

叠加 §2–§4 的机器精度内证,npol=2 EDT(链 + born)达到生产认证。
验收与取证脚本存于 qe-edt 仓库 `post/`:`doubling_accept.py`、`parallel_test.py`、
`analyze_split.py`、`svdouble2.py`。

## 6. 对 SOC 生产的直接后果

- SOC(lspinorb=T)链一律使用 **EDT 自产 born edmat**(`do_full_block=.true.,
  born_only=.true.`),禁止投喂 EDI-direct 的 noncolin/SOC edmat。
- 顺带发现旧 SOC O$_\mathrm{S}$ 超胞是**未弛豫几何**(O 替位 $z$ 差 1.0 bohr),
  与非 SOC 生产线(弛豫)不一致——四只 SOC 超胞(洁净、O$_\mathrm{S}$、V$_\mathrm{S}$、
  Se$_\mathrm{S}$)已按战役几何重算排队。
- 50 vs 100 Ry 收敛检查:SOC 劈裂本身 50 Ry 已收敛(K-VB 149.0 meV、K-CB 3.0 meV,
  差 <0.1 meV),但活性窗绝对能量漂移 max 16 meV → SOC 栈维持 100 Ry。

## 7. 第二幕:按对仪器化转储与接口真凶

把 EDI 装上 `EDI_DBG_PAIR` 探针(对指定 $(k_i,k_f)$ 对转储装配前的
$m_\mathrm{loc}$ 与 $m_\mathrm{nl}$ 矩阵),与裁判/已认证的 EDT KB 三方对账:

| 对账(对 $(k_1,k_{29})$,80 元素) | 结果 |
|---|---|
| EDI $m_\mathrm{loc}$ vs 独立裁判 | $6.9\times10^{-12}$ |
| EDI $m_\mathrm{nl}$ vs EDT KB(已认证) | $9.5\times10^{-14}$ |
| 但整装 .bin vs EDT-born 全矩阵 | max 1.68 |

数值全对、摆放全错 → 罪在装配/写盘。7 种轴假设一次扫描:

$$\max\bigl|E_\mathrm{EDI} - T_\mathrm{EDT}^{T(\mathrm{bands})}\bigr| = 8.7\times10^{-14}$$

——块内带轴纯转置(无 $k$ 交换、无共轭)。修复为 EDI 写入端转置 + 格式版本
20260809(edi_v8d);**所有 20260702 版的 EDI-direct .bin 都携带带转置**,须重
生成或读时转置(conv.py 产的同版本号文件是 EDT 约定,不受影响)。同批修复:
B8($q_\mathrm{cryst}$ 改用 NSCF 表示,防 $e^{i\Delta G\cdot r}$ 假相位)与
full_q $\sigma_1$-only。修复后同波函数 EDI-vs-EDT 全矩阵闭环与标量跨码
unit test 见验证表(运行中补录)。

**总教训**:跨码对不上时,第一步不是审引擎,而是**按对仪器化转储把
"引擎 vs 接口"分离**——两个各自正确的引擎完全可以隔着一个转置互相"定罪"。

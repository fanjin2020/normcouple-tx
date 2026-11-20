# 论文3：审稿补实验一键程序

## 服务器一条命令

将 `paper3_revision_oneclick.py` 放到服务器工作目录，执行：

```bash
python paper3_revision_oneclick.py --setup --project-root /home/research/git/repository/paper1_public_evidence --output /home/research/git/repository/paper3_revision_results_v1 --workers 8
```

只需传一个 `.py`，不需要ZIP。它会展开校验过的源码，建立私有Python环境，并安装缺少的R软件包。**不下载研究数据，不修改旧稿、旧结果或AI声明。** 需要已有Python3.10+（建议3.11/3.12）、R4.4+及Linux编译工具链。首次安装软件需要联网，与重新下载研究数据不同。双P40不参与计算，默认8个CPU进程。

默认200次重复，原来的12个内部参考-目标单元与30个外部单元。外部仍是已有GSE174367，不增加数据集。首次安装和正式模拟可能耗时较长，软件测试时间不能用来推算正式运行时间。

## 2026-10-13 汇总修复：已完成实验不要重跑

版本 1.0.2 修复 CSV 把 `null` 识别成缺失值的问题，并统一条件指标的供体等权及配对 MC 区间。新增每个计划条件、方法、重复和基因家族完整性检查。原采样、随机种子、归一化、去卷积和 HC3 模型均未改变。

**已经完成 results_v2 的用户，只运行下面的小型修复脚本，不运行完整模拟入口、不加 `--setup`：**

```bash
/home/research/git/repository/paper1_public_evidence/paper3_revision_source_5f1ccb4fd07c/.venv/bin/python \
  paper3_repair_summary_oneclick.py \
  --source /home/research/git/repository/paper3_revision_results_v2 \
  --output /home/research/git/repository/paper3_final_results_fixed
```

`--source` 必须是服务器原始完整结果目录，不能是轻量导出目录。脚本读取并校验保存的 chunk 结果，恢复全部条件；经计划、真值和基因家族核验后恢复旧逐基因表的空 null 标签。不会重新采样、调用 R、下载安装、删除或覆盖原结果，也不改原 protocol/complete 签名。目标目录必须不存在。

看到 `REPAIR FINISHED` 后，将 `paper3_final_results_fixed` 文件夹直接传回，不用 ZIP。默认不复制逐次大表、计数矩阵、R 环境、缓存或拟合数组；保留最终汇总、全部逐基因结果、机制记录、图表、失败原因和修复审计。完整原始结果继续留在服务器。

`SUMMARY_PROVENANCE.json` 区分原实验版本与修复汇总版本；复制的 `protocol.json` 是原实验记录，不表示用新版重新运行实验。源级和配对区间均重新计算，不能沿用旧条件指标区间。MC 算法及缺失处理写入 `SUMMARY_METHODS.md`；有供体完全无有效值时，该条件指标标记为无法估计，不静默删除供体。

以后从头运行的新实验使用更新后的 `paper3_revision_oneclick.py`，需要新的输出目录。不要用新版完整模拟入口继续旧版本结果，也不要编辑旧签名绕过校验。

## 2026-10-12 第二次修复：分组整数类型（历史记录）

新版1.0.1修复真实输入的pandas分组编号为 `numpy.int64`、冻结随机种子生成器无法JSON序列化的问题。在任务准备入口统一转为Python整数；不转成字符串，不修改冻结引擎或随机种子算法。回归测试覆盖内部／外部任务、pandas实际分组类型，以及与原生整数任务的抽样和原引擎数值一致性。

如果已完成R安装并在 `int64 is not JSON serializable` 处停止，上传最新单文件脚本后直接使用原私有Python运行，**不加 `--setup`，无需重新安装依赖**：

```bash
/home/research/git/repository/paper1_public_evidence/paper3_revision_source_5f1ccb4fd07c/.venv/bin/python \
  paper3_revision_oneclick.py \
  --project-root /home/research/git/repository/paper1_public_evidence \
  --output /home/research/git/repository/paper3_revision_results_v2 \
  --r-library /home/research/git/repository/paper3_revision_results_v1/R_library \
  --workers 8
```

此时旧输出已有协议记录，必须使用新的 `results_v2` 目录以免混用不同源码签名。保留旧目录里的R库及日志。数据仍复用现有缓存，只重新校验和读取，不重新下载。`RUN` 表示批次开始，不表示完成；单元完成会显示 `DONE`，完整运行结束显示 `FINISHED`。

## 2026-10-12 第一次修复：R安装

针对服务器日志中的两项失败，安装层已修改，实验算法与参数未改：

- `--setup` 检查 CMake；缺失或版本过低时，在私有Python环境安装 `cmake==3.31.10` 二进制wheel，并把其目录传给R的PATH。不需要sudo。
- 不再把 `BiocManager` 作为安装前置条件，也不再访问其在线 `config.yaml`。通过R自带的 `install.packages` 从明确兼容的Bioconductor版本仓库解析并安装依赖：R4.4→3.20、R4.5→3.22、R4.6→3.23。拒绝未知R版本，HTTPS验证保持开启。
- 已安装成功的依赖保留；失败时打印最后100行R日志，并保留上次日志副本。此修改不保证服务器上所有网络或编译问题均已排除。
- 实际安装版本仍记入运行记录。这是兼容的软件仓库选择，不代表所有R包都锁定为本地测试的版本。

对于日志中原来展开到 `paper3_revision_source_5f1ccb4fd07c`、尚未开始实验的这次运行，上传更新后的单文件脚本到原服务器目录，然后执行：

```bash
PAPER3_ENV_DIR=/home/research/git/repository/paper1_public_evidence/paper3_revision_source_5f1ccb4fd07c/.venv \
python paper3_revision_oneclick.py --setup \
  --project-root /home/research/git/repository/paper1_public_evidence \
  --output /home/research/git/repository/paper3_revision_results_v1 \
  --workers 8
```

这里指定旧Python环境，是为了跳过已经完成的大型Python包下载。不要删除旧源码目录、`.venv`或结果目录。此次在创建实验协议前失败，因此可沿用原输出目录；若之后已开始正式实验，再更换代码时仍需遵守原有协议校验规则。

依据：[官方R/Bioconductor对应表](https://bioconductor.org/about/release-announcements/)、[CMake wheel](https://pypi.org/project/cmake/3.31.10/)。

## 新增实验

1. 伪计数0.1、1、10，CPM倍率10^4、10^6、10^8；联合缩放等价性与纯对数正值测试；同参数下配对比较观察／扰动前固定分母。
2. 官方edgeR TMM与DESeq2 `poscounts`尺寸因子，不使用改变基因标签估计因子。
3. 弱效应0.1／0.25、原效应0.5、期望新增总量平衡的双向效应、低／高参考表达基因集。
4. 官方MEAD1.0.2点估计进入相同下游回归；原生输出及参考产量转换均保留。
5. FWL核验、完整失败记录、条件／交付覆盖率、供体异质性、配对MC误差与断点校验。

**MEAD边界：** 调用未修改的 `MEAD_getX`、`MEAD_est`，不是仿制算法。但本实验只检验其点估计作为组成协变量的效果，下游仍用HC3；不宣称已传播MEAD比例协方差，也不是完整默认预处理基准。详细定义见 `PROTOCOL.md`。

## 路径和断点

- 中断后原命令重跑。完成批次通过SHA-256核验后跳过；代码、参数、包版本改变时需要新输出目录。
- 复用已知LIBD／SEA-AD预处理计数和映射；扫描项目及工作目录中的GSE174367缓存。文件缺失或哈希不符会停止，不悄悄换数据。
- 可用 `--libd-dir /path --seaad-dir /path --gene-map /path/map.csv --external-cache /path/cache` 指定输入。
- `--stage2-results /path/to/old/results` 提供旧结果目录以定位缓存和额外核对；旧结果不是必需输入。
- `--smoke --output paper3_revision_smoke`：真实数据少量供体、3次重复，只用于软件检查。
- `--synthetic --output paper3_revision_software_test --workers 1`：人工数据端到端检查，会实际调用R方法，不是科研证据。
- `--preflight`：只核对输入和依赖。`--report-only`：只汇总已完成且校验通过的结果。
- 软件失败看 `RUN_STATUS.json`、`install_R.log`、`R_check.log` 或 `units/*/chunk_*/bridge.log`。算法失败保留为不可用结果；依赖等基础设施错误直接停止。
- 镜像可通过 `PAPER3_CRAN`、`PAPER3_BIO_MIRROR` 指定；默认Bioconductor使用官方列出的Posit镜像。

## 发回哪些结果

完整保留 `paper3_revision_results_v1` 文件夹即可，无需自行压成ZIP。至少发回 `RUN_STATUS.json`、`protocol.json`、`tables/`、`reference_inputs/` 和 `units/*/chunk_*/r_controls.json`。

- `tables/method_status.csv`、`failure_reasons.csv`：先看可用性，再比较表现。
- `tables/source_summary.csv`：来源分层均值及固定供体库条件下的MC误差。
- `tables/paired_summary.csv`：配对方法差，以及同一伪计数／倍率下观察减固定分母的差。
- `tables/donor_summary.csv`：供体差异；`per_gene_summary.csv`：偏差、经验SD、估计SE、覆盖率、可用性。
- `tables/mechanism_replicates.csv`：注入前后预测与未改变基因系数的变化。
- `units/*/chunk_*/normalization_factors.csv`：实际使用的分母。
- `fit_arrays.npy`：可审计中间数组，不是投稿系统附件；含义由同目录 `fit_index.csv` 和源码 `ARRAY_KEYS` 给出。

## 公开复现输入

新实验不再强制依赖旧结果CSV。可独立导出字节不变的预处理文件与哈希清单：

```bash
python paper3_revision_oneclick.py --project-root /home/research/git/repository/paper1_public_evidence --export-inputs /home/research/git/repository/paper3_verified_inputs --output /home/research/git/repository/paper3_export_receipt
```

另一台机器用 `--input-bundle /path/paper3_verified_inputs` 即可。**本地导出不等于公开存储；仍需核对原数据许可并提供真实可访问地址。** 本程序不补造缺失的原始预处理历史。

## 源码分工

`revision_core.py`：输入和计数干预；`r_bridge.R`：官方R方法；`revision_metrics.py`：推断指标；`revision_driver.py`：调度、断点与汇总；`test_revision.py`：单元测试；`vendor/`：原引擎与固定版本源码。来源、许可见 `THIRD_PARTY.md`。

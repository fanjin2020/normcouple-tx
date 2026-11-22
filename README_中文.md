# NormCouple-Tx：转录组归一化耦合诊断框架

本地软件版本：2.3.0。CBC 新增实验及绘图统一由 `bash run_cbc.sh` 启动。[一步命令、输入与结果](docs/CBC_CN.md)。原六项核查入口为 `bash run_evidence.sh`。[原证据核查说明](docs/EVIDENCE.md)

公开仓库：[fanjin2020/normcouple-tx](https://github.com/fanjin2020/normcouple-tx)。

## 已有数据复用

CBC 增强使用首稿服务器记录的 Python 3.12.13 和 [requirements.txt](requirements.txt) 中的核心计算库版本。先运行 `python reproduce.py check-enhancement-env` 核对实际版本。已归档的实验按原记录使用 [requirements-archived.txt](requirements-archived.txt)。程序不会自动安装或改动环境。`openpyxl` 只用于可选的 Excel 导出。

下面的 `/data/...` 是示例路径，必须换成服务器上的实际目录：

```bash
python reproduce.py prepare-inputs --output ../paper3_inputs \
  --cache-root /data/paper3_prepared_inputs --cache-root /data/public_bulk_cache
python reproduce.py check-inputs --inputs ../paper3_inputs
```

已有文件逐个核对SHA256后导入，原缓存不修改。不加 `--allow-downloads` 就不会联网下载。LIBD、SEA-AD计数库和基因映射也可从固定的作者原始文件重建，7个输出文件已实测与论文输入逐字节一致：

```bash
python reproduce.py rebuild-inputs --source-dir /data/public_raw \
  --output ../paper3_inputs --r-library /data/R_library --allow-downloads
python reproduce.py prepare-inputs --output ../paper3_inputs --allow-downloads
```

源文件地址、哈希、类别映射、抽样种子均随包登记。已有原始下载可放在source-dir中复用；第二条命令补齐GEO输入。

## 一次运行

```bash
python reproduce.py run --inputs ../paper3_inputs --work ../paper3_run \
  --r-library /data/R_library --workers 36
```

默认运行九个阶段：归一化、机制与功效、标志基因扰动网格、随机面板、耦合、共同支持、标志基因疾病关联、跨估计器模拟、双队列临床对照。最终汇总读取本次计算结果。跨估计器比较涵盖 dtangle、MuSiC 和线性 SVR，软件版本见 docs/INPUTS.md。

```bash
python reproduce.py status --work ../paper3_run
python reproduce.py run --inputs ../paper3_inputs --work ../paper3_run \
  --r-library /data/R_library --workers 36 --resume
```

`run_server.sh` 是同一入口的shell包装，使用当前Python，也可设置 `PAPER3_PYTHON` 指定现有解释器。完整计算在Linux运行，本地Windows可做数据导出、绘图和软件测试。

## 不重算，只复现投稿图表

```bash
python reproduce.py export --published --formats pdf png --workbook --output ../paper3_submission_v222
```

这条命令不下载数据、不重新拟合模型。它导出9张主图、11张补图和含 Guide、S1–S21 的 Excel。PDF 用作独立图件，300 dpi PNG 可直接嵌入 Word。仅需要补表和 CSV 时，用 `python reproduce.py tables --published --workbook --output ../paper3_tables`；两种导出方式的图号和数据来源一致。

`export --work ../paper3_run` 从计算结果生成 S1–S19；加 `--evidence ../paper3_review/evidence` 可从原记录重建 S20 和 S21。最小离线示例：`python reproduce.py example --output ../normcouple_example`。

## 最终输出

- 存档导出含21张补充表。S20给出单基因、五基因和全家族干预的实测剂量与阈值跨越数；S21给出按供体等权汇总的组织块成像误差。
- 另有11个绘图或诊断数据文件，其中1个为gzip压缩；共32个数值文件。
- 9张主图和11张补图，默认仅PDF。完整图9使用S13和S20，新补图S11使用S21。没有实测证据的计算导出只画图9的有限剂量结果，不生成S11。
- 图号与CSV对应表、表题/列说明、来源及SHA256清单。

表格与绘图共用同一份CSV，不再重复输出同值的“表格版”和“绘图版”。仅指定 `--workbook` 时生成Excel；仅指定 `--formats pdf tiff` 等选项时生成其他图像格式。日志、矩阵、checkpoint在工作目录保留供续跑使用，不放进最终结果目录，也不装进GitHub包。

```bash
python reproduce.py verify
python -m unittest discover -s tests -v
python reproduce.py package --output ../normcouple-tx.zip
```

源码统一保存在 `pipeline/`、`kernels/`、`publication/`。科学模块按内容去重，运行时还原内部导入结构，保持算法、随机数与原校验身份。详细要求见英文README和 `docs/`。

自有数据可用同一入口的`sensitivity`子命令，输入固定结局、设计和多组组成估计，输出系数、标准误、权重距离与BH状态变化。它描述指定设计下的敏感性，不把阈值跨越解释为真实临床假阳性。完整示例见`docs/fixed_response_audit.md`。

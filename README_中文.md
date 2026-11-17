# NormCouple-Tx：转录组归一化耦合诊断框架

版本：2.1.0。统一入口为 `reproduce.py`，不需要分别运行修复、增强或汇总脚本。

公开仓库：[fanjin2020/normcouple-tx](https://github.com/fanjin2020/normcouple-tx)。

## 已有数据复用

先激活已有计算环境。统计重跑要求 `requirements.txt` 的固定版本；程序不会自动安装或改动环境。`openpyxl` 只用于可选的Excel导出。

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

默认依次完成8个阶段：归一化比较、机制/功效分析、marker扰动网格、随机marker面板、耦合与真实bulk分析、共同支持及多步临床分析、marker疾病关联、dtangle与MuSiC跨估计器对照。前置结果由同一次运行生成，最终汇总读取新计算结果，不会用随包存档的论文数字替代缺失阶段。MuSiC固定版本要求见docs/INPUTS.md，可用--music-library追加已有R库。

```bash
python reproduce.py status --work ../paper3_run
python reproduce.py run --inputs ../paper3_inputs --work ../paper3_run \
  --r-library /data/R_library --workers 36 --resume
```

`run_server.sh` 是同一入口的shell包装，使用当前Python，也可设置 `PAPER3_PYTHON` 指定现有解释器。完整计算在Linux运行，本地Windows可做数据导出、绘图和软件测试。

## 不重算，只复现投稿图表

```bash
python reproduce.py export --published --output ../paper3_published
```

这里明确使用存档论文数据，输出清单标记为 `published_data_replay`，不等同于从原始计数重算。已经完成的新计算使用 `export --work ../paper3_run`。

## 最终输出

- 16个补充表CSV：TableS1至TableS16。
- 11个额外绘图/诊断数据文件，其中1个gzip压缩；共27个数值文件。
- 9张主图和8张补图，默认仅PDF。
- 图号与CSV对应表、表题/列说明、来源及SHA256清单。

表格与绘图共用同一份CSV，不再重复输出同值的“表格版”和“绘图版”。仅指定 `--workbook` 时生成Excel；仅指定 `--formats pdf tiff` 等选项时生成其他图像格式。日志、矩阵、checkpoint在工作目录保留供续跑使用，不放进最终结果目录，也不装进GitHub包。

```bash
python reproduce.py verify
python -m unittest discover -s tests -v
python reproduce.py package --output ../paper3-normalization-coupling_v2.1.0.zip
```

源码统一保存在 `pipeline/`、`kernels/`、`publication/`。科学模块按内容去重，运行时还原内部导入结构，保持算法、随机数与原校验身份。详细要求见英文README和 `docs/`。

自有数据可用同一入口的`sensitivity`子命令，输入固定结局、设计和多组组成估计，输出系数、标准误、权重距离与BH状态变化。它描述指定设计下的敏感性，不把阈值跨越解释为真实临床假阳性。完整示例见`docs/fixed_response_audit.md`。

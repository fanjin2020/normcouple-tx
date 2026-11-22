# 论文3 CBC 一步运行

在解压后的代码目录执行。CBC 使用首稿服务器记录的 Python 3.12.13 和
[`requirements.txt`](../requirements.txt) 中的 NumPy、pandas 和 SciPy 版本。
h5py 和 Matplotlib 按服务器数据读取和绘图环境固定；首稿仿真记录没有登记这两项。
历史结果复跑用 [`requirements-archived.txt`](../requirements-archived.txt)。

```bash
python reproduce.py check-enhancement-env
bash run_cbc.sh --pseudobulk ../pseudobulk.npz \
  --imaging-cache ../cbc_cache/libd --bulk-cache ../cbc_cache/geo \
  --output ../cbc_results --workers 36 --allow-downloads
```

把论文3原 `paper3_marker_anchor_results_v1/prepared/pseudobulk.npz` 放在与代码目录同级，
或者把参数改成服务器已有文件路径。已有 LIBD 与 GSE125583 缓存可直接填路径。
省略 `--allow-downloads` 时不会联网；启用后只补齐缺失的校验匹配公开输入。
代码不安装依赖、不改动原数据和原结果。版本不匹配时，在读取数据或创建结果目录前停止。
可用 `PAPER3_PYTHON` 指定符合要求的现有 Python。Windows 下把 `bash run_cbc.sh` 换成
`python reproduce.py cbc-upgrade` 即可。

断点续跑时使用同一条命令、同一结果目录，加 `--resume`。每批通过校验后才计入进度。
默认运行新增的配对建库、锚点与覆盖验证，并重绘原 20 张图，不重跑旧模拟。

只回传 `cbc_results/final/`：六个新增结果 CSV、Fig. S12、原有图表的重绘数据和
两份核验记录。不要回传隐藏的工作批次目录。默认只导出 PDF；Word 插图可加
`--formats pdf png`。详细来源、判断标准和文件对应关系见 [CBC.md](CBC.md)。

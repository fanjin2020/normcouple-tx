# 论文3机制增强一键分析 1.0.0

本轮只增加机制验证、marker-only对照、MEAD实际列空间敏感性及GSE174367真实bulk应用入口。
不修改当前投稿文件，不预写实验结论，不扩展新队列，不打ZIP包。

## 一步执行

将单文件入口放在服务器已有数据目录后执行：

```bash
python paper3_mechanism_oneclick.py --setup --project-root /home/research/git/repository/paper1_public_evidence --output /home/research/git/repository/paper3_mechanism_results_v1 --workers 12
```

需要Python 3.10以上、R 4.4以上及可用的编译工具链。优先复用服务器
`/home/research/git/repository/paper3_revision_results_v1/R_library`；也可用`--r-library`指定。
Python依赖在独立环境中安装，不覆盖原环境。源码解包到有内容哈希的独立目录。
公共文件优先复用本地。内部LIBD/SEA-AD输入沿用原路径并校验原哈希，缺失时不偷换输入。
公共外部单核矩阵缺失时才下载GEO已知文件，不下载FASTQ、不绕过受控访问。

## 真正执行什么

1. 核验公式的有限变化、导数、纯log行位移及相同a/c性质。
2. 先在每个内部来源的第一个固定目标供体跑首批数值试验。只检验候选预测有限、配对不变性成立，不按FDR好坏筛选候选。
3. 通过数值门槛后执行原42个单元、80个mixture、200次重复。原随机种子和测试基因不变。内部保留原7种效应设置，外部保留null/up05；各加两种预定义marker受扰动情形。
4. 比较observed、fixed、TMM、marker matched、marker尺度对照、仅换分母对照和oracle。MEAD lambda=1的原生比例及参考yield转换比例，在原始非marker扰动条件下同时使用full-six及adaptive规则。不是联合传播MEAD比例协方差。
5. 保存每次拟合的比例矩阵、残差化暴露及基因拟合数组，避免以后只能看到汇总表。输出FDP/FDR估计、覆盖率、delivered coverage、功效、比例RMSE、可用性、秩、杠杆值和暴露R2。
6. 输出供体等权汇总、配对MC区间、供体差异和失败原因。MC区间条件于固定供体银行，不能写成供体总体置信区间。

首批试验不是性能筛选：不要求marker-only优于fixed/TMM，不要求100%下游可用。
`--pilot-only`会在首批后停下；同一完整命令再次运行可继续。参数/依赖/源码变动须换输出目录。
每批打印单元、重复起点和耗时。已完成批次有校验记录，原命令可断点续跑。
计算使用CPU，不会因为有P40而自动获得GPU加速。默认12进程，每进程BLAS单线程。

## 真实bulk的重要边界

已实际读取GEO公开的`GSE174367_bulkRNA_processed.rda.gz`：只有`normExpr.reg`和`targets`，
是18,849基因、90个样本的处理后表达，而不是原始计数；元数据为44 AD、46 Control，均为FC。
不能把它指数变换或四舍五入为“原始计数”。代码会下载该文件和单核元数据，检查SampleID、
供体、脑区、年龄/性别/疾病标签，输出明确匹配和排除记录；处理后表达关联标为探索性结果。

要完成所要求的observed/TMM/marker-only真实bulk比较，还需要可信的原始bulk基因计数。
脚本不会自动发起大体积FASTQ重建或登录受控库。原文数据入口为
https://www.synapse.org/Synapse:syn22079621 ，其可访问性和许可需由实际账号确认。
提供一个CSV（首列gene_id，Ensembl基因ID；其余列精确SampleID）以及
项目根目录的`paper3_bulk_raw_manifest.json`：

```json
{
  "accession": "GSE174367",
  "measurement": "raw_gene_counts",
  "counts_path": "/absolute/path/to/real_bulk_counts.csv",
  "sha256": "替换为该文件真实SHA256",
  "provenance": "填写实际来源、资源版本与计数生成流程"
}
```

这只是输入格式示例，未把占位符当成已确认数据。原始计数缺失时，整体状态为
`PARTIAL_REAL_BULK_COUNTS_OR_ACCESS_NEEDED`，模拟仍可完成，退出码2表示有待办而非模拟失败。
获得原始计数后用同一命令加`--only bulk`继续，不重跑模拟。
原始计数应用固定FC、AD对Control、年龄/性别/PMI/RIN/测序组完整病例；用同一个TMM
log2-CPM+1表达结局比较不同组成输入，避免混合改变结局尺度。BH基因家族排除所有参考marker，
各方法一致；真实数据不计算或宣称真实FDR。已有单核的恢复核比例仅作一致性对照，非组织组成真值。
真实疾病分析没有干预前分母，因此不伪造fixed分支。

## 输出与归档

- `final_results/`：带回本机的紧凑结果，CSV、图、状态、协议与机制说明，不含庞大参考缓存。
- `audit_tables/`：全部逐次诊断和配对差值，可能较大，留在服务器；供体及基因汇总已在`final_results/`。
- `units/`：断点和逐次中间量；`fit_index.csv`标明NPY第一维；`axes.json`标明基因和字段。
- `real_bulk/`：公开数据核对、供体配对及疾病关联结果，单独STATUS.json。
- `RUN_STATUS.json`：先读它，不能只看有无图或命令退出。
- 相邻的`paper3_mechanism_results_v1_prepared_inputs/`：默认按哈希复制准备输入，不放入最终结果目录。
  可用`--input-bundle`在另一台机器运行。空间不足可加`--no-input-export`，此时输入交付未完成。
- 各资源公开发布及权限审查仍需人工完成；本地导出不等于公开DOI。

只重绘：在解包的源码目录执行`python extension_figures.py /path/to/final_results`。
软件示例：`python paper3_mechanism_oneclick.py --setup --synthetic --output /tmp/paper3_extension_fixture --workers 1`。
人工合成示例永远标为`SOFTWARE_FIXTURE_ONLY`，不能当论文证据。

本轮尚未使用服务器真实输入完整跑完；本地软件测试及真实GEO结构核对记录单独提供。
机制推导是可审阅的草稿，不能用一次数值自检代替数学审阅或科学验证。

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22047484.svg)](https://doi.org/10.5281/zenodo.22047484)

v1.0.0归档版本：https://doi.org/10.5281/zenodo.22047485

# 血红蛋白与谵妄双向孟德尔随机化：公开代码

本仓库收录血红蛋白浓度与谵妄双向两样本孟德尔随机化研究的冻结分析代码、经核验的汇总结果、图表及复现说明。

## 公开状态

- 公开版本：`v1.0.0`
- GitHub仓库：https://github.com/RIan-5801/hemoglobin-delirium-bidirectional-MR
- Zenodo版本DOI：https://doi.org/10.5281/zenodo.22047485
- Zenodo Concept DOI：https://doi.org/10.5281/zenodo.22047484
- 许可证：MIT

## 已完成的核验

- 排除了3个原始GWAS文件、清洗数据、LD参考文件和PLINK二进制程序。
- 排除了论文写作文件、文献下载、旧失败版本及本地运行历史。
- 保留了各分析阶段最终采用的冻结脚本。
- 使用方向元数据修正后的最终结果矩阵。
- 独立核对了置信区间、主要P值、前向OR换算、异质性P值和SNP数量。
- 保留了Chen前向MR-PRESSO在冻结配置下未完成的真实状态。
- 已纳入完整的`renv.lock`，记录R 4.6.1及89个R包。
- 本机环境恢复成功，`R/validate_frozen_results.R`的9项检查均通过。

## 数据与复现边界

本仓库不重新分发GWAS汇总统计数据、LD参考数据或本地可执行程序。相关数据来源、标识符及预期文件名见`docs/data_sources.md`。

部分冻结脚本保留原始项目路径，用于维持分析审计来源，因此本仓库不宣称能够在缺少原始数据和本地资源的情况下，一键重跑全部GWAS处理流程。

详细方法、代码映射及核验结果见英文`README.md`和`docs/`目录。

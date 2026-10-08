# MatterGPT Carbon Workflow

使用 MatterGPT 实现训练碳结构模型、条件生成和 Oganov–Valle 指纹分析的初步研究资料库。

## 当前状态

已提供训练/验证数据、模型权重及配套配置、训练与生成启动脚本、两份五组分析脚本以及历史分析汇总表。尚未完成端到端运行验证，历史结果未重新计算。

仍待补充：原实验 SLICES 提交号、依赖环境、硬度标签单位及来源、数据划分方法、解码和筛选流程、十份分析输入及 C32 参考结构。

## 内容

| 目录 | 内容 |
|---|---|
| scripts/ | 训练、生成启动脚本及五组分析代码 |
| data/ | train_data_c.csv（186 条）和 val_data_c.csv（50 条） |
| model/ | 8-31.pt、8-31.ini、8-31_vocab.json |
| results/ | 已保存的五组内部分析汇总表，非本次复算结果 |
| docs/ | 运行说明 |

## 训练和生成

模型实现来自 [SLICES / MatterGPT](https://github.com/xiaohang007/SLICES)，本仓库不重复分发其完整实现。训练数据列依次为 SLICES、Vickers Hardness、crystal_system。

训练参数：run_name=8-31，batch_size=8，max_epochs=100，n_embd=128，n_layer=4，n_head=4，learning_rate=3e-4。

生成目标：40、50、60、70、80、90、100。gen_size=5000，batch_size=8。按当前上游实现，每个目标采样 5000 次，总计尝试 35000 条；有效条数及重复情况取决于生成结果。本地历史代码版本仍待核对，不将重复调用自动解释为不同随机种子实验。

详细步骤见 [运行说明](docs/usage.md)。

## 五组分析

输入文件编号两两合并：(0000,1000)、(2000,3000)、(4000,5000)、(6000,7000)、(8000,9000)。默认指纹参数为 Rmax=15 Å、dR=0.05 Å、sigma=0.075 Å。

距离 D=0.5×(1−cosine)，展示相似度 S=1−D。余弦距离为零表示非零指纹方向相同，不证明两个晶体严格等价。结构相似性、重复结构判定和对数据库的新颖性检查应分别说明。

## 来源与许可

上游实现及其许可证见 [SLICES](https://github.com/xiaohang007/SLICES)（当前根许可证为 LGPL 2.1，实际使用版本待固定）。本研究使用原实现并调整运行参数。部分分析代码使用 AI 辅助编写。

本初步版本尚未为新增代码、数据和模型指定再使用许可证；公开可见不等于授予任意使用许可。第三方材料仍受其原有条款约束。后续将补齐来源和适用 LICENSE。

参考文献：

- [MatterGPT: A Generative Transformer for Multi-Property Inverse Design of Solid-State Materials](https://arxiv.org/abs/2408.07608).
- Oganov, A. R. and Valle, M. How to quantify energy landscapes of solids. J. Chem. Phys. 130, 104504 (2009). https://doi.org/10.1063/1.3079326

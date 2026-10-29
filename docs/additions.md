# 补充资料

新增十份分析输入、两幅内部比较结果图（PNG/PDF）、逐结构 C32 比较结果和 C32 绘图脚本。均来自已有实验文件，本次未重新计算。绘图脚本的运行入口和依赖输入仍需核对。

## 文件名映射

- `resultsC0000 .csv` → `data/resultsC0000.csv`
- `resultsC1000 .csv` → `data/resultsC1000.csv`
- `resultsC2000.csv` → `data/resultsC2000.csv`
- `resultsC3000.csv` → `data/resultsC3000.csv`
- `resultsC4000.csv` → `data/resultsC4000.csv`
- `resultsC5000.csv` → `data/resultsC5000.csv`
- `resultsC6000.csv` → `data/resultsC6000.csv`
- `resultsC7000.csv` → `data/resultsC7000.csv`
- `resultsC8000.csv` → `data/resultsC8000.csv`
- `resultsC9000.csv` → `data/resultsC9000.csv`

## 仍需提供

1. C32-C2-m-no-opt.vasp 参考结构及来源。
2. 原实验 SLICES 提交号（如可获取）；已确认使用未修改的上游代码，无需重复提供源码。
3. 原实验环境版本可后续补充；安装先参考 SLICES 官方说明，不将其当作已核实的原实验环境。
4. Vickers Hardness 单位、标签来源、是否归一化；训练/验证划分方法及增强顺序。
5. 从生成 SLICES 到 resultsC 文件的解码、优化、筛选命令；各次生成与文件编号的对应关系、随机种子处理。
6. 许可证选择、版权署名和数据/模型的可再分发条款。

重复生成不自动证明使用不同随机种子。下载文件是否与服务器最终版本完全一致，仍需作者确认。

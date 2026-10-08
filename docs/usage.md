# 运行说明

这是初步版本，尚未验证安装环境或完整流程。请优先使用原实验中运行成功的 SLICES 版本和环境。

## 训练

在本仓库根目录的 Linux / Bash 环境执行：

```bash
bash scripts/train_carbon.sh /path/to/SLICES
```

替换路径为 SLICES checkout。脚本使用本仓库的训练和验证 CSV，进入上游 MatterGPT/1_train_generate 后调用 train.py。若上游 model/ 中已经存在 8-31 同名文件会停止，避免覆盖。

## 使用已有模型生成

将本仓库 model/ 中三个文件复制到上游 MatterGPT/1_train_generate/model/，先确认不存在需保留的同名文件。

```bash
bash scripts/generate_carbon.sh /path/to/SLICES C10000.csv
```

输出写入本仓库 results/generated/。重复执行时更换输出文件名。脚本未额外设置随机种子，保留用户提供的调用方式。

模型配置：vocab_size=168，block_size=628，num_props=1，n_layer=4，n_head=4，n_embd=128，sym_dim=7；JSON 词表包含 168 个 token。仅核对配置与词表，未加载权重进行推理。

## 五组分析

准备 resultsC0000.csv、resultsC1000.csv，依次至 resultsC9000.csv。脚本读取 CSV 中的 POSCAR 文本，因此不是直接读取生成的 SLICES 字符串。解码和筛选步骤尚待补充。

```bash
python scripts/five_group_oganov_similarity.py /path/to/analysis_inputs
python scripts/compare_five_groups_to_C32_oganov.py /path/to/analysis_inputs --reference /path/to/C32-C2-m-no-opt.vasp
```

依赖包括 numpy、pandas、matplotlib、scipy 和 pymatgen，具体版本待核对。两份程序会在输入目录写结果，建议使用单独工作目录，避免覆盖历史结果。分析输入和参考结构当前未包含在仓库中。

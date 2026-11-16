# ProtPhenoAge

**中文** | [English](README.md)

> 基于血浆蛋白组学的蛋白表型年龄（ProtPhenoAge）计算管道
> vbeta.20261116（2026-11-16）

ProtPhenoAge 使用 Olink 血浆蛋白数据预测表型年龄（PhenoAge），对应论文发表于 *Journal of Advanced Research*（2026）。本管道基于 Python，提供三个功能：

1. **UKB 模型复现**：使用 UK Biobank 训练好的蛋白表型年龄模型，或在 UKB 数据上复现训练
2. **独立数据集蛋白训练**：在自己队列的蛋白数据上训练模型（BorutaShap 特征筛选 + Optuna 调参 + XGBoost）
3. **独立数据集蛋白验证**：加载训练好的模型，对新的蛋白矩阵输出 ProtPhenoAge 预测

`model/` 中提供的 UKB 模型可直接加载用于预测，无需训练数据。

## 蛋白矩阵预处理

1. **长表转宽表**：按 `protein_id` 展开 `result`
2. **去重与重复测量**：按 `eid` 去重；对重复测量者，逐蛋白取首次非缺失值
3. **蛋白 QC**：剔除质控不合格蛋白（UKB 中为 `GLIPR1`）；剔除缺失比例 > 20% 的蛋白（保留非缺失比例 > 80%）
4. **样本 QC**：仅保留白人（可选）；排除基线癌症（可选）
5. **缺失值填补**：KNN（k = 10，scale = TRUE，加权平均）
6. **标准化**：对每个蛋白做基于秩的逆正态变换（RINT）
7. **表型年龄计算**：若在自己队列中训练，需先在训练集计算 PhenoAge，参考 [BioAge](https://github.com/dayoonkwon/BioAge/)

> 由于 Olink 检测的蛋白量为相对定量，每个数据集最好在内部重新训练；条件不允许时可直接使用 UKB 训练好的模型参数/模型文件。

## 仓库结构

| 路径 | 说明 |
|---|---|
| `ProtPhenoAge4github.qmd` | 管道规范文档（Quarto） |
| `script/ProtPhenoAge4UKB.py` | 推理脚本：加载 UKB 模型，对任意蛋白矩阵预测 |
| `model/` | 训练好的 UKB 模型文件（`ProtPhenoAge_UKB.json`）及元数据 |
| `data/` | 训练/测试数据（UKB 原始数据受使用协议限制，不随仓库分发） |

## 环境依赖

```bash
pip install pandas numpy scipy scikit-learn xgboost
# 独立训练默认流程另需:
pip install optuna BorutaShap
```

## 快速开始

### 1. 使用 UKB 模型预测新矩阵

```bash
python script/ProtPhenoAge4UKB.py \
  --input data/test_prot_1000.csv \
  --model model/ProtPhenoAge_UKB.json \
  --output output/test_prot_1000_protphenoage.csv
```

- `--input` 必填；`--model` 默认 `model/ProtPhenoAge_UKB.json`；`--output` 默认 `<输入文件名>_protphenoage.csv`（与输入同目录）
- 输入含 ID 列时输出该 ID，否则输出 `row`（1 起始行号）
- `--keep-all`：输出完整输入矩阵并追加 `protphenoage` 列
- 模型含 185 个蛋白特征，按列名对齐；**建议输入完整的 2,920 蛋白矩阵**，只包含这 185 个模型蛋白也可以运行；多余列自动忽略，缺少的模型蛋白按缺失值处理并给出警告

也可在 Python 中直接加载模型文件预测：

```python
import xgboost as xgb

booster = xgb.Booster(model_file='model/ProtPhenoAge_UKB.json')
features = booster.feature_names
pred = booster.predict(xgb.DMatrix(df.reindex(columns=features),
                                   feature_names=features))
```

> 输入蛋白值须与训练数据处于同一预处理尺度（每个队列内部完成缺失值插补 + RINT）。Olink 为相对定量，跨队列直接迁移可能存在系统偏差，条件允许时建议在自身队列内部重新训练。

### 2. 独立队列训练新模型

- 正在开发

## 输入数据格式

输入文件统一为 **CSV**：英文逗号 `,` 分隔、UTF-8 编码、第一行为列名（表头）。矩阵遵循「一行为一个样本，一列为一个变量」的约定。

### 蛋白矩阵（预测 `--input` / 训练 `--x`）

- **行（row）**：每行一个样本（参与者）；行序即样本顺序
- **行名（首列）**：行索引列，表头通常为空（`""`，pandas 读入为 `Unnamed: 0`）；也可以是样本 ID 列（列名如 `eid`、`f.eid`、`sample_id`，脚本自动识别并从特征中剔除）；该列可省略（如 `data/test_prot_1000.csv`）
- **列（column）**：除首列外，每个蛋白一列，列名为 Olink 蛋白符号（如 `A1BG`、`GDF15`、`IL6`）。模型按**列名**对齐蛋白，与列顺序无关；**建议提供完整的 2,920 蛋白矩阵**，仅包含 185 个模型蛋白也可运行；多余列自动忽略，缺少的模型蛋白按缺失值处理并给出警告
- **单元格**：蛋白丰度数值（Olink NPX），须与训练数据同尺度，即已完成质控、缺失值插补与 RINT 预处理；缺失可留空或写 `NA`
- **分隔符**：英文逗号 `,`

示例（取自 `data/x_train.csv`，实际共 2,920 个蛋白列）：

```csv
,A1BG,AAMDC,AARSD1,ABCA2,...,ZPR1
31180,0.0678,-0.95205,-0.91035,-0.1899,...
```

### 标签文件（训练/评估 `--y`）

- **行（row）**：与蛋白矩阵行一一对应。无公共 ID 时按行顺序对齐（两者行数必须一致）；有公共 ID 时按 ID 匹配
- **列（column）**：
  - 首列：行索引（表头为空，可省略）
  - `f.eid`：样本/参与者编号
  - `phenoage`：表型年龄标签（数值，单位：年），需事先按 [BioAge](https://github.com/dayoonkwon/BioAge/) 方法计算
- **分隔符**：英文逗号 `,`

示例（`data/y_train.csv`）：

```csv
,f.eid,phenoage
31180,5184712,57.7314167660788
```

### 预测输出

输出为 `<输入文件名>_protphenoage.csv`（CSV，逗号分隔，含表头）：

- 输入含 ID 列时首列为该 ID，否则为 `row`（1 起始行号）
- `protphenoage`：预测的蛋白表型年龄（年）

示例：

```csv
row,protphenoage
1,41.804337
2,56.852463
```

## 引用

如果本管道对你的研究有帮助，请引用：

> Wang Y, Sun Y, Yang F, et al. **ProtPhenoAge: Integrating plasma proteomics to predict Aging-Related disease Risks**. *Journal of Advanced Research*, 2026. doi:10.1016/j.jare.2026.05.022

[论文链接（ScienceDirect）](https://www.sciencedirect.com/science/article/pii/S2090123226004170)

```bibtex
@article{wang2026protophenoage,
  title   = {ProtPhenoAge: Integrating plasma proteomics to predict Aging-Related disease Risks},
  author  = {Wang, Yifan and Sun, Yue and Yang, Fan and Li, Musu and Qi, Tianchen and Lu, Zixuan and Wang, Qian and Bu, Qingyin and Sun, Longhu and Wo, Hongmei and Zhao, Yang and Yi, Honggang and Dai, Juncheng},
  journal = {Journal of Advanced Research},
  year    = {2026},
  doi     = {10.1016/j.jare.2026.05.022}
}
```

## 致谢

- PhenoAge 计算参考：[BioAge](https://github.com/dayoonkwon/BioAge/)
- 特征筛选：[BorutaShap](https://github.com/Ekeany/Boruta-Shap)
- 调参与模型：Optuna、XGBoost

## 许可协议

本项工作（论文与本仓库内容）采用 [CC BY-NC-ND 4.0](https://creativecommons.org/licenses/by-nc-nd/4.0/) 许可协议（© 2026 The Authors）。分享转载需注明出处；不得用于商业用途，不得修改后分发。如需商业使用或其他授权，请联系作者。

## 联系方式

王一帆 — wangyfly123@163.com

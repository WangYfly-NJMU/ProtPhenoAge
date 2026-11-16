#!/usr/bin/env python3
"""使用训练好的 ProtPhenoAge UKB 模型，对新蛋白矩阵进行推断。

输入: 蛋白矩阵 CSV（可为 x_test 格式、test_prot_1000 格式或任意宽表）
输出: 预测结果 CSV（默认紧凑格式; --keep-all 输出完整矩阵 + 预测列）

运行示例:
    python script/ProtPhenoAge4UKB.py --input data/test_prot_1000.csv
    python script/ProtPhenoAge4UKB.py --input F:/UKB/data/new.csv --output F:/UKB/out.csv
"""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import xgboost as xgb

PROTEINS = [
    'ANGPTL4', 'RET', 'AFP', 'PPY', 'XCL1', 'CDHR2', 'S100G', 'EFEMP1', 'MMP9', 'ACTA2',
    'AGER', 'COL4A1', 'GP2', 'ADA', 'IGF1R', 'PCDH9', 'CD1C', 'FOLR2', 'LTBP2', 'SPINK6',
    'NTF3', 'PROK1', 'OSM', 'MEPE', 'CA14', 'CD93', 'LMOD1', 'CHRDL1', 'ADAMTS15', 'CILP',
    'MOG', 'ITGB5', 'TNFRSF11B', 'MMP1', 'CTSV', 'CCDC80', 'CXCL9', 'COL15A1', 'CAPS', 'L1CAM',
    'XG', 'ENG', 'KLK14', 'GAS6', 'SCG3', 'CELSR2', 'ANGPT2', 'CXCL12', 'CD8A', 'LRG1', 'ROBO1',
    'SLITRK6', 'ACRV1', 'LTA', 'GALNT2', 'TNXB', 'CRTAC1', 'CXCL8', 'NADK', 'RLN2', 'CHI3L1',
    'CLEC4A', 'GDF15', 'IL6', 'ASAH2', 'PRELP', 'PTPRR', 'ITGAV', 'PSG1', 'ADM', 'KLK7', 'TCOF1',
    'FBLN2', 'CXCL14', 'PKD1', 'FST', 'CTSO', 'CRYBB2', 'KLK3', 'INHBB', 'IL12A_IL12B', 'TTR',
    'S100A12', 'CLEC3B', 'PODXL2', 'GDF2', 'DPP4', 'SCARF2', 'TNFRSF10C', 'CSF3', 'TSPAN1',
    'WIF1', 'MATN3', 'NEFL', 'PRND', 'DSC2', 'SERPINA3', 'LAMP2', 'CLEC5A', 'LILRA5', 'SUSD5',
    'CCL25', 'CD248', 'ADAMTS13', 'LILRB4', 'BCAN', 'ENPP2', 'KLK4', 'CLEC4D', 'CHIT1', 'ERBB4',
    'LECT2', 'SPINT1', 'HSD11B1', 'FASLG', 'WFDC2', 'GHRL', 'CLEC14A', 'APCS', 'GFAP', 'HBQ1',
    'CDON', 'CTHRC1', 'FSHB', 'KIT', 'NTproBNP', 'FABP1', 'NOTCH3', 'NCAN', 'CNTN5', 'AGRN',
    'CCN5', 'ZNRD2', 'CEACAM16', 'PLAT', 'C9', 'FGA', 'GFRAL', 'CRH', 'NTF4', 'CDCP1', 'WNT9A',
    'KIRREL2', 'DSG4', 'ANXA1', 'FCN1', 'FUT3_FUT5', 'PLA2G7', 'CA4', 'CX3CL1', 'ZNF830',
    'LRIG1', 'ELN', 'C1QL2', 'MLN', 'FGF21', 'AGRP', 'CDH3', 'GIP', 'IGDCC4', 'PTH1R', 'HYAL1',
    'NCF2', 'GAST', 'CR2', 'BOC', 'DPEP2', 'MET', 'EDA2R', 'AHSP', 'PAEP', 'SERPINC1', 'CXCL17',
    'IL17D', 'TCN2', 'MMP12', 'SHBG', 'CSF1', 'PTPRN2', 'KEL', 'TGFA', 'HAVCR1', 'GZMA', 'KLK8',
    'CGA',
]

ID_CANDIDATES = {'eid', 'f.eid', 'id', 'iid', 'fid', 'sample', 'sample_id', 'subject_id'}


def parse_args():
    p = argparse.ArgumentParser(description='ProtPhenoAge UKB 模型推断：输入蛋白矩阵，输出表型年龄预测')
    p.add_argument('--input', type=str, required=True, help='输入蛋白矩阵 CSV')
    p.add_argument('--output', type=str, default=None,
                   help='输出 CSV 路径（默认: <输入文件名>_protphenoage.csv，与输入同目录）')
    p.add_argument('--model', type=str, default=None,
                   help='模型 JSON 路径（默认: <项目>/model/ProtPhenoAge_UKB.json）')
    p.add_argument('--keep-all', action='store_true',
                   help='输出完整输入矩阵并在末尾追加 protphenoage 列（默认只输出 ID + 预测）')
    return p.parse_args()


def detect_id_column(columns):
    for c in columns:
        if str(c).strip().lower() in ID_CANDIDATES:
            return c
    return None


def main():
    args = parse_args()
    input_path = Path(args.input)
    model_path = (Path(args.model) if args.model
                  else Path(__file__).resolve().parent.parent / 'model' / 'ProtPhenoAge_UKB.json')
    output_path = (Path(args.output) if args.output
                   else input_path.with_name(f'{input_path.stem}_protphenoage.csv'))

    if not input_path.is_file():
        raise SystemExit(f'[错误] 输入文件不存在: {input_path}')
    if not model_path.is_file():
        raise SystemExit(f'[错误] 模型文件不存在: {model_path}')

    booster = xgb.Booster(model_file=str(model_path))
    feature_names = booster.feature_names or PROTEINS
    print(f'模型: {model_path}')
    meta_path = model_path.with_name(model_path.stem + '.metadata.json')
    if meta_path.is_file():
        meta = json.loads(meta_path.read_text(encoding='utf-8'))
        print(f"  训练时间: {meta.get('created_at', '?')} | "
              f"xgboost: {meta.get('environment', {}).get('xgboost', '?')} | "
              f"树数: {meta.get('n_boost_round_final', '?')}")

    header = pd.read_csv(input_path, nrows=0)
    columns = [c for c in header.columns if not str(c).startswith('Unnamed')]
    available = [p for p in feature_names if p in columns]
    missing = [p for p in feature_names if p not in columns]
    if not available:
        raise SystemExit(f'[错误] 输入中未找到任何模型所需蛋白（应为 {len(feature_names)} 个，'
                         f'如 {feature_names[:5]} ...），请检查文件格式')
    if missing:
        print(f'[警告] 输入缺少 {len(missing)}/{len(feature_names)} 个模型蛋白，'
              f'将按缺失值处理（XGBoost 走默认分支）: '
              f'{missing[:20]}{" ..." if len(missing) > 20 else ""}')

    id_col = detect_id_column(columns)
    extra_cols = [c for c in ('phenoage',) if c in columns]
    usecols = columns if args.keep_all else list(dict.fromkeys(
        ([id_col] if id_col else []) + available + extra_cols))
    df = pd.read_csv(input_path, usecols=usecols)

    X = df.reindex(columns=feature_names).apply(pd.to_numeric, errors='coerce')
    n_na = int(X.isna().sum().sum())
    if n_na:
        print(f'[警告] 特征矩阵含 {n_na} 个缺失值（缺失蛋白或非数值），XGBoost 将按缺失处理')

    pred = booster.predict(xgb.DMatrix(X, feature_names=feature_names))

    if args.keep_all:
        out = df.copy()
        out['protphenoage'] = pred
    else:
        ids = df[id_col].to_numpy() if id_col else np.arange(1, len(df) + 1)
        out = pd.DataFrame({id_col if id_col else 'row': ids})
        if extra_cols:
            out['phenoage'] = df['phenoage'].to_numpy()
        out['protphenoage'] = pred
    if extra_cols:
        out['protPA_accel'] = out['phenoage'] - out['protphenoage']

    output_path.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(output_path, index=False)
    print(f'预测完成: {len(out)} 行 x {len(feature_names)} 蛋白 -> {output_path}')
    print('protphenoage: min=%.2f, mean=%.2f, sd=%.2f, max=%.2f'
          % (float(np.min(pred)), float(np.mean(pred)), float(np.std(pred)), float(np.max(pred))))


if __name__ == '__main__':
    main()

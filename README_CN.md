# 脑小血管病EEG-fNIRS分析

[English](README.md) | **中文**

基于MATLAB的EEG-fNIRS同步数据处理、特征提取及神经血管耦合分析完整流水线，用于脑小血管病（CSVD）亚型的脑功能评估。

## 项目概述

本项目构建了一套完整的多模态神经影像分析框架，系统比较三种脑小血管病亚型——血管病危险因素相关SVD、CADASIL（伴皮层下梗死和白质脑病的常染色体显性遗传性脑动脉病）和CAA（脑淀粉样血管病）——与健康对照（HC）在多种认知运动任务中的神经电生理与血流动力学特征差异，并从神经血管耦合角度揭示各亚型特异性的脑功能损害模式。

### 任务范式

| 任务 | 范式描述 | EEG特征 | fNIRS感兴趣区 |
|------|---------|---------|-------------|
| **听觉Oddball** | 标准/偏差音调检测 | P300事件相关电位、Delta/Theta频段功率 | 前额叶、顶叶、颞叶 |
| **手指运动** | 左/右手运动执行 | Mu/Beta事件相关去同步化（ERD）、侧化指数 | 双侧中央区 |
| **语义分类** | 视觉词汇分类判断 | Alpha/Beta相对功率、ERD | 额区、颞区 |

## 仓库结构

```
├── eeg/                                # EEG数据与处理脚本
│   ├── audi/                           # 听觉Oddball任务
│   │   ├── audi_EEG_ERP_Processing/
│   │   │   └── scripts/               # STEP0–STEP14 MATLAB脚本
│   │   └── sub-*/                      # 各被试原始数据（.edf）
│   ├── motor/                          # 手指运动任务
│   │   ├── motor_EEG_Processing/
│   │   │   └── scripts/               # STEP0–STEP7 MATLAB脚本
│   │   └── sub-*/
│   └── cate/                           # 语义分类任务
│       ├── cate_EEG_Processing/
│       │   └── scripts/               # STEP0–STEP7 MATLAB脚本
│       └── sub-*/
│
├── fnirs/                              # fNIRS数据与处理脚本
│   ├── audi/
│   │   ├── audi_fnirs_processing/     # STEP1–STEP6 MATLAB脚本
│   │   └── sub-*/                      # 各被试数据（.snirf）
│   ├── motor/
│   │   ├── motor_fnirs_processing/
│   │   └── sub-*/
│   ├── cate/
│   │   ├── cate_fnirs_processing/
│   │   └── sub-*/
│   └── nirs_to_snirf.m               # 数据格式转换工具
│
└── eeg_fnirs_processing/              # EEG-fNIRS耦合分析
    ├── audi_eeg_fnirs_processing/     # STEP1–STEP4（完整耦合分析流水线）
    ├── motor_eeg_fnirs_processing/    # STEP1–STEP2（Block级相关 + EEG-informed GLM）
    └── cate_eeg_fnirs_processing/     # STEP1–STEP2（Block级相关 + EEG-informed GLM）
```

## 处理流水线

### 1. EEG处理

基于EEGLAB/ERPLAB框架的顺序化MATLAB脚本：

| 步骤 | 脚本 | 功能描述 |
|------|------|---------|
| 0 | `STEP0_marker.m` | 从trigger通道提取事件标记，校正LCD显示器延迟 |
| 1 | `STEP1_Import_Raw_EEG_Shift_DS_Reref_Hpfilt.m` | 导入原始EEG（1000 Hz）→ 降采样至256 Hz → 重参考至双侧乳突均值（A1/A2）→ 0.1 Hz高通滤波 |
| 2 | `STEP2_ICA_Prep.m` | ICA预处理准备（通道插值、坏段剔除） |
| 3 | `STEP3_Run_ICA.m` | 扩展Infomax ICA分解（60导）结合ICLabel自动成分分类 |
| 4 | `STEP4_Remove_ICA_Components.m` | 人工审核并移除眼电、肌电等伪迹成分 |
| 5 | `STEP5_Elist_Bin_Epoch.m` | 创建事件列表、按bin分类、截取分段（−200至2000 ms） |
| 6 | `STEP6_Artifact_Rejection.m` | 伪迹检测：简单电压阈值法 + 移动窗口峰-峰值法 |
| 7+ | 任务特异性分析 | **Oddball**：`STEP7_Average_ERPs.m` → ERP叠加平均与成分测量；**运动**：`STEP7_ERD.m` → 时频ERD分析（小波）；**语义**：`STEP7_power.m` → 频谱功率提取 |

### 2. fNIRS处理

基于Homer3框架的顺序化MATLAB脚本：

| 步骤 | 脚本 | 功能描述 |
|------|------|---------|
| 1 | `STEP1_preprocessing.m` | 加载SNIRF → 光强转光密度（OD）→ 通道质量评估（SNR ≥ 2, SCI ≥ 0.6）→ 运动伪迹校正（TDDR + 小波去噪, IQR = 1.5）→ 带通滤波（0.01–0.2 Hz, Butterworth）→ 修正Beer-Lambert定律（PPF = 6.0）→ 输出ΔHbO/ΔHbR |
| 2 | `STEP2_qc.m` | 质量控制可视化：原始光强 → OD → 校正后OD → HbO/HbR时间序列；生成 `QC_summary.csv` |
| 3 | `STEP3_BlockAverage.m` | Block平均分析：识别任务block → 基线校正（−5至0 s）→ 提取block平均血流动力学响应 |
| 4 | `STEP4_GLM.m` | 个体水平GLM：刺激序列与gamma HRF（k=6, θ=1, 持续30 s）卷积 → 设计矩阵含DCT漂移回归量（截止周期256 s）→ OLS/AR-IRLS拟合 → 输出β系数和t统计量 |
| 5 | `STEP5_group_block.m` | 组水平Block分析：单样本t检验 + FDR多重比较校正（α = 0.05） |
| 6 | `STEP6_group_GLM.m` | 组水平GLM分析 |

### 3. EEG-fNIRS耦合分析

多层次神经血管耦合分析：

| 步骤 | 脚本 | 分析方法 |
|------|------|---------|
| 1 | `STEP1_correlation.m` | **Pearson相关分析**：EEG特征（P300波幅、频段功率、ERD）与fNIRS血氧信号（ΔHbO/ΔHbR）的相关性。Oddball任务为试次级分析，运动/语义任务为Block级分析 |
| 2 | `STEP2_EEG_informed_fnirs.m` | **EEG-informed GLM**：以EEG特征（z-score标准化）调制HRF作为参数化回归量纳入fNIRS GLM，量化神经电活动对血流动力学变异的解释贡献 |
| 3 | `STEP3_timelag_correlation_and_coherence.m` | **互相关与相干分析**：EEG频段功率的Hilbert包络 → 插值至fNIRS采样率 → 互相关求最优时滞 → 幅度平方相干性分析（仅Oddball任务） |
| 4 | `STEP4_dynamic_coupling.m` | **动态耦合建模**：ARX系统辨识（na=nb=2, 延迟2.5 s）→ 状态空间建模 → 自适应Kalman滤波（遗忘因子0.98）估计时变耦合强度（仅Oddball任务） |

## 核心参数

### EEG参数

| 参数 | 设定值 |
|------|-------|
| 采样率（原始 → 处理后） | 1000 Hz → 256 Hz |
| 参考电极 | 双侧乳突均值（A1/A2） |
| 高通滤波 | 0.1 Hz |
| 低通滤波（ERP分析） | 20 Hz |
| 分段时间窗 | −200至2000 ms |
| 基线校正 | −200至0 ms |
| ICA算法 | 扩展Infomax |
| P300时间窗 | 300–500 ms（Cz/Pz电极） |

### fNIRS参数

| 参数 | 设定值 |
|------|-------|
| 通道质量标准 | SNR ≥ 2, SCI ≥ 0.6 |
| 运动伪迹校正 | TDDR + 小波去噪（IQR = 1.5） |
| 带通滤波 | 0.01–0.2 Hz（Butterworth） |
| Beer-Lambert PPF | 6.0 |
| HRF模型 | Gamma函数（k = 6, θ = 1, 持续30 s） |
| GLM漂移去除 | DCT高通（截止周期256 s） |
| Block基线窗 | −5至0 s |

## 数据格式

| 模态 | 原始格式 | 处理后格式 |
|------|---------|----------|
| EEG | `.edf`（欧洲数据格式） | `.set`（EEGLAB格式） |
| fNIRS | `.snirf`（标准近红外光谱格式） | `.mat`（MATLAB格式） |
| 行为数据 | `.txt`（E-Prime导出） | `.csv` |
| 统计结果 | — | `.csv`、`.mat`、`.xlsx` |

## 环境依赖

- **MATLAB** R2020b 或更高版本
- **EEGLAB**（v2021+），需安装插件：
  - ERPLAB（事件相关电位分析）
  - ICLabel（ICA成分自动分类）
  - clean_rawdata（数据清洗）
- **Homer3**（fNIRS处理工具箱）
- **Signal Processing Toolbox**（MATLAB信号处理工具箱）
- **Statistics and Machine Learning Toolbox**（MATLAB统计工具箱）

## 使用说明

各处理模块采用顺序化的STEP脚本设计，脚本自动识别 `sub-*` 格式的被试目录。

```matlab
% 示例：运行听觉Oddball任务的EEG预处理
cd eeg/audi/audi_EEG_ERP_Processing/scripts/
STEP0_marker
STEP1_Import_Raw_EEG_Shift_DS_Reref_Hpfilt
STEP2_ICA_Prep
STEP3_Run_ICA
STEP4_Remove_ICA_Components   % 需要人工审核ICA成分
STEP5_Elist_Bin_Epoch
STEP6_Artifact_Rejection
STEP7_Average_ERPs

% 示例：运行fNIRS预处理
cd fnirs/audi/audi_fnirs_processing/
STEP1_preprocessing
STEP2_qc
STEP3_BlockAvarage
STEP4_GLM

% 示例：运行EEG-fNIRS耦合分析
cd eeg_fnirs_processing/audi_eeg_fnirs_processing/
STEP1_correlation
STEP2_EEG_informed_fnirs
STEP3_timelag_correlation_and_coherence
STEP4_dynamic_coupling
```

> **注意：** 各步骤须按顺序执行。EEG处理中的 `STEP4_Remove_ICA_Components` 需要人工视觉审核ICA成分。

## 分析流程图

```
原始EEG数据（.edf）                  原始fNIRS数据（.snirf）
       │                                    │
       ▼                                    ▼
┌──────────────────┐             ┌───────────────────────┐
│   EEG预处理       │             │   fNIRS预处理          │
│   STEP 0–6       │             │   STEP 1–4            │
│   ·降采样/重参考   │             │   ·通道质量控制         │
│   ·ICA伪迹去除    │             │   ·运动伪迹校正         │
│   ·坏段标记       │             │   ·带通滤波            │
│                  │             │   ·Beer-Lambert转换    │
└────────┬─────────┘             └───────────┬───────────┘
         │                                   │
         ▼                                   ▼
┌──────────────────┐             ┌───────────────────────┐
│   特征提取        │             │   特征提取              │
│   ·P300（Oddball）│             │   ·Block平均ΔHbO/HbR  │
│   ·ERD（运动）    │             │   ·GLM β系数           │
│   ·频谱功率（语义）│             │   ·达峰时间             │
└────────┬─────────┘             └───────────┬───────────┘
         │                                   │
         └───────────┬───────────────────────┘
                     ▼
     ┌────────────────────────────────────┐
     │   EEG-fNIRS神经血管耦合分析         │
     │   ·Block/试次级Pearson相关          │
     │   ·EEG-informed GLM               │
     │   ·互相关与频域相干性               │
     │   ·ARX/Kalman动态耦合建模          │
     └──────────────┬─────────────────────┘
                    ▼
     ┌────────────────────────────────────┐
     │   组水平统计分析                    │
     │   ·Welch t检验（vs HC）            │
     │   ·Benjamini-Hochberg FDR校正      │
     └────────────────────────────────────┘
```


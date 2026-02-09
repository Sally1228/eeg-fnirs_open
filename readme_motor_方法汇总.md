# Motor Task 方法学：EEG-fNIRS 神经血管偶联分析

## 1. 实验设计

### 1.1 任务范式

本研究采用单手握拳运动任务(motor task)来研究感觉运动皮层的神经血管偶联功能。任务包括两种条件：

- **左手握拳条件 (Left hand)**：持续握拳30秒，随后休息17秒 (Marker = 1)
- **右手握拳条件 (Right hand)**：持续握拳30秒，随后休息17秒 (Marker = 2)

每个被试完成8个block，其中4个block为左手握拳，4个block为右手握拳，左右手条件随机交替呈现。每个block开始时有event marker标记，单次握拳动作无额外标记。

**任务设计理由**：
- 30秒的握拳时长足够引发稳定的皮层激活和血流动力学响应
- 17秒的休息期允许血氧信号返回基线，避免相邻block间的信号重叠
- 左右手分别激活对侧运动皮层，便于对侧-同侧对比分析

### 1.2 数据采集

**EEG数据采集**：
- 采样率：1000 Hz (原始采样)
- 电极系统：64导EEG系统，按照10-10国际系统放置
- 参考电极：双耳乳突电极(A1/A2)
- 目标ROI：左中央区(C3, C1, C5, CP3, FC3)，右中央区(C4, C2, C6, FC4, CP4)

**fNIRS数据采集**：
- 采样率：11 Hz
- 波长：3波长系统(用于改善血红蛋白浓度估计)
- 光极布局：覆盖双侧感觉运动皮层
- 目标ROI：左中央区(Source-Detector pairs: 14-10, 10-10, 10-14, 14-14, 14-18, 18-18, 18-14)；右中央区(Source-Detector pairs: 11-11, 15-11, 15-15, 11-15, 19-15, 19-19, 15-19)

---

## 2. EEG数据处理流程

### 2.1 Event Marker提取 (STEP0)

**目的**：从trigger通道中提取block起始的event marker。

**方法**：
- 使用通道38、39、40作为trigger通道
- 对每个通道进行二值化处理(k-means聚类，重复5次)
- 将trigger信号编码为事件码(event code = weights * binary_bits，weights = [1, 2, 4])
- 检测event code的上升沿(从0变为非0)作为block起始时刻
- 将event时刻写入EEG.event结构

**参数**：
- `zero_tol = 1e-6`：零值容忍度
- `min_sep_factor = 5`：最小信号分离系数，用于判断高低电平差异

**输出**：
- EEG.set文件，包含标准化的event结构
- marker_events.csv和marker_counts.csv(可选)，用于质量检查

### 2.2 预处理 (STEP1)

**目的**：降采样、重参考、高通滤波和陷波滤波。

**降采样 (Downsampling)**：
- 从1000 Hz降采样到256 Hz
- 自动应用抗混叠低通滤波器(Nyquist频率 = 128 Hz)
- **理由**：减少计算负担，256 Hz足够捕捉运动相关的神经振荡(μ/β频段 < 30 Hz)

**重参考 (Re-referencing)**：
- 参考方式：平均耳电极参考(A1和A2的平均)
- **理由**：耳乳突电极远离运动皮层，可以减少局部肌电干扰，并提供稳定的参考

**高通滤波 (High-pass filtering)**：
- 滤波器类型：Butterworth，2阶，非因果(双向)
- 截止频率：0.1 Hz (半幅度衰减点)
- 滚降：12 dB/octave
- **理由**：去除直流偏移和极低频漂移，保留运动相关频段(μ/β频段 ≥ 8 Hz)

**陷波滤波 (Notch filtering)**：
- 频率：50 Hz ± 1 Hz (49-51 Hz)
- **理由**：去除工频干扰(工频干扰影响比较大，干扰后续质量检查，先滤掉50hz)

**输出**：
- `*_eeg_ds_reref_hpfilt_notch.set`

### 2.3 ICA准备 (STEP2)

**目的**：去除特别嘈杂的连续数据片段，为ICA提供干净的数据。

**方法**：
- 使用滑动窗口连续伪迹检测算法(`pop_continuousartdet`)
- **参数**（从Excel文件`ICA_Prep_Values_motor.xlsx`中读取，每个被试可独立设置）：
  - `AmpthValue`：幅度阈值(μV)，通常设为100-150 μV
  - `WindowValue`：滑动窗口大小(ms)，通常设为200-500 ms
  - `StepValue`：窗口滑动步长(ms)，通常设为50-100 ms
- 检测到的伪迹段被标记并从后续ICA计算中排除

**输出**：

- `*_eeg_ds_reref_hpfilt_ica_pre2.set`

### 2.4 独立成分分析 (STEP3)

**目的**：分解EEG信号为独立成分，主要去除眼动。

**方法**：
- ICA算法：Extended Infomax (`pop_runica`)
- 输入通道：1-60通道(不包括耳电极通道61-62)
- **自动分类**：使用ICLabel插件自动分类ICA成分为：Brain, Muscle, Eye, Heart, Line Noise, Channel Noise, Other
- 将ICA权重从准备好的数据集(`ica_pre2`)转移到完整的连续数据集(`hpfilt_notch`)

**输出**：
- `*_eeg_ds_reref_hpfilt_ica_weighted.set`：包含ICA权重的完整数据集
- `*_motor_ICA_Weights.png`：ICA成分地形图(包含ICLabel分类标签)

### 2.5 ICA成分移除 (STEP4)

**目的**：移除与眼动、肌电等伪迹相关的ICA成分。

**方法**：
- 从Excel文件(`ICA_Components_motor.xlsx`)中读取每个被试需要移除的成分编号【★这个文件需要根据上一步ICA输出的图自己整理】
- 使用`pop_subcomp`从数据中减去选定的ICA成分

**输出**：
- `*_eeg_ds_reref_hpfilt_ica_corr.set`：ICA校正后的连续数据

### 2.6 Event List和Epoching (STEP5)

**目的**：创建event列表、分配bin、分epoch。

**方法**：

**Event List创建**：
- 使用`pop_creabasiceventlist`创建所有event的记录
- 保存到`*_motor_Eventlist.txt`

**Bin分配 (Binlister)**：
- 使用`pop_binlister`根据bin描述文件(`BDF_motor.txt`)将event分配到bin
- Bin定义：
  - Bin 1：左手握拳(Marker = 1)
  - Bin 2：右手握拳(Marker = 2)
- 结果保存到`*_motor_Eventlist_Bins.txt`

**Epoching**：
- Epoch窗口：**-2000到40000 ms**(相对于block起始marker)
- **不进行baseline校正**(`'none'`)
- 理由：
  - -2000ms允许捕捉准备电位(readiness potential)
  - 40000ms覆盖整个30秒握拳期 + 部分休息期
  - 在ERD分析中再进行baseline校正，保留原始数据的灵活性

**输出**：
- `*_eeg_ds_reref_hpfilt_ica_corr_elist_bins_epoch.set`

### 2.7 伪迹剔除 (STEP6)

**目的**：插值坏通道，识别并标记包含伪迹的epoch。

**坏通道插值 (Channel Interpolation)**：
- 从Excel文件(`Interpolate_Channels_motor.xlsx`)读取每个被试的坏通道列表
- 使用球面插值法(`spherical`)进行插值
- 排除耳电极(61, 62)不参与插值

**伪迹检测**：

1. **简单电压阈值法 (Simple Voltage Threshold, SVT)**：
   - 检测epoch中电压超出阈值的样本
   - 参数(从`AR_Parameters_for_SVT_CRAP_motor.xlsx`读取)：
   - Flag：[1 2]表示同时标记artifact和user flag
   
2. **移动窗口峰峰值法 (Moving Window Peak-to-Peak, MW)**：
   - 在滑动窗口内计算峰峰值，检测是否超过阈值
   - 参数(从`AR_Parameters_for_MW_CRAP_motor.xlsx`读取)：
   - Flag：[1 3]

**注意**：epoch被标记但**不被删除**，在后续平均时使用`'Exclude epochs marked during artifact detection'`选项排除这些epoch。

**输出**：
- `*_eeg_ds_reref_hpfilt_ica_corr_elist_bins_epoch_interp_ar.set`

### 2.8 运动相关去同步化分析 (ERD Analysis, STEP7)

**目的**：计算μ和β频段的事件相关去同步化(Event-Related Desynchronization, ERD)，比较对侧vs同侧ROI。

**ERD定义**：
- ERD = [(任务期功率 - 基线期功率) / 基线期功率] × 100%
- 负值表示去同步化(功率降低)，正值表示同步化(功率增加)
- 运动执行通常引起对侧感觉运动皮层μ/β频段的ERD

**ROI定义**：
- 左侧ROI：C3, C1, C5, CP3, FC3
- 右侧ROI：C4, C2, C6, CP4,FC4

**对侧-同侧映射**：
- 左手握拳 → 对侧ROI = 右侧，同侧ROI = 左侧
- 右手握拳 → 对侧ROI = 左侧，同侧ROI = 右侧

**时间窗口**：
- Baseline窗口：**-2000到0 ms**
- 分析窗口：**-2000到40000 ms**(整个epoch)
- 任务窗口：**0到30000 ms**(用于计算平均ERD值)

**频段定义**：
- μ频段：**8-12 Hz**
- β频段：**13-30 Hz**

**时频分析参数 (`newtimef`)**：
- 频率范围：**4-40 Hz**(完整的时频分析范围)
- 小波周期数：**[3, 0.5]**
  - 3个周期用于低频(提高频率分辨率)
  - 0.5个周期用于高频(提高时间分辨率)
- 时间点数：**200**个时间点(在分析窗口内均匀分布)
- Padding比率：**2**(用于FFT)
- Baseline校正方式：dB单位(10*log10(功率))，相对于baseline窗口

**ERD计算步骤**：

1. **Good trials筛选**：
   - 排除所有EEG.reject字段中标记的trial
   - 排除包含'boundary' event的trial

2. **逐trial时频分析**：
   - 对每个good trial，提取ROI通道的平均信号
   - 使用`newtimef`计算时频谱(event-related spectral perturbation, ERSP)
   - ERSP返回dB单位(相对于baseline)

3. **频段平均**：
   - 在μ频段(8-12 Hz)和β频段(13-30 Hz)内对频率维度取平均
   - 得到μ-ERD和β-ERD的时间过程(单位：dB)

4. **dB转百分比**：
   - ERD(%) = (10^(ERD_dB/10) - 1) × 100

5. **跨trials平均**：
   - 将相同条件、相同ROI的trials平均，得到个体水平的ERD时间过程

**统计汇总**：

对每个被试、每个条件、每个频段计算：

1. **ROI平均ERD**：
   - 在任务窗口(0-30000 ms)内对时间维度取平均
   - 结果保存到`STEP7_ERD_Summary.xlsx` (Sheet: ROI)
   - 包含：Subject, Condition, ROI, Band, MeanERD_Pct, GoodTrials, TaskStart_ms, TaskEnd_ms

2. **对侧-同侧对比**：
   - 对侧平均ERD：ContraMeanERD_Pct
   - 同侧平均ERD：IpsiMeanERD_Pct
   - 对侧-同侧差值：ContraMinusIpsi_Pct = ContraMeanERD - IpsiMeanERD
   - 侧化指数(Lateralization Index, LI)：
     - LI = (ContraMag - IpsiMag) / (ContraMag + IpsiMag)
     - 其中Mag = -ERD(取绝对值，因为ERD为负值)
   - 结果保存到`STEP7_ERD_Summary.xlsx` (Sheet: ContraIpsi)

**可视化输出**：

1. **ERD时间过程图**：
   - μ和β频段分别绘制
   - 每个条件显示对侧ROI vs 同侧ROI
   - 包含95%置信区间(CI = 1.96 × SEM)
   - 文件：`ERD_Mu_Timecourse.png`, `ERD_Beta_Timecourse.png`

2. **对侧-同侧差值时间过程**：
   - μ和β频段分别绘制
   - 显示ContraERD - IpsiERD的时间演变
   - 文件：`ERD_Mu_ContraMinusIpsi.png`, `ERD_Beta_ContraMinusIpsi.png`

**输出文件**：
- `STEP7_ERD_Summary.xlsx`：ROI和对侧-同侧统计汇总
- `STEP7_ERD_Timecourses.mat`：完整的ERD时间序列数据(用于后续EEG-fNIRS偶联分析)
- PNG图像：组平均的时间过程图

**理由**：
- μ/β频段ERD是运动皮层激活的经典电生理指标
- 对侧皮层激活显著强于同侧，提供侧化性验证
- 时频分析比传统滤波-希尔伯特变换提供更好的时频权衡
- dB单位便于比较，百分比单位便于解释

---

## 3. fNIRS数据处理流程

### 3.1 预处理 (STEP1)

**目的**：从原始SNIRF文件进行全流程预处理，得到血红蛋白浓度时间序列。

#### 3.1.1 数据加载

**方法**：
- 使用Homer3的`SnirfClass`加载.snirf文件
- 提取：
  - `intensity`：原始光强时间序列(nT × nMeas)
  - `probe`：光极布局(源、探测器位置、波长)
  - `time`：时间向量(秒)

**波长检测**：
- 从`measurementList`中推断波长数(nW)
- 本研究使用3波长系统(nW = 3)

#### 3.1.2 通道质量控制 (Channel QC)

**目的**：识别并排除信噪比低、源探距离不当或头皮耦合差的通道。

**非正强度处理**：
- 检测光强 ≤ 0的样本(对数转换无效)
- 如果某通道的非正强度比例 > 5%，则跳过该run
- 否则将非正强度的通道标记为不活跃

**信噪比(SNR)检测**：
- SNR = mean(intensity) / std(intensity)
- 阈值：**SNR ≥ 2**
- 低于阈值的通道被标记为低SNR

**源探距离(SD)检测**：
- 从probe的2D位置计算欧氏距离
- 阈值：**0-45 mm**
- 超出范围的通道被排除

**头皮耦合指数(Scalp Coupling Index, SCI)**：

SCI评估光极与头皮的接触质量，基于心脏搏动信号的耦合强度。

**原理**：
- 良好的光极-头皮耦合应在心搏频段(0.5-2.5 Hz)显示高相关
- 对同一源-探测器对的两个波长信号，心搏引起的光强变化应高度相关

**计算步骤**：
1. 对两个波长的光强信号进行Z-score标准化
2. 使用3阶Butterworth带通滤波器提取心搏频段(0.5-2.5 Hz)
3. 在滑动窗口(10秒)内计算两个波长信号的Pearson相关系数
4. 取相关系数的中位数作为该源-探测器对的SCI值
5. **阈值**：SCI ≥ 0.6

**通道活跃性判断**：
- 对于每个源-探测器对，所有波长(3个)必须同时满足SNR、SD和SCI条件
- 只有当一个源-探测器对的所有波长都通过质量控制时，该对的测量才被保留(`mlActAuto = 1`)

**自动dRange估计**：
- 不同fNIRS设备的光强scaling不同，使用自动dRange更稳健
- 对每个波长，取光强均值的1%和99%分位数
- dRange = [P1 × 0.1, P99 × 10]

**输出**：
- `mlActAuto`：通道活跃性mask (1=keep, 0=reject)

#### 3.1.3 光密度转换

**方法**：
- 使用`hmrR_Intensity2OD`将光强转换为光密度(Optical Density, OD)
- OD = -log(I / I0)，其中I0是参考光强(通常为首个样本或平均值)

#### 3.1.4 运动伪迹检测

**目的**：识别头动引起的突然光强变化。

**方法**：
- 使用`hmrR_MotionArtifactByChannel`在OD信号上检测运动伪迹
- **参数**：
  - `tMotion = 0.5` s：运动artifact的最小持续时间
  - `tMask = 1.0` s：检测到artifact后mask的时长
  - `STDEVthresh = 20.0`：标准差阈值(倍数)
  - `AMPthresh = 0.5`：OD幅度阈值

**算法原理**：
- 计算滑动窗口内OD的标准差和幅度变化
- 当标准差或幅度超过阈值时标记为运动伪迹
- 生成两个mask：
  - `tIncAuto`：时间点mask (nT × 1)，1=干净，0=artifact
  - `tIncAutoCh`：逐通道时间点mask (nT × nMeas)

**输出**：
- `tIncAuto`, `tIncAutoCh`：运动伪迹mask
- QC指标：`motionFrac` = 被标记为运动的时间比例

#### 3.1.5 运动校正

**两步法**：结合TDDR和小波变换进行运动校正。

**1. TDDR (Temporal Derivative Distribution Repair)**：

**原理**：
- 运动伪迹在时域导数分布中表现为异常的尖峰
- TDDR通过修正导数分布的尾部来减轻运动伪迹

**方法**：
- 在OD信号上应用TDDR (`TDDR.m`)
- 采样率：自动从时间向量计算(fs = 1/median(diff(t)))
- 仅对有限值通道进行处理

**2. Wavelet运动校正**：

**原理**：
- 运动伪迹在小波域表现为高幅度系数
- 通过识别和衰减异常小波系数来校正运动伪迹

**方法**：
- 使用`hmrR_MotionCorrectWavelet`进行小波校正
- **参数**：
  - `iqrWavelet = 1.5`：IQR系数阈值，用于识别异常值

**理由**：
- TDDR在时域处理，适合连续性伪迹
- Wavelet在频域处理，适合瞬时尖峰
- 两者结合可以有效去除各种类型的运动伪迹

#### 3.1.6 带通滤波

**目的**：去除生理干扰(呼吸、Mayer波、心跳)和极低频漂移，保留任务相关的血流动力学信号。

**方法**：
- 使用`hmrR_BandpassFilt`对OD信号进行带通滤波
- **参数**：
  - 高通：**0.01 Hz**
  - 低通：**0.20 Hz**

**理由**：
- 血流动力学响应的典型频率 < 0.2 Hz
- 0.01 Hz高通去除基线漂移
- 0.20 Hz低通去除心跳(~1 Hz)、呼吸(~0.3 Hz)和高频噪声
- 保留任务相关的缓慢血流动力学变化

#### 3.1.7 血红蛋白浓度计算

**目的**：使用修正的Beer-Lambert定律(Modified Beer-Lambert Law, MBLL)将OD转换为血红蛋白浓度。

**方法**：
- 使用`hmrR_OD2Conc`进行转换
- **参数**：
  - PPF (Partial Pathlength Factor)：**6.0** (对所有3个波长)

**MBLL公式**：
- ΔC = [ε × PPF × L]^(-1) × ΔOD
- 其中：
  - ΔC：血红蛋白浓度变化矩阵[ΔHbO, ΔHbR, ΔHbT]
  - ε：消光系数矩阵(波长 × 血红蛋白类型)
  - PPF：部分路径长度因子
  - L：源-探测器距离

**PPF选择理由**：
- PPF = 6.0是成人头部的典型值
- 不同组织和年龄的PPF略有差异，但6.0是常用起始值

**输出**：
- `dc`：血红蛋白浓度时间序列，包含HbO、HbR、HbT通道
- `dc.measurementList`：包含`dataTypeLabel` ("HbO", "HbR", "HbT")

#### 3.1.8 质量控制指标

**输出指标**：
- `nMeasTotal`：总测量通道数
- `nMeasActive`：通过QC的活跃通道数
- `activeMeasFrac`：活跃通道比例
- `motionFrac`：运动伪迹时间比例
- `nFailSNRMeas/Pair`：未通过SNR的通道/源探对数
- `nFailSDMeas/Pair`：未通过SD的通道/源探对数
- `nFailSCIMeas/Pair`：未通过SCI的通道/源探对数

**QC flags**：
- `flag_lowActive`：活跃通道 < 60%
- `flag_highMotion`：运动时间 > 20%

**输出文件**：
- `*_task-motor_fnirs_preproc.mat`：包含dc, dod, qc, params, mlActAuto, tIncAuto, tIncAutoCh, stim
- `QC_summary.csv`：所有被试的QC指标汇总

### 3.2 质量可视化 (STEP2)

**目的**：可视化预处理的每个步骤，便于质量检查。

**输出图像**：

1. **Overview图**：
   - 原始光强(log10尺度，按波长平均)
   - 原始OD(按波长平均)
   - 处理后OD(TDDR+Wavelet+Filter，按波长平均)
   - 血红蛋白浓度(HbO, HbR, HbT跨通道平均)

2. **QC Masks图**：
   - 通道活跃性mask (`mlActAuto`)
   - 时间运动mask (`tIncAuto`)
   - 逐通道运动mask (`tIncAutoCh`, heatmap)
   - 运动检测flags

3. **SCI Summary图**：
   - 每个源-探测器对的SCI中位数(bar plot)
   - SCI阈值线(0.6)

**输出目录**：
- `sub-XXX/figures/`

### 3.3 Block平均分析 (STEP3)

**目的**：对每个block进行epoch并计算HbO/HbR的平均时间过程。

#### 3.3.1 Block检测

**方法**：
- 从`stim(1)`(左手)和`stim(2)`(右手)中提取event onset时刻
- 对每个条件，根据onset排序
- 使用onset间隔检测block边界：
  - 如果两个连续onset的间隔 > `blockGapSec`(20秒)，则认为是不同block的起始
  - **理由**：block内trial间隔为47秒(30s任务 + 17s休息)，block间隔通常更长

**Block参数估计**：
- Block持续时间 = 最后一个onset - 第一个onset + 30秒
- 取所有block持续时间的中位数作为后续分析的`blockDurationSecUsed`

#### 3.3.2 ROI通道选择

**ROI定义**：
- 左中央区：Source-Detector pairs = [14-10, 10-10, 10-14, 14-14, 14-18, 18-18, 18-14]
- 右中央区：Source-Detector pairs = [11-11, 15-11, 15-15, 11-15, 19-15, 19-19, 15-19]

**对侧ROI映射**：
- 左手握拳 → 使用右中央区ROI(对侧激活)
- 右手握拳 → 使用左中央区ROI(对侧激活)

**QC mask应用**：
- 使用`mlActAuto`排除预处理中标记的坏通道
- 通过源-探测器对匹配，将dod的QC mask映射到dc(血红蛋白浓度)
- 只有通过QC的通道参与ROI平均

**ROI平均**：
- 对ROI内所有HbO通道取平均 → ROI HbO时间序列
- 对ROI内所有HbR通道取平均 → ROI HbR时间序列

#### 3.3.3 Block Epoching和平均

**Epoch窗口**：
- 起始：`baselineWindowSec[1]` = **-5秒**(相对于block onset)
- 结束：`blockDurationSecUsed + postBlockSec` = **block持续时间 + 10秒**
- 例如，如果block持续30秒，则epoch为 [-5, 40]秒

**Baseline校正**：
- Baseline窗口：**-5到0秒**
- 方法：每个epoch的信号减去baseline窗口的平均值
- HbO_corrected = HbO - mean(HbO[baseline])
- HbR_corrected = HbR - mean(HbR[baseline])

**跨block平均**：
- 对同一条件的所有block，计算平均时间过程
- 计算标准误(SEM)和95%置信区间(CI = 1.96 × SEM)

#### 3.3.4 Block均值计算

**任务窗口**：**0到blockDurationSec**(例如0-30秒)

对每个block，计算任务窗口内的平均HbO和HbR值：
- `hboBlockMean = mean(HbO[0:30s])`
- `hbrBlockMean = mean(HbR[0:30s])`

#### 3.3.5 输出

**MAT文件** (`*_task-motor_blockavg.mat`)：
- `results.conditions(c).epochTime`：epoch时间轴
- `results.conditions(c).hboAvg/hbrAvg`：平均时间过程
- `results.conditions(c).hboSem/hbrSem`：标准误
- `results.conditions(c).hboCi/hbrCi`：95% CI
- `results.conditions(c).hboBlockMeans/hbrBlockMeans`：每个block的均值

**CSV文件** (`*_task-motor_blockavg_table.csv`)：
- 列：subjId, condition, roi, chromophore, blockIndex, blockMean, nChannels

**图像**：
1. `*_step-blockavg_qc_blockboundaries.png`：Block检测质量检查(onset序列、inter-block gaps)
2. `*_step-blockavg_qc_channels.png`：ROI通道数(HbO/HbR)
3. `*_step-blockavg_signal.png`：ROI平均时间序列(全程)
4. `*_step-blockavg_result.png`：Block平均时间过程(带CI)
5. `*_step-blockavg_trend.png`：Block均值趋势(检查疲劳效应或漂移)

### 3.4 一般线性模型分析 (GLM, STEP4)

**目的**：使用GLM估计任务相关的血流动力学响应，控制漂移和自相关。

#### 3.4.1 HRF建模

**Gamma函数HRF**：
- h(t) = [t^(shape-1) × exp(-t/scale)] / [Γ(shape) × scale^shape]
- **参数**：
  - shape = **6**
  - scale = **1**
  - 持续时间 = **30秒**

**理由**：
- Gamma函数可以模拟血流动力学响应的上升和下降
- 峰值时刻 ≈ (shape-1) × scale = 5秒，符合典型HRF
- 30秒持续时间足够长以捕捉完整的血流动力学响应

#### 3.4.2 任务回归量构建

**左手回归量**：
1. 从`stim(1)`提取所有左手event (onset, duration, amplitude)
2. 构建方波函数`uLeft`：在每个event的[onset, onset+duration]区间内，幅度 = amplitude
3. 与HRF卷积：`xLeft = conv(uLeft, hrf)`
4. 裁剪到原始时间长度：`xLeft = xLeft[1:nT]`

**右手回归量**：
- 同样方法构建`xRight`

**理由**：
- 卷积将神经事件转换为预期的血流动力学响应
- Duration = 30秒对应持续握拳期

#### 3.4.3 设计矩阵构建

**基础回归量**：
- `xLeft`：左手任务回归量
- `xRight`：右手任务回归量
- `intercept`：常数项(截距)

**漂移建模** (使用DCT)：
- 方法：离散余弦变换(Discrete Cosine Transform)
- **参数**：
  - `dctCutoffSec = 256`秒：高通截止周期
  - DCT basis数量：K = floor(2 × totalDuration / dctCutoffSec)

**DCT basis**：
- dct_k(n) = cos(π × (n + 0.5) × k / nT)，k = 1, 2, ..., K
- 每个basis代表一个低频正弦趋势

**理由**：
- DCT是fMRI和fNIRS分析中标准的漂移处理方法
- 等效于高通滤波，保留 > 1/(2×dctCutoffSec) Hz的频率
- 比多项式漂移更稳健，避免边缘效应

**最终设计矩阵X**：
- X = [xLeft, xRight, intercept, dct_1, dct_2, ..., dct_K]
- 列数(回归量数) p = 2(任务) + 1(截距) + K(漂移)

#### 3.4.4 ROI和通道筛选

**ROI筛选**：
- 仅保留左中央区和右中央区ROI内的通道
- 使用源-探测器对匹配

**对侧ROI-only策略**：
- 如果启用`useContraRoiOnly = true`：
  - 左手任务回归量只与右中央区ROI的通道回归
  - 右手任务回归量只与左中央区ROI的通道回归
- **理由**：对侧皮层激活最强，减少多重比较，提高统计效能

**QC mask**：
- 排除`mlActAuto = 0`的通道

**运动mask**：
- 使用`tIncAuto`和`tIncAutoCh`排除被标记为运动伪迹的时间点

#### 3.4.5 GLM估计 (AR-IRLS)

**AR-IRLS方法**：
- AR：Auto-Regressive (自回归)，建模时间序列的自相关
- IRLS：Iteratively Reweighted Least Squares (迭代加权最小二乘)，稳健回归

**算法**：
1. 估计AR模型：对残差拟合AR(Pmax)模型，估计自相关结构
2. 预白化(Pre-whitening)：使用AR模型对设计矩阵X和数据y进行白化
3. 稳健拟合：使用IRLS对白化后的数据拟合GLM，降低异常值的权重
4. 迭代收敛：重复步骤1-3直到收敛

**参数**：
- `arIrlsOrderSec = 1`秒
- `Pmax = round(arIrlsOrderSec / dt)`
  - 例如，dt = 1/11 ≈ 0.091秒，则Pmax ≈ 11

**理由**：
- fNIRS数据存在显著的时间自相关(生理噪声、仪器漂移)
- 忽略自相关会导致统计推断偏差(膨胀的t值)
- AR-IRLS比OLS更稳健，是fNIRS分析的金标准方法

**估计输出**：
- `beta`：回归系数 (p × nMeas)
- `SE`：标准误 (p × nMeas)
- `t`：t统计量 (p × nMeas)
- `p`：p值 (p × nMeas)
- `dfe`：有效自由度

#### 3.4.6 输出

**逐通道结果表** (`*_task-motor_glm_table.csv`)：
- 列：subjId, task, measIndex, src, det, chromophore, roi, term, beta, SE, t, p
- 每行代表一个通道、一个回归量的GLM结果

**ROI汇总表** (`*_task-motor_glm_roi_table.csv`)：
- 对每个ROI、每个血红蛋白类型、每个回归量：
  - `meanBeta`：ROI内通道的平均beta
  - `sdBeta`：标准差
  - `seBeta`：标准误 = sdBeta / sqrt(nChannels)

**MAT文件** (`*_task-motor_glm.mat`)：
- 包含beta, SE, t, p, dof, regNames, X, chrom, src, det, roiLabel等

**HRF图** (`*_step-glm_hrf.png`)：
- 上panel：Gamma HRF曲线
- 下panel：卷积后的左/右手回归量

---

## 4. EEG-fNIRS偶联分析

### 4.1 Block水平相关分析 (STEP1)

**目的**：评估EEG ERD与fNIRS血流动力学响应在block水平的相关性。

#### 4.1.1 数据准备

**EEG特征提取**：
- 从epoched EEG数据(STEP6输出)中提取每个block的ERD
- 对每个good trial：
  1. 提取对侧ROI通道的平均信号
  2. 计算baseline窗口([-2000, 0] ms)和任务窗口([0, 30000] ms)的功率谱密度(PSD)
     - 使用FFT：`Pxx = |fft(signal)|^2 / (fs × N)`
     - 双边谱转单边谱(2×)
  3. 在μ(8-12 Hz)和β(13-30 Hz)频段内积分PSD
     - `mu_power = trapz(f[μ], Pxx[μ])`
     - `beta_power = trapz(f[β], Pxx[β])`
  4. 计算ERD：
     - `mu_ERD = (mu_task / mu_baseline - 1) × 100`
     - `beta_ERD = (beta_task / beta_baseline - 1) × 100`

**fNIRS特征提取**：
- 从STEP3的block average表中读取每个block的HbO和HbR均值
- Block均值 = 任务窗口([0, 30000] ms)内的平均浓度

#### 4.1.2 Block对齐

**EEG block排序**：
- 根据event latency对同一条件的所有epoch排序
- 分配block index：1, 2, 3, ...

**fNIRS block排序**：
- 根据onset时刻对同一条件的所有event排序
- Block index从STEP3的table中读取

**对齐策略**：
- 使用block index对齐EEG和fNIRS
- 只保留两种模态都存在的block (intersection)
- 如果EEG某个block被artifact rejection剔除，则该block不参与相关分析

#### 4.1.3 相关分析

**被试内相关**：
- 对每个被试、每个条件(左/右/全部)：
  - Pearson相关：corr(EEG_MuERD, FNIRS_HbO)
  - Pearson相关：corr(EEG_MuERD, FNIRS_HbR)
  - Pearson相关：corr(EEG_BetaERD, FNIRS_HbO)
  - Pearson相关：corr(EEG_BetaERD, FNIRS_HbR)
- 最小block数要求：`minBlocksForCorr = 3`
- p值计算：t检验，dof = n - 2

**组水平相关**：
- 跨所有被试，pooling所有block
- 相同的4组相关分析
- 可以评估群体水平的神经血管偶联强度

#### 4.1.4 输出

**Block特征表** (`motor_eeg_fnirs_block_features.csv`)：
- 列：Subject, Condition, BlockIndex, EEG_MuERD_Pct, EEG_BetaERD_Pct, FNIRS_HbO_BlockMean, FNIRS_HbR_BlockMean
- 每行代表一个对齐的block

**对齐表** (`motor_eeg_fnirs_block_alignment.csv`)：
- 列：Subject, Condition, EEG_TotalBlocks, EEG_GoodBlocks, FNIRS_Blocks, MatchedBlocks
- 汇总每个被试的block对齐情况

**相关表** (`motor_eeg_fnirs_block_correlation.csv`)：
- 列：Subject, Analysis(Left/Right/All), EEGBand(Mu/Beta), FNIRSChrom(HbO/HbR), N, R, P, Detail
- 包含被试内和组水平的相关结果

**MAT文件** (`motor_eeg_fnirs_block_results.mat`)：
- 包含所有表格和参数

**理由**：
- Block水平相关评估任务驱动的神经血管偶联
- 对侧ROI映射确保对比同源皮层区域
- μ和β频段分别反映不同的运动皮层动态
- HbO通常比HbR显示更强的任务相关变化

### 4.2 EEG-informed fNIRS GLM (STEP2)

**目的**：将EEG ERD的时间过程作为回归量，评估EEG能否解释fNIRS信号的额外方差。

#### 4.2.1 EEG ERD时间过程提取

**对每个被试、每个条件、每个block**：

1. **时频分析**：
   - 对每个good trial，提取对侧ROI通道的平均信号
   - 使用`newtimef`计算时频谱(ERSP)
   - 参数与STEP7相同：
     - Baseline窗口：[-2000, 0] ms
     - 分析窗口：[-2000, 30000] ms
     - 频率范围：4-40 Hz
     - 小波周期：[3, 0.5]
     - 时间点数：200

2. **频段平均**：
   - μ频段(8-12 Hz)平均 → μ-ERD时间过程(dB)
   - β频段(13-30 Hz)平均 → β-ERD时间过程(dB)

3. **dB转百分比**：
   - ERD(%) = (10^(ERD_dB/10) - 1) × 100

4. **Block平均**：
   - 对同一条件、同一block index的所有good trials平均
   - 得到每个block的μ-ERD和β-ERD时间过程

**时间对齐**：
- EEG时间轴：相对于block onset的时间(ms)
- fNIRS时间轴：绝对时间(秒)
- 对齐：fNIRS时间 = block_onset + EEG时间

#### 4.2.2 EEG回归量构建

**对每个fNIRS时间点t**：

1. **确定对应的block和相对时间**：
   - 找到包含t的block (onset ≤ t < onset + duration)
   - 计算相对时间：rel_time = t - onset

2. **插值ERD值**：
   - 从该block的ERD时间过程中，在rel_time处插值(线性插值)
   - 得到该时间点的μ-ERD和β-ERD值

3. **累加到回归量**：
   - `uEegMuLeft[t] += mu_ERD_left_block[rel_time]`
   - `uEegBetaLeft[t] += beta_ERD_left_block[rel_time]`
   - 对左手和右手分别构建

4. **HRF卷积**：
   - `xEegMuLeft = conv(uEegMuLeft, hrf)`
   - 同样处理MuRight, BetaLeft, BetaRight

5. **Z-score标准化** (如果启用)：
   - 对每个EEG回归量，在非零时间点上计算均值和标准差
   - Z-score = (x - mean) / std
   - **理由**：使EEG回归量的尺度与任务回归量可比

6. **正交化** (可选)：
   - 将EEG回归量相对于任务回归量正交化
   - xEeg_orth = xEeg - X_task × (X_task \ xEeg)
   - **理由**：评估EEG的独特贡献，去除与任务回归量共享的方差

#### 4.2.3 缺失EEG block的处理

**问题**：
- 如果某个fNIRS block对应的EEG block由于artifact rejection缺失
- 该block的时间段内无EEG ERD信息

**策略** (如果启用`excludeMissingEegBlocks = true`)：
- 记录缺失EEG block的fNIRS时间段 → `missingBlockMask`
- 在GLM拟合时，排除这些时间点(tMask = tMask & ~missingBlockMask)
- **理由**：避免将缺失EEG信息的时间段归因于EEG回归量

#### 4.2.4 GLM设计矩阵

**完整设计矩阵X**：
- **任务回归量** (如果启用)：xLeft, xRight
- **EEG回归量** (如果启用)：
  - xEegMuLeft, xEegMuRight (如果includeMu = true)
  - xEegBetaLeft, xEegBetaRight (如果includeBeta = true)
- **截距**：ones(nT, 1)
- **漂移**：DCT basis (与STEP4相同)

**回归量命名**：
- 'left', 'right', 'eeg_mu_left', 'eeg_mu_right', 'eeg_beta_left', 'eeg_beta_right', 'intercept', 'dct01', 'dct02', ...

#### 4.2.5 GLM拟合

**方法**：
- 与STEP4相同，使用AR-IRLS
- 对每个fNIRS通道：
  1. 构建时间mask：排除运动伪迹、missing EEG blocks、非有限值
  2. AR(Pmax)预白化
  3. IRLS稳健拟合
  4. 得到beta, SE, t, p

**对侧ROI-only策略**：
- 左手任务和EEG回归量只与右中央区ROI回归
- 右手任务和EEG回归量只与左中央区ROI回归

#### 4.2.6 结果解释

**任务回归量的beta**：
- 反映基本的任务相关血流动力学响应
- 与STEP4的GLM结果类似(如果没有EEG回归量)

**EEG回归量的beta**：
- **显著的正beta**：
  - EEG ERD越负(去同步化越强)，fNIRS HbO越高
  - 表明神经电活动抑制与血流增加正相关
  - 符合神经血管偶联的预期：运动执行 → μ/β ERD → 局部代谢需求 → 血流增加
- **显著的负beta**：
  - EEG ERD越负，fNIRS HbR越低
  - 氧合血红蛋白增加，脱氧血红蛋白减少
- **不显著**：
  - EEG不能解释fNIRS的额外方差
  - 可能原因：神经血管解耦、空间不匹配、时间延迟

**模型比较**：
- 比较仅任务回归量的GLM vs 任务+EEG回归量的GLM
- 评估R²增量：EEG回归量解释的额外方差

#### 4.2.7 输出

**逐通道结果表** (`*_task-motor_eeg_informed_glm_table.csv`)：
- 列：subjId, task, measIndex, src, det, chromophore, roi, term, beta, SE, t, p

**ROI汇总表** (`*_task-motor_eeg_informed_glm_roi_table.csv`)：
- 每个ROI、每个血红蛋白类型、每个回归量的meanBeta, sdBeta, seBeta

**对齐表** (`*_task-motor_eeg_fnirs_alignment.csv`)：
- 列：Subject, Condition, BlockIndex, FnirsOnsetSec, FnirsDurSec, EegBlockAvailable, EegGoodTrials, EegMuMean, EegBetaMean
- 记录每个fNIRS block是否有对应的EEG数据

**汇总表** (`STEP2_eeg_informed_fnirs_summary.csv`)：
- 列：Subject, EEG_LeftBlocks, EEG_LeftGood, EEG_RightBlocks, EEG_RightGood, FNIRS_LeftBlocks, FNIRS_RightBlocks, Matched_Left, Matched_Right, MissingTimePct

**回归量图** (`*_task-motor_eeg_informed_regressors.png`)：
- 任务回归量时间过程
- EEG μ回归量时间过程
- EEG β回归量时间过程

**MAT文件** (`*_task-motor_eeg_informed_glm.mat`)：
- 包含beta, SE, t, p, dof, regNames, X, 以及EEG相关信息

**理由**：
- EEG-informed GLM直接评估神经血管偶联的时间动态
- 控制漂移和自相关，提供稳健的统计推断
- 可以评估不同频段(μ vs β)对血流动力学的不同贡献
- 对侧ROI策略确保空间匹配

---



## 代码文件列表

### EEG处理代码
```
eeg/motor/motor_EEG_Processing/scripts/
├── STEP0_marker.m                           # Marker提取
├── STEP1_Import_Raw_EEG_Shift_DS_Reref_Hpfilt.m  # 预处理
├── STEP2_ICA_Prep.m                        # ICA准备
├── STEP3_Run_ICA.m                         # ICA运行
├── STEP4_Remove_ICA_Components.m           # ICA成分移除
├── STEP5_Elist_Bin_Epoch.m                 # Epoching
├── STEP6_Artifact_Rejection.m              # 伪迹剔除
└── STEP7_ERD.m                             # ERD分析
```

### fNIRS处理代码
```
fnirs/motor/motor_fnirs_processing/
├── STEP1_preprocessing.m                    # 预处理
├── STEP2_qc.m                              # QC可视化
├── STEP3_BlockAvarage.m                    # Block平均
├── STEP4_GLM.m                             # GLM分析
├── STEP5_group_block.m                     # 组水平block平均
└── STEP6_group_GLM.m                       # 组水平GLM
```

### EEG-fNIRS偶联分析代码
```
eeg_fnirs_processing/motor_eeg_fnirs_processing/
├── STEP1_motor_block_correlation.m         # Block相关分析
└── STEP2_eeg_informed_fnirs.m              # EEG-informed GLM
```


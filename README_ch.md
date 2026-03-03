# Complex Geometry

<font color=yellow>
<strong>
<em>
复合物（蛋白+小分子）局部几何特征的提取
</em>
</strong>
</font>

## 1. 库选择策略

| 功能 | 使用库 | 选择原因 |
|------|--------|----------|
| 结构加载 | Biotite | 速度更快 |
| 质心计算 | Biotite | C实现，高效 |
| SASA计算 | Biotite | 优化的Shrake-Rupley算法 |
| CellList邻域搜索 | Biotite | 高效的空间索引 |
| 原子选择 | MDAnalysis | 更灵活的选择语法 |
| 接触计算 | MDAnalysis | 成熟的选择器支持 |
| 氢键分析 | MDAnalysis | 更好的原子类型识别 |
| 轨迹支持 | MDAnalysis | 原生支持轨迹文件 |

---

## 2. 公式详解

本文特征计算方法基于以下研究：

> **Predicting the impacts of mutations on protein-ligand binding affinity based on molecular dynamics simulations and machine learning methods**

### 2.1 Protein-ligand Closeness (cn)

**公式**：
$$c_n = -\frac{1}{|BSR|} \sum_{R_k \in BSR} \| \bar{X}_{R_k} - \bar{X}_{LIG} \|$$

其中：
- $BSR$ = 结合位点残基集合
- $\bar{X}_{R_k}$ = 残基$R_k$的几何中心
- $\bar{X}_{LIG}$ = 配体的几何中心

**实现（使用Biotite）**：
```python
# Biotite: 快速质心计算
ligand_center = bstruc.centroid(self.ligand_atoms)

total_dist = 0.0
for res_id, res_name, res_atoms in self.binding_residues:
    res_center = bstruc.centroid(res_atoms)
    total_dist += np.linalg.norm(res_center - ligand_center)

cn = -total_dist / len(self.binding_residues)
```

**物理解释**：负的平均距离，结合越紧密值越大（越正）。

---

### 2.2 Solvent Accessible Surface Area (SASA)

**公式**：
$$SASA = \sum_{i \in P} SASA(i)$$

其中$P$是蛋白质原子集合。

**实现（使用Biotite）**：
```python
# Biotite: 优化的SASA计算
sasa_vals = bstruc.sasa(self.protein_atoms, probe_radius=1.4)
return float(np.sum(sasa_vals))
```

**物理解释**：蛋白质暴露于溶剂的表面积，结合位点SASA越小通常意味着更紧密的结合。

---

### 2.3 Protein-ligand Orientation (ot)

**公式**：
$$o_t = \frac{1}{|BSR|} \sum_{R_k \in BSR} \left| \frac{\pi}{2} - \angle LBB_k \right|$$

几何解释：
```
        R_k (残基)
         /
        /
       B (结合位点中心) ---- L (配体)
```

**实现（使用Biotite）**：
```python
# 计算结合位点中心
all_coords = np.vstack([r[2].coord for r in self.binding_residues])
bs_center = np.mean(all_coords, axis=0)

# 计算每个残基的角度偏差
for res_id, res_name, res_atoms in self.binding_residues:
    res_center = bstruc.centroid(res_atoms)
    vec_res = res_center - bs_center
    vec_lig = ligand_center - bs_center
    
    cos_angle = np.dot(vec_res, vec_lig) / (norm(vec_res) * norm(vec_lig))
    angle = np.arccos(cos_angle)
    deviation = abs(np.pi/2 - angle)
```

**物理解释**：配体相对于结合位点的空间取向，值越大表示取向越偏离初始反应位置。

---

### 2.4 Protein-ligand Contacts (ct)

**公式**：
$$c_t = \sum_{i \in P} \mathbb{I}[\min_{j \in L} d(i,j) < t]$$

**实现（使用MDAnalysis）**：
```python
# MDAnalysis: 灵活的距离计算
for lig_atom in self.ligand_sel:
    nearby = self.protein_sel.select_atoms(
        f'point {x} {y} {z} {cutoff}'
    )
    contacts += len(nearby)
```

**扩展：接触类型细分**
```python
# 疏水接触 vs 极性接触
if elem == 'C':
    hydrophobic += 1  # 碳原子
elif elem in ['N', 'O', 'S']:
    polar += 1  # 氮、氧、硫
```

**物理解释**：
- 疏水接触：范德华相互作用
- 极性接触：可能形成氢键或盐桥

---

### 2.5 Interfacial Hydrogen Bonds (hb)

**公式**：
$$h_b = n_{D \in P, A \in L} + n_{D \in L, A \in P}$$

**氢键判定标准**：
- 距离 < 3.5 Å
- 角度 > 135°（供体-H-受体）

**实现（使用MDAnalysis）**：
```python
# 获取供体和受体
donors = bs_sel.select_atoms('name N')  # 氮原子作为供体
acceptors = bs_sel.select_atoms('name O')  # 氧原子作为受体

# 检查蛋白供体 -> 配体受体
for d in donors:
    for a in lig_acceptors:
        dist = np.linalg.norm(d.position - a.position)
        if dist < distance_cutoff:
            hbonds += 1
```

**物理解释**：氢键是特异性结合的关键作用力。

---

### 2.6 距离特征

| 特征 | 公式 | 实现 |
|------|------|------|
| **质心距离** | $\|\bar{X}_{P} - \bar{X}_{L}\|$ | Biotite `mass_center()` |
| **几何中心距离** | $\|C_P - C_L\|$ | Biotite `centroid()` |
| **最小距离** | $\min d(i,j)$ | Biotite `CellList` |

---

## 3. 代码架构

### 3.1 类结构

```python
class UnifiedProteinLigandGeometry:
    """"""
    
    def __init__(self, pdb_file, chain_id, ligand_resname, binding_cutoff):
        # 
        self.u = mda.Universe(pdb_file)      # MDAnalysis
        self.atoms = load_structure(pdb_file) # Biotite
        
        # 初始化选择
        self._setup_selections()
        self._find_binding_site()
    
    # ===== Biotite 功能 =====
    def calculate_closeness(self) -> float: ...
    def calculate_sasa(self) -> float: ...
    def calculate_orientation(self) -> float: ...
    def calculate_distance_features(self) -> Dict[str, float]: ...
    
    # ===== MDAnalysis 功能 =====
    def calculate_contacts(self, cutoff) -> Dict[str, int]: ...
    def calculate_hydrogen_bonds(self, distance, angle) -> Dict[str, int]: ...
    
    # ===== 统一接口 =====
    def calculate_all_features(self) -> Dict[str, float]: ...
```

### 3.2 结合位点识别

```python
def _find_binding_site(self):
    """使用Biotite识别结合位点"""
    ligand_center = bstruc.centroid(self.ligand_atoms)
    
    for res in protein:
        res_center = bstruc.centroid(res)
        dist = np.linalg.norm(res_center - ligand_center)
        
        if dist < self.binding_cutoff:  # 默认10Å
            self.binding_residues.append(res)
```

---

## 4. 示例

### 4.1 安装

```bash
# 克隆仓库
git clone git@github.com:MilesYyh/complex_geometry.git
cd complex-geometry/complex_geometry
pip install numpy biotite mdanalysis -i https://pypi.mirrors.ustc.edu.cn/simple
pip install -e .
```

### 4.2 示例运行

`example/` 文件夹包含PDB结构 `1ZZT.pdb` (二氢叶酸还原酶与底物DP9) 及示例输出。

```bash
# 运行示例
cd example
python -c "
from complex_geometry import ProteinLigandGeometry
geom = ProteinLigandGeometry('1ZZT.pdb', ligand_resname='DP9')
features = geom.calculate_all_features()
print(features)
"

# 或使用命令行
cpx-calc 1ZZT.pdb -l DP9
```

### 4.3 输出结果

**PDB: 1ZZT (二氢叶酸还原酶) + DP9 (底物)**

```
closeness                     : -7.7819
sasa                          : 19381.1641
sasa_binding_site             : 1190.0000
orientation                   : 0.3126

contacts_total                : 5
contacts_hydrophobic          : 0
contacts_polar                : 5

hydrogen_bonds                : 0

com_distance                  : 18.6272
cog_distance                 : 18.6731
min_distance                  : 0.0000

n_protein_atoms               : 3220
n_ligand_atoms                : 46
n_binding_residues           : 15
```

#### 特征解释

| 特征 | 值 | 解释 |
|------|-----|------|
| closeness | -7.78 | 负值表示配体与结合位点平均距离约7.8Å |
| sasa | 19381 Å² | 蛋白质总表面积 |
| sasa_bs | 1190 Å² | 结合位点表面积 |
| orientation | 0.31 | 较小值表示取向较理想 |
| contacts | 5 | 5对原子接触 |
| polar | 5 | 极性接触 |

---

## 5.扩展应用

### 5.1 突变体特征差异

```python
def calculate_mutation_delta(features_wt, features_mut):
    """计算突变前后的特征差异"""
    delta = {}
    for key in features_wt:
        delta[f'delta_{key}'] = features_mut[key] - features_wt[key]
    return delta
```

### 5.2 MD轨迹时间序列分析

```python
u = mda.Universe('structure.pdb', 'trajectory.xtc')

for ts in u.trajectory:
    # 每一帧计算特征
    features = calculate_all_features()
    time_series.append(features)
```

---

## 7. 附录：完整特征列表

| 类别 | 特征名 | 类型 | 描述 |
|------|--------|------|------|
| 几何 | closeness | float | 蛋白质-配体的结合程度 |
| 几何 | orientation | float | 蛋白质-配体取向 |
| 几何 | cog_distance | float | 几何中心距离 |
| 几何 | com_distance | float | 质心距离 |
| 几何 | min_distance | float | 最小原子距离 |
| 表面积 | sasa | float | 蛋白质总SASA |
| 表面积 | sasa_binding_site | float | 结合位点SASA |
| 接触 | contacts_total | int | 总接触数 |
| 接触 | contacts_hydrophobic | int | 疏水接触数 |
| 接触 | contacts_polar | int | 极性接触数 |
| 氢键 | hydrogen_bonds | int | 氢键数 |
| 信息 | n_protein_atoms | int | 蛋白原子数 |
| 信息 | n_ligand_atoms | int | 配体原子数 |
| 信息 | n_binding_residues | int | 结合位点残基数 |

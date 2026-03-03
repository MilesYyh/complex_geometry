# 发布指南

将 complex-geometry 发布到 PyPI。

## 1. 包结构

```
complex_geometry/
├── pyproject.toml          # 包配置文件
├── README.md               # 项目说明
├── LICENSE                 # MIT许可证
├── src/
│   └── complex_geometry/
│       ├── __init__.py     # 包入口
│       ├── core.py         # 核心功能
│       └── cli.py          # 命令行工具
└── tests/                  # 测试目录
```

## 2. 本地开发安装

```bash
# 克隆仓库
git clone https://github.com/yourusername/complex-geometry.git
cd complex-geometry

# 开发模式安装（可编辑）
pip install -e .

# 安装所有依赖
pip install -e ".[all]"
```

## 3. 使用方式

### 作为 Python 库

```python
from complex_geometry import ProteinLigandGeometry

# 不指定 chain_id，自动检测第一个链
geom = ProteinLigandGeometry('protein.pdb')

# 指定 chain_id
geom = ProteinLigandGeometry('protein.pdb', chain_id='A')

# 指定配体名称
geom = ProteinLigandGeometry('protein.pdb', ligand_resname='LIG')

# 调整结合位点阈值
geom = ProteinLigandGeometry('protein.pdb', binding_cutoff=12.0)

features = geom.calculate_all_features()
print(features)
```

### 作为命令行工具

```bash
# 自动检测链
cpx-calc protein.pdb

# 指定链
cpx-calc protein.pdb -c A

# 指定配体
cpx-calc protein.pdb -l LIG

# 指定特征
cpx-calc protein.pdb -l LIG --features closeness sasa contacts

# 输出到文件
cpx-calc protein.pdb -l LIG -o results.json
```

## 4. 发布到 PyPI

### 步骤 1: 准备账户

```bash
# 注册 PyPI 账户
https://pypi.org/account/register/

# 创建 API token
https://pypi.org/manage/account/#api-tokens
```

### 步骤 2: 配置发布工具

```bash
# 安装发布工具
pip install build twine

# 创建 ~/.pypirc 文件
cat > ~/.pypirc << EOF
[pypi]
username = __token__
password = pypi-Agznde...

[testpypi]
username = __token__
password = pypi-Agznde...
EOF
```

### 步骤 3: 构建和发布

```bash
# 1. 清理旧构建
rm -rf build/ dist/ *.egg-info

# 2. 构建源码包和 wheel
python -m build

# 3. 检查构建产物
twine check dist/*

# 4. 发布到 PyPI
twine upload dist/*

# 或者先发布到 TestPyPI 测试
twine upload --repository testpypi dist/*
```

### 步骤 4: 版本标签（推荐）

```bash
# 创建 git 标签
git tag -a v0.1.0 -m "Release v0.1.0"

# 推送到远程
git push origin v0.1.0
```

## 5. 发布到 Bioconda（可选）

```bash
# Fork https://github.com/bioconda/bioconda-recipes
# 修改 recipes/complex_geometry/meta.yaml
# 提交 PR
```

# Paper2Perturb

**[SimuCella](https://mindflow-cn.github.io/simucella/) 的数据整理智能体。**

Paper2Perturb 从科学论文和公共单细胞数据仓库中提取扰动实验信息，核验论文证据，并生成标准化 h5ad 数据。项目以可复用的 Codex/Claude Code skills 组织完整流程。

[English](README.md)

## 工作流程

```text
论文 PDF
  -> Markdown 与图表
  -> 扰动实验元数据与论文证据
  -> 公共 scRNA-seq 数据下载
  -> 细胞类型注释
  -> 标准化 h5ad
  -> 元数据与 h5ad 一致性校验
```

当前主要面向人类单细胞小分子药物扰动研究。符合条件的数据写入 `result.xlsx`，不适合主数据集但仍有整理价值的实验写入 `result_excluded.xlsx`。

## Skills

| Skill | 功能 |
|---|---|
| `extract-perturbation` | 提取论文、药物、细胞、剂量、时间、靶基因、调控方向和证据。 |
| `build-h5ad` | 下载公共 scRNA-seq 数据并生成标准化 h5ad 与 JSON。 |
| `validate-metadata` | 同时执行字段规则检查和论文证据核验。 |
| `validate-h5ad` | 检查 h5ad、`test_case.json` 和表格元数据的一致性。 |

`validate-metadata` 同时包含确定性字段校验和论文证据核验。这是一次元数据审计的两个阶段，因此不拆分成两个 skill。

## 快速开始

需要 Python 3.10+，以及 Codex 或 Claude Code。PDF 转换需要 MinerU API key，图表理解需要 Qwen API key。

```bash
git clone https://github.com/mindflow-cn/paper2perturb.git
cd paper2perturb
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
cp .env.example .env
```

在 `.env` 中配置 MinerU 和 Qwen API：

```bash
MINERU_API_KEY=your-mineru-jwt-token
Qwen_API_KEY=your-dashscope-api-key
```

`MINERU_API_KEY` 用于自动将 PDF 转换为 Markdown，`Qwen_API_KEY` 用于提取过程中理解论文图表。也可以通过环境变量 `DASHSCOPE_API_KEY` 配置 Qwen；该变量的优先级高于 `.env` 中的 `Qwen_API_KEY`。

将全部 skill 软链接到工作项目：

```bash
./scripts/install-skills.sh /path/to/working-project
```

该命令会同时安装到 `.agents/skills/` 和 `.claude/skills/`。安装后可以直接请求 agent：

```text
使用 $extract-perturbation 处理 papers/example.pdf。
使用 $build-h5ad 处理 result.xlsx 中的 GSE139129。
使用 $validate-metadata 对照本地论文核验 PMID 34591417。
```

Agent 客户端也可以根据自然语言请求自动触发相应 skill。

每个 skill 的完整输入、输出和执行规则见对应目录中的 `SKILL.md`。

## 仓库结构

```text
paper2perturb/
├── skills/                  # Agent skills 及其实现
├── scripts/                 # 项目维护与安装工具
├── README.md
├── README_CN.md
├── requirements.txt
└── .env.example
```

论文、下载的表达矩阵、h5ad 文件、结果表格、日志和 API 凭据等运行时内容不会纳入版本控制。

## 参与贡献

每个 skill 应保持自包含：执行说明写入 `SKILL.md`，可复用代码放在 `scripts/`，需要按需读取的详细领域规则放在 `references/`。提交修改前请运行 `AGENTS.md` 中列出的校验命令。基础项目校验命令为：

```bash
python3 scripts/validate_project.py
```

## 许可证

Apache License 2.0，详见 [LICENSE](LICENSE)。

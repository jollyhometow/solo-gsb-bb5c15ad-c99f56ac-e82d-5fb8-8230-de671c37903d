# MSPC — 质谱同位素分布预计算 API

基于 **Python + FastAPI + SQLite** 的质谱预计算服务。客户端提交含元素计数的分子式（可带同位素标记），服务按各元素天然同位素丰度将分布表示为"质量-强度多项式"，逐元素卷积得到完整同位素包络，合并质量差小于 0.0001 Da 的峰，按基峰归一化后返回前 N 个峰及分子单同位素质量。

## 快速启动（单条命令）

```bash
docker compose up --build
```

启动过程会自动完成依赖安装（镜像构建时）与 SQLite 数据库初始化（服务启动时建表），无需任何手工迁移。

## 访问地址

| 地址 | 说明 |
|---|---|
| http://localhost:8000 | 服务根地址（端口 8000） |
| http://localhost:8000/docs | Swagger 交互式文档 |
| http://localhost:8000/health | 健康检查 |

## API 用法

### `POST /api/v1/envelope`

请求体：

```json
{
  "formula": "C6H12O6",
  "top_n": 5
}
```

- `formula`（必填）：分子式，如 `C6H12O6`；同位素标记用方括号表示，如 `C6H6[13C]2`（2 个纯 ¹³C 原子）、`C2H5[15N]O[18O]`。
- `top_n`（可选，默认 20，最大 200）：返回强度最高的前 N 个峰。

响应：

```json
{
  "formula": "C6H12O6",
  "composition": {"C": 6, "H": 12, "O": 6},
  "labels": {},
  "monoisotopic_mass": 180.063388,
  "base_peak_intensity": 100.0,
  "peak_count": 5,
  "peaks": [
    {"mass": 180.063388, "intensity": 100.0},
    {"mass": 181.066743, "intensity": 6.489437}
  ],
  "cached": false
}
```

- `monoisotopic_mass`：单同位素质量（各元素取天然丰度最高的同位素，标记原子取标记同位素质量）。
- `peaks`：基峰强度归一化为 100，按质量升序返回前 N 个峰。
- `cached`：结果是否来自 SQLite 缓存（相同 `formula` + `top_n` 的重复请求直接命中缓存）。

### `POST /api/v1/enrichment-fit`

由观测峰列表网格搜索估计单一同位素标签的富集度（enrichment）。

请求体：

```json
{
  "formula": "C6H12O6",
  "label": "13C",
  "peaks": [
    {"mass": 180.063388, "intensity": 100.0},
    {"mass": 181.062131, "intensity": 26.5},
    {"mass": 182.065485, "intensity": 3.4}
  ],
  "step": 1.0,
  "tolerance": 0.01
}
```

- `formula`（必填）：分子式，规则同 `/api/v1/envelope`。
- `label`（必填）：待估计的单个同位素标签，格式为质量数 + 元素符号，如 `13C`、`15N`、`2H`。
- `peaks`（必填，1～200 个）：观测峰，质量必须**严格递增**；强度在拟合前会归一化到基峰 100。
- `step`（可选，默认 1.0，范围 0.5～10）：丰度搜索步长（百分点）。候选丰度为 0%、step、2·step、…、100%（两端点始终包含）。
- `tolerance`（可选，默认 0.01，必须为正）：观测峰与预测峰配对的质量容差（Da）。

对每个候选丰度，将被标记元素的同位素分布替换为 `(1-f)·天然分布 + f·纯标记同位素`，复用包络卷积生成归一化包络（基峰 = 100）；每个观测峰与容差内**最近**的预测质量峰配对，未匹配的观测峰预测强度按 0 计，以强度均方误差（MSE）最小的丰度为最佳；同分时取较低丰度。

响应：

```json
{
  "formula": "C6H12O6",
  "label": "13C",
  "best_abundance_percent": 20.0,
  "mse": 0.0012,
  "search_points": 101,
  "pairs": [
    {
      "observed_mass": 180.063388,
      "observed_intensity": 100.0,
      "predicted_mass": 180.063388,
      "predicted_intensity": 100.0,
      "mass_error": 0.0
    },
    {
      "observed_mass": 185.5,
      "observed_intensity": 1.2,
      "predicted_mass": null,
      "predicted_intensity": 0.0,
      "mass_error": null
    }
  ]
}
```

- `best_abundance_percent`：最佳丰度（百分数）。
- `mse`：最佳丰度下的强度均方误差（对全部观测峰取平均）。
- `search_points`：实际搜索的候选丰度点数。
- `pairs`：最佳丰度下的逐峰配对；未匹配峰的 `predicted_mass`/`mass_error` 为 `null`，`predicted_intensity` 为 0。`observed_intensity` 为归一化（基峰 = 100）后的值。

错误定位（均返回 422）：

- 非法标签：`{"detail": {"error": "invalid_label", "message": "...", "position": 0}}`，`position` 为标签串内的 0 基位置。
- 空观测、非递增或重复质量：`{"detail": {"error": "invalid_observed_peaks", "message": "...", "index": 1}}`，`index` 为出错峰的 0 基下标（整体性问题如空列表时为 `null`）。
- 非法分子式：同 `/api/v1/envelope` 的 `invalid_formula`。
- `step`/`tolerance`/峰数越界等参数校验失败返回 FastAPI 默认的 422（`detail` 中含 `loc` 定位）。

### 错误定位

非法分子式或标记返回 HTTP 422，并指出出错字符的 **0 基位置**：

```json
{
  "detail": {
    "error": "invalid_formula",
    "message": "Element 'C' has no known isotope with mass number 20; check the label mass",
    "position": 2
  }
}
```

可检测的错误包括：未知元素符号、元素计数为 0、未闭合的 `[`、标记格式错误、未知标记质量数、孤立数字、非法字符等。

### 其他端点

- `GET /api/v1/elements`：列出全部支持的元素及其天然同位素质量/丰度。
- `GET /health`：健康检查。

## 配置（环境变量）

均可在 `docker-compose.yml` 的 `environment` 段覆盖：

| 变量 | 默认值 | 说明 |
|---|---|---|
| `MSPC_DB_PATH` | 容器内 `/data/mspc.db`，本地 `./mspc.db` | SQLite 数据库文件路径 |
| `MSPC_MERGE_TOLERANCE` | `0.0001` | 峰合并质量容差（Da） |
| `MSPC_PRUNE_THRESHOLD` | `1e-12` | 卷积过程中按相对强度剪枝的阈值 |
| `MSPC_TOP_N_MAX` | `200` | `top_n` 上限 |
| `MSPC_TOP_N_DEFAULT` | `20` | `top_n` 默认值 |

## 本地开发（不用 Docker）

```bash
pip install -r requirements-dev.txt
uvicorn app.main:app --reload          # 服务: http://127.0.0.1:8000
pytest                                 # 运行测试
```

## 实现说明

- **多项式卷积**：每种元素的同位素模式表示为多项式 `Σ aᵢ·x^mᵢ`（丰度 × 质量的幂），分子包络为各元素多项式（按原子数取幂，使用快速幂）之积；多项式乘法即 (质量, 强度) 峰列表的卷积。纯同位素标记贡献固定质量位移。
- **富集度拟合**：候选丰度 `f` 将被标记元素的同位素分布替换为 `(1-f)·天然分布 + f·δ(标记同位素)`，复用同一套卷积代码（`isotope_envelope` 接受可选的同位素表覆盖）生成候选包络；在 0%～100% 网格上按强度 MSE 搜索最优，同分取较低丰度。
- **峰合并**：排序后相邻质量差 < 0.0001 Da 的峰融合为一个，质量取强度加权平均，强度取和。
- **单同位素质量**：按 IUPAC 惯例取各元素天然丰度最高的同位素质量之和。
- **同位素数据**：`app/isotopes.py` 内置约 50 种常见元素的 AME/CIAAW 同位素质量与天然丰度，可按需扩充。
- **持久化**：SQLite 存计算结果缓存（`envelope_cache`）与请求日志（`request_log`），schema 随服务启动自动创建。

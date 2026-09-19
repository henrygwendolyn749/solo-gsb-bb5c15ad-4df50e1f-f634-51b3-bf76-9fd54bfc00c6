# PCR 引物二聚体筛查 API

基于 **Python + FastAPI + SQLite** 的引物二聚体（primer dimer）批量筛查服务。
调用方提交一批带名称的 IUPAC DNA 引物及两个阈值（连续互补阈值、3′端连续互补阈值），
服务枚举每对引物**自二聚体**与**两两异二聚体的全部错位对齐**，按 IUPAC 简并碱基
**集合语义**判断碱基互补，返回每对的最危险对齐、最长连续互补长度、涉及任一 3′端的
连续长度及中文风险原因；同分时按固定偏移规则取唯一结果。批次与结果持久化到
SQLite，支持事后查询。对已保存批次还可通过 `POST /api/batches/{id}/pool-plan`
基于异二聚体结果构建冲突图，按 DSatur 风格确定性贪心把引物分配到若干引物池；
通过 `GET /api/batches/{base_id}/compare/{target_id}` 按引物名称构造稳定配对键，
对比两个批次的配对风险升降（escalated / deescalated / introduced / resolved）；
通过 `POST /api/batches/{id}/replacement-ranking` 对某条待替换引物给出
2～50 条候选序列，纯内存模拟替换后的相关配对，按 high 配对数 / 最大 3′ 连续
互补 / 最大连续互补三个越小越优的指标做非支配分层（Pareto rank）。

---

## 一条命令启动

```bash
docker compose up --build -d
```

启动后访问：

| 用途 | 地址 |
| --- | --- |
| 交互式 API 文档（Swagger UI） | http://localhost:8000/docs |
| ReDoc 文档 | http://localhost:8000/redoc |
| 健康检查 | http://localhost:8000/health |
| 批次接口 | http://localhost:8000/api/batches |
| 引物分池 | `POST http://localhost:8000/api/batches/{id}/pool-plan` |
| 替换候选排序 | `POST http://localhost:8000/api/batches/{id}/replacement-ranking` |
| 批次风险对比 | `GET http://localhost:8000/api/batches/{base_id}/compare/{target_id}` |

> 监听端口 **8000**。SQLite 数据库文件保存在命名卷 `primer-data` 内
> （容器路径 `/data/primer_dimer.db`）。镜像构建期安装依赖、容器启动期自动建库
> （`docker-entrypoint.sh` → `python -m app.cli initdb` → `uvicorn`）。

停止与查看日志：

```bash
docker compose logs -f
docker compose down            # 保留数据卷
docker compose down -v         # 同时删除数据库卷
```

---

## 必要配置（环境变量，均有默认值，不配也能直接启动）

| 变量 | 默认值 | 说明 |
| --- | --- | --- |
| `PRIMER_DB_PATH` | `/data/primer_dimer.db` | SQLite 数据库文件路径 |
| `PRIMER_HOST` | `0.0.0.0` | 监听地址 |
| `PRIMER_PORT` | `8000` | 监听端口（容器内；改宿主机端口请编辑 `docker-compose.yml` 的 `ports`） |
| `PRIMER_MAX_PRIMERS` | `500` | 单批次引物数量上限 |
| `PRIMER_MAX_LENGTH` | `200` | 单条引物长度上限 |
| `PRIMER_DB_RESET` | `0` | 设为 `1` 时启动清空并重建表（仅用于重置数据） |

---

## API 概览

### `POST /api/batches` — 提交批次（同步计算后返回完整结果）

```bash
curl -s -X POST http://localhost:8000/api/batches \
  -H 'Content-Type: application/json' \
  -d '{
    "primers": [
      {"name": "F", "sequence": "GGGGGAAAAA"},
      {"name": "R", "sequence": "GGGGGTTTTT"}
    ],
    "continuous_threshold": 4,
    "three_prime_threshold": 3
  }'
```

请求字段：

| 字段 | 说明 |
| --- | --- |
| `primers[].name` | 引物名称，非空，批次内唯一 |
| `primers[].sequence` | 5′→3′ IUPAC DNA 序列，接受 `A T C G R Y S W K M B D H V N`（大小写、首尾空白自动规范化；含 `-`、空格等非 IUPAC 字符返回 422） |
| `continuous_threshold` | 连续互补长度阈值（≥1），达到即提示双链结合稳定 |
| `three_prime_threshold` | 3′端连续互补长度阈值（≥1），达到即提示聚合酶延伸风险，风险等级 `high` |

返回（节选）：

```json
{
  "id": "f98cdd27-...",
  "primer_count": 2,
  "pair_count": 3,
  "high_risk_count": 1, "medium_risk_count": 0, "low_risk_count": 2,
  "pairs": [
    {
      "pair_index": 2,
      "name_top": "F", "name_bottom": "R", "is_self_dimer": false,
      "alignments_checked": 19,
      "best_shift": 5,
      "overlap_len": 10,
      "continuous_complementary_length": 5,
      "three_prime_run_length": 5,
      "top3_run_length": 5, "bottom3_run_length": 5,
      "involved_three_prime_ends": ["F:3′", "R:3′"],
      "risk_level": "high",
      "risk_reasons": [
        "3′端连续互补达到阈值，存在聚合酶延伸风险：top 3′端连续互补 5 bp，bottom 3′端连续互补 5 bp",
        "最长连续互补 5 bp 达到阈值 4 bp，双链结合稳定"
      ],
      "alignment": {
        "top":    "GGGGGAAAAA     ",
        "match":  "     |||||     ",
        "bottom": "     TTTTTGGGGG"
      }
    }
  ]
}
```

### `GET /api/batches/{id}` — 批次结果查询

* `GET /api/batches/{id}` 返回该批次全部配对；
* `GET /api/batches/{id}?risky_only=true` 只返回 `high`/`medium` 风险配对；
* 不存在的批次返回 `404`。

### `GET /api/batches?limit=50&offset=0` — 批次列表（按创建时间倒序）

### `POST /api/batches/{id}/pool-plan` — 对已保存批次做引物分池

对**已保存批次**的异二聚体筛查结果构建冲突图并分池：

* **自二聚体忽略**，只取异二聚体配对；风险等级命中 `isolated_risk_levels`
  （`high` / `medium` / `low`，可多选、自动去重）的异二聚体两端引物之间连一条冲突边；
* 贪心选择下一个引物，依次比较 **① 已见邻池数**（相邻引物已占据的不同池数，
  DSatur saturation）、**② 冲突度**（冲突图中的边数）、**③ 批次位置**（靠前者优先）；
* 贪心选择放置池，依次比较 **① 新增冲突数**（放入后新产生的同池冲突边数，少者优先）、
  **② 池内人数**（少者优先）、**③ 池号**（小者优先）；
* 容量允许（`pool_count × pool_capacity ≥ 引物数`）时**全部引物必被分配**；
  容量总和不足返回 `422`；批次不存在返回 `404`。全部比较均为全序，相同输入结果稳定可复现。

```bash
curl -s -X POST http://localhost:8000/api/batches/f98cdd27-.../pool-plan \
  -H 'Content-Type: application/json' \
  -d '{
    "pool_count": 3,
    "pool_capacity": 8,
    "isolated_risk_levels": ["high", "medium"]
  }'
```

请求字段：

| 字段 | 说明 |
| --- | --- |
| `pool_count` | 池数量，≥1 |
| `pool_capacity` | 每池容量（每池最多引物数），≥1 |
| `isolated_risk_levels` | 需隔离的异二聚体风险等级列表，取值 `high`/`medium`/`low`，至少 1 个；重复值去重 |

返回示例：

```json
{
  "batch_id": "f98cdd27-...",
  "pool_count": 3,
  "pool_capacity": 8,
  "isolated_risk_levels": ["high", "medium"],
  "primer_count": 3,
  "assigned_count": 3,
  "total_conflicts": 0,
  "pools": [
    {
      "pool": 0,
      "capacity": 8,
      "members": ["R"],
      "member_positions": [1],
      "conflict_count": 0
    }
  ],
  "conflicts": [
    {
      "pool": 0,
      "position_a": 1, "position_b": 4,
      "name_a": "R", "name_b": "P5",
      "risk_level": "high"
    }
  ]
}
```

* `pools[].members` / `member_positions` 为各池成员，按批次位置升序；
  `conflict_count` 为该池**内部**同池冲突边数；
* `conflicts` 平铺全部池的同池冲突边（按池号、两端批次位置升序），
  长度恒等于 `total_conflicts`，也等于各池 `conflict_count` 之和；
* 容量允许但冲突图无法完全隔离时（如池数少于色数），算法仍分配全部引物，
  不可避免的冲突通过 `conflicts` 如实返回。

### `GET /api/batches/{base_id}/compare/{target_id}` — 批次风险对比

按引物**名称**构造稳定配对键，对比 base → target 两个已保存批次的配对风险：

* **配对键**：自二聚体保留单名（如 `F`）；异二聚体忽略批次内提交顺序，
  两名按字典序排列（如 `F|R`，与提交时谁先谁后无关）；
* **风险等级全序** `low < medium < high`：两批次共有的配对，等级上升进入
  `escalated`、下降进入 `deescalated`、不变不计入明细（计入
  `summary.unchanged_pair_count`）；仅存在于 target 的配对进入 `introduced`，
  仅存在于 base 的进入 `resolved`；
* **同名引物序列变化**仍按同一配对键比较，并在明细中标记
  `sequence_changed: true`（配对中任一同名引物序列不同即为 true）；
* 四类明细各自按配对名确定性排序，相同输入结果稳定可复现；
  任一批次不存在返回 `404`。

```bash
curl -s http://localhost:8000/api/batches/f98cdd27-.../compare/3fa9c1e2-...
```

返回示例：

```json
{
  "summary": {
    "base_batch_id": "f98cdd27-...",
    "target_batch_id": "3fa9c1e2-...",
    "base_pair_count": 15,
    "target_pair_count": 15,
    "common_pair_count": 10,
    "unchanged_pair_count": 7,
    "escalated_count": 1,
    "deescalated_count": 2,
    "introduced_count": 5,
    "resolved_count": 5
  },
  "escalated": [
    {
      "pair_key": "A1|A2",
      "name_a": "A1", "name_b": "A2", "is_self_dimer": false,
      "sequence_changed": true,
      "base_risk_level": "medium", "target_risk_level": "high"
    }
  ],
  "deescalated": [
    {
      "pair_key": "D1|D2",
      "name_a": "D1", "name_b": "D2", "is_self_dimer": false,
      "sequence_changed": true,
      "base_risk_level": "high", "target_risk_level": "low"
    }
  ],
  "introduced": [
    {
      "pair_key": "New",
      "name_a": "New", "name_b": "New", "is_self_dimer": true,
      "sequence_changed": false,
      "base_risk_level": null, "target_risk_level": "low"
    }
  ],
  "resolved": [
    {
      "pair_key": "A1|Gone",
      "name_a": "A1", "name_b": "Gone", "is_self_dimer": false,
      "sequence_changed": false,
      "base_risk_level": "high", "target_risk_level": null
    }
  ]
}
```

* `introduced` 明细的 `base_risk_level` 为 `null`，`resolved` 明细的
  `target_risk_level` 为 `null`；
* `summary.common_pair_count = unchanged + escalated + deescalated`，
  `introduced_count = target_pair_count - common_pair_count`，
  `resolved_count = base_pair_count - common_pair_count`。

### `POST /api/batches/{id}/replacement-ranking` — 替换引物候选的内存模拟与非支配分层

对**已保存批次**中的一条待替换引物，提交 2～50 条带**唯一标签**的候选 IUPAC 序列，
服务对每个候选做**纯内存**模拟：用候选序列替换目标引物，复用批次筛查同款二聚体
算法（`analyze_pair`）与**批次原阈值**，只计算

* 候选自身的**自二聚体** 1 个；
* 候选与批次中**其余**引物（不含目标）的异二聚体 N−1 个。

**不写库、不持久化，也不重算与候选无关的原批次配对**（其余引物彼此之间的配对
保持原样，无需参与评估）。异二聚体的 top/bottom 方向沿用“批次位置靠前者为 top”
的规则，因此结果与真实新建批次完全一致。

按批次原阈值汇总三个**越小越优**的指标：

| 指标 | 含义 |
| --- | --- |
| `high_risk_count` | 候选相关配对中风险等级为 `high` 的数量 |
| `max_three_prime_run` | 各配对 3′ 端连续互补长度（`terminal_run`）的最大值 |
| `max_continuous_complementary` | 各配对最长连续互补长度的最大值 |

**非支配分层**（Pareto 排序）：候选 A 支配候选 B，当且仅当 A 三个指标均不劣于 B
且至少一个严格更优；指标完全相同或互有胜负（tradeoff）的候选同层。第 1 层为
不被任何候选支配者，逐层剥离依次为第 2、3……层。返回的 `ranked_candidates`
按 `(rank, input_order)` 升序——**同层保持请求中的输入顺序**，结果稳定可复现。

```bash
curl -s -X POST http://localhost:8000/api/batches/f98cdd27-.../replacement-ranking \
  -H 'Content-Type: application/json' \
  -d '{
    "target_name": "F",
    "candidates": [
      {"label": "good", "sequence": "GGGGG"},
      {"label": "t-run", "sequence": "TTTT"},
      {"label": "pal", "sequence": "ACGT"}
    ]
  }'
```

请求字段：

| 字段 | 说明 |
| --- | --- |
| `target_name` | 批次内待替换引物名称（首尾空白自动规范化；不存在返回 `422`） |
| `candidates` | 2～50 条候选，每条含唯一 `label`（非空、≤100 字符）与 `sequence`（IUPAC 校验规则同批次提交） |
| `candidates[].label` | 候选标签，同一请求内必须唯一，重复返回 `422` |

批次不存在返回 `404`；目标引物不存在、候选标签重复、候选数越界、空标签或
非法 IUPAC 序列均返回 `422`。

返回示例（节选）：

```json
{
  "batch_id": "f98cdd27-...",
  "target_name": "F",
  "continuous_threshold": 4,
  "three_prime_threshold": 3,
  "candidate_count": 3,
  "ranked_candidates": [
    {
      "rank": 1,
      "input_order": 0,
      "label": "good",
      "sequence": "GGGGG",
      "evaluated_pair_count": 3,
      "high_risk_count": 0,
      "max_three_prime_run": 0,
      "max_continuous_complementary": 0,
      "high_risk_pairs": [],
      "max_three_prime_pairs": [ { "...PairOut，同批次筛查字段；另附 label/other_name..." } ],
      "max_continuous_pairs": []
    }
  ]
}
```

每个候选结果中：

* `evaluated_pair_count` 恒等于 `1 +（批次引物数 − 1）`；
* `high_risk_pairs` 为达到 high 的配对；
* `max_three_prime_pairs` / `max_continuous_pairs` 为**达到对应最大值**的全部配对
  （最大值为 0 时即全部配对，便于定位“零风险由哪些配对构成”）；
* 配对对象复用批次筛查的 `PairOut` 全部字段，另附：
  `label`（候选标签，配对中的候选一侧）与 `other_name`（自二聚体为标签本身，
  异二聚体为批次原引物名）；异二聚体的 `name_top`/`name_bottom` 保持批次位置方向；
* 每个候选的 `pair_index` 独立从 0 起，计算顺序为自二聚体在前，
  异二聚体按其余引物的批次位置升序。

## 配对结果字段含义

| 字段 | 含义 |
| --- | --- |
| `is_self_dimer` | 是否为自二聚体（同一引物与自身对齐） |
| `alignments_checked` | 实际枚举的错位数量，恒为 `len(top)+len(bottom)-1` |
| `best_shift` | 最危险对齐的带符号错位（见下） |
| `overlap_len` | 该对齐重叠碱基数 |
| `continuous_complementary_length` | 重叠区内最长连续互补长度 |
| `top3_run_length` / `bottom3_run_length` | 触及对应 3′端的连续互补长度（端未参与则为 0） |
| `three_prime_run_length` | 两端的最大值（`max(top3, bottom3)`） |
| `involved_three_prime_ends` | 3′端连续互补长度 > 0 的端名称 |
| `risk_level` | `high`（3′端达阈值）/ `medium`（仅连续互补达阈值）/ `low` |
| `risk_reasons` | 中文风险原因列表 |
| `alignment` | 等宽三行可视化；中间行 `\|`=互补，`.`=不互补，空格=未重叠 |

## 错位 `shift` 定义与同分规则

两条引物按反向平行对齐，top 固定 5′→3′，bottom 按 3′→5′ 打印。
`shift`（记为 d）是 **bottom 3′端碱基所在的全局列坐标**，top 的 5′端 T0 位于列 0：

下列三行与 API 的 `alignment` 输出格式一致（等宽，空格为未重叠，`|` 为互补）：

```
d = 0（齐平）        d > 0（bottom 3′伸入 top）   d < 0（bottom 3′悬空于左侧）
T0 T1 T2             T0 T1 T2                        T0 T1 T2
|  |  |              .  |  |  |                  |  |  .
   B2 B1 B0             B2 B1 B0              B2 B1 B0
（渲染为连续字符行，例如 d=0 时 top="TTT"、bottom="BBB"）
```

* `d = 0`：bottom 3′端碱基与 top 5′端 T0 同列；
* `d > 0`：bottom 3′端向 top 3′方向伸入 top 内部；
* `d < 0`：bottom 3′端悬空于 top 5′端左侧。

枚举范围为 `-(b-1) <= d <= a-1`，共 `a+b-1` 个对齐，全部检查，不做启发式剪枝。

**最危险对齐**按以下键取最大者：

1. `max(top3_run, bottom3_run)`（聚合酶延伸风险，最重要）；
2. 最长连续互补长度（双链稳定性）；
3. 重叠碱基数。

**同分固定偏移规则**（保证结果唯一、可复现）：风险键完全相同时依次取

1. `|d|` 最小者（悬空端最短、最稳定构象）；
2. `d` 为正者优先，其次 0，最后负值；
3. 再按 `d` 升序兜底（全序，无并列残留）。

## IUPAC 简并碱基互补规则

每个简并码映射为碱基集合（如 `R={A,G}`、`Y={C,T}`、`N={A,C,G,T}`）。
同列两碱基“可能互补”当且仅当一个碱基集合的 Watson-Crick 互补集合与另一个
集合**交集非空**。因此 `R` 与 `Y` 互补、`N` 与任意碱基互补、`S={G,C}` 与
`W={A,T}` 不互补。只要存在至少一种简并实例化方式可形成 WC 配对即判为互补。

---

## 本地开发（不使用 Docker）

```bash
pip install -r requirements.txt
export PRIMER_DB_PATH=./data/primer_dimer.db
python -m app.cli initdb
uvicorn app.main:app --host 0.0.0.0 --port 8000
```

运行测试：

```bash
pip install pytest
python -m pytest -q
```

## 项目结构

```
app/
  config.py     # 环境变量配置
  iupac.py      # IUPAC 碱基集合与互补矩阵
  dimer.py      # 错位枚举、连续互补/3′端长度、风险判定、同分规则
  pooling.py    # 异二聚体冲突图构建与 DSatur 风格确定性贪心分池
  compare.py    # 批次间稳定配对键构造与风险升降分类（escalated 等四类）
  ranking.py    # 三指标非支配分层（Pareto rank），同层按输入顺序稳定排序
  schemas.py    # Pydantic 模型与输入校验
  db.py         # SQLite 建表、写入、查询
  services.py   # 批次计算编排（N 自二聚体 + N(N-1)/2 异二聚体）、分池、对比与替换排序编排
  main.py       # FastAPI 路由
  cli.py        # initdb 命令
tests/          # pytest 算法与 API 测试
Dockerfile
docker-compose.yml
docker-entrypoint.sh
requirements.txt
```

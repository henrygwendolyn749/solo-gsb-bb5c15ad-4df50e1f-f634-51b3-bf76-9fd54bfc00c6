"""替换引物候选排序接口（POST /api/batches/{id}/replacement-ranking）端到端测试。"""

from __future__ import annotations

import os

import pytest


@pytest.fixture(scope="module")
def client(tmp_path_factory):
    tmp_dir = tmp_path_factory.mktemp("ranking_db")
    os.environ["PRIMER_DB_PATH"] = str(tmp_dir / "test.db")

    from app import db, main

    db.init_db()
    from fastapi.testclient import TestClient

    with TestClient(main.app) as c:
        yield c


PAYLOAD = {
    "primers": [
        {"name": "F", "sequence": "GGGGGAAAAA"},   # pos 0（目标）
        {"name": "R", "sequence": "GGGGGTTTTT"},   # pos 1
        {"name": "X", "sequence": "AAAAAAAAAA"},   # pos 2
    ],
    "continuous_threshold": 4,
    "three_prime_threshold": 3,
}


@pytest.fixture(scope="module")
def batch_id(client):
    return client.post("/api/batches", json=PAYLOAD).json()["id"]


CANDIDATES = [
    {"label": "good", "sequence": "GGGGG"},    # 与 R/X 无互补，自二聚体低风险
    {"label": "t-run", "sequence": "TTTT"},    # 与 X 形成 high
    {"label": "pal", "sequence": "ACGT"},      # 自二聚体 high
]


def _rank(client, batch_id, target="F", candidates=None, **extra):
    body = {"target_name": target, "candidates": candidates or CANDIDATES, **extra}
    return client.post(f"/api/batches/{batch_id}/replacement-ranking", json=body)


# ---------- 正常分层 ----------

def test_ranking_fronts_and_ordering(client, batch_id):
    data = _rank(client, batch_id).json()

    assert data["batch_id"] == batch_id
    assert data["target_name"] == "F"
    assert data["continuous_threshold"] == 4
    assert data["three_prime_threshold"] == 3
    assert data["candidate_count"] == 3

    by_label = {c["label"]: c for c in data["ranked_candidates"]}

    # good 不被支配 -> 第 1 层；t-run / pal 均为 1 high、4/4 指标 -> 第 2 层
    assert by_label["good"]["rank"] == 1
    assert by_label["t-run"]["rank"] == 2
    assert by_label["pal"]["rank"] == 2

    good = by_label["good"]
    assert good["high_risk_count"] == 0
    assert good["max_three_prime_run"] == 0
    assert good["max_continuous_complementary"] == 0

    for cand in data["ranked_candidates"]:
        # 3 引物批次：1 自二聚体 + 2 异二聚体
        assert cand["evaluated_pair_count"] == 3

    # 整体顺序：按 (rank, input_order)
    labels = [c["label"] for c in data["ranked_candidates"]]
    assert labels == ["good", "t-run", "pal"]


def test_metrics_match_dimer_engine_independently(client, batch_id):
    """用现有 analyze_pair 独立重算，校验三个汇总指标与达成最大值的配对。"""
    from app.dimer import analyze_pair

    batch = client.get(f"/api/batches/{batch_id}").json()
    others = {p["name"]: p["sequence"] for p in batch["primers"] if p["name"] != "F"}

    data = _rank(client, batch_id).json()
    for cand in data["ranked_candidates"]:
        seq = cand["sequence"]
        pairs = [
            analyze_pair(
                cand["label"], seq, cand["label"], seq, 4, 3, is_self_dimer=True
            )
        ]
        for name, oseq in others.items():
            pairs.append(
                analyze_pair(cand["label"], seq, name, oseq, 4, 3)
            )

        assert cand["high_risk_count"] == sum(
            1 for p in pairs if p.risk_level == "high"
        )
        max_term = max(p.terminal_run for p in pairs)
        max_cont = max(p.max_complementary_run for p in pairs)
        assert cand["max_three_prime_run"] == max_term
        assert cand["max_continuous_complementary"] == max_cont

        assert {p["other_name"] for p in cand["max_three_prime_pairs"]} == {
            (p.name_top if p.name_bottom == cand["label"] else p.name_bottom)
            if not p.is_self_dimer
            else cand["label"]
            for p in pairs
            if p.terminal_run == max_term
        }
        assert {p["other_name"] for p in cand["max_continuous_pairs"]} == {
            (p.name_top if p.name_bottom == cand["label"] else p.name_bottom)
            if not p.is_self_dimer
            else cand["label"]
            for p in pairs
            if p.max_complementary_run == max_cont
        }
        assert {p["other_name"] for p in cand["high_risk_pairs"]} == {
            (p.name_top if p.name_bottom == cand["label"] else p.name_bottom)
            for p in pairs
            if p.risk_level == "high"
        }
        for pdct in cand["high_risk_pairs"]:
            assert pdct["label"] == cand["label"]
            assert not (pdct["is_self_dimer"] and pdct["other_name"] != cand["label"])


def test_pair_direction_preserves_batch_position(client, batch_id):
    """目标位于 pos 0：异二聚体 top 为候选标签；自二聚体两端均为标签。"""
    data = _rank(client, batch_id).json()
    trun = next(c for c in data["ranked_candidates"] if c["label"] == "t-run")
    high = trun["high_risk_pairs"]
    assert len(high) == 1
    p = high[0]
    assert p["is_self_dimer"] is False
    assert p["name_top"] == "t-run"
    assert p["name_bottom"] == "X"
    assert p["other_name"] == "X"
    assert p["continuous_complementary_length"] == 4
    assert p["three_prime_run_length"] == 4
    # 每个候选的 pair_index 独立从 0 起
    assert p["pair_index"] == 2

    pal = next(c for c in data["ranked_candidates"] if c["label"] == "pal")
    sp = pal["high_risk_pairs"][0]
    assert sp["is_self_dimer"] is True
    assert sp["name_top"] == sp["name_bottom"] == "pal"
    assert sp["other_name"] == "pal"
    assert sp["pair_index"] == 0


def test_single_primer_batch_evaluates_self_dimer_only(client):
    created = client.post(
        "/api/batches",
        json={
            "primers": [{"name": "solo", "sequence": "AAAATTTT"}],
            "continuous_threshold": 3,
            "three_prime_threshold": 3,
        },
    ).json()
    data = _rank(
        client,
        created["id"],
        target="solo",
        candidates=[
            {"label": "aa", "sequence": "AAAAAAAA"},
            {"label": "at", "sequence": "AAAATTTT"},
        ],
    ).json()

    by = {c["label"]: c for c in data["ranked_candidates"]}
    assert by["aa"]["evaluated_pair_count"] == 1
    assert by["at"]["evaluated_pair_count"] == 1
    assert by["aa"]["rank"] == 1
    assert by["at"]["rank"] == 2
    # 最大值为 0 时返回全部（唯一）配对
    zero_pairs = by["aa"]["max_three_prime_pairs"]
    assert len(zero_pairs) == 1
    assert zero_pairs[0]["other_name"] == "aa"
    # at 自二聚体 8bp 全互补 -> high
    assert by["at"]["high_risk_count"] == 1
    assert by["at"]["max_three_prime_run"] == 8


def test_same_metrics_keep_input_order_and_share_rank(client, batch_id):
    """指标完全相同的候选同层，且同层按输入顺序稳定排序。"""
    candidates = [
        {"label": "a1", "sequence": "TTTT"},
        {"label": "a2", "sequence": "TTTT"},
        {"label": "g", "sequence": "GGGGG"},
    ]
    data = _rank(client, batch_id, candidates=candidates).json()
    ranked = data["ranked_candidates"]
    assert [c["label"] for c in ranked] == ["g", "a1", "a2"]
    assert ranked[1]["rank"] == ranked[2]["rank"] == 2
    assert ranked[1]["input_order"] == 0 and ranked[2]["input_order"] == 1

    # 调换输入顺序：层级不变，同层顺序随输入变化
    swapped = [candidates[2], candidates[0], candidates[1]]
    data2 = _rank(client, batch_id, candidates=swapped).json()
    by = {c["label"]: c["rank"] for c in data2["ranked_candidates"]}
    assert by == {"g": 1, "a1": 2, "a2": 2}
    assert [c["label"] for c in data2["ranked_candidates"]] == ["g", "a1", "a2"]


def test_does_not_write_database(client, batch_id):
    before = client.get(f"/api/batches/{batch_id}").json()
    _rank(client, batch_id)
    after = client.get(f"/api/batches/{batch_id}").json()
    # 批次、引物与配对结果均不被修改
    assert after["pair_count"] == before["pair_count"] == 6
    assert after["primer_count"] == 3
    assert after["high_risk_count"] == before["high_risk_count"]
    assert [p["name"] for p in after["primers"]] == ["F", "R", "X"]
    listing = client.get("/api/batches").json()
    assert all(b["id"] != batch_id or b["pair_count"] == 6 for b in listing["batches"])


# ---------- 错误处理 ----------

def test_missing_batch_returns_404(client):
    resp = _rank(client, "no-such-batch")
    assert resp.status_code == 404


def test_missing_target_returns_422(client, batch_id):
    resp = _rank(client, batch_id, target="GHOST")
    assert resp.status_code == 422


def test_duplicate_labels_return_422(client, batch_id):
    resp = _rank(
        client,
        batch_id,
        candidates=[
            {"label": "x", "sequence": "AAAA"},
            {"label": "x", "sequence": "CCCC"},
        ],
    )
    assert resp.status_code == 422


def test_candidate_count_bounds(client, batch_id):
    one = _rank(client, batch_id, candidates=[{"label": "a", "sequence": "AAAA"}])
    assert one.status_code == 422

    fifty_one = [
        {"label": f"c{i:02d}", "sequence": "AAAA"} for i in range(51)
    ]
    too_many = _rank(client, batch_id, candidates=fifty_one)
    assert too_many.status_code == 422

    # 边界：2 条与 50 条合法
    two = _rank(
        client, batch_id,
        candidates=[{"label": "a", "sequence": "AAAA"},
                    {"label": "b", "sequence": "CCCC"}],
    )
    assert two.status_code == 200, two.text

    fifty = [
        {"label": f"c{i:02d}", "sequence": "AAAA"} for i in range(50)
    ]
    ok = _rank(client, batch_id, candidates=fifty)
    assert ok.status_code == 200, ok.text
    assert ok.json()["candidate_count"] == 50


def test_blank_label_and_bad_iupac_return_422(client, batch_id):
    blank = _rank(
        client, batch_id,
        candidates=[{"label": "  ", "sequence": "AAAA"},
                    {"label": "b", "sequence": "CCCC"}],
    )
    assert blank.status_code == 422

    bad_seq = _rank(
        client, batch_id,
        candidates=[{"label": "a", "sequence": "ACGX"},
                    {"label": "b", "sequence": "CCCC"}],
    )
    assert bad_seq.status_code == 422


def test_target_name_trimmed(client, batch_id):
    """target_name 首尾空白规范化后仍可匹配。"""
    resp = _rank(client, batch_id, target="  F  ")
    assert resp.status_code == 200
    assert resp.json()["target_name"] == "F"

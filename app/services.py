"""批次计算编排：枚举全部自二聚体与两两异二聚体并落库。"""

from __future__ import annotations

import json
import uuid

from . import db
from .compare import PairDiff, compare_pair_risks
from .dimer import PairResult, analyze_pair
from .pooling import build_pool_plan
from .ranking import non_dominated_ranks


def _pair_to_row(pair_index: int, r: PairResult) -> tuple:
    return (
        pair_index,
        r.name_top,
        r.name_bottom,
        1 if r.is_self_dimer else 0,
        r.alignments_checked,
        r.best_shift,
        r.overlap_len,
        r.max_complementary_run,
        r.top3_run,
        r.bottom3_run,
        r.terminal_run,
        json.dumps(list(r.involved_three_prime_ends), ensure_ascii=False),
        r.risk_level,
        json.dumps(list(r.risk_reasons), ensure_ascii=False),
        r.top_line,
        r.match_line,
        r.bottom_line,
    )


def run_batch(
    *,
    primers: list[tuple[str, str]],
    continuous_threshold: int,
    three_prime_threshold: int,
) -> str:
    """计算一个批次。

    primers 为 (name, normalized_sequence) 列表。
    枚举内容：N 个自二聚体 + N*(N-1)/2 个无序异二聚体。
    """
    batch_id = str(uuid.uuid4())
    created_at = db.utc_now_iso()

    results: list[PairResult] = []
    n = len(primers)

    # 自二聚体
    for name, seq in primers:
        results.append(
            analyze_pair(
                name, seq, name, seq,
                continuous_threshold, three_prime_threshold,
                is_self_dimer=True,
            )
        )
    # 两两异二聚体（无序对；top 固定为提交顺序中靠前者）
    for i in range(n):
        for j in range(i + 1, n):
            name_i, seq_i = primers[i]
            name_j, seq_j = primers[j]
            results.append(
                analyze_pair(
                    name_i, seq_i, name_j, seq_j,
                    continuous_threshold, three_prime_threshold,
                    is_self_dimer=False,
                )
            )

    risk_counts = {"high": 0, "medium": 0, "low": 0}
    for r in results:
        risk_counts[r.risk_level] += 1

    pair_rows = [
        _pair_to_row(idx, r) for idx, r in enumerate(results)
    ]
    primer_rows = [
        (pos, name, seq) for pos, (name, seq) in enumerate(primers)
    ]

    with db.get_conn() as conn:
        db.insert_batch(
            conn,
            batch_id=batch_id,
            created_at=created_at,
            continuous_threshold=continuous_threshold,
            three_prime_threshold=three_prime_threshold,
            primer_count=n,
            pair_count=len(results),
            risk_counts=risk_counts,
            primers=primer_rows,
            pair_rows=pair_rows,
        )
    return batch_id


def plan_batch_pools(
    *,
    batch_id: str,
    pool_count: int,
    pool_capacity: int,
    isolated_risk_levels: list[str],
) -> dict:
    """基于已保存批次的异二聚体结果构建冲突图并确定性分池。

    批次不存在返回 None（由路由层转 404）；容量总和不足抛出
    ValueError（由路由层转 422）。
    """
    with db.get_conn() as conn:
        batch = db.get_batch(conn, batch_id)
        if batch is None:
            return None
        primers = db.get_primers(conn, batch_id)
        pairs = db.get_pairs(conn, batch_id, risky_only=False)

    n = len(primers)
    if pool_count * pool_capacity < n:
        raise ValueError(
            f"池容量总和 {pool_count}×{pool_capacity}="
            f"{pool_count * pool_capacity} 不足以容纳 {n} 条引物"
        )

    # 批次内引物名称唯一（入库时已保证）
    position_by_name = {row["name"]: row["position"] for row in primers}
    names_by_position = {pos: name for name, pos in position_by_name.items()}

    # 仅取异二聚体（忽略自二聚体），且风险等级在隔离集合内
    levels = set(isolated_risk_levels)
    edges: list[tuple[int, int, str]] = []
    for row in pairs:
        if row["is_self_dimer"] or row["risk_level"] not in levels:
            continue
        edges.append(
            (
                position_by_name[row["name_top"]],
                position_by_name[row["name_bottom"]],
                row["risk_level"],
            )
        )

    plan = build_pool_plan(
        primer_count=n,
        edges=edges,
        pool_count=pool_count,
        pool_capacity=pool_capacity,
    )

    # 每个池的冲突边数
    conflict_count_by_pool = [0] * pool_count
    conflicts_out: list[dict] = []
    for pool_id, edge in plan.conflicts:
        conflict_count_by_pool[pool_id] += 1
        conflicts_out.append(
            {
                "pool": pool_id,
                "position_a": edge.position_a,
                "position_b": edge.position_b,
                "name_a": names_by_position[edge.position_a],
                "name_b": names_by_position[edge.position_b],
                "risk_level": edge.risk_level,
            }
        )

    pools_out = [
        {
            "pool": pool_id,
            "capacity": pool_capacity,
            "members": [names_by_position[pos] for pos in members],
            "member_positions": members,
            "conflict_count": conflict_count_by_pool[pool_id],
        }
        for pool_id, members in enumerate(plan.pool_members)
    ]

    return {
        "batch_id": batch_id,
        "pool_count": pool_count,
        "pool_capacity": pool_capacity,
        "isolated_risk_levels": isolated_risk_levels,
        "primer_count": n,
        "assigned_count": sum(len(m) for m in plan.pool_members),
        "total_conflicts": plan.total_conflicts,
        "pools": pools_out,
        "conflicts": conflicts_out,
    }


def _pair_result_to_dict(pair_index: int, r: PairResult) -> dict:
    """把内存中的 PairResult 转为与 PairOut 一致的字典（不落库）。"""
    return {
        "pair_index": pair_index,
        "name_top": r.name_top,
        "name_bottom": r.name_bottom,
        "is_self_dimer": r.is_self_dimer,
        "alignments_checked": r.alignments_checked,
        "best_shift": r.best_shift,
        "overlap_len": r.overlap_len,
        "continuous_complementary_length": r.max_complementary_run,
        "three_prime_run_length": r.terminal_run,
        "top3_run_length": r.top3_run,
        "bottom3_run_length": r.bottom3_run,
        "involved_three_prime_ends": list(r.involved_three_prime_ends),
        "risk_level": r.risk_level,
        "risk_reasons": list(r.risk_reasons),
        "alignment": {
            "top": r.top_line,
            "match": r.match_line,
            "bottom": r.bottom_line,
        },
    }


def rank_batch_replacements(
    *,
    batch_id: str,
    target_name: str,
    candidates: list[tuple[str, str]],
) -> dict | None:
    """对已保存批次做“替换某引物”的内存模拟并按多指标非支配分层。

    批次不存在返回 None（路由层转 404）；目标引物名不存在抛出
    ValueError（路由层转 422）。每个候选仅计算其自二聚体以及与批次中
    **其余**引物的异二聚体，复用 analyze_pair 与批次原阈值，不写库，
    也不重算与该候选无关的原批次配对。

    指标（越小越优）：high 配对数、最大 3′ 连续互补、最大连续互补。
    """
    with db.get_conn() as conn:
        batch = db.get_batch(conn, batch_id)
        if batch is None:
            return None
        primers = db.get_primers(conn, batch_id)

    cont_thr = batch["continuous_threshold"]
    three_thr = batch["three_prime_threshold"]

    target_idx = next(
        (i for i, row in enumerate(primers) if row["name"] == target_name),
        None,
    )
    if target_idx is None:
        raise ValueError(f"批次中不存在名为 {target_name!r} 的引物")

    others = [
        (i, row["name"], row["sequence"])
        for i, row in enumerate(primers)
        if i != target_idx
    ]

    evaluated: list[dict] = []
    for input_order, (label, seq) in enumerate(candidates):
        pairs: list[PairResult] = []
        # 自二聚体（候选标签即引物名）
        pairs.append(
            analyze_pair(
                label, seq, label, seq,
                cont_thr, three_thr,
                is_self_dimer=True,
            )
        )
        # 与其余引物的异二聚体：保持批次位置决定的 top/bottom 方向
        for other_idx, other_name, other_seq in others:
            if target_idx < other_idx:
                top_name, top_seq = label, seq
                bottom_name, bottom_seq = other_name, other_seq
            else:
                top_name, top_seq = other_name, other_seq
                bottom_name, bottom_seq = label, seq
            pairs.append(
                analyze_pair(
                    top_name, top_seq, bottom_name, bottom_seq,
                    cont_thr, three_thr,
                    is_self_dimer=False,
                )
            )

        high_count = sum(1 for r in pairs if r.risk_level == "high")
        max_terminal = max(r.terminal_run for r in pairs)
        max_cont = max(r.max_complementary_run for r in pairs)

        pair_dicts = [_pair_result_to_dict(i, r) for i, r in enumerate(pairs)]
        # 配对中“候选一侧”的引物名统一为标签
        for pdct, r in zip(pair_dicts, pairs):
            pdct["label"] = label
            pdct["other_name"] = (
                label
                if r.is_self_dimer
                else (r.name_bottom if r.name_top == label else r.name_top)
            )

        evaluated.append(
            {
                "input_order": input_order,
                "label": label,
                "sequence": seq,
                "evaluated_pair_count": len(pairs),
                "high_risk_count": high_count,
                "max_three_prime_run": max_terminal,
                "max_continuous_complementary": max_cont,
                "high_risk_pairs": [
                    pdct for pdct, r in zip(pair_dicts, pairs)
                    if r.risk_level == "high"
                ],
                "max_three_prime_pairs": [
                    pdct for pdct, r in zip(pair_dicts, pairs)
                    if r.terminal_run == max_terminal
                ],
                "max_continuous_pairs": [
                    pdct for pdct, r in zip(pair_dicts, pairs)
                    if r.max_complementary_run == max_cont
                ],
            }
        )

    # 非支配分层：三个指标越小越优
    vectors = [
        (
            item["high_risk_count"],
            item["max_three_prime_run"],
            item["max_continuous_complementary"],
        )
        for item in evaluated
    ]
    ranks = non_dominated_ranks(vectors)
    for item, rank in zip(evaluated, ranks):
        item["rank"] = rank

    # 同层按输入顺序稳定排序
    ranked = sorted(evaluated, key=lambda item: (item["rank"], item["input_order"]))

    return {
        "batch_id": batch_id,
        "target_name": target_name,
        "continuous_threshold": cont_thr,
        "three_prime_threshold": three_thr,
        "candidate_count": len(candidates),
        "ranked_candidates": ranked,
    }


def _diff_to_out(d: PairDiff) -> dict:
    return {
        "pair_key": d.pair_key,
        "name_a": d.name_a,
        "name_b": d.name_b,
        "is_self_dimer": d.is_self_dimer,
        "sequence_changed": d.sequence_changed,
        "base_risk_level": d.base_risk_level,
        "target_risk_level": d.target_risk_level,
    }


def compare_batches(*, base_batch_id: str, target_batch_id: str) -> dict | None:
    """对比两个已保存批次的配对风险（base → target）。

    按引物名称构造稳定配对键（自二聚体保留单名，异二聚体忽略批次内
    提交顺序），比较风险等级 low<medium<high，输出 escalated /
    deescalated / introduced / resolved 四类明细与汇总。任一批次
    不存在返回 None（由路由层转 404）。
    """
    with db.get_conn() as conn:
        base = db.get_batch(conn, base_batch_id)
        target = db.get_batch(conn, target_batch_id)
        if base is None or target is None:
            return None
        base_primers = db.get_primers(conn, base_batch_id)
        target_primers = db.get_primers(conn, target_batch_id)
        base_pairs = db.get_pairs(conn, base_batch_id, risky_only=False)
        target_pairs = db.get_pairs(conn, target_batch_id, risky_only=False)

    result = compare_pair_risks(
        base_pairs=[
            (r["name_top"], r["name_bottom"], bool(r["is_self_dimer"]), r["risk_level"])
            for r in base_pairs
        ],
        target_pairs=[
            (r["name_top"], r["name_bottom"], bool(r["is_self_dimer"]), r["risk_level"])
            for r in target_pairs
        ],
        base_sequences={r["name"]: r["sequence"] for r in base_primers},
        target_sequences={r["name"]: r["sequence"] for r in target_primers},
    )

    return {
        "summary": {
            "base_batch_id": base_batch_id,
            "target_batch_id": target_batch_id,
            "base_pair_count": len(base_pairs),
            "target_pair_count": len(target_pairs),
            "common_pair_count": result.common_pair_count,
            "unchanged_pair_count": result.unchanged_pair_count,
            "escalated_count": len(result.escalated),
            "deescalated_count": len(result.deescalated),
            "introduced_count": len(result.introduced),
            "resolved_count": len(result.resolved),
        },
        "escalated": [_diff_to_out(d) for d in result.escalated],
        "deescalated": [_diff_to_out(d) for d in result.deescalated],
        "introduced": [_diff_to_out(d) for d in result.introduced],
        "resolved": [_diff_to_out(d) for d in result.resolved],
    }

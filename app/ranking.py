"""多指标非支配分层（non-dominated sorting）。

每个候选由若干“越小越优”的指标向量描述。分层规则：

* 候选 A **支配**候选 B 当且仅当 A 在所有指标上不劣于 B（``<=``），
  且至少一个指标严格更优（``<``）；指标完全相同的候选互不支配；
* 第一层（rank=1）为当前未分层集合中不被任何其他候选支配者；
* 逐层剥离，rank 从 1 递增。

同层顺序由调用方按输入顺序稳定给出（本模块返回 rank，不做重排）。
候选数上限 50，直接采用 O(M²) 朴素枚举，结果确定可复现。
"""

from __future__ import annotations

from collections.abc import Sequence


def dominates(a: Sequence[int], b: Sequence[int]) -> bool:
    """指标向量 a 是否支配 b：逐维不劣且至少一维严格更优。"""
    return all(x <= y for x, y in zip(a, b)) and any(x < y for x, y in zip(a, b))


def non_dominated_ranks(values: Sequence[Sequence[int]]) -> list[int]:
    """返回每个候选的层级（1 起）。

    指标完全相同的候选落入同一层；层级分配与输入顺序无关，仅由指标决定。
    """
    remaining = set(range(len(values)))
    ranks = [0] * len(values)
    rank = 1
    while remaining:
        front = [
            i
            for i in remaining
            if not any(
                j != i and dominates(values[j], values[i])
                for j in remaining
            )
        ]
        for i in front:
            ranks[i] = rank
        remaining.difference_update(front)
        rank += 1
    return ranks

"""非支配分层算法单测：越小越优、同层并列、逐层剥离。"""

from __future__ import annotations

from app.ranking import dominates, non_dominated_ranks


def test_dominates():
    assert dominates((0, 1, 2), (1, 2, 3))
    # 至少一维严格更优
    assert dominates((0, 1, 2), (0, 1, 3))
    # 相等不构成支配
    assert not dominates((1, 2, 3), (1, 2, 3))
    # 各有胜负（tradeoff）互不支配
    assert not dominates((0, 2, 0), (1, 0, 1))
    assert not dominates((1, 0, 1), (0, 2, 0))


def test_ranks_single_and_two_fronts():
    assert non_dominated_ranks([(0, 0, 0), (1, 1, 1)]) == [1, 2]
    assert non_dominated_ranks([(2, 2, 2), (0, 0, 0), (1, 1, 1)]) == [3, 1, 2]


def test_equal_metrics_share_front():
    assert non_dominated_ranks([(1, 2, 3), (1, 2, 3)]) == [1, 1]
    # 两个等价最优 + 一个被支配者
    assert non_dominated_ranks([(1, 2, 3), (2, 3, 4), (1, 2, 3)]) == [1, 2, 1]


def test_tradeoff_non_comparable_share_front():
    # 三个候选在不同指标上占优，两两不可比 -> 同层
    assert non_dominated_ranks(
        [(0, 5, 5), (5, 0, 5), (5, 5, 0)]
    ) == [1, 1, 1]


def test_rank_depends_only_on_values_not_order():
    values = [(2, 2, 2), (0, 0, 0), (1, 1, 1)]
    ranks = non_dominated_ranks(values)
    # 打乱输入后，按指标查找的层级保持一致
    perm = [2, 0, 1]
    shuffled = [values[i] for i in perm]
    shuffled_ranks = non_dominated_ranks(shuffled)
    for idx, rank in zip(perm, shuffled_ranks):
        assert rank == ranks[idx]

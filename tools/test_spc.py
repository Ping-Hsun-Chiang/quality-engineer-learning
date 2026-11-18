"""用筆記中的手算範例驗證 spc.py。執行：python -m pytest tools -q"""

import numpy as np
import pytest

from spc import capability, grade_ca, grade_cpk, nelson_rules, xbar_r_limits


def test_xbar_r_limits_match_hand_calculation():
    # 兩組 n = 5，X̄ 分別為 10.02、10.04，R 皆為 0.10 → X̿ = 10.03，R̄ = 0.10
    data = np.array([
        [9.97, 10.07, 10.02, 10.02, 10.02],
        [9.99, 10.09, 10.04, 10.04, 10.04],
    ])
    lim = xbar_r_limits(data)
    assert lim.xbarbar == pytest.approx(10.03)
    assert lim.rbar == pytest.approx(0.10)
    assert lim.x_ucl == pytest.approx(10.0877)
    assert lim.x_lcl == pytest.approx(9.9723)
    assert lim.r_ucl == pytest.approx(0.2114)
    assert lim.r_lcl == 0
    assert lim.sigma_within == pytest.approx(0.10 / 2.326)


def test_capability_matches_hand_calculation():
    # 規格 10.00 ± 0.15，平均 10.03，σ = 0.10 / 2.326
    data = np.array([[10.02, 10.04]])
    cap = capability(data, lsl=9.85, usl=10.15, sigma_within=0.10 / 2.326)
    assert cap.ca == pytest.approx(0.20)
    assert cap.cp == pytest.approx(1.163, abs=1e-3)
    assert cap.cpk == pytest.approx(0.930, abs=1e-3)
    assert cap.cpk == pytest.approx(cap.cp * (1 - abs(cap.ca)))
    assert grade_ca(cap.ca) == "B"
    assert grade_cpk(cap.cpk) == "C"


def test_rule_1_point_beyond_3_sigma():
    x = np.array([0.0, 0.5, -0.5, 3.5, 0.0])
    assert 1 in nelson_rules(x, 0, 1)[3]


def test_rule_2_nine_points_same_side():
    x = np.array([0.5, 0.3, 0.6, 0.2, 0.4, 0.5, 0.3, 0.6, 0.4])
    hits = nelson_rules(x, 0, 1)
    assert 2 in hits[8]
    assert all(2 not in hits.get(i, []) for i in range(8))


def test_rule_3_six_points_trending():
    x = np.array([-1.0, -0.6, -0.2, 0.1, 0.5, 0.9])
    assert 3 in nelson_rules(x, 0, 1)[5]


def test_rule_7_fifteen_points_hugging_center():
    x = np.tile([0.1, -0.1], 8)[:15]
    assert 7 in nelson_rules(x, 0, 1)[14]


def test_stable_random_process_has_few_alarms():
    rng = np.random.default_rng(0)
    x = rng.normal(0, 1, 25)
    hits = nelson_rules(x, 0, 1)
    assert 1 not in sum(hits.values(), [])

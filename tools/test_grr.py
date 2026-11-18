"""用筆記中的手算範例驗證 grr.py。執行：python -m pytest tools -q"""

import numpy as np
import pytest

from grr import analyze, grr_from_stats, make_demo_data


def test_matches_hand_calculation():
    # 3 人 × 10 零件 × 3 次，R̿ = 0.020，X̄diff = 0.015，Rp = 0.40，公差 0.60
    res = grr_from_stats(rbarbar=0.020, xdiff=0.015, rp=0.40, n=10, r=3, k=3, tolerance=0.60)
    assert res.ev == pytest.approx(0.01182, abs=1e-5)
    assert res.av == pytest.approx(0.00754, abs=1e-5)
    assert res.grr == pytest.approx(0.01402, abs=1e-5)
    assert res.pv == pytest.approx(0.12584, abs=1e-5)
    assert res.tv == pytest.approx(0.12662, abs=1e-5)
    assert res.pct(res.grr) == pytest.approx(0.111, abs=1e-3)
    assert res.ndc == 12
    assert res.pct_tol == pytest.approx(0.140, abs=1e-3)
    assert res.verdict == "有條件接受"


def test_variances_add_in_squares_not_linearly():
    res = grr_from_stats(rbarbar=0.020, xdiff=0.015, rp=0.40, n=10, r=3, k=3)
    pct_ev, pct_av, pct_grr = res.pct(res.ev), res.pct(res.av), res.pct(res.grr)
    assert pct_ev + pct_av > pct_grr
    assert pct_grr == pytest.approx(np.hypot(pct_ev, pct_av))


def test_av_is_zero_when_negative_under_root():
    # 量測人員之間幾乎沒有差異時，根號內為負，AV 取 0
    res = grr_from_stats(rbarbar=0.05, xdiff=0.001, rp=0.40, n=10, r=3, k=3)
    assert res.av == 0
    assert res.grr == pytest.approx(res.ev)


def test_perfect_gauge_has_zero_grr():
    true = np.linspace(9.8, 10.2, 10)
    data = np.broadcast_to(true[None, :, None], (3, 10, 3)).copy()
    res, _ = analyze(data)
    assert res.ev == 0 and res.av == 0 and res.grr == 0


def test_raw_data_summary_matches_definition():
    data, _ = make_demo_data()
    res, s = analyze(data)
    ranges = data.max(axis=2) - data.min(axis=2)
    assert s["rbarbar"] == pytest.approx(ranges.mean())
    assert s["rp"] == pytest.approx(np.ptp(data.mean(axis=(0, 2))))
    assert res.pct(res.tv) == pytest.approx(1.0)

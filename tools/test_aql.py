"""用筆記中的查表範例驗證 aql.py。執行：python -m pytest tools -q"""

import pytest

from aql import accept_probability, code_letter, defect_rate_at, lookup


@pytest.mark.parametrize("lot, level, expected", [
    (8, "II", "A"), (9, "II", "B"), (500, "II", "H"), (501, "II", "J"),
    (2000, "II", "K"), (2000, "I", "H"), (2000, "III", "L"), (2000, "S-4", "G"),
    (600000, "II", "Q"), (600000, "III", "R"),
])
def test_code_letter(lot, level, expected):
    assert code_letter(lot, level) == expected


@pytest.mark.parametrize("aql, n, ac, re, used", [
    ("0.65", 125, 2, 3, "K"),   # 筆記範例：主要缺點
    ("1.5", 125, 5, 6, "K"),    # 筆記範例：次要缺點
    ("0.25", 200, 1, 2, "L"),   # ↓ 往下到 L
    ("0.15", 80, 0, 1, "J"),    # ↑ 往上到 J
    ("6.5", 125, 14, 15, "K"),
])
def test_lookup_lot_2000(aql, n, ac, re, used):
    plan = lookup(2000, aql)
    assert (plan.n, plan.ac, plan.re, plan.used_code) == (n, ac, re, used)


def test_up_arrow_in_lower_right_of_table():
    # 代字 R、AQL 1.0 落在表格右下的 ↑ 區域，應使用 Q 列的 21/22
    plan = lookup(600000, "1.0", "III")
    assert (plan.code, plan.used_code, plan.n, plan.ac, plan.re) == ("R", "Q", 1250, 21, 22)


def test_full_inspection_when_sample_exceeds_lot():
    plan = lookup(10, "0.10")
    assert plan.full_inspection
    assert plan.n == 10 and plan.ac == 0


def test_unsupported_aql_raises():
    with pytest.raises(ValueError):
        lookup(2000, "0.33")


def test_oc_curve_values_from_notes():
    assert float(accept_probability(125, 2, 0.0065)) == pytest.approx(0.951, abs=1e-3)
    assert float(accept_probability(125, 2, 0.02)) == pytest.approx(0.543, abs=1e-3)
    assert defect_rate_at(125, 2, 0.10) == pytest.approx(0.042, abs=1e-3)


def test_c0_plan_accepts_only_zero_defects():
    assert float(accept_probability(50, 0, 0.0)) == 1.0
    assert float(accept_probability(50, 0, 0.01)) == pytest.approx(0.99**50)

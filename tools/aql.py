"""AQL 抽樣計畫查詢（ISO 2859-1 正常檢驗單次抽樣）與 OC 曲線。

用法：
    python tools/aql.py --demo
    python tools/aql.py --lot 2000 --aql 0.65 1.5
    python tools/aql.py --lot 2000 --aql 0.65 --level I --defects 2
    python tools/aql.py --lot 2000 --aql 0.65 --compare 50,0 200,1
"""

from __future__ import annotations

import argparse
import sys
from dataclasses import dataclass
from math import comb
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

LETTERS = "ABCDEFGHJKLMNPQR"
SAMPLE_SIZE = dict(zip(LETTERS, [2, 3, 5, 8, 13, 20, 32, 50, 80, 125, 200, 315, 500, 800, 1250, 2000]))

# 表一：樣本代字（批量上限, S-1, S-2, S-3, S-4, I, II, III）
LEVELS = ["S-1", "S-2", "S-3", "S-4", "I", "II", "III"]
CODE_TABLE = [
    (8, "AAAAAAB"),
    (15, "AAAAABC"),
    (25, "AABBBCD"),
    (50, "ABBCCDE"),
    (90, "BBCCCEF"),
    (150, "BBCDDFG"),
    (280, "BCDEEGH"),
    (500, "BCDEFHJ"),
    (1200, "CCEFGJK"),
    (3200, "CDEGHKL"),
    (10000, "CDFGJLM"),
    (35000, "CDFHKMN"),
    (150000, "DEGJLNP"),
    (500000, "DEGJMPQ"),
    (float("inf"), "DEHKNQR"),
]

# 表二的優先 AQL 數列；Ac/Re 沿對角線（代字序號 + AQL 序號）排列
AQL_SERIES = ["0.010", "0.015", "0.025", "0.040", "0.065", "0.10", "0.15", "0.25",
              "0.40", "0.65", "1.0", "1.5", "2.5", "4.0", "6.5"]
SUPPORTED_AQL = AQL_SERIES[2:]  # 0.025 ～ 6.5
DIAGONAL = {14: (0, 1), 17: (1, 2), 18: (2, 3), 19: (3, 4), 20: (5, 6),
            21: (7, 8), 22: (10, 11), 23: (14, 15), 24: (21, 22)}


@dataclass
class Plan:
    lot: int
    level: str
    aql: str
    code: str  # 由批量查到的代字
    used_code: str  # 依箭頭調整後實際使用的代字
    n: int
    ac: int
    re: int
    full_inspection: bool

    @property
    def arrow(self) -> str:
        if self.used_code == self.code:
            return ""
        return "↓" if LETTERS.index(self.used_code) > LETTERS.index(self.code) else "↑"


def code_letter(lot: int, level: str = "II") -> str:
    if lot < 2:
        raise ValueError("批量至少為 2")
    col = LEVELS.index(level)
    for upper, codes in CODE_TABLE:
        if lot <= upper:
            return codes[col]
    raise AssertionError("unreachable")


def normalize_aql(aql: str | float) -> str:
    value = float(aql)
    for s in AQL_SERIES:
        if abs(float(s) - value) < 1e-9:
            if s not in SUPPORTED_AQL:
                break
            return s
    raise ValueError(f"AQL {aql} 不在支援範圍：{', '.join(SUPPORTED_AQL)}")


def lookup(lot: int, aql: str | float, level: str = "II") -> Plan:
    """查詢正常檢驗單次抽樣計畫，依箭頭規則調整樣本數。"""
    aql = normalize_aql(aql)
    code = code_letter(lot, level)
    c, a = LETTERS.index(code), AQL_SERIES.index(aql)
    d = c + a
    if d < 14 or d == 16:  # ↓：往下到第一個有數字的格子
        c = 14 - a if d < 14 else c + 1
    elif d == 15 or d > 24:  # ↑：往上到第一個有數字的格子
        c = c - 1 if d == 15 else 24 - a
    used = LETTERS[c]
    ac, re = DIAGONAL[c + a]
    n = SAMPLE_SIZE[used]
    full = n >= lot
    if full:
        n = lot
    return Plan(lot, level, aql, code, used, n, ac, re, full)


def accept_probability(n: int, ac: int, p: float | np.ndarray) -> float | np.ndarray:
    """二項分布：樣本 n 中不良數 ≤ Ac 的機率。"""
    p = np.asarray(p, dtype=float)
    return sum(comb(n, d) * p**d * (1 - p) ** (n - d) for d in range(ac + 1))


def defect_rate_at(n: int, ac: int, pa: float) -> float:
    """找出允收機率等於 pa 的批不良率（二分法）。"""
    lo, hi = 0.0, 1.0
    for _ in range(60):
        mid = (lo + hi) / 2
        if accept_probability(n, ac, mid) > pa:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def describe(plan: Plan) -> str:
    lines = [f"批量 {plan.lot:,}，檢驗水準 {plan.level}，AQL {plan.aql}"]
    if plan.arrow:
        lines.append(f"  代字 {plan.code} 的格子為「{plan.arrow}」→ 改用代字 {plan.used_code}")
    else:
        lines.append(f"  代字 {plan.code}")
    if plan.full_inspection:
        lines.append(f"  樣本數 ≥ 批量 → 全檢 {plan.n} 個，Ac {plan.ac}、Re {plan.re}")
    else:
        lines.append(f"  樣本數 n = {plan.n}，Ac {plan.ac}、Re {plan.re}")
        p95 = defect_rate_at(plan.n, plan.ac, 0.95)
        p10 = defect_rate_at(plan.n, plan.ac, 0.10)
        pa_aql = float(accept_probability(plan.n, plan.ac, float(plan.aql) / 100))
        lines.append(f"  不良率 = AQL 時允收機率 {pa_aql:.1%}；允收 95% 的不良率 {p95:.2%}；"
                     f"允收 10% 的不良率（LTPD）{p10:.2%}")
    return "\n".join(lines)


def plot_oc(curves: list[tuple[str, int, int]], aql_pct: float | None, out: Path) -> None:
    plt.rcParams["font.sans-serif"] = ["Microsoft JhengHei", "PingFang TC", "Noto Sans CJK TC", "DejaVu Sans"]
    plt.rcParams["axes.unicode_minus"] = False
    colors = ["#264653", "#e76f51", "#2a9d8f", "#e9c46a", "#8ecae6"]

    p_max = max(defect_rate_at(n, ac, 0.02) for _, n, ac in curves)
    p = np.linspace(0, min(p_max * 1.1, 1), 400)

    fig, ax = plt.subplots(figsize=(10, 6))
    for i, ((label, n, ac), color) in enumerate(zip(curves, colors)):
        ax.plot(p * 100, accept_probability(n, ac, p) * 100, color=color, lw=2,
                label=f"{label}（n = {n}，Ac = {ac}）")
        p10 = defect_rate_at(n, ac, 0.10)
        ax.plot(p10 * 100, 10, "o", color=color, ms=6)
        ax.annotate(f"LTPD {p10:.1%}", (p10 * 100, 10), textcoords="offset points",
                    xytext=(6, 8 + 14 * i), fontsize=8, color=color)
    if aql_pct is not None:
        ax.axvline(aql_pct, color="#999", ls=":", lw=1.2)
        ax.text(aql_pct, 102, f"AQL {aql_pct}%", ha="center", fontsize=9, color="#666")
        for (_, n, ac), color in zip(curves, colors):
            pa = float(accept_probability(n, ac, aql_pct / 100)) * 100
            ax.plot(aql_pct, pa, "s", color=color, ms=6)
            ax.annotate(f"{pa:.1f}%", (aql_pct, pa), textcoords="offset points",
                        xytext=(6, -12), fontsize=8, color=color)
    ax.axhline(95, color="#ccc", lw=0.8, ls="--")
    ax.axhline(10, color="#ccc", lw=0.8, ls="--")
    ax.text(p[-1] * 100, 95, "Pa = 95%（α = 5%）", ha="right", va="bottom", fontsize=8, color="#888")
    ax.text(p[-1] * 100, 10, "Pa = 10%（β = 10%）", ha="right", va="bottom", fontsize=8, color="#888")
    ax.set_xlim(0, p[-1] * 100)
    ax.set_ylim(0, 108)
    ax.set_xlabel("批不良率 p（%）")
    ax.set_ylabel("允收機率 Pa（%）")
    ax.set_title("OC 曲線：不同抽樣計畫的鑑別力比較")
    ax.legend(loc="center right", bbox_to_anchor=(1, 0.6), fontsize=9)
    ax.grid(alpha=0.25)

    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=120, bbox_inches="tight")
    plt.close(fig)


def parse_pair(text: str) -> tuple[int, int]:
    n, ac = text.split(",")
    return int(n), int(ac)


def main() -> None:
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")

    parser = argparse.ArgumentParser(description="AQL 抽樣計畫查詢與 OC 曲線")
    parser.add_argument("--lot", type=int, help="批量 N")
    parser.add_argument("--aql", nargs="+", help="一個或多個 AQL，例如 0.65 1.5")
    parser.add_argument("--level", default="II", choices=LEVELS, help="檢驗水準（預設 II）")
    parser.add_argument("--defects", type=int, nargs="+", help="樣本中的不良數（依 --aql 順序），用來判定允收或拒收")
    parser.add_argument("--compare", nargs="*", default=[], help="額外比較的抽樣計畫，格式 n,Ac，例如 50,0")
    parser.add_argument("--out", type=Path, default=Path("tools/output/oc_curve.png"), help="輸出圖檔路徑")
    parser.add_argument("--demo", action="store_true", help="批量 2,000、AQL 0.65，並與 C=0 及其他計畫比較")
    args = parser.parse_args()

    if args.demo:
        args.lot, args.aql, args.compare = 2000, ["0.65", "1.5"], ["50,0", "200,1"]
    if not (args.lot and args.aql):
        parser.error("請提供 --lot 與 --aql，或使用 --demo")

    plans = [lookup(args.lot, a, args.level) for a in args.aql]
    for plan in plans:
        print(describe(plan))

    if args.defects:
        if len(args.defects) != len(plans):
            parser.error("--defects 的數量要和 --aql 相同")
        print("\n判定")
        verdicts = []
        for plan, d in zip(plans, args.defects):
            ok = d <= plan.ac
            verdicts.append(ok)
            print(f"  AQL {plan.aql}：不良 {d} 個（Ac {plan.ac}、Re {plan.re}）→ {'允收' if ok else '拒收'}")
        print(f"  整批：{'允收' if all(verdicts) else '拒收'}")

    curves = [(f"查表結果 AQL {p.aql}", p.n, p.ac) for p in plans[:1] if not p.full_inspection]
    curves += [("C=0 抽樣" if ac == 0 else "對照計畫", n, ac) for n, ac in map(parse_pair, args.compare)]
    if curves:
        plot_oc(curves, float(plans[0].aql), args.out)
        print(f"\nOC 曲線已輸出：{args.out}")


if __name__ == "__main__":
    main()

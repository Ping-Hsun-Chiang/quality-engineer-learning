"""X̄-R 管制圖、Nelson 判異規則與製程能力分析。

用法：
    python tools/spc.py --demo
    python tools/spc.py --csv tools/data/shaft_od.csv --lsl 9.85 --usl 10.15 --baseline 25

CSV 格式：每一列是一組（子群），每一欄是一個量測值，第一列為欄位名稱。
"""

from __future__ import annotations

import argparse
import sys
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# 管制圖係數（n = 2～10）
CONSTANTS = {
    #  n: (A2,    D3,    D4,    d2)
    2: (1.880, 0.000, 3.267, 1.128),
    3: (1.023, 0.000, 2.574, 1.693),
    4: (0.729, 0.000, 2.282, 2.059),
    5: (0.577, 0.000, 2.114, 2.326),
    6: (0.483, 0.000, 2.004, 2.534),
    7: (0.419, 0.076, 1.924, 2.704),
    8: (0.373, 0.136, 1.864, 2.847),
    9: (0.337, 0.184, 1.816, 2.970),
    10: (0.308, 0.223, 1.777, 3.078),
}

RULE_TEXT = {
    1: "1 點超出 3σ",
    2: "連續 9 點在中心線同一側",
    3: "連續 6 點持續上升或下降",
    4: "連續 14 點上下交替",
    5: "連續 3 點中有 2 點在同側 2σ 以外",
    6: "連續 5 點中有 4 點在同側 1σ 以外",
    7: "連續 15 點在 1σ 以內",
    8: "連續 8 點在 1σ 以外（兩側皆有）",
}


@dataclass
class XbarRLimits:
    n: int
    xbarbar: float
    rbar: float
    x_ucl: float
    x_lcl: float
    r_ucl: float
    r_lcl: float
    sigma_within: float  # 組內標準差 R̄ / d2

    @property
    def sigma_xbar(self) -> float:
        """X̄ 的標準差，用來劃分 A/B/C 區。"""
        return (self.x_ucl - self.xbarbar) / 3


@dataclass
class Capability:
    mean: float
    sigma_within: float
    sigma_overall: float
    ca: float
    cp: float
    cpk: float
    pp: float
    ppk: float


def xbar_r_limits(data: np.ndarray) -> XbarRLimits:
    """由子群數據（列 = 組，欄 = 量測值）計算 X̄-R 管制界限。"""
    n = data.shape[1]
    if n not in CONSTANTS:
        raise ValueError(f"每組樣本數 n = {n}，僅支援 2～10")
    a2, d3, d4, d2 = CONSTANTS[n]
    xbar = data.mean(axis=1)
    r = data.max(axis=1) - data.min(axis=1)
    xbarbar, rbar = xbar.mean(), r.mean()
    return XbarRLimits(
        n=n,
        xbarbar=xbarbar,
        rbar=rbar,
        x_ucl=xbarbar + a2 * rbar,
        x_lcl=xbarbar - a2 * rbar,
        r_ucl=d4 * rbar,
        r_lcl=d3 * rbar,
        sigma_within=rbar / d2,
    )


def nelson_rules(x: np.ndarray, center: float, sigma: float) -> dict[int, list[int]]:
    """回傳 {點的索引: [觸犯的規則編號]}，點的索引為觸發規則的最後一點。"""
    z = (np.asarray(x) - center) / sigma
    hits: dict[int, list[int]] = {}

    def flag(i: int, rule: int) -> None:
        hits.setdefault(i, [])
        if rule not in hits[i]:
            hits[i].append(rule)

    for i in range(len(z)):
        # 規則 1：1 點超出 3σ
        if abs(z[i]) > 3:
            flag(i, 1)
        # 規則 2：連續 9 點在同一側
        if i >= 8:
            w = z[i - 8 : i + 1]
            if np.all(w > 0) or np.all(w < 0):
                flag(i, 2)
        # 規則 3：連續 6 點持續上升或下降
        if i >= 5:
            d = np.diff(x[i - 5 : i + 1])
            if np.all(d > 0) or np.all(d < 0):
                flag(i, 3)
        # 規則 4：連續 14 點上下交替
        if i >= 13:
            d = np.sign(np.diff(x[i - 13 : i + 1]))
            if np.all(d != 0) and np.all(d[1:] * d[:-1] < 0):
                flag(i, 4)
        # 規則 5：3 點中有 2 點在同側 2σ 以外
        if i >= 2:
            w = z[i - 2 : i + 1]
            if (w > 2).sum() >= 2 or (w < -2).sum() >= 2:
                flag(i, 5)
        # 規則 6：5 點中有 4 點在同側 1σ 以外
        if i >= 4:
            w = z[i - 4 : i + 1]
            if (w > 1).sum() >= 4 or (w < -1).sum() >= 4:
                flag(i, 6)
        # 規則 7：連續 15 點在 1σ 以內
        if i >= 14 and np.all(np.abs(z[i - 14 : i + 1]) < 1):
            flag(i, 7)
        # 規則 8：連續 8 點在 1σ 以外，且兩側皆有
        if i >= 7:
            w = z[i - 7 : i + 1]
            if np.all(np.abs(w) > 1) and (w > 0).any() and (w < 0).any():
                flag(i, 8)
    return hits


def capability(data: np.ndarray, lsl: float, usl: float, sigma_within: float) -> Capability:
    """計算 Ca、Cp、Cpk（組內 σ）與 Pp、Ppk（整體 σ）。"""
    values = data.ravel()
    mean = values.mean()
    sigma_overall = values.std(ddof=1)
    m, t = (usl + lsl) / 2, usl - lsl

    def cp_cpk(sigma: float) -> tuple[float, float]:
        cp = t / (6 * sigma)
        cpk = min(usl - mean, mean - lsl) / (3 * sigma)
        return cp, cpk

    cp, cpk = cp_cpk(sigma_within)
    pp, ppk = cp_cpk(sigma_overall)
    return Capability(
        mean=mean,
        sigma_within=sigma_within,
        sigma_overall=sigma_overall,
        ca=(mean - m) / (t / 2),
        cp=cp,
        cpk=cpk,
        pp=pp,
        ppk=ppk,
    )


def grade_cpk(cpk: float) -> str:
    for threshold, grade in [(1.67, "A+"), (1.33, "A"), (1.00, "B"), (0.67, "C")]:
        if cpk >= threshold:
            return grade
    return "D"


def grade_ca(ca: float) -> str:
    for threshold, grade in [(0.125, "A"), (0.25, "B"), (0.50, "C")]:
        if abs(ca) <= threshold:
            return grade
    return "D"


def make_demo_data(seed: int = 7) -> np.ndarray:
    """模擬軸外徑（規格 10.00 ± 0.15）：前 25 組穩定，後 10 組因刀具磨損逐漸變大。"""
    rng = np.random.default_rng(seed)
    stable = rng.normal(10.03, 0.043, size=(25, 5))
    drift = 10.03 + 0.012 * np.arange(1, 11)[:, None]
    worn = rng.normal(drift, 0.043, size=(10, 5))
    return np.round(np.vstack([stable, worn]), 3)


def plot(
    data: np.ndarray,
    lim: XbarRLimits,
    hits: dict[int, list[int]],
    cap: Capability,
    lsl: float,
    usl: float,
    baseline: int,
    out: Path,
) -> None:
    plt.rcParams["font.sans-serif"] = ["Microsoft JhengHei", "PingFang TC", "Noto Sans CJK TC", "DejaVu Sans"]
    plt.rcParams["axes.unicode_minus"] = False

    xbar = data.mean(axis=1)
    r = data.max(axis=1) - data.min(axis=1)
    idx = np.arange(1, len(xbar) + 1)

    fig = plt.figure(figsize=(13, 8))
    gs = fig.add_gridspec(2, 3, width_ratios=[2.2, 2.2, 1.6], hspace=0.35, wspace=0.3)
    ax_x = fig.add_subplot(gs[0, :2])
    ax_r = fig.add_subplot(gs[1, :2], sharex=ax_x)
    ax_h = fig.add_subplot(gs[:, 2])

    # X̄ 圖
    s = lim.sigma_xbar
    for k in (1, 2):
        ax_x.axhline(lim.xbarbar + k * s, color="#bbb", lw=0.8, ls=":")
        ax_x.axhline(lim.xbarbar - k * s, color="#bbb", lw=0.8, ls=":")
    ax_x.axhline(lim.xbarbar, color="#2a9d8f", lw=1.2, label=f"CL = {lim.xbarbar:.4f}")
    ax_x.axhline(lim.x_ucl, color="#e76f51", lw=1.2, ls="--", label=f"UCL = {lim.x_ucl:.4f}")
    ax_x.axhline(lim.x_lcl, color="#e76f51", lw=1.2, ls="--", label=f"LCL = {lim.x_lcl:.4f}")
    ax_x.plot(idx, xbar, marker="o", ms=4, color="#264653", lw=1)
    for i, rules in hits.items():
        ax_x.plot(idx[i], xbar[i], "o", ms=8, mfc="none", mec="#d62828", mew=1.8)
        ax_x.annotate(",".join(map(str, sorted(rules))), (idx[i], xbar[i]),
                      textcoords="offset points", xytext=(0, 8), ha="center", fontsize=8, color="#d62828")
    if baseline < len(xbar):
        for ax in (ax_x, ax_r):
            ax.axvline(baseline + 0.5, color="#888", lw=1)
        ax_x.text(baseline + 0.7, lim.x_lcl, "監控階段 →", va="bottom", fontsize=9, color="#666")
        ax_x.text(baseline + 0.3, lim.x_lcl, "← 基準階段", va="bottom", ha="right", fontsize=9, color="#666")
    ax_x.set_title(r"$\bar{X}$ 管制圖（紅圈 = 異常，數字 = 觸犯的規則）")
    ax_x.legend(loc="upper left", fontsize=8)

    # R 圖
    ax_r.axhline(lim.rbar, color="#2a9d8f", lw=1.2, label=f"CL = {lim.rbar:.4f}")
    ax_r.axhline(lim.r_ucl, color="#e76f51", lw=1.2, ls="--", label=f"UCL = {lim.r_ucl:.4f}")
    ax_r.axhline(lim.r_lcl, color="#e76f51", lw=1.2, ls="--", label=f"LCL = {lim.r_lcl:.4f}")
    ax_r.plot(idx, r, marker="o", ms=4, color="#264653", lw=1)
    out_r = (r > lim.r_ucl) | (r < lim.r_lcl)
    ax_r.plot(idx[out_r], r[out_r], "o", ms=8, mfc="none", mec="#d62828", mew=1.8)
    ax_r.set_title("R 管制圖")
    ax_r.set_xlabel("組別")
    ax_r.legend(loc="upper left", fontsize=8)

    # 直方圖（基準期數據）
    values = data[:baseline].ravel()
    ax_h.hist(values, bins=12, color="#8ecae6", edgecolor="white", density=True)
    xs = np.linspace(min(values.min(), lsl) - 0.02, max(values.max(), usl) + 0.02, 300)
    pdf = np.exp(-0.5 * ((xs - cap.mean) / cap.sigma_within) ** 2) / (cap.sigma_within * np.sqrt(2 * np.pi))
    ax_h.plot(xs, pdf, color="#264653", lw=1.2)
    for v, name in [(lsl, "LSL"), (usl, "USL")]:
        ax_h.axvline(v, color="#d62828", lw=1.5)
        ax_h.text(v, ax_h.get_ylim()[1] * 0.95, name, ha="center", color="#d62828", fontsize=9)
    ax_h.axvline((lsl + usl) / 2, color="#999", lw=1, ls=":")
    ax_h.set_title(f"直方圖（基準期 {baseline} 組）")
    ax_h.set_yticks([])
    text = (
        f"平均 = {cap.mean:.4f}\n"
        f"σ組內 = {cap.sigma_within:.4f}\n"
        f"σ整體 = {cap.sigma_overall:.4f}\n"
        f"Ca  = {cap.ca:+.1%}（{grade_ca(cap.ca)}）\n"
        f"Cp  = {cap.cp:.2f}\n"
        f"Cpk = {cap.cpk:.2f}（{grade_cpk(cap.cpk)}）\n"
        f"Pp  = {cap.pp:.2f}\n"
        f"Ppk = {cap.ppk:.2f}"
    )
    ax_h.text(0.5, -0.07, text, transform=ax_h.transAxes, va="top", ha="center", fontsize=10,
              linespacing=1.5, bbox=dict(boxstyle="round", fc="white", ec="#ccc"))

    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=120, bbox_inches="tight")
    plt.close(fig)


def report(lim: XbarRLimits, hits: dict[int, list[int]], cap: Capability, baseline: int, total: int) -> str:
    lines = [
        f"管制界限（以第 1～{baseline} 組計算，n = {lim.n}）",
        f"  X̄ 圖：CL = {lim.xbarbar:.4f}  UCL = {lim.x_ucl:.4f}  LCL = {lim.x_lcl:.4f}",
        f"  R 圖 ：CL = {lim.rbar:.4f}  UCL = {lim.r_ucl:.4f}  LCL = {lim.r_lcl:.4f}",
        "",
        f"製程能力（基準期）",
        f"  Ca = {cap.ca:+.1%}（{grade_ca(cap.ca)}）  Cp = {cap.cp:.2f}  Cpk = {cap.cpk:.2f}（{grade_cpk(cap.cpk)}）"
        f"  Pp = {cap.pp:.2f}  Ppk = {cap.ppk:.2f}",
        "",
        f"判異結果（共 {total} 組）",
    ]
    if not hits:
        lines.append("  未發現異常")
    for i in sorted(hits):
        rules = "；".join(f"規則 {k}：{RULE_TEXT[k]}" for k in sorted(hits[i]))
        lines.append(f"  第 {i + 1:>2} 組  {rules}")
    return "\n".join(lines)


def main() -> None:
    # Windows 終端機預設 cp950，無法顯示 X̄ 等符號
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")

    parser = argparse.ArgumentParser(description="X̄-R 管制圖與製程能力分析")
    parser.add_argument("--csv", type=Path, help="子群數據 CSV（每列一組）")
    parser.add_argument("--lsl", type=float, help="規格下限")
    parser.add_argument("--usl", type=float, help="規格上限")
    parser.add_argument("--baseline", type=int, help="用前幾組計算管制界限（預設：全部）")
    parser.add_argument("--out", type=Path, default=Path("tools/output/xbar_r.png"), help="輸出圖檔路徑")
    parser.add_argument("--demo", action="store_true", help="使用內建的軸外徑範例數據")
    args = parser.parse_args()

    if args.demo:
        data = make_demo_data()
        lsl, usl, baseline = 9.85, 10.15, 25
    else:
        if not (args.csv and args.lsl is not None and args.usl is not None):
            parser.error("請提供 --csv、--lsl、--usl，或使用 --demo")
        data = np.loadtxt(args.csv, delimiter=",", skiprows=1)
        lsl, usl = args.lsl, args.usl
        baseline = args.baseline or len(data)

    lim = xbar_r_limits(data[:baseline])
    xbar = data.mean(axis=1)
    hits = nelson_rules(xbar, lim.xbarbar, lim.sigma_xbar)
    cap = capability(data[:baseline], lsl, usl, lim.sigma_within)

    print(report(lim, hits, cap, baseline, len(data)))
    plot(data, lim, hits, cap, lsl, usl, baseline, args.out)
    print(f"\n圖檔已輸出：{args.out}")


if __name__ == "__main__":
    main()

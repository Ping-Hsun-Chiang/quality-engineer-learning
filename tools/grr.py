"""GR&R 量具重複性與再現性分析（AIAG 平均值與全距法）。

用法：
    python tools/grr.py --demo
    python tools/grr.py --csv tools/data/grr_od.csv --tol 0.60

CSV 格式（長表）：appraiser,part,trial,value，第一列為欄位名稱。
"""

from __future__ import annotations

import argparse
import csv
import math
import sys
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# AIAG MSA 第四版係數
K1 = {2: 0.8862, 3: 0.5908}  # 依重複次數 r
K2 = {2: 0.7071, 3: 0.5231}  # 依量測人數 k
K3 = {2: 0.7071, 3: 0.5231, 4: 0.4467, 5: 0.4030, 6: 0.3742,
      7: 0.3534, 8: 0.3375, 9: 0.3249, 10: 0.3146}  # 依零件數 n
D4 = {2: 3.267, 3: 2.574}  # R 圖上限係數（子群大小 = r）
A2 = {2: 1.880, 3: 1.023}  # X̄ 圖係數（子群大小 = r）


@dataclass
class GRRResult:
    ev: float
    av: float
    grr: float
    pv: float
    tv: float
    ndc: int
    tolerance: float | None = None

    def pct(self, value: float) -> float:
        return value / self.tv

    @property
    def pct_tol(self) -> float | None:
        return None if self.tolerance is None else 6 * self.grr / self.tolerance

    @property
    def verdict(self) -> str:
        p = self.pct(self.grr)
        if p < 0.10:
            return "可接受"
        if p <= 0.30:
            return "有條件接受"
        return "不可接受"


def grr_from_stats(rbarbar: float, xdiff: float, rp: float, n: int, r: int, k: int,
                   tolerance: float | None = None) -> GRRResult:
    """由 R̿、X̄diff、Rp 計算 GR&R。"""
    ev = rbarbar * K1[r]
    av_sq = (xdiff * K2[k]) ** 2 - ev**2 / (n * r)
    av = math.sqrt(av_sq) if av_sq > 0 else 0.0
    grr = math.hypot(ev, av)
    pv = rp * K3[n]
    tv = math.hypot(grr, pv)
    ndc = int(1.41 * pv / grr) if grr > 0 else 0
    return GRRResult(ev, av, grr, pv, tv, ndc, tolerance)


def summarize(data: np.ndarray) -> dict:
    """data 形狀為 (k 人, n 零件, r 次)。"""
    ranges = data.max(axis=2) - data.min(axis=2)  # (k, n)
    cell_means = data.mean(axis=2)  # (k, n)
    appraiser_means = data.mean(axis=(1, 2))
    part_means = data.mean(axis=(0, 2))
    return {
        "ranges": ranges,
        "cell_means": cell_means,
        "rbar_by_appraiser": ranges.mean(axis=1),
        "rbarbar": ranges.mean(),
        "appraiser_means": appraiser_means,
        "xdiff": appraiser_means.max() - appraiser_means.min(),
        "part_means": part_means,
        "rp": part_means.max() - part_means.min(),
        "xbarbar": data.mean(),
    }


def analyze(data: np.ndarray, tolerance: float | None = None) -> tuple[GRRResult, dict]:
    k, n, r = data.shape
    s = summarize(data)
    return grr_from_stats(s["rbarbar"], s["xdiff"], s["rp"], n, r, k, tolerance), s


def load_csv(path: Path) -> tuple[np.ndarray, list[str]]:
    with open(path, newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    appraisers = sorted({row["appraiser"] for row in rows})
    parts = sorted({int(row["part"]) for row in rows})
    trials = sorted({int(row["trial"]) for row in rows})
    data = np.full((len(appraisers), len(parts), len(trials)), np.nan)
    for row in rows:
        data[appraisers.index(row["appraiser"]), parts.index(int(row["part"])),
             trials.index(int(row["trial"]))] = float(row["value"])
    if np.isnan(data).any():
        raise ValueError("數據不完整：每位量測人員都要量測每個零件相同次數")
    return data, appraisers


def make_demo_data(seed: int = 3) -> tuple[np.ndarray, list[str]]:
    """10 個零件（外徑約 9.8～10.2），3 位量測人員各量 3 次。

    C 員量法偏小、且手法較不穩定。
    """
    rng = np.random.default_rng(seed)
    true = np.array([9.82, 9.90, 9.95, 9.98, 10.00, 10.02, 10.05, 10.09, 10.13, 10.19])
    bias = np.array([0.000, 0.004, -0.008])
    noise = np.array([0.008, 0.008, 0.013])
    data = true[None, :, None] + bias[:, None, None] + rng.normal(0, 1, (3, 10, 3)) * noise[:, None, None]
    return np.round(data, 3), ["A", "B", "C"]


def save_csv(data: np.ndarray, appraisers: list[str], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["appraiser", "part", "trial", "value"])
        for i, a in enumerate(appraisers):
            for j in range(data.shape[1]):
                for t in range(data.shape[2]):
                    w.writerow([a, j + 1, t + 1, f"{data[i, j, t]:.3f}"])


def report(res: GRRResult, s: dict, appraisers: list[str], shape: tuple[int, int, int]) -> str:
    k, n, r = shape
    r_ucl = D4[r] * s["rbarbar"]
    x_half = A2[r] * s["rbarbar"]
    outside = np.abs(s["cell_means"] - s["xbarbar"]) > x_half
    lines = [
        f"研究設定：{k} 位量測人員 × {n} 個零件 × {r} 次",
        f"  R̿ = {s['rbarbar']:.5f}   X̄diff = {s['xdiff']:.5f}   Rp = {s['rp']:.5f}",
        "  各量測人員：" + "  ".join(
            f"{a}（平均 {m:.4f}，R̄ {rb:.4f}）"
            for a, m, rb in zip(appraisers, s["appraiser_means"], s["rbar_by_appraiser"])
        ),
        "",
        f"{'來源':<8}{'標準差':>10}{'%總變異':>10}",
        f"{'EV 重複性':<6}{res.ev:>11.5f}{res.pct(res.ev):>10.1%}",
        f"{'AV 再現性':<6}{res.av:>11.5f}{res.pct(res.av):>10.1%}",
        f"{'GRR':<8}{res.grr:>11.5f}{res.pct(res.grr):>10.1%}",
        f"{'PV 零件':<7}{res.pv:>11.5f}{res.pct(res.pv):>10.1%}",
        f"{'TV 總變異':<6}{res.tv:>11.5f}",
        "",
        f"%GRR = {res.pct(res.grr):.1%} → {res.verdict}",
        f"ndc  = {res.ndc} → {'可接受' if res.ndc >= 5 else '不可接受'}",
    ]
    if res.pct_tol is not None:
        lines.append(f"%公差 = {res.pct_tol:.1%}（公差 {res.tolerance}）")
    lines += [
        f"主要變異來源：{'量具（EV > AV）' if res.ev > res.av else '人員／方法（AV > EV）'}",
        "",
        f"R 圖：UCL = {r_ucl:.4f}，超出 {int((s['ranges'] > r_ucl).sum())} 點（期望 0 點）",
        f"X̄ 圖：界限 = {s['xbarbar']:.4f} ± {x_half:.4f}，超出 {outside.mean():.0%}（期望 > 50%）",
    ]
    return "\n".join(lines)


def plot(data: np.ndarray, res: GRRResult, s: dict, appraisers: list[str], out: Path) -> None:
    plt.rcParams["font.sans-serif"] = ["Microsoft JhengHei", "PingFang TC", "Noto Sans CJK TC", "DejaVu Sans"]
    plt.rcParams["axes.unicode_minus"] = False
    k, n, r = data.shape
    colors = ["#264653", "#e76f51", "#2a9d8f"]

    fig, axes = plt.subplots(2, 2, figsize=(13, 8.5))
    fig.subplots_adjust(hspace=0.4, wspace=0.25)

    # 1. 變異來源
    ax = axes[0, 0]
    labels = ["EV\n重複性", "AV\n再現性", "GRR", "PV\n零件"]
    vals = [res.pct(v) * 100 for v in (res.ev, res.av, res.grr, res.pv)]
    bars = ax.bar(labels, vals, color=["#8ecae6", "#8ecae6", "#e76f51", "#2a9d8f"])
    ax.bar_label(bars, fmt="%.1f%%", fontsize=9)
    ax.axhline(10, color="#2a9d8f", ls="--", lw=1)
    ax.axhline(30, color="#d62828", ls="--", lw=1)
    ax.text(-0.45, 10, "10% 可接受", va="bottom", ha="left", fontsize=8, color="#2a9d8f")
    ax.text(-0.45, 30, "30% 不可接受", va="bottom", ha="left", fontsize=8, color="#d62828")
    ax.set_ylim(0, 115)
    ax.set_ylabel("% 總變異")
    title = f"變異來源（%GRR = {res.pct(res.grr):.1%}，ndc = {res.ndc}）"
    ax.set_title(title)

    # 2. R 圖（依量測人員）
    ax = axes[0, 1]
    r_ucl = D4[r] * s["rbarbar"]
    for i, a in enumerate(appraisers):
        x = np.arange(n) + i * (n + 1)
        ax.plot(x, s["ranges"][i], marker="o", ms=4, color=colors[i % 3], lw=1, label=a)
    ax.axhline(s["rbarbar"], color="#2a9d8f", lw=1)
    ax.axhline(r_ucl, color="#e76f51", ls="--", lw=1.2)
    ax.text(k * (n + 1) - 1.5, r_ucl, "UCL", va="bottom", ha="right", fontsize=8, color="#e76f51")
    ax.set_xticks([i * (n + 1) + (n - 1) / 2 for i in range(k)], appraisers)
    ax.set_title("R 圖：每人重複量測的一致性（應全部在 UCL 內）")
    ax.set_ylim(0, max(r_ucl, s["ranges"].max()) * 1.25)

    # 3. X̄ 圖（依量測人員）
    ax = axes[1, 0]
    half = A2[r] * s["rbarbar"]
    for i, a in enumerate(appraisers):
        x = np.arange(n) + i * (n + 1)
        ax.plot(x, s["cell_means"][i], marker="o", ms=4, color=colors[i % 3], lw=1)
    ax.axhline(s["xbarbar"], color="#2a9d8f", lw=1)
    ax.axhspan(s["xbarbar"] - half, s["xbarbar"] + half, color="#e76f51", alpha=0.15, label="量測雜訊範圍")
    ax.set_xticks([i * (n + 1) + (n - 1) / 2 for i in range(k)], appraisers)
    outside = (np.abs(s["cell_means"] - s["xbarbar"]) > half).mean()
    ax.set_title(rf"$\bar{{X}}$ 圖：零件差異 vs 量測雜訊（{outside:.0%} 超出，應 > 50%）")
    ax.legend(loc="upper left", fontsize=8)

    # 4. 零件 × 量測人員
    ax = axes[1, 1]
    parts = np.arange(1, n + 1)
    for i, a in enumerate(appraisers):
        ax.plot(parts, s["cell_means"][i], marker="o", ms=4, color=colors[i % 3], lw=1.2, label=f"{a} 員")
    ax.set_xticks(parts)
    ax.set_xlabel("零件編號")
    ax.set_title("各量測人員的零件平均（線應幾乎重疊）")
    ax.legend(fontsize=8)

    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=120, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")

    parser = argparse.ArgumentParser(description="GR&R 平均值與全距法")
    parser.add_argument("--csv", type=Path, help="長表 CSV：appraiser,part,trial,value")
    parser.add_argument("--tol", type=float, help="公差寬度（USL − LSL），用來計算 %公差")
    parser.add_argument("--out", type=Path, default=Path("tools/output/grr.png"), help="輸出圖檔路徑")
    parser.add_argument("--demo", action="store_true", help="使用內建範例數據")
    args = parser.parse_args()

    if args.demo:
        data, appraisers = make_demo_data()
        tol = args.tol or 0.60
    elif args.csv:
        data, appraisers = load_csv(args.csv)
        tol = args.tol
    else:
        parser.error("請提供 --csv，或使用 --demo")

    k, n, r = data.shape
    if r not in K1 or k not in K2 or n not in K3:
        sys.exit(f"不支援的設定：k = {k}（需 2～3）、n = {n}（需 2～10）、r = {r}（需 2～3）")

    res, s = analyze(data, tol)
    print(report(res, s, appraisers, data.shape))
    plot(data, res, s, appraisers, args.out)
    print(f"\n圖檔已輸出：{args.out}")


if __name__ == "__main__":
    main()

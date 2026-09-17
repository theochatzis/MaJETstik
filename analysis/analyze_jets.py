#!/usr/bin/env python3
"""Stage 5 -- read the NanoAOD file, print a summary, make plots.

This stage does no simulation; it is the example analysis.  Everything here is
ordinary uproot + numpy + matplotlib, and it is deliberately written the way a
first real analysis looks, so it can be copied and modified.

What it measures
----------------
* **Jet pT spectrum** -- steeply falling, as QCD radiation always is.
  Reconstructed and generator-level jets are overlaid; the gap between them at
  low pT is where detector thresholds and resolution bite.
* **Jet mass, raw and groomed** -- soft-drop removes soft wide-angle radiation,
  so the groomed distribution is pulled towards lower masses and is much less
  sensitive to the underlying event.
* **N-subjettiness ratio tau21 = tau2/tau1** -- small values indicate a jet
  with two distinct prongs (a boosted W/Z/H); ordinary quark and gluon jets
  populate the high-tau21 region.  With R = 0.4 jets almost everything is
  one-pronged, which is itself the lesson: substructure needs large-R jets.
* **Jet energy response** pT(reco)/pT(gen) for jets matched within DeltaR < R/2
  -- centred slightly below 1 because some of the jet's energy escapes the
  jet cone and because neutrinos are invisible.
* **Reconstruction efficiency vs generator jet pT** -- the fraction of truth
  jets that have a reconstructed partner, the standard turn-on curve.
* **Z -> e+e- invariant mass** -- the calibration candle: a peak at the PDG Z
  mass of 91.19 GeV, with a low-mass tail from final-state radiation.

Input
-----
``data/<run_name>_nano.root`` from stage 4.

Output
------
PNG figures in ``data/plots/`` and a summary table on the console and in
``data/logs/analyze_jets.log``.

Run standalone
--------------
    source setup_env.sh
    python3 analysis/analyze_jets.py --run-name smoke

How to modify
-------------
Each plot is a small self-contained function taking the loaded arrays; add
another one and call it from :func:`analyze`.  To work interactively instead,
see ``examples/analysis_tutorial.ipynb``.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Dict

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from majetstik.config import PipelineConfig
from majetstik.logging_setup import banner, get_logger
from majetstik.provenance import build_provenance, merge_provenance

# ---------------------------------------------------------------------------
# Plot styling.  Three categorical colours, used in a fixed order and never
# cycled; validated for colour-vision deficiency as a set.  Every figure with
# more than one series carries a legend, so colour is never the only cue.
# ---------------------------------------------------------------------------
SERIES = ("#2a78d6", "#eb6834", "#1baf7a")   # blue, orange, aqua
SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_MUTED = "#52514e"
GRID = "#d8d7d2"
PDG_Z_MASS = 91.1876                          # GeV, PDG 2024


def _style():
    """Apply a recessive, print-friendly style to matplotlib."""
    import matplotlib
    matplotlib.use("Agg")                     # no display needed
    import matplotlib.pyplot as plt

    plt.rcParams.update({
        "figure.figsize": (6.4, 4.2),
        "figure.facecolor": SURFACE,
        "axes.facecolor": SURFACE,
        "savefig.facecolor": SURFACE,
        "axes.edgecolor": GRID,
        "axes.labelcolor": INK,
        "axes.titlesize": 12,
        "axes.titleweight": "semibold",
        "axes.titlecolor": INK,
        "axes.titlelocation": "left",
        # Leave room above the axes for the one-line subtitle added by _finish.
        "axes.titlepad": 22,
        "axes.grid": True,
        "grid.color": GRID,
        "grid.linewidth": 0.8,
        "grid.alpha": 0.7,
        "xtick.color": INK_MUTED,
        "ytick.color": INK_MUTED,
        "text.color": INK,
        "legend.frameon": False,
        "font.size": 10,
        "lines.linewidth": 2.0,
    })
    return plt


def _finish(fig, ax, path: Path, subtitle: str = "") -> Path:
    """Common axis cleanup and save."""
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    if subtitle:
        # Sits in the gap opened by axes.titlepad, below the left-aligned title.
        ax.text(0.0, 1.01, subtitle, transform=ax.transAxes,
                fontsize=9, color=INK_MUTED, ha="left", va="bottom")
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=150)
    return path


def load(cfg: PipelineConfig):
    """Load the NanoAOD tree into awkward arrays."""
    import uproot
    if not cfg.nanoaod_path.is_file():
        raise FileNotFoundError(
            f"{cfg.nanoaod_path} not found -- run stage 4 first:\n"
            f"    python3 analysis/write_nanoaod.py --run-name {cfg.run_name}")
    with uproot.open(cfg.nanoaod_path) as f:
        return f["Events"].arrays(), str(f["ReadMe"])


def hadronic_jet_mask(events):
    """Jets that are genuine hadronic jets, not an isolated lepton.

    Particle flow hands the jet algorithm *everything* it reconstructed,
    electrons included, so a Z -> ee event really does produce two jets that
    are electrons.  Stage 4 flagged them; a hadronic analysis removes them.
    This is exactly what an experiment's "overlap removal" does.
    """
    return events["Jet_isLepton"] == 0


# ---------------------------------------------------------------------------
# individual figures
# ---------------------------------------------------------------------------
def plot_jet_pt(events, cfg, plt) -> Path:
    import awkward as ak
    reco = np.asarray(ak.flatten(events["Jet_pt"][hadronic_jet_mask(events)]))
    gen = np.asarray(ak.flatten(events["GenJet_pt"]))
    if reco.size == 0:
        return None

    hi = max(np.percentile(reco, 99.5), cfg.jet_pt_min * 2)
    bins = np.linspace(cfg.jet_pt_min, hi, 40)
    fig, ax = plt.subplots()
    ax.hist(gen, bins=bins, histtype="step", color=SERIES[1],
            label=f"Generator jets ({gen.size})")
    ax.hist(reco, bins=bins, histtype="step", color=SERIES[0],
            label=f"Reconstructed jets ({reco.size})")
    ax.set_yscale("log")
    ax.set_xlabel(r"jet $p_T$  [GeV]")
    ax.set_ylabel("jets / bin")
    ax.set_title("Jet transverse-momentum spectrum")
    ax.legend()
    return _finish(fig, ax, cfg.plot_dir / "jet_pt.png",
                   f"anti-$k_T$ R={cfg.jet_r}, {cfg.run_name}")


def plot_jet_mass(events, cfg, plt) -> Path:
    import awkward as ak
    mask = hadronic_jet_mask(events)
    raw = np.asarray(ak.flatten(events["Jet_mass"][mask]))
    sd = np.asarray(ak.flatten(events["Jet_softdrop_mass"][mask]))
    sd = sd[sd >= 0]                                   # -1 = grooming failed
    if raw.size == 0:
        return None

    bins = np.linspace(0, max(np.percentile(raw, 99), 10.0), 40)
    fig, ax = plt.subplots()
    ax.hist(raw, bins=bins, histtype="step", color=SERIES[0],
            label="ungroomed mass")
    ax.hist(sd, bins=bins, histtype="step", color=SERIES[1],
            label=f"soft-drop mass ($z_{{cut}}$={cfg.softdrop_zcut}, "
                  rf"$\beta$={cfg.softdrop_beta:g})")
    ax.set_xlabel(r"jet mass  [GeV]")
    ax.set_ylabel("jets / bin")
    ax.set_title("Jet mass before and after grooming")
    ax.legend()
    return _finish(fig, ax, cfg.plot_dir / "jet_mass.png",
                   "grooming removes soft, wide-angle radiation")


def plot_nsubjettiness(events, cfg, plt) -> Path:
    import awkward as ak
    mask = hadronic_jet_mask(events)
    tau1 = np.asarray(ak.flatten(events["Jet_tau1"][mask]))
    tau2 = np.asarray(ak.flatten(events["Jet_tau2"][mask]))
    good = (tau1 > 0) & (tau2 >= 0)
    if good.sum() == 0:
        return None
    ratio = tau2[good] / tau1[good]

    fig, ax = plt.subplots()
    ax.hist(ratio, bins=np.linspace(0, 1.2, 40), histtype="step",
            color=SERIES[0])
    ax.set_xlabel(r"$\tau_{21} = \tau_2/\tau_1$")
    ax.set_ylabel("jets / bin")
    ax.set_title("N-subjettiness ratio")
    return _finish(fig, ax, cfg.plot_dir / "jet_tau21.png",
                   "low values = two-pronged jet; most R=0.4 jets are one-pronged")


def plot_response(events, cfg, plt) -> Dict[str, float]:
    """Jet energy response pT(reco)/pT(gen) for matched jets."""
    import awkward as ak
    mask = hadronic_jet_mask(events)
    jet_pt = events["Jet_pt"][mask]
    gidx = events["Jet_genJetIdx"][mask]
    gen_pt = events["GenJet_pt"]

    ratios = []
    for i in range(len(jet_pt)):
        gi = np.asarray(gidx[i])
        if gi.size == 0:
            continue
        g = np.asarray(gen_pt[i])
        sel = (gi >= 0)
        if sel.sum() == 0 or g.size == 0:
            continue
        ratios.append(np.asarray(jet_pt[i])[sel] / g[gi[sel]])
    if not ratios:
        return None, {}
    ratio = np.concatenate(ratios)

    fig, ax = plt.subplots()
    ax.hist(ratio, bins=np.linspace(0.3, 1.7, 40), histtype="step",
            color=SERIES[0])
    ax.axvline(1.0, color=INK_MUTED, linewidth=1.2, linestyle="--")
    stats = {"mean": float(np.mean(ratio)), "std": float(np.std(ratio)),
             "median": float(np.median(ratio)), "n": int(ratio.size)}
    ax.text(0.03, 0.95,
            f"mean {stats['mean']:.3f}\nRMS  {stats['std']:.3f}\n"
            f"{stats['n']} matched jets",
            transform=ax.transAxes, va="top", ha="left",
            fontsize=9, color=INK_MUTED)
    ax.set_xlabel(r"$p_T^{\rm reco} / p_T^{\rm gen}$")
    ax.set_ylabel("jets / bin")
    ax.set_title("Jet energy response")
    path = _finish(fig, ax, cfg.plot_dir / "jet_response.png",
                   f"matched within $\\Delta R$ < {cfg.jet_r/2:g}")
    return path, stats


def plot_efficiency(events, cfg, plt):
    """Fraction of generator jets with a reconstructed partner, vs gen pT.

    This is the standard "turn-on curve": efficiency rises from zero below the
    reconstruction threshold to a plateau once jets are comfortably above it.

    Two details decide whether the number means anything:

    * **Match against *all* reconstructed jets, not just the hadronic ones.**
      The generator jet collection contains jets built from the two electrons
      just as the reconstructed collection does.  Vetoing lepton-like jets in
      the numerator while leaving their partners in the denominator would
      depress the efficiency by ~40% for no physical reason.
    * **Restrict the denominator to the detector acceptance.**  Delphes
      clusters generator jets over all pseudorapidity, including forward beam
      remnants that no detector covers; counting those as inefficiency would
      again say nothing about reconstruction.
    """
    import awkward as ak
    gidx = events["Jet_genJetIdx"]
    gen_pt = events["GenJet_pt"]
    gen_eta = events["GenJet_eta"]

    all_gen, matched_gen = [], []
    for i in range(len(gen_pt)):
        g = np.asarray(gen_pt[i])
        if g.size == 0:
            continue
        in_acceptance = np.abs(np.asarray(gen_eta[i])) < cfg.jet_eta_max
        all_gen.append(g[in_acceptance])

        gi = np.asarray(gidx[i])
        hit = np.unique(gi[gi >= 0]) if gi.size else np.array([], dtype=int)
        hit = hit[(hit < g.size)]
        hit = hit[in_acceptance[hit]]
        matched_gen.append(g[hit])
    if not all_gen:
        return None, {}
    all_gen = np.concatenate(all_gen)
    matched_gen = np.concatenate(matched_gen) if matched_gen else np.array([])
    if all_gen.size == 0:
        return None, {}

    bins = np.linspace(cfg.jet_pt_min, max(np.percentile(all_gen, 99), 100), 12)
    n_all, _ = np.histogram(all_gen, bins=bins)
    n_match, _ = np.histogram(matched_gen, bins=bins)
    centres = 0.5 * (bins[1:] + bins[:-1])
    with np.errstate(divide="ignore", invalid="ignore"):
        eff = np.where(n_all > 0, n_match / n_all, np.nan)
        # Binomial (approximate) uncertainty on an efficiency.
        err = np.where(n_all > 0,
                       np.sqrt(np.clip(eff * (1 - eff), 0, None) / np.maximum(n_all, 1)),
                       np.nan)

    fig, ax = plt.subplots()
    ax.errorbar(centres, eff, yerr=err, fmt="o", markersize=5,
                color=SERIES[0], ecolor=SERIES[0], capsize=3, linewidth=1.5)
    ax.axhline(1.0, color=INK_MUTED, linewidth=1.0, linestyle="--")
    ax.set_ylim(0, 1.15)
    ax.set_xlabel(r"generator jet $p_T$  [GeV]")
    ax.set_ylabel("reconstruction efficiency")
    ax.set_title("Jet reconstruction efficiency")
    overall = float(matched_gen.size / all_gen.size)
    path = _finish(fig, ax, cfg.plot_dir / "jet_efficiency.png",
                   f"generator jets with $|\\eta|$ < {cfg.jet_eta_max:g}; "
                   f"overall {overall:.1%}")
    return path, {"overall_efficiency": overall,
                  "n_gen_jets": int(all_gen.size),
                  "n_matched": int(matched_gen.size)}


def plot_dielectron_mass(events, cfg, plt):
    """Invariant mass of the two highest-pT opposite-charge electrons."""
    import awkward as ak
    n_e = ak.num(events["Electron_pt"])
    sel = n_e >= 2
    if int(ak.sum(sel)) == 0:
        return None, {}
    pt = events["Electron_pt"][sel][:, :2]
    eta = events["Electron_eta"][sel][:, :2]
    phi = events["Electron_phi"][sel][:, :2]
    q = events["Electron_charge"][sel][:, :2]

    px, py = pt * np.cos(phi), pt * np.sin(phi)
    pz = pt * np.sinh(eta)
    # Electrons are effectively massless at these energies (m_e = 0.511 MeV).
    e = np.sqrt(px ** 2 + py ** 2 + pz ** 2)
    m2 = ((e[:, 0] + e[:, 1]) ** 2 - (px[:, 0] + px[:, 1]) ** 2
          - (py[:, 0] + py[:, 1]) ** 2 - (pz[:, 0] + pz[:, 1]) ** 2)
    mass = np.asarray(np.sqrt(np.maximum(m2, 0)))
    opposite = np.asarray(q[:, 0] * q[:, 1] < 0)
    mass = mass[opposite]
    if mass.size == 0:
        return None, {}

    fig, ax = plt.subplots()
    ax.hist(mass, bins=np.linspace(60, 120, 40), histtype="step",
            color=SERIES[0])
    ax.axvline(PDG_Z_MASS, color=SERIES[1], linewidth=1.5, linestyle="--")
    ax.annotate(f"PDG $m_Z$ = {PDG_Z_MASS:.2f} GeV",
                xy=(PDG_Z_MASS, ax.get_ylim()[1] * 0.92),
                xytext=(4, 0), textcoords="offset points",
                fontsize=9, color=SERIES[1], va="top")
    ax.set_xlabel(r"$m(e^+e^-)$  [GeV]")
    ax.set_ylabel("events / bin")
    ax.set_title("Dielectron invariant mass")
    path = _finish(fig, ax, cfg.plot_dir / "dielectron_mass.png",
                   f"{mass.size} opposite-charge pairs")
    return path, {"n_pairs": int(mass.size), "mean_mass": float(np.mean(mass)),
                  "median_mass": float(np.median(mass))}


# ---------------------------------------------------------------------------
def summarize(events, log) -> Dict[str, float]:
    """Print the numbers a physicist checks first."""
    import awkward as ak

    n_events = len(events["event"])
    mask = hadronic_jet_mask(events)
    n_jets = int(ak.sum(ak.num(events["Jet_pt"])))
    n_had = int(ak.sum(ak.num(events["Jet_pt"][mask])))
    jet_pt = ak.flatten(events["Jet_pt"][mask])

    rows = [
        ("events", f"{n_events}"),
        ("jets (all)", f"{n_jets}  ({n_jets / max(n_events,1):.2f} / event)"),
        ("jets (hadronic, Jet_isLepton==0)",
         f"{n_had}  ({n_had / max(n_events,1):.2f} / event)"),
        ("electrons", f"{int(ak.sum(ak.num(events['Electron_pt'])))}"),
        ("muons", f"{int(ak.sum(ak.num(events['Muon_pt'])))}"),
        ("b-tagged jets",
         f"{int(ak.sum(ak.flatten(events['Jet_btag']) == 1))}"),
        ("mean MET", f"{float(ak.mean(events['MET_pt'])):.1f} GeV"),
    ]
    if len(jet_pt):
        pts = np.asarray(jet_pt)
        rows += [
            ("jet pT  mean / median", f"{pts.mean():.1f} / {np.median(pts):.1f} GeV"),
            ("jet pT  max", f"{pts.max():.1f} GeV"),
            ("mean constituents/jet",
             f"{float(ak.mean(ak.flatten(events['Jet_nConstituents'][mask]))):.1f}"),
        ]

    log.info("-" * 60)
    log.info("SUMMARY")
    log.info("-" * 60)
    for label, value in rows:
        log.info("  %-34s %s", label, value)
    log.info("-" * 60)
    return {"n_events": n_events, "n_jets": n_jets, "n_hadronic_jets": n_had}


def analyze(cfg: PipelineConfig, log=None) -> dict:
    log = log or get_logger("analyze_jets", cfg.log_dir)
    banner(log, "STAGE 5/5  Analysis and plots")

    events, readme = load(cfg)
    log.info("input        : %s", cfg.nanoaod_path)
    stats = summarize(events, log)

    if not cfg.make_plots:
        log.info("plotting disabled (make_plots=False)")
        return stats

    plt = _style()
    cfg.plot_dir.mkdir(parents=True, exist_ok=True)
    made = []

    for fn in (plot_jet_pt, plot_jet_mass, plot_nsubjettiness):
        try:
            p = fn(events, cfg, plt)
            if p:
                made.append(p)
        except Exception as exc:                       # a bad plot must not
            log.warning("%s failed: %s", fn.__name__, exc)   # kill the stage

    for fn in (plot_response, plot_efficiency, plot_dielectron_mass):
        try:
            p, extra = fn(events, cfg, plt)
            if p:
                made.append(p)
            stats.update({f"{fn.__name__}.{k}": v for k, v in extra.items()})
        except Exception as exc:
            log.warning("%s failed: %s", fn.__name__, exc)

    for p in made:
        log.info("plot         : %s", p)
    log.info("wrote %d figures to %s", len(made), cfg.plot_dir)

    if "plot_response.mean" in stats:
        log.info("jet energy response : %.3f +- %.3f",
                 stats["plot_response.mean"], stats["plot_response.std"])
    if "plot_efficiency.overall_efficiency" in stats:
        log.info("jet reco efficiency : %.1f%%",
                 100 * stats["plot_efficiency.overall_efficiency"])
    if "plot_dielectron_mass.median_mass" in stats:
        log.info("median m(ee)        : %.2f GeV  (PDG m_Z = %.2f)",
                 stats["plot_dielectron_mass.median_mass"], PDG_Z_MASS)

    stats["plots"] = [str(p) for p in made]
    merge_provenance(cfg.provenance_path,
                     build_provenance(cfg, "analyze_jets", **stats))
    return stats


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description="Stage 5: analyse the NanoAOD file and make plots.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("--config", help="JSON config file")
    ap.add_argument("--run-name", help="label for input files")
    ap.add_argument("--no-plots", action="store_true", help="print the summary only")
    ap.add_argument("--print-readme", action="store_true",
                    help="dump the branch documentation stored in the ROOT file")
    args = ap.parse_args(argv)

    cfg = PipelineConfig.from_file(args.config) if args.config else PipelineConfig()
    if args.run_name: cfg.run_name = args.run_name
    if args.no_plots: cfg.make_plots = False
    cfg.validate()
    cfg.ensure_dirs()

    if args.print_readme:
        print(load(cfg)[1])
        return 0

    analyze(cfg)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

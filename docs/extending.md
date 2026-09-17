# Extending MaJETstik

Concrete recipes, ordered from "change one flag" to "add a new stage". Each one
says what to edit and what you should see afterwards.

---

## 1. Change the physics process

The generator is driven entirely by a Pythia command card, so a new process is
a new text file.

```bash
cp generators/cards/zee_jets.cmnd generators/cards/ttbar.cmnd
```

Edit the hard-process block:

```
! top pair production
Top:gg2ttbar = on
Top:qqbar2ttbar = on
```

Run it:

```bash
python3 main.py --pythia-card generators/cards/ttbar.cmnd \
                --run-name ttbar --events 2000
```

You should see b-tagged jets appear (`Jet_btag == 1`) and higher jet
multiplicity. Other useful starting points, all Pythia setting names:

| Process | Settings |
|---------|----------|
| QCD dijets | `HardQCD:all = on`, `PhaseSpace:pTHatMin = 200.` |
| $W \to \ell\nu$ | `WeakSingleBoson:ffbar2W = on` (gives real MET) |
| Higgs, gluon fusion | `HiggsSM:gg2H = on` |
| Boosted $W$ from top | top pairs with `PhaseSpace:pTHatMin = 400.` |

The full settings index is at <https://pythia.org/latest-manual/Welcome.html>.

---

## 2. Change the detector

The stock Delphes cards cover most experiments:

```bash
ls $DELPHES_DIR/cards
python3 main.py --delphes-card $DELPHES_DIR/cards/delphes_card_ATLAS.tcl \
                --run-name atlas
```

To modify a detector rather than swap it, edit
`simulation/cards/majetstik_cms.tcl`. Some worthwhile experiments:

- **Magnetic field**: `set Bz` in `module ParticlePropagator`. Lowering it from
  3.8 T degrades low-$p_T$ track momentum resolution — watch the jet energy
  response in `data/plots/jet_response.png` widen.
- **Calorimeter resolution**: the `set ResolutionFormula` blocks in `ECal` and
  `HCal`, written as $a/\sqrt{E} \oplus b$. Inflating the stochastic term
  broadens the jet mass distribution.
- **Tracking efficiency**: `ChargedHadronTrackingEfficiency`. Dropping it
  removes charged candidates from particle flow and biases jets low.

---

## 3. Change the b-tagging working point

Delphes parameterises b-tagging as an efficiency versus $p_T$, $\eta$ and true
flavour. In `module BTagging BTagging` of the card, each `add EfficiencyFormula
{flavour} {formula}` block gives one:

```tcl
# light jets faking a b (the mistag rate)
add EfficiencyFormula {0} {0.01 + 0.000038*pt}
# true c jets
add EfficiencyFormula {4} {0.25*tanh(0.018*pt)*(1/(1+0.0013*pt))}
# true b jets
add EfficiencyFormula {5} {0.85*tanh(0.0025*pt)*(25.0/(1+0.063*pt))}
```

A tighter working point means lower `{5}` and much lower `{0}`. Re-run stages
2–5 (`python3 main.py --from 2`) and compare the b-tagged jet count in the
stage-5 summary.

---

## 4. Add a substructure observable

The substructure kernel is the C++ string `_KERNEL_SOURCE` in
`majetstik/fastjet.py`, JIT-compiled by ROOT at import. Adding an observable
takes four edits. As an example, energy correlation functions:

1. **Include and compute** inside `Clusterer::cluster`:

   ```cpp
   #include "fastjet/contrib/EnergyCorrelator.hh"
   ...
   fastjet::contrib::EnergyCorrelatorD2 d2(nsub_beta_);
   out.ecfD2.push_back(cons.size() >= 2 ? d2(jet) : -1.f);
   ```

2. **Add the field** to `struct ClusterResult`:
   `std::vector<float> ecfD2;`

3. **Declare it in Python** by adding `"ecfD2": np.float32` to
   `JetClusterer.FIELDS` — the wrapper copies every field listed there.

4. **Propagate it** through `analysis/cluster_jets.py` (add `"ecfD2"` to the
   `per_jet` dict), `analysis/write_nanoaod.py` (add `Jet_ecfD2` to `columns`)
   and `BRANCH_DOC` (so the embedded ReadMe stays honest — stage 4 warns if a
   branch is undocumented).

The available contribs are listed by
`ls $(fastjet-config --prefix)/include/fastjet/contrib`.

---

## 5. Cluster generator-level particles instead of detector output

To see exactly what the detector does to a jet, cluster the same event twice.
In `analysis/cluster_jets.py`, `read_pf_candidates` reads the `EFlow*`
collections; the truth equivalent is the `Particle` collection with
`Particle.Status == 1` (and neutrinos, PDG 12/14/16, removed). Cluster both,
match them in $\Delta R$, and the ratio of the two jet $p_T$ values is the jet
energy response — already plotted in stage 5 using Delphes' own `GenJet`
collection.

---

## 6. Add pileup

Real LHC collisions come ~60 at a time. Delphes models this by overlaying
minimum-bias events:

1. Generate a minimum-bias sample (`SoftQCD:all = on`) and convert it to
   Delphes' pileup format with `root2pileup`.
2. Use a pileup-enabled card, e.g. `delphes_card_CMS_PileUp.tcl`, pointing
   `set PileUpFile` at that file and `set MeanPileUp` at the desired value.

This is where grooming earns its keep: compare `Jet_mass` with
`Jet_softdrop_mass` before and after adding pileup and the groomed mass will be
far more stable.

---

## 7. Speed it up

Generation is the slow stage (~40 events/s), and it is trivially parallel
because the seed fully determines the output:

```bash
for i in $(seq 1 8); do
  python3 main.py --run-name batch$i --seed $((1000 + i)) \
                  --events 5000 --to 4 &
done
wait
```

Each job writes its own `data/batch<i>_nano.root`; `uproot.concatenate` or
`uproot.iterate` reads them as one sample. **Use a different seed per job** —
identical seeds produce identical events, which silently multiplies your
statistics by nothing.

---

## 8. Add a sixth stage

Stages are plain scripts with one entry-point function, registered in the
`STAGES` dict of `main.py`:

```python
STAGES = {
    ...
    6: ("My selection", "analysis.my_selection", "select"),
}
```

The function takes a `PipelineConfig` and returns a dict. Follow the pattern of
`analysis/analyze_jets.py`: derive paths from `cfg`, log through
`get_logger`, and record what you did with `merge_provenance` so the run stays
reproducible.

---

## Where things live

| I want to change... | Edit |
|---------------------|------|
| the physics process | `generators/cards/*.cmnd` |
| the detector | `simulation/cards/majetstik_cms.tcl` |
| jet algorithm, radius, thresholds | CLI flags, or `majetstik/config.py` defaults |
| substructure observables | `_KERNEL_SOURCE` in `majetstik/fastjet.py` |
| output branches | `analysis/write_nanoaod.py` (+ `BRANCH_DOC`) |
| plots | `analysis/analyze_jets.py` |
| pipeline stages | `STAGES` in `main.py` |

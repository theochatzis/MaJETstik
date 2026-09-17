# Output format

The final product of the pipeline is `data/<run_name>_nano.root`: one flat
`Events` tree modelled on CMS NanoAOD, plus two documentation objects.

```
<run_name>_nano.root
├── Events      TTree      one entry per collision event
├── ReadMe      TObjString the branch guide below, embedded in the file
└── Metadata    TObjString full provenance as JSON (seed, cards, versions)
```

The `ReadMe` and `Metadata` objects travel *inside* the file, so a ROOT file
that has been copied, renamed or emailed still says exactly how it was made.

```python
import uproot
with uproot.open("data/zee_jets_nano.root") as f:
    print(f["ReadMe"])        # branch documentation
    print(f["Metadata"])      # JSON provenance
    events = f["Events"]
```

## The one rule of the format

Branches are named `<Collection>_<variable>`. **Within an event, every branch
sharing a prefix has the same length**, given by the `n<Collection>` branch.
So `Jet_pt[3]`, `Jet_eta[3]` and `Jet_btag[3]` all describe the same jet, and
there are `nJet` of them. Nothing else is needed to read the file — no custom
classes, no dictionaries, no experiment software.

Types are `float` for kinematics, `int` for counts, indices and identifiers.
Angles are in radians, momenta, energies and masses in GeV.

---

## Branches

### Event level

| Branch | Description |
|--------|-------------|
| `run` | Run number. Always 1; present so the file looks like real NanoAOD. |
| `event` | Event number, 0-based, in generation order. |
| `nJet` | Number of reconstructed jets in this event. |
| `nElectron` | Number of reconstructed electrons. |
| `nMuon` | Number of reconstructed muons. |
| `nGenJet` | Number of generator-level (truth) jets. |
| `nPFCand` | Number of particle-flow candidates. |
| `nParticle` | Number of stored truth particles. |

### Jets  (`nJet` entries per event)

| Branch | Description |
|--------|-------------|
| `Jet_pt` | Jet transverse momentum [GeV]. Anti-kT, R as recorded in Metadata. |
| `Jet_eta` | Jet pseudorapidity eta = -ln(tan(theta/2)). |
| `Jet_phi` | Jet azimuthal angle [rad], in (-pi, pi]. |
| `Jet_mass` | Jet invariant mass [GeV], from the sum of constituent four-vectors. |
| `Jet_nConstituents` | Number of particle-flow candidates clustered into the jet. |
| `Jet_btag` | 1 if the matched Delphes jet is b-tagged, else 0; -1 if no Delphes jet was found within DeltaR < 0.2. Delphes applies a parameterised efficiency (~70% for true b jets, ~1% for light). |
| `Jet_flavor` | PDG id of the parton the matched Delphes jet came from (5 = b, 4 = c, 1-3 = light, 21 = gluon, 0 = undefined); -1 if unmatched. |
| `Jet_tau1` | 1-subjettiness: how well the jet is described by a single subjet. Small for an ordinary quark or gluon jet. -1 if undefined (too few constituents). |
| `Jet_tau2` | 2-subjettiness. tau2/tau1 is the standard two-prong (W/Z/H) tagger. |
| `Jet_tau3` | 3-subjettiness. tau3/tau2 is the standard three-prong (top) tagger. |
| `Jet_softdrop_mass` | Jet mass [GeV] after soft-drop grooming, which recursively removes soft wide-angle radiation. -1 if grooming rejected the jet. |
| `Jet_softdrop_pt` | Transverse momentum [GeV] of the groomed jet. |
| `Jet_isLepton` | 1 if an isolated reconstructed electron or muon lies within DeltaR < 0.2 of the jet axis, i.e. this 'jet' is really a lepton that particle flow also handed to the jet algorithm. Require Jet_isLepton == 0 for hadronic jets. |
| `Jet_genJetIdx` | Index into the GenJet_* arrays of the closest generator jet within DeltaR < 0.2, or -1 if there is none. Use it to measure jet energy response and reconstruction efficiency. |

### Jet constituent map

| Branch | Description |
|--------|-------------|
| `JetConstituent_jetIdx` | Flat jet->constituent map, part 1: for each (jet, constituent) pair, the index into Jet_*. Two flat arrays are used instead of a nested array because a NanoAOD tree stores only one level of variable length. |
| `JetConstituent_pfIdx` | Flat jet->constituent map, part 2: the index into the PFCand_* arrays of that constituent. Example: the constituents of jet j are PFCand_pt[JetConstituent_pfIdx[JetConstituent_jetIdx == j]]. |

### Particle-flow candidates  (`nPFCand` entries per event)

| Branch | Description |
|--------|-------------|
| `PFCand_pt` | Particle-flow candidate transverse momentum [GeV]. |
| `PFCand_eta` | Particle-flow candidate pseudorapidity. |
| `PFCand_phi` | Particle-flow candidate azimuthal angle [rad]. |
| `PFCand_mass` | Particle-flow candidate mass [GeV]; 0 for calorimeter towers. |
| `PFCand_charge` | Electric charge in units of e; 0 for neutral candidates. |
| `PFCand_pdgId` | Species of the candidate: the track's PDG id for charged candidates, 22 for photon towers, 130 for neutral-hadron towers. Neutral towers carry no true species information. |

### Electrons  (`nElectron` entries per event)

| Branch | Description |
|--------|-------------|
| `Electron_pt` | Electron transverse momentum [GeV], after detector smearing. |
| `Electron_eta` | Electron pseudorapidity. |
| `Electron_phi` | Electron azimuthal angle [rad]. |
| `Electron_charge` | Electron charge in units of e (-1 for e-, +1 for e+). |
| `Electron_iso` | Relative isolation: scalar sum of the pT of other particles in a cone around the electron, divided by the electron pT. Prompt electrons from a Z have small values (< 0.1); electrons inside jets have large ones. |

### Muons  (`nMuon` entries per event)

| Branch | Description |
|--------|-------------|
| `Muon_pt` | Muon transverse momentum [GeV], after detector smearing. |
| `Muon_eta` | Muon pseudorapidity. |
| `Muon_phi` | Muon azimuthal angle [rad]. |
| `Muon_charge` | Muon charge in units of e. |
| `Muon_iso` | Relative isolation, defined as for Electron_iso. |

### Missing transverse energy  (one value per event)

| Branch | Description |
|--------|-------------|
| `MET_pt` | Missing transverse momentum [GeV]: the magnitude of the negative vector sum of all reconstructed objects. Large when neutrinos (or new invisible particles) escape. |
| `MET_phi` | Azimuthal direction [rad] of the missing transverse momentum. |
| `GenMET_pt` | Truth-level missing transverse momentum [GeV], from neutrinos only. |
| `GenMET_phi` | Truth-level missing transverse momentum direction [rad]. |

### Generator-level jets  (`nGenJet` entries per event)

| Branch | Description |
|--------|-------------|
| `GenJet_pt` | Generator-level jet pT [GeV]. Clustered from the stable truth particles (status 1, neutrinos removed) with exactly the same jet definition as Jet_*, so the two may be compared directly whatever the jet radius is set to. |
| `GenJet_eta` | Generator-level jet pseudorapidity. |
| `GenJet_phi` | Generator-level jet azimuthal angle [rad]. |
| `GenJet_mass` | Generator-level jet mass [GeV]. |

### Truth particles  (`nParticle` entries per event)

| Branch | Description |
|--------|-------------|
| `Particle_pt` | Truth particle transverse momentum [GeV]. |
| `Particle_eta` | Truth particle pseudorapidity. |
| `Particle_phi` | Truth particle azimuthal angle [rad]. |
| `Particle_mass` | Truth particle mass [GeV]. |
| `Particle_pdgId` | PDG identification code (11 = e-, 13 = mu-, 22 = gamma, 211 = pi+, 23 = Z, 2212 = proton, 21 = gluon, ...). |
| `Particle_status` | HepMC status code: 1 = stable final-state particle, 2 = decayed, 4 = beam particle, 21-29 = hard process. |
---

## Working with the jet → constituent map

A ROOT tree stores only **one** level of variable length per branch, so a
genuinely nested "list of constituents per jet per event" cannot be a single
branch. NanoAOD's solution — used by CMS's own PFNano and copied here — is to
flatten the relation into two parallel arrays, one entry per (jet, constituent)
pair:

```
JetConstituent_jetIdx : 0 0 0 0 1 1 1 2 2 2 2 2 ...   which jet
JetConstituent_pfIdx  : 4 9 2 7 1 3 8 0 5 6 11 12 ...  which PF candidate
```

To recover the constituents of jet `j` in one event:

```python
import numpy as np, uproot

with uproot.open("data/zee_jets_nano.root") as f:
    arrays = f["Events"].arrays(
        ["JetConstituent_jetIdx", "JetConstituent_pfIdx",
         "PFCand_pt", "PFCand_eta", "PFCand_phi", "PFCand_pdgId"],
        entry_stop=1)

jet_idx = np.asarray(arrays["JetConstituent_jetIdx"][0])
pf_idx  = np.asarray(arrays["JetConstituent_pfIdx"][0])
pf_pt   = np.asarray(arrays["PFCand_pt"][0])

j = 0
constituents = pf_idx[jet_idx == j]          # indices into the PFCand_* arrays
print(f"jet {j}: {len(constituents)} constituents, "
      f"scalar sum pT = {pf_pt[constituents].sum():.1f} GeV")
```

`Jet_nConstituents[j]` equals `(jet_idx == j).sum()`, so it can be used as a
cross-check or to size buffers.

---

## Two things that surprise people

**1. Some "jets" are electrons.** Particle flow reconstructs electrons and then
hands them to the jet algorithm along with everything else, so a
$Z \to e^+e^-$ event genuinely produces two jets that *are* the electrons. This
is not a bug; a real experiment removes them with "overlap removal". Stage 4
flags them, so a hadronic analysis should start with

```python
hadronic = events["Jet_isLepton"] == 0
jet_pt   = events["Jet_pt"][hadronic]
```

**2. `Jet_btag` and `Jet_flavor` can be −1.** FastJet clusters four-vectors and
knows nothing about flavour; those two variables are borrowed from the Delphes
jet matched within $\Delta R < R/2$. When there is no match — most often
because Delphes' own overlap removal already discarded that jet as an electron
— the value is −1, meaning "not available", which is distinct from
`Jet_btag == 0`, meaning "not b-tagged".

---

## Sizes

Roughly, for the default Z+jets configuration:

| Contents | per event |
|----------|-----------|
| Full file as shipped | ~5–6 kB |
| Without `Particle_*` (`--no-truth`) | ~2 kB |
| Without `PFCand_*` and the constituent map | ~0.5 kB |

The truth record dominates. Shrink it with `truth_pt_min` and
`truth_status_filter` in the config, or drop it with `--no-truth`.

---

## Reading it without Python

The point of the format is that nothing special is needed:

```cpp
// ROOT / C++
TFile f("data/zee_jets_nano.root");
TTree *t = (TTree*)f.Get("Events");
t->Draw("Jet_pt", "Jet_isLepton==0");
```

```python
# ROOT RDataFrame
import ROOT
df = ROOT.RDataFrame("Events", "data/zee_jets_nano.root")
h = df.Filter("nJet > 0").Histo1D("Jet_pt")
```

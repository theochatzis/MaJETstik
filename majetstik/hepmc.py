"""Convert a Pythia 8 event record into a valid HepMC3 file.

Why this module exists
----------------------
Delphes does not talk to Pythia directly in this setup: it reads an event file.
HepMC3 is the standard interchange format for that -- a text description of
every particle in the event and the vertices connecting them.

The C++ Pythia distribution ships a ready-made writer (``Pythia8::Pythia8ToHepMC``),
but the *Python* bindings in the LCG view do not expose it, so we build the
HepMC3 graph ourselves with ``pyHepMC3``.  That is only ~40 lines, and it makes
the event-record structure explicit, which is instructive in its own right.

The two subtleties
------------------
1. **One end vertex per particle.**  In HepMC3 a particle may have at most one
   decay ("end") vertex.  Pythia's ``motherList()`` can report several mothers
   for one particle (e.g. the two incoming partons of a hard scatter), and the
   same mother can appear in several daughters' mother lists.  If you blindly
   create a new vertex per distinct mother set, ``add_particle_in`` silently
   *reassigns* a mother's end vertex and leaves earlier vertices dangling, and
   the resulting file fails to parse.  We therefore reuse a mother's existing
   end vertex whenever it already has one.

2. **Parents must be written before children.**  HepMC3's ASCII writer encodes
   a simple production vertex inline, as the id of the parent particle, and the
   reader resolves those ids as it streams the file.  Pythia's initial-state
   shower entries (status -41/-42, backward evolution) legitimately list
   mothers that sit *later* in the event record, so writing particles in plain
   index order produces forward references that the reader cannot resolve --
   it then reports "too few vertices were parsed" and drops the event.
   :func:`topological_order` fixes this by emitting mothers first.

Both problems produce a silently *empty* Delphes output rather than a crash,
which is exactly the kind of bug this docstring is meant to save you from.
"""

from __future__ import annotations

from collections import deque
from pathlib import Path
from typing import Dict, List, Sequence, Tuple


def topological_order(event) -> Tuple[List[int], Dict[int, List[int]]]:
    """Order Pythia indices so that every mother precedes all of its daughters.

    Parameters
    ----------
    event : pythia8.Event
        The event record, ``event[0]`` being the overall system pseudo-particle
        which we skip.

    Returns
    -------
    order : list of int
        Pythia indices 1..N-1, topologically sorted (Kahn's algorithm).
    mothers : dict
        index -> list of its mother indices, already filtered to the valid
        range, so callers do not have to repeat that filtering.

    Notes
    -----
    If the mother relation ever contained a cycle (it should not), the
    remaining particles are appended in index order rather than dropped, so the
    function is total.
    """
    n = event.size()
    mothers = {i: [m for m in event[i].motherList() if 1 <= m < n]
               for i in range(1, n)}

    indegree = {i: len(ms) for i, ms in mothers.items()}
    children: Dict[int, List[int]] = {i: [] for i in mothers}
    for i, ms in mothers.items():
        for m in ms:
            children[m].append(i)

    # Seed with the beam particles (no mothers).  Sorting keeps the output
    # deterministic, which matters for reproducibility.
    queue = deque(sorted(i for i, d in indegree.items() if d == 0))
    order: List[int] = []
    while queue:
        i = queue.popleft()
        order.append(i)
        for c in children[i]:
            indegree[c] -= 1
            if indegree[c] == 0:
                queue.append(c)

    if len(order) < len(indegree):                # cycle guard, should not happen
        seen = set(order)
        order.extend(i for i in sorted(indegree) if i not in seen)
    return order, mothers


def event_to_hepmc3(pythia_event, event_number: int, hepmc3_module):
    """Build one :class:`HepMC3.GenEvent` from a Pythia event record.

    Units are GeV and mm, matching both Pythia's internal units and what
    Delphes expects.

    The HepMC status code comes from ``Particle.statusHepMC()``, which maps
    Pythia's internal status onto the HepMC convention:
    ``1`` = final-state, ``2`` = decayed, ``4`` = beam particle.  Delphes uses
    status 1 to decide what enters the detector, and the rest to reconstruct
    the truth record.
    """
    hm = hepmc3_module
    ev = hm.GenEvent(hm.Units.GEV, hm.Units.MM)
    ev.set_event_number(event_number)

    order, mothers = topological_order(pythia_event)

    # 1. Create every GenParticle up front (cheap, keeps the next loop clear).
    gen: Dict[int, object] = {}
    for i in order:
        pa = pythia_event[i]
        gp = hm.GenParticle(
            hm.FourVector(pa.px(), pa.py(), pa.pz(), pa.e()),
            pa.id(), pa.statusHepMC())
        gp.set_generated_mass(pa.m())
        gen[i] = gp

    # 2. Wire them into vertices, mothers first (see module docstring).
    for i in order:
        pa = pythia_event[i]
        mums = mothers[i]

        if not mums:
            # Beam particle: HepMC3 represents it as a particle with no
            # production vertex and status 4.
            ev.add_particle(gen[i])
            continue

        vertex = None
        for m in mums:
            if gen[m].end_vertex():
                vertex = gen[m].end_vertex()      # reuse; never reassign
                break
        if vertex is None:
            # Production point of the daughter == decay point of the mother.
            vertex = hm.GenVertex(hm.FourVector(pa.xProd(), pa.yProd(),
                                                pa.zProd(), pa.tProd()))
            ev.add_vertex(vertex)
        for m in mums:
            if not gen[m].end_vertex():
                vertex.add_particle_in(gen[m])
        vertex.add_particle_out(gen[i])

    return ev


class HepMC3Writer:
    """Small context-manager wrapper around ``HepMC3.WriterAscii``.

    Example
    -------
    >>> with HepMC3Writer("events.hepmc") as writer:     # doctest: +SKIP
    ...     for i in range(100):
    ...         if pythia.next():
    ...             writer.write(pythia.event, i)
    """

    def __init__(self, path: str | Path):
        from pyHepMC3 import HepMC3 as hm

        self._hm = hm
        self.path = Path(path)
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self._writer = hm.WriterAscii(str(self.path))
        self.n_written = 0

    def write(self, pythia_event, event_number: int) -> None:
        ev = event_to_hepmc3(pythia_event, event_number, self._hm)
        self._writer.write_event(ev)
        self.n_written += 1

    def close(self) -> None:
        # Closing flushes the END_EVENT_LISTING trailer; without it the file is
        # truncated and Delphes reads zero events.
        self._writer.close()

    def __enter__(self) -> "HepMC3Writer":
        return self

    def __exit__(self, *exc) -> None:
        self.close()


def count_events(path: str | Path) -> int:
    """Count events in a HepMC3 ASCII file by reading it back.

    Used as a sanity check after generation: a file that Delphes will silently
    treat as empty is detected here instead.
    """
    from pyHepMC3 import HepMC3 as hm

    reader = hm.ReaderAscii(str(path))
    n = 0
    while True:
        ev = hm.GenEvent()
        if not reader.read_event(ev) or reader.failed():
            break
        n += 1
    return n

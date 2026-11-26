#!/usr/bin/env python3
"""SnapGene primer write defaults.

Whenever primers are added or renamed on a .dna file, use this helper so the
Primers wrapper includes HybridizationParams and each BindingSite has
Component.hybridizedRange (same as location). Without these, SnapGene applies
stricter built-in thresholds and primers are hidden on the map until the user
manually sets "minimum continuous match" to 8.

Always write with sgffp.SgffWriter(..., preserve=True). Never raw-byte-replace
inside a SnapGene package.
"""
from __future__ import annotations
from typing import Optional

HYBRIDIZATION_PARAMS = {
    "minContinuousMatchLen": "8",
    "allowMismatch": "1",
    "minMeltingTemperature": "40",
    "showAdditionalFivePrimeMatches": "1",
    "minimumFivePrimeAnnealing": "15",
}


def _location(start: int, end: int) -> str:
    """sgffp BindingSite: 0-based start, exclusive end → SnapGene 1-based inclusive location."""
    if end <= start:
        return f"{start}-{start}"
    return f"{start}-{end - 1}"


def ensure_binding_site_component(bs) -> None:
    loc = _location(bs.start, bs.end)
    bases = bs.annealed_bases or ""
    bs.extras = dict(bs.extras or {})
    bs.extras["Component"] = {"hybridizedRange": loc, "bases": bases}


def ensure_primer_display_defaults(sgff) -> None:
    """Fill HybridizationParams + Component on existing primers. Call before write-back."""
    primers = sgff.primers
    _ = primers.items  # load
    wrap = primers._wrapper_extras
    wrap["HybridizationParams"] = dict(HYBRIDIZATION_PARAMS)
    wrap.setdefault("nextValidID", str(len(primers.items)))
    for p in primers.items:
        for bs in p.binding_sites:
            ensure_binding_site_component(bs)
    primers._sync()


def add_primer_with_defaults(
    sgff,
    name: str,
    sequence: str,
    *,
    bind_start: Optional[int] = None,
    bind_end: Optional[int] = None,
    bind_strand: str = "+",
    annealed_bases: Optional[str] = None,
    melting_temperature: Optional[float] = 60.0,
    description: Optional[str] = None,
):
    """Add one primer and immediately apply display defaults (including a simplified site)."""
    from sgffp.models.primer import SgffPrimer, SgffBindingSite

    sites = []
    if bind_start is not None:
        end = bind_end if bind_end is not None else bind_start + len(sequence)
        ann = annealed_bases if annealed_bases is not None else sequence
        for sim in (False, True):
            sites.append(
                SgffBindingSite(
                    start=bind_start,
                    end=end,
                    bound_strand=bind_strand,
                    annealed_bases=ann,
                    melting_temperature=melting_temperature,
                    simplified=sim,
                )
            )
    extras = {"recentID": str(len(sgff.primers.items))}
    if description:
        extras["description"] = description
    sgff.primers.add(SgffPrimer(name=name, sequence=sequence, binding_sites=sites, extras=extras))
    ensure_primer_display_defaults(sgff)
    return sgff

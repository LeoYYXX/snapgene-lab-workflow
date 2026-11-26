#!/usr/bin/env python3
"""Build the fictional example project (FASTA, SnapGene .dna, tables).

All sequences are short and synthetic. This script never reads a live lab tree.
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path
from textwrap import wrap

SKILL_ROOT = Path(__file__).resolve().parents[1]
REPO_ROOT = SKILL_ROOT.parents[1]
EX = REPO_ROOT / "examples"

COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def rc(seq: str) -> str:
    return seq.translate(COMPLEMENT)[::-1]


def wrap_fa(name: str, seq: str, comment: str) -> str:
    body = "\n".join(wrap(seq, 70, break_on_hyphens=False))
    return f">{name} {comment}\n{body}\n"


def translate(dna: str) -> str:
    from Bio.Seq import Seq

    return str(Seq(dna).translate(to_stop=False))


def write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")
    print("wrote", path.relative_to(REPO_ROOT), f"({path.stat().st_size} bytes)")


def add_feat(sg, name, typ, start, end, color, **kw):
    """start/end are 0-based half-open."""
    sg.add_feature(name, typ, start, end, strand="+", color=color, **kw)


def try_write_dna(path: Path, seq: str, topology: str, *, document_name: str, organism: str, notes_html: str, features, primers=(), accession=""):
    try:
        from sgffp import SgffObject, SgffWriter

        scripts = str(Path(__file__).resolve().parent)
        if scripts not in sys.path:
            sys.path.insert(0, scripts)
        from sgffp_primer_write import add_primer_with_defaults, ensure_primer_display_defaults
    except ImportError as exc:
        print("sgffp not available, skip", path.name, exc)
        return

    sg = SgffObject.new(seq, topology=topology)
    for feat in features:
        add_feat(sg, **feat)
    for pr in primers:
        add_primer_with_defaults(sg, **pr)
    ensure_primer_display_defaults(sg)
    notes = sg.notes
    notes.set("Description", notes_html)
    notes.set("Organism", organism)
    notes.set("Type", "Synthetic" if topology == "circular" or document_name.startswith("EV") or document_name.startswith("p") else "Natural")
    if accession:
        notes.set("AccessionNumber", accession)
    try:
        notes.set("DocumentName", document_name)
    except Exception:
        pass
    path.parent.mkdir(parents=True, exist_ok=True)
    SgffWriter(str(path)).write(sg, preserve=False)
    print("wrote", path.relative_to(REPO_ROOT), f"({path.stat().st_size} bytes)")


# ── colors (lab-style convention, generic) ────────────────────────────────
GREEN, ORANGE, GRAY, BLUE = "#008000", "#ff9900", "#808080", "#3366ff"
PURPLE, RED, TEAL = "#9900cc", "#cc0000", "#008080"

# ── EXA1 / EXA1L genomic ──────────────────────────────────────────────────
# Coordinates are 0-based half-open.
# upstream 60, 5'UTR 30, CDS1 60, intron 30, CDS2 93 (includes TAA), 3'UTR 30, downstream 30.
# CDS = 153 bp = 51 codons (50 aa + stop).

EXA1_UP = "GCTAGCCTAGGTTACCGATCCTAGGTTACATCGCTAGCCTAGGTTACCGATCCTAGGTTA"  # 60
EXA1_UTR5 = "TCTCAGGTAACTTAGGCCTTGAACCTTGCA"  # 30, divergent from paralog
EXA1_CDS1 = (
    "ATGAAGCTGGCTATCGGTGCTCCATTCGCTAAGATCGAGGCTCGTCTGAACGGTGCTATC"
)  # 60, 20 codons
EXA1_IVS = "GTAAGTAACCATGGATCCGGTACCTTGCAG"  # 30, gt…ag
# CDS2 includes an in-frame XhoI (CTCGAG = Leu-Glu) and a terminal TAA.
EXA1_CDS2 = (
    "CCAGCTGCTTACAAGGAGATCGCTAACCTGCTCGAGGACTTCATCGAGAAGGCTCCAGGT"
    "AACCTGGCTATCGACGGTAAGATCGAGCTGTAA"
)  # 93, in-frame XhoI CTCGAG
EXA1_UTR3 = "GGAATTCCTGTAAGCCATTGGCCTAAGTCA"  # 30
EXA1_DOWN = "TGCATGCATGCCTAGGATCCTAGGATGCAT"  # 30

# EXA1L: same CDS length, ~90% CDS identity, completely different UTRs.
EXA1L_UP = "CCGATTGGAACCTTACCGGATGGAACCTTACCCGATTGGAACCTTACCGGATGGAACCTA"  # 60
EXA1L_UTR5 = "AGACCGTTGGCATTCAGGAACTGGAACCTA"  # 30
# Silent / conservative changes vs EXA1; keep enough identity that a CDS primer would co-amplify.
EXA1L_CDS1 = (
    "ATGAAGCTGGCCATCGGTGCCCCATTCGCTAAGATTGAGGCTCGTCTCAACGGTGCTATC"
)  # 60 (3 mismatches vs EXA1)
EXA1L_IVS = "GTAAGTAACCTTGGATCCGGTACCTTCCAG"  # 30
EXA1L_CDS2 = (
    "CCAGCCGCTTACAAGGAGATCGCCAACCTGCTCGAGGATTTCATCGAGAAGGCCCCAGGT"
    "AACCTGGCCATCGACGGTAAGATTGAGCTGTAA"
)  # 93
EXA1L_UTR3 = "CCTTGGAAACAGTCGATTCCAAGGTACAAG"  # 30
EXA1L_DOWN = "GCTTGCTTGCCCTAGGATCCTAGGCTTGCT"  # 30


def assemble_gene(*parts: str) -> str:
    seq = "".join(parts)
    assert all(ch in "ACGT" for ch in seq), set(seq)
    return seq


EXA1 = assemble_gene(EXA1_UP, EXA1_UTR5, EXA1_CDS1, EXA1_IVS, EXA1_CDS2, EXA1_UTR3, EXA1_DOWN)
EXA1L = assemble_gene(EXA1L_UP, EXA1L_UTR5, EXA1L_CDS1, EXA1L_IVS, EXA1L_CDS2, EXA1L_UTR3, EXA1L_DOWN)
assert len(EXA1) == len(EXA1L) == 333
EXA1_CDS = EXA1_CDS1 + EXA1_CDS2
EXA1L_CDS = EXA1L_CDS1 + EXA1L_CDS2
assert len(EXA1_CDS) == len(EXA1L_CDS) == 153
assert EXA1_CDS.startswith("ATG") and EXA1_CDS.endswith("TAA")
assert "CTCGAG" in EXA1_CDS and "CTCGAG" in EXA1L_CDS
assert "*" not in translate(EXA1_CDS)[:-1] and translate(EXA1_CDS).endswith("*")
assert "*" not in translate(EXA1L_CDS)[:-1] and translate(EXA1L_CDS).endswith("*")
EXA1_CDNA = EXA1_UTR5 + EXA1_CDS + EXA1_UTR3
EXA1L_CDNA = EXA1L_UTR5 + EXA1L_CDS + EXA1L_UTR3
EXA1_CDS_NOSTOP = EXA1_CDS[:-3]
EXA1_PROMOTER_INSERT = EXA1_UP[10:] + EXA1_UTR5  # F core may start inward; ends before gene ATG

# Gene-feature coordinates (same layout for both paralogs)
GENE_FEATURES = [
    dict(name="upstream", typ="misc_feature", start=0, end=60, color=GRAY),
    dict(name="5'UTR", typ="5'UTR", start=60, end=90, color=GREEN),
    dict(name="CDS1", typ="CDS", start=90, end=150, color=ORANGE),
    dict(name="intron1", typ="intron", start=150, end=180, color=GRAY),
    dict(name="CDS2", typ="CDS", start=180, end=273, color=ORANGE),
    dict(name="3'UTR", typ="3'UTR", start=273, end=303, color=BLUE),
    dict(name="downstream", typ="misc_feature", start=303, end=333, color=GRAY),
]


# ── EV-FLAG (circular, synthetic) ─────────────────────────────────────────
# Layout (0-based):
#   0-80     synthetic ori-like region
#   80-170   synthetic bacterial marker CDS (NOT a real resistance gene)
#   170-226  P-synth promoter
#   226-246  F arm window ending with XhoI (CTCGAG)
#   246-248  stuffer TT (dropped by XhoI+SpeI / IF)
#   248-268  SpeI (ACTAGT) + linker GGTGGA + start of tag arm
#   then 3xFLAG-like tag + TAA + terminator + filler
#
# F arm (20 nt, plus strand, ends with XhoI): taken from the actual backbone.
# R arm (20 nt, as synthesized): RC of (ACTAGTGGTGGA + first 8 nt of tag).

EV_FLAG_ORI = "CGTACGTACGTACGTAGCTAGCTAGCTACGTACGTACGTACGTAGCTAGCTAGCTACGTACGTACGTACGTAGCTAGC"  # 80
EV_FLAG_MARKER = (
    "ATGCTGACCGCTCGTATCGCTGCTCGTCTGAACGCTCCAGGTCTGAAGGCTATCGACGCT"
    "CGTCTGCGTGCTATCGCTTAA"
)  # 90, 30 aa + stop, fictional marker
EV_FLAG_PROM = "GCTAGCTAGCTAGCTCCTAGGATCCTATATAAGGCTAGCCTAGGATCCGCTA"  # 52
# 14 nt prefix + XhoI = 20 nt F arm
EV_FLAG_PRE_XHO = "TTGGTACCGAATTC"  # 14
EV_FLAG_XHO = "CTCGAG"
EV_FLAG_STUFFER = "TT"
EV_FLAG_SPE = "ACTAGT"
EV_FLAG_LINKER = "GGTGGA"
# Standard 3xFLAG AA tag (DYKDHDGDYKDHDIDYKDDDDK) — a public tag sequence, not lab data.
EV_FLAG_TAG = "GACTACAAAGACCATGACGGTGATTATAAAGATCATGACATCGATTACAAGGATGACGATGACAAG"
EV_FLAG_STOP = "TAA"
EV_FLAG_TERM = "GCGGCCGCACTAGTCCATGGGCTAGCGAATTC"  # 32, synthetic
EV_FLAG_FILL = "CATTGCATTGCATTGCATTGCATTGCATTGCATTGCATTGCATTGCATTGCATTG"  # 54

EV_FLAG = assemble_gene(
    EV_FLAG_ORI,
    EV_FLAG_MARKER,
    EV_FLAG_PROM,
    EV_FLAG_PRE_XHO,
    EV_FLAG_XHO,
    EV_FLAG_STUFFER,
    EV_FLAG_SPE,
    EV_FLAG_LINKER,
    EV_FLAG_TAG,
    EV_FLAG_STOP,
    EV_FLAG_TERM,
    EV_FLAG_FILL,
)

# Coordinates
c = 0

def span(n):
    global c
    a = c
    c += n
    return a, c

s_ori = span(len(EV_FLAG_ORI))
s_mark = span(len(EV_FLAG_MARKER))
s_prom = span(len(EV_FLAG_PROM))
s_pre = span(len(EV_FLAG_PRE_XHO))
s_xho = span(len(EV_FLAG_XHO))
s_stu = span(len(EV_FLAG_STUFFER))
s_spe = span(len(EV_FLAG_SPE))
s_lnk = span(len(EV_FLAG_LINKER))
s_tag = span(len(EV_FLAG_TAG))
s_stp = span(len(EV_FLAG_STOP))
s_term = span(len(EV_FLAG_TERM))
s_fill = span(len(EV_FLAG_FILL))
assert c == len(EV_FLAG)
F_ARM_FLAG = EV_FLAG[s_pre[0] : s_xho[1]]  # 20 nt ending in CTCGAG
assert F_ARM_FLAG.endswith("CTCGAG") and len(F_ARM_FLAG) == 20
DOWN_ARM_FLAG = EV_FLAG[s_spe[0] : s_spe[0] + 20]  # ACTAGTGGTGGA + 8 nt of tag
assert DOWN_ARM_FLAG.startswith("ACTAGT")
R_ARM_FLAG = rc(DOWN_ARM_FLAG)
assert EV_FLAG.count(F_ARM_FLAG) == 1 and EV_FLAG.count(DOWN_ARM_FLAG) == 1

POE_INSERT = EXA1_CDS_NOSTOP
POE_PRODUCT = EV_FLAG[: s_xho[1]] + POE_INSERT + EV_FLAG[s_spe[0] :]
assert len(POE_PRODUCT) == len(EV_FLAG) - len(EV_FLAG_STUFFER) + len(POE_INSERT)


# ── EV-GUS (circular, transcriptional reporter) ───────────────────────────
# EcoRI … stuffer … NcoI (CCATGG, reporter ATG is nt 2–4) … synthetic reporter CDS
# Promoter insert ends BEFORE the gene's own ATG. Reverse primer has no ±TAA.

EV_GUS_ORI = "TACGTACGTACGTAGCTACGTACGTACGTAGCTACGTACGTACGTAGCTACGTACGTACGTAGCTACGTACGTACGTA"  # 80
EV_GUS_MARKER = (
    "ATGGCTCGTAAGCTGATCGCTCGTCCAGGTCTGAACGCTATCGCTCGTGCTCTG"
    "AAGGACATCGCTAAGCGTTAA"
)  # 84
EV_GUS_PROM = "TATATAAGGCTAGCCTAGGTCCGGATCGCTAGCTCCTAGGAAT"  # 43, ends just before EcoRI prefix
EV_GUS_PRE_ECO = "ACGTTGCTAG"  # 10
EV_GUS_ECO = "GAATTC"
EV_GUS_STUFFER = "GCTGAGGCTTAATATGATTA"  # 20, fictional cassette, not a live-lab sequence
EV_GUS_NCO = "CCATGG"  # reporter ATG is ATG
# Short synthetic reporter (NOT uidA / real GUS). 90 bp = 29 aa + stop, starts at ATG of NcoI.
EV_GUS_REP_AFTER_ATG = (
    "CTGGATCCAGCTATCGGTGCTCCATTCGCTAAGATCGAGGCTCGTCTGAACGGTGCTATC"
    "CCAGCTGCTTACAAGAGATAA"
)  # 84 after ATG; total reporter CDS = ATG + 84 = 87? 
# NcoI is CCATGG, so ATG is already the start. Reporter remainder follows GG (the last 2 of CCATGG are the first 2 of ATG... 
# CCATGG = C + ATG + G. The first amino acid is M (ATG), next codon starts at the G plus two more bases.
# We'll define the reporter CDS as starting at the ATG inside NcoI: ATG + G + remainder.
# Let me set remainder so total CDS is 90 bp including stop.

# After NcoI CCATGG, the last G of the site is already the first base of codon 2
# (ATG GCT ...). The tail is the remainder after ATGG.
EV_GUS_REP_TAIL = (
    "CTGGATCCAGCTATCGTGCTCCATTCGCTAAGATCGAGGCTCGTCTGAACGGTGCTATC"
    "CCAGCTGCTTACAAGGAATAA"
)  # 80 nt; reporter CDS = ATGG + tail = 84 nt
# ATG is 3; tail length must be %3==0. len of string above: I'll compute at runtime.

EV_GUS_TERM = "GCGGCCGCCTGCAGGCTAGC"  # 20
EV_GUS_FILL = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"  # 66

# Build after validating tail
_rep_cds = "ATGG" + EV_GUS_REP_TAIL  # NcoI contributes ATGG; start codon is ATG at index 2 of CCATGG
assert len(_rep_cds) % 3 == 0, len(_rep_cds)
assert _rep_cds.startswith("ATGG") and _rep_cds.endswith("TAA")
assert "*" not in translate(_rep_cds)[:-1]
assert translate(_rep_cds).startswith("MAG")

EV_GUS = assemble_gene(
    EV_GUS_ORI,
    EV_GUS_MARKER,
    EV_GUS_PROM,
    EV_GUS_PRE_ECO,
    EV_GUS_ECO,
    EV_GUS_STUFFER,
    EV_GUS_NCO,
    EV_GUS_REP_TAIL,  # starts with G... continuing from CCATGG
    EV_GUS_TERM,
    EV_GUS_FILL,
)

c = 0
g_ori = span(len(EV_GUS_ORI))
g_mark = span(len(EV_GUS_MARKER))
g_prom = span(len(EV_GUS_PROM))
g_pre = span(len(EV_GUS_PRE_ECO))
g_eco = span(len(EV_GUS_ECO))
g_stu = span(len(EV_GUS_STUFFER))
g_nco = span(len(EV_GUS_NCO))
g_tail = span(len(EV_GUS_REP_TAIL))
g_term = span(len(EV_GUS_TERM))
g_fill = span(len(EV_GUS_FILL))
assert c == len(EV_GUS), (c, len(EV_GUS))

# F arm 19 nt ending with EcoRI; R arm is RC of (C of NcoI + ATGG + first of tail)
# Downstream plus-strand arm starts at the first C of CCATGG so the reporter ATG is arm nt 2-4.
F_ARM_GUS = EV_GUS[g_eco[1] - 19 : g_eco[1]]
assert F_ARM_GUS.endswith("GAATTC") and len(F_ARM_GUS) == 19
DOWN_ARM_GUS = EV_GUS[g_nco[0] : g_nco[0] + 18]  # CCATGG + 12 nt
R_ARM_GUS = rc(DOWN_ARM_GUS)
assert EV_GUS.count(F_ARM_GUS) == 1 and EV_GUS.count(DOWN_ARM_GUS) == 1

# Promoter insert: EXA1 nt 10–90 (upstream trimmed + full 5'UTR), ends before ATG.
PPRO_INSERT = EXA1_PROMOTER_INSERT
assert not PPRO_INSERT.endswith("ATG")
PPRO_PRODUCT = EV_GUS[: g_eco[1]] + PPRO_INSERT + EV_GUS[g_nco[0] :]

# UTR primer cores (gene-specific, on cDNA). Chosen from divergent UTRs.
UTR_EXA1_F_CORE = EXA1_UTR5[:22]  # 22 nt at start of 5'UTR
UTR_EXA1_R_CORE = rc(EXA1_UTR3[-22:])
UTR_EXA1L_F_CORE = EXA1L_UTR5[:22]
UTR_EXA1L_R_CORE = rc(EXA1L_UTR3[-22:])

# IF cores: F core = start of CDS (may include ATG); R core for -TAA is RC of the last 22 nt of CDS minus stop.
IF_EXA1_F_CORE = EXA1_CDS[:22]
IF_EXA1_R_CORE = rc(EXA1_CDS_NOSTOP[-22:])
IF_EXA1_F = F_ARM_FLAG + IF_EXA1_F_CORE
IF_EXA1_R = R_ARM_FLAG + IF_EXA1_R_CORE

# Promoter-GUS cores: F starts at insert 5' (which is EXA1_UP[10:]); R ends at end of 5'UTR (before ATG).
PROM_F_CORE = PPRO_INSERT[:22]
PROM_R_CORE = rc(PPRO_INSERT[-22:])
PROM_F = F_ARM_GUS + PROM_F_CORE
PROM_R = R_ARM_GUS + PROM_R_CORE


def gene_notes(symbol, accession, function, region):
    return (
        "<html><body>"
        f"symbol: {symbol}<br/>"
        f"accession: {accession}<br/>"
        "organism: Example plant<br/>"
        f"function: {function}<br/>"
        "annotation_versions: identical<br/>"
        "isoform: .1<br/>"
        f"(longest CDS {len(EXA1_CDS)} bp)<br/>"
        f"region: {region}<br/>"
        "created: 2026-11-26 16:00 +0800 (Asia/Shanghai)"
        "</body></html>"
    )


def main():
    # Sanity: CDS identity
    ident = sum(a == b for a, b in zip(EXA1_CDS, EXA1L_CDS)) / len(EXA1_CDS)
    utr_ident = sum(a == b for a, b in zip(EXA1_UTR5 + EXA1_UTR3, EXA1L_UTR5 + EXA1L_UTR3)) / 60
    print(f"CDS identity EXA1 vs EXA1L: {ident:.1%}; UTR identity: {utr_ident:.1%}")
    assert ident >= 0.85, ident
    assert utr_ident < 0.4, utr_ident

    # FASTA
    write(EX / "sequences" / "EXA1_gene.fa", wrap_fa(
        "EXA1_gene", EXA1,
        "FICTIONAL Example plant gene EXA1 (EXApl.01G000100). Genomic, plus strand. NOT a real gene.",
    ))
    write(EX / "sequences" / "EXA1L_gene.fa", wrap_fa(
        "EXA1L_gene", EXA1L,
        "FICTIONAL paralog EXA1L (EXApl.03G000200). Divergent UTRs, similar CDS. NOT a real gene.",
    ))
    write(EX / "sequences" / "EXA1_cdna.fa", wrap_fa(
        "EXA1_cdna", EXA1_CDNA, "FICTIONAL spliced cDNA logic: 5'UTR + CDS + 3'UTR.",
    ))
    write(EX / "sequences" / "EXA1L_cdna.fa", wrap_fa(
        "EXA1L_cdna", EXA1L_CDNA, "FICTIONAL spliced cDNA of the EXA1L paralog.",
    ))
    write(EX / "sequences" / "EXA1_cds.fa", wrap_fa(
        "EXA1_cds", EXA1_CDS, "FICTIONAL EXA1 CDS including the native TAA stop.",
    ))
    write(EX / "sequences" / "EXA1_cds_nostop.fa", wrap_fa(
        "EXA1_cds_nostop", EXA1_CDS_NOSTOP,
        "FICTIONAL EXA1 CDS with the native stop removed (for a C-terminal tag, -TAA primer).",
    ))
    write(EX / "sequences" / "EXA1_promoter_insert.fa", wrap_fa(
        "EXA1_promoter_insert", PPRO_INSERT,
        "FICTIONAL promoter insert: upstream[10:60] + 5'UTR, ends before the gene ATG. No ±TAA.",
    ))
    write(EX / "sequences" / "EV-FLAG.fa", wrap_fa(
        "EV-FLAG", EV_FLAG,
        f"FICTIONAL empty backbone, circular, {len(EV_FLAG)} bp. C-terminal FLAG-like tag after SpeI. NOT a real vector.",
    ))
    write(EX / "sequences" / "EV-GUS.fa", wrap_fa(
        "EV-GUS", EV_GUS,
        f"FICTIONAL empty backbone, circular, {len(EV_GUS)} bp. Transcriptional reporter after NcoI. Reporter CDS is synthetic, not uidA.",
    ))
    write(EX / "sequences" / "pOE001_EV-FLAG_EXA1.fa", wrap_fa(
        "pOE001_EV-FLAG_EXA1", POE_PRODUCT,
        "FICTIONAL in-silico In-Fusion product: EXA1 CDS (−stop) into EV-FLAG. Not a wet-lab clone.",
    ))
    write(EX / "sequences" / "pPRO001_pro-EXA1-GUS.fa", wrap_fa(
        "pPRO001_pro-EXA1-GUS", PPRO_PRODUCT,
        "FICTIONAL in-silico In-Fusion product: EXA1 promoter+5'UTR into EV-GUS. Transcriptional fusion.",
    ))

    # Feature table
    feat_rows = [
        ["file", "name", "type", "start_0based", "end_0based", "color"],
    ]
    for f in GENE_FEATURES:
        feat_rows.append(["DNA001_EXA1_EXApl.01G000100", f["name"], f["typ"], f["start"], f["end"], f["color"]])
    write(EX / "annotations" / "EXA1_features.tsv", "".join("\t".join(map(str, r)) + "\n" for r in feat_rows))

    schemes = f"""# FICTIONAL In-Fusion schemes. Arms are taken from the example backbones.
# Paths are relative to this file.
schemes:
  EVFLAG-EXA1:
    description: "C-terminal FLAG fusion. Reverse primer is -TAA (native stop removed)."
    backbone: ../sequences/EV-FLAG.fa
    insert: ../sequences/EXA1_cds_nostop.fa
    f_arm: "{F_ARM_FLAG}"
    r_arm: "{R_ARM_FLAG}"
    translate_from: "ATG"
    expect_motif: "DYKD"
    keep_native_stop: false
    linear: false
    f_compensation_nt: 0
    r_compensation_nt: 0
    reverse_primer_suffix: "-TAA"
  EVGUS-proEXA1:
    description: "Transcriptional promoter-reporter fusion. No ±TAA. Reporter starts at its own ATG."
    backbone: ../sequences/EV-GUS.fa
    insert: ../sequences/EXA1_promoter_insert.fa
    f_arm: "{F_ARM_GUS}"
    r_arm: "{R_ARM_GUS}"
    translate_from: "{EV_GUS[g_nco[0]+1:g_nco[0]+1+30]}"  # reporter ATG region; overridden below
    expect_motif: null
    keep_native_stop: true
    linear: false
    f_compensation_nt: 0
    r_compensation_nt: 0
    reverse_primer_suffix: none
"""
    # The translate_from for GUS should be the reporter ATG on the PRODUCT, which is
    # CCATGG → start at ATG (unique if we search inside reporter). We'll set ATG and
    # let the checker pick the unique ATG in the "insert or reporter" — actually the
    # promoter insert may have no ATG (we designed it that way). Reporter ATG should
    # be unique if insert has none.
    assert "ATG" not in PPRO_INSERT
    schemes = f"""# FICTIONAL In-Fusion schemes. Arms are taken from the example backbones.
# Paths are relative to this file.
schemes:
  EVFLAG-EXA1:
    description: "C-terminal FLAG fusion. Reverse primer is -TAA (native stop removed)."
    backbone: ../sequences/EV-FLAG.fa
    insert: ../sequences/EXA1_cds_nostop.fa
    f_arm: "{F_ARM_FLAG}"
    r_arm: "{R_ARM_FLAG}"
    translate_from: "ATG"
    expect_motif: "DYKD"
    keep_native_stop: false
    linear: false
    f_compensation_nt: 0
    r_compensation_nt: 0
    reverse_primer_suffix: "-TAA"
  EVGUS-proEXA1:
    description: "Transcriptional promoter-reporter fusion. No ±TAA. Reporter starts at its own ATG."
    backbone: ../sequences/EV-GUS.fa
    insert: ../sequences/EXA1_promoter_insert.fa
    f_arm: "{F_ARM_GUS}"
    r_arm: "{R_ARM_GUS}"
    translate_from: "ATG"
    expect_motif: null
    keep_native_stop: true
    linear: false
    orf_start: {g_eco[1] + len(PPRO_INSERT) + 2}
    f_compensation_nt: 0
    r_compensation_nt: 0
    reverse_primer_suffix: none
"""
    write(EX / "annotations" / "if_schemes.yaml", schemes)
    import json
    # orf_start is 0-based on the rebuilt product
    poe_orf = None  # unique insert ATG is auto-detected
    ppro_orf = g_eco[1] + len(PPRO_INSERT) + 2  # CCATGG: ATG starts at the A (index +2)
    schemes_json = {
        "schemes": {
            "EVFLAG-EXA1": {
                "description": "C-terminal FLAG fusion. Reverse primer is -TAA (native stop removed).",
                "backbone": "../sequences/EV-FLAG.fa",
                "insert": "../sequences/EXA1_cds_nostop.fa",
                "f_arm": F_ARM_FLAG,
                "r_arm": R_ARM_FLAG,
                "translate_from": "ATG",
                "expect_motif": "DYKD",
                "keep_native_stop": False,
                "linear": False,
                "f_compensation_nt": 0,
                "r_compensation_nt": 0,
                "reverse_primer_suffix": "-TAA",
            },
            "EVGUS-proEXA1": {
                "description": "Transcriptional promoter-reporter fusion. No ±TAA. Reporter starts at its own ATG.",
                "backbone": "../sequences/EV-GUS.fa",
                "insert": "../sequences/EXA1_promoter_insert.fa",
                "f_arm": F_ARM_GUS,
                "r_arm": R_ARM_GUS,
                "translate_from": "ATG",
                "expect_motif": None,
                "keep_native_stop": True,
                "linear": False,
                "orf_start": ppro_orf,
                "f_compensation_nt": 0,
                "r_compensation_nt": 0,
                "reverse_primer_suffix": "none",
            },
        }
    }
    (EX / "annotations" / "if_schemes.json").write_text(
        json.dumps(schemes_json, indent=2) + "\n", encoding="utf-8"
    )
    print("wrote", (EX / "annotations" / "if_schemes.json").relative_to(REPO_ROOT))


    # primers.csv
    primers = [
        {
            "id": "PID-000001",
            "name": "UTR-EXA1-F",
            "sequence": UTR_EXA1_F_CORE,
            "direction": "F",
            "function": "UTR",
            "target": "EXA1",
            "status": "draft",
            "synth_status": "unsynthesized",
            "vendor_id": "",
            "legacy_name": "UTR-EXA1-F",
            "notes": "cDNA UTR-anchored clone; core only; FICTIONAL",
        },
        {
            "id": "PID-000002",
            "name": "UTR-EXA1-R",
            "sequence": UTR_EXA1_R_CORE,
            "direction": "R",
            "function": "UTR",
            "target": "EXA1",
            "status": "draft",
            "synth_status": "unsynthesized",
            "vendor_id": "",
            "legacy_name": "UTR-EXA1-R",
            "notes": "cDNA UTR-anchored clone; core only; FICTIONAL",
        },
        {
            "id": "PID-000003",
            "name": "IF-EVFLAG-EXA1-F",
            "sequence": IF_EXA1_F,
            "direction": "F",
            "function": "IF",
            "target": "EXA1",
            "status": "to_order",
            "synth_status": "unsynthesized",
            "vendor_id": "",
            "legacy_name": "IF-EVFLAG-EXA1-F",
            "notes": "In-Fusion into EV-FLAG; Tm on core only; FICTIONAL",
        },
        {
            "id": "PID-000004",
            "name": "IF-EVFLAG-EXA1-R-TAA",
            "sequence": IF_EXA1_R,
            "direction": "R",
            "function": "IF",
            "target": "EXA1",
            "status": "to_order",
            "synth_status": "unsynthesized",
            "vendor_id": "",
            "legacy_name": "IF-EVFLAG-EXA1-R-TAA",
            "notes": "C-terminal tag: native stop removed. Suffix -TAA. FICTIONAL",
        },
        {
            "id": "PID-000005",
            "name": "IF-EVGUS-proEXA1-F",
            "sequence": PROM_F,
            "direction": "F",
            "function": "IF",
            "target": "proEXA1",
            "status": "draft",
            "synth_status": "unsynthesized",
            "vendor_id": "",
            "legacy_name": "IF-EVGUS-proEXA1-F",
            "notes": "Promoter-reporter IF; no ±TAA; FICTIONAL",
        },
        {
            "id": "PID-000006",
            "name": "IF-EVGUS-proEXA1-R",
            "sequence": PROM_R,
            "direction": "R",
            "function": "IF",
            "target": "proEXA1",
            "status": "draft",
            "synth_status": "unsynthesized",
            "vendor_id": "",
            "legacy_name": "IF-EVGUS-proEXA1-R",
            "notes": "Core ends before the gene ATG. Transcriptional fusion. FICTIONAL",
        },
        {
            "id": "PID-000007",
            "name": "Pr0001-UTR-EXA1L-F",
            "sequence": UTR_EXA1L_F_CORE,
            "direction": "F",
            "function": "UTR",
            "target": "EXA1L",
            "status": "ordered",
            "synth_status": "ordered",
            "vendor_id": "Pr0001",
            "legacy_name": "UTR-EXA1L-F",
            "notes": "Example of a functional name after the order was actually sent. FICTIONAL",
        },
        {
            "id": "PID-000008",
            "name": "Pr0002-UTR-EXA1L-R",
            "sequence": UTR_EXA1L_R_CORE,
            "direction": "R",
            "function": "UTR",
            "target": "EXA1L",
            "status": "ordered",
            "synth_status": "ordered",
            "vendor_id": "Pr0002",
            "legacy_name": "UTR-EXA1L-R",
            "notes": "Example of a functional name after the order was actually sent. FICTIONAL",
        },
    ]
    fields = list(primers[0].keys())
    ppath = EX / "primers" / "primers.csv"
    with ppath.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fields)
        w.writeheader()
        w.writerows(primers)
    print("wrote", ppath.relative_to(REPO_ROOT))

    # Order template + draft (CSV; also a simple xlsx)
    header_lines = [
        "Example Vendor — primer synthesis order (FICTIONAL template)",
        "Researcher,Jane Researcher",
        "Email,researcher@example.org",
        "Phone,00000000000",
        "Institute,Example University",
        "Group,Example Lab",
        "Date,2026-11-26",
        "",
    ]
    cols = ["ID", "Primer name", "Sequence 5'to3'", "Purity", "OD per tube", "Tubes", "Delivery form", "Notes"]
    # template is empty rows
    tmpl = "".join(l + "\n" for l in header_lines) + ",".join(cols) + "\n"
    write(EX / "orders" / "templates" / "primer_order_template.csv", tmpl)

    draft_rows = [
        ["1", "IF-EVFLAG-EXA1-F", IF_EXA1_F, "OPC", "2", "1", "", ""],
        ["2", "IF-EVFLAG-EXA1-R-TAA", IF_EXA1_R, "OPC", "2", "1", "", ""],
    ]
    draft = "".join(l + "\n" for l in header_lines) + ",".join(cols) + "\n"
    for r in draft_rows:
        draft += ",".join(r) + "\n"
    write(EX / "orders" / "drafts" / "primer_order_example.csv", draft)

    try:
        from openpyxl import Workbook
        from openpyxl.styles import Font

        def fill_sheet(wb_path, rows):
            wb = Workbook()
            ws = wb.active
            ws.title = "order"
            ws["A1"] = "Example Vendor — primer synthesis order (FICTIONAL)"
            ws["A1"].font = Font(bold=True)
            meta = [
                ("Researcher", "Jane Researcher"),
                ("Email", "researcher@example.org"),
                ("Phone", "00000000000"),
                ("Institute", "Example University"),
                ("Group", "Example Lab"),
                ("Date", "2026-11-26"),
            ]
            for i, (k, v) in enumerate(meta, start=2):
                ws.cell(i, 1, k)
                ws.cell(i, 2, v)
            start = 9
            for j, col in enumerate(cols, start=1):
                cell = ws.cell(start, j, col)
                cell.font = Font(bold=True)
            for i, row in enumerate(rows, start=start + 1):
                for j, val in enumerate(row, start=1):
                    ws.cell(i, j, val)
            wb.save(wb_path)

        fill_sheet(EX / "orders" / "templates" / "primer_order_template.xlsx", [])
        fill_sheet(EX / "orders" / "drafts" / "primer_order_example.xlsx", draft_rows)
        print("wrote xlsx order template and draft")
    except ImportError:
        print("openpyxl not available; CSV drafts only")

    # SnapGene .dna (fictional, built from scratch — not exported from any lab library)
    notes1 = gene_notes("EXA1", "EXApl.01G000100", "example kinase-like (fictional)", "chr1:1000-1332 (+)")
    notes2 = gene_notes("EXA1L", "EXApl.03G000200", "paralog of EXA1 (fictional)", "chr3:5000-5332 (+)")

    # Primer bind coords on genomic file: UTR-F binds 5'UTR start (60-82), UTR-R on 3'UTR end minus strand.
    try_write_dna(
        EX / "plasmids" / "genes" / "DNA001_EXA1_EXApl.01G000100.dna",
        EXA1, "linear",
        document_name="DNA001_EXA1_EXApl.01G000100",
        organism="Example plant",
        accession="EXApl.01G000100",
        notes_html=notes1,
        features=GENE_FEATURES,
        primers=[
            dict(name="UTR-EXA1-F", sequence=UTR_EXA1_F_CORE, bind_start=60, bind_end=82, bind_strand="+", annealed_bases=UTR_EXA1_F_CORE),
            dict(name="UTR-EXA1-R", sequence=UTR_EXA1_R_CORE, bind_start=281, bind_end=303, bind_strand="-", annealed_bases=UTR_EXA1_R_CORE),
        ],
    )
    try_write_dna(
        EX / "plasmids" / "genes" / "DNA002_EXA1L_EXApl.03G000200.dna",
        EXA1L, "linear",
        document_name="DNA002_EXA1L_EXApl.03G000200",
        organism="Example plant",
        accession="EXApl.03G000200",
        notes_html=notes2,
        features=GENE_FEATURES,
        primers=[
            dict(name="Pr0001-UTR-EXA1L-F", sequence=UTR_EXA1L_F_CORE, bind_start=60, bind_end=82, bind_strand="+", annealed_bases=UTR_EXA1L_F_CORE),
            dict(name="Pr0002-UTR-EXA1L-R", sequence=UTR_EXA1L_R_CORE, bind_start=281, bind_end=303, bind_strand="-", annealed_bases=UTR_EXA1L_R_CORE),
        ],
    )
    try_write_dna(
        EX / "plasmids" / "empty" / "EV-FLAG.dna",
        EV_FLAG, "circular",
        document_name="EV-FLAG",
        organism="synthetic",
        notes_html="<html><body>name: EV-FLAG<br/>purpose: C-terminal FLAG-like fusion<br/>chemistry: In-Fusion<br/>linearization: XhoI+SpeI<br/>compensation: 0/0<br/>created: 2026-11-26 (Asia/Shanghai)<br/>note: FICTIONAL backbone, not a real vector.</body></html>",
        features=[
            dict(name="synthetic_ori", typ="rep_origin", start=s_ori[0], end=s_ori[1], color=GRAY),
            dict(name="marker_example", typ="CDS", start=s_mark[0], end=s_mark[1], color=ORANGE),
            dict(name="P-synth", typ="promoter", start=s_prom[0], end=s_prom[1], color=GREEN),
            dict(name="XhoI", typ="misc_feature", start=s_xho[0], end=s_xho[1], color=TEAL),
            dict(name="SpeI", typ="misc_feature", start=s_spe[0], end=s_spe[1], color=TEAL),
            dict(name="3xFLAG_like", typ="CDS", start=s_tag[0], end=s_stp[1], color=PURPLE),
            dict(name="terminator", typ="terminator", start=s_term[0], end=s_term[1], color=GRAY),
        ],
    )
    try_write_dna(
        EX / "plasmids" / "empty" / "EV-GUS.dna",
        EV_GUS, "circular",
        document_name="EV-GUS",
        organism="synthetic",
        notes_html="<html><body>name: EV-GUS<br/>purpose: promoter-reporter transcriptional fusion<br/>chemistry: In-Fusion<br/>linearization: EcoRI+NcoI<br/>compensation: 0/0<br/>note: reporter starts at its own ATG (inside NcoI). FICTIONAL backbone. Reporter is not uidA.</body></html>",
        features=[
            dict(name="synthetic_ori", typ="rep_origin", start=g_ori[0], end=g_ori[1], color=GRAY),
            dict(name="marker_example", typ="CDS", start=g_mark[0], end=g_mark[1], color=ORANGE),
            dict(name="P-min", typ="promoter", start=g_prom[0], end=g_prom[1], color=GREEN),
            dict(name="EcoRI", typ="misc_feature", start=g_eco[0], end=g_eco[1], color=TEAL),
            dict(name="stuffer", typ="misc_feature", start=g_stu[0], end=g_stu[1], color=GRAY),
            dict(name="NcoI", typ="misc_feature", start=g_nco[0], end=g_nco[1], color=TEAL),
            dict(name="reporter_example", typ="CDS", start=g_nco[0] + 2, end=g_tail[1], color=ORANGE),
            dict(name="terminator", typ="terminator", start=g_term[0], end=g_term[1], color=GRAY),
        ],
    )
    # Construct features (approximate)
    poe_insert_start = s_xho[1]
    poe_insert_end = s_xho[1] + len(POE_INSERT)
    poe_tag_start = poe_insert_end + (s_tag[0] - s_spe[0])
    poe_tag_end = poe_tag_start + len(EV_FLAG_TAG) + 3
    try_write_dna(
        EX / "plasmids" / "constructs" / "pOE001_EV-FLAG_EXA1.dna",
        POE_PRODUCT, "circular",
        document_name="pOE001_EV-FLAG_EXA1",
        organism="synthetic",
        notes_html="<html><body>name: pOE001_EV-FLAG_EXA1<br/>type: in-silico design file (not a wet-lab clone or glycerol stock)<br/>insert: EXA1 CDS minus stop<br/>backbone: EV-FLAG<br/>frame: M-[EXA1]-TSGG-DYKDHDG… (verify with infusion_frame_check.py)<br/>wet-lab clone: unknown / unrecorded</body></html>",
        features=[
            dict(name="P-synth", typ="promoter", start=s_prom[0], end=s_prom[1], color=GREEN),
            dict(name="EXA1", typ="CDS", start=poe_insert_start, end=poe_insert_end, color=ORANGE),
            dict(name="3xFLAG_like", typ="CDS", start=poe_tag_start, end=poe_tag_end, color=PURPLE),
        ],
        primers=[
            dict(name="IF-EVFLAG-EXA1-F", sequence=IF_EXA1_F, bind_start=poe_insert_start, bind_end=poe_insert_start + 22, bind_strand="+", annealed_bases=IF_EXA1_F_CORE),
        ],
    )
    ppro_ins_start = g_eco[1]
    ppro_ins_end = g_eco[1] + len(PPRO_INSERT)
    ppro_rep_start = ppro_ins_end + 2  # CCATGG: ATG starts at the A
    ppro_rep_end = ppro_rep_start + len("ATG" + EV_GUS_REP_TAIL)
    try_write_dna(
        EX / "plasmids" / "constructs" / "pPRO001_pro-EXA1-GUS.dna",
        PPRO_PRODUCT, "circular",
        document_name="pPRO001_pro-EXA1-GUS",
        organism="synthetic",
        notes_html="<html><body>name: pPRO001_pro-EXA1-GUS<br/>type: in-silico design file (not a wet-lab clone)<br/>insert: EXA1 promoter + 5'UTR, ends before gene ATG<br/>backbone: EV-GUS<br/>fusion: transcriptional; reporter uses its own ATG<br/>wet-lab clone: unknown / unrecorded</body></html>",
        features=[
            dict(name="upstream", typ="misc_feature", start=ppro_ins_start, end=ppro_ins_start + len(EXA1_UP[10:]), color=GRAY),
            dict(name="5'UTR", typ="5'UTR", start=ppro_ins_start + len(EXA1_UP[10:]), end=ppro_ins_end, color=GREEN),
            dict(name="reporter_example", typ="CDS", start=ppro_rep_start, end=ppro_rep_end, color=ORANGE),
        ],
    )

    print("EXA1", len(EXA1), "EV-FLAG", len(EV_FLAG), "EV-GUS", len(EV_GUS))
    print("F_ARM_FLAG", F_ARM_FLAG, "R_ARM_FLAG", R_ARM_FLAG)
    print("F_ARM_GUS", F_ARM_GUS, "R_ARM_GUS", R_ARM_GUS)
    print("POE", len(POE_PRODUCT), "PPRO", len(PPRO_PRODUCT))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Rebuild an In-Fusion product from the actual backbone sequence and check tag frame.

Homology arms must be substrings of the backbone (F arm on the plus strand,
R arm as the reverse-complement of the plus-strand downstream arm). Recombination
rebuilds those arm sequences. Do not compute phase from sticky-end cut offsets
(a common off-by-one that can introduce a premature stop at a rebuilt SpeI).

If the junction or the downstream tag frame cannot be verified, this script
exits non-zero and says so. It never invents a translation.

Examples:
  python3 infusion_frame_check.py --scheme ../../examples/annotations/if_schemes.yaml --id EVFLAG-EXA1
  python3 infusion_frame_check.py --backbone ev.fa --insert insert.fa --f-arm CTCGAG... --r-arm ...
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def rc(seq: str) -> str:
    return seq.translate(COMPLEMENT)[::-1]


def read_fasta(path: Path) -> tuple[str, str]:
    name = path.stem
    chunks: list[str] = []
    with path.open(encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                name = line[1:].split()[0]
                continue
            chunks.append(line.replace(" ", "").replace("\t", ""))
    seq = "".join(chunks).upper()
    if not seq:
        sys.exit(f"empty sequence: {path}")
    if any(ch not in "ACGTN" for ch in seq):
        bad = sorted({ch for ch in seq if ch not in "ACGTN"})
        sys.exit(f"non-ACGTN bases in {path}: {bad}")
    return name, seq


def translate(dna: str) -> str:
    try:
        from Bio.Seq import Seq

        return str(Seq(dna).translate(to_stop=False))
    except ImportError:
        table = {
            "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
            "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
            "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
            "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
            "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
            "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
            "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
            "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
            "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
            "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
            "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
            "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
            "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
            "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
            "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
            "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
        }
        aa = []
        for i in range(0, len(dna) - len(dna) % 3, 3):
            aa.append(table.get(dna[i : i + 3], "X"))
        return "".join(aa)


def find_unique(haystack: str, needle: str, label: str) -> int:
    idx = haystack.find(needle)
    if idx < 0:
        sys.exit(f"cannot verify: {label} not found on the backbone. Stop; do not invent a junction.")
    if haystack.find(needle, idx + 1) >= 0:
        sys.exit(f"cannot verify: {label} is not unique on the backbone. Stop; do not invent a junction.")
    return idx


def rebuild(backbone: str, f_arm: str, r_arm: str, insert: str, circular: bool = True) -> str:
    """F arm is plus-strand; R arm is the synthesized reverse primer 5'→3' (RC of downstream arm)."""
    f_arm = f_arm.upper()
    r_arm = r_arm.upper()
    insert = insert.upper()
    down = rc(r_arm)
    f_at = find_unique(backbone, f_arm, "F homology arm")
    d_at = find_unique(backbone, down, "plus-strand sequence complementary to the R arm")
    f_end = f_at + len(f_arm)
    if circular and d_at < f_end:
        # wrap: keep from d_at? no: product is backbone[:f_end] wait that's wrong
        # Layout on a circle: ...[F arm][dropped stuffer][downstream arm]...
        # If d_at < f_end, the window wraps around the origin. Unusual for MCS.
        sys.exit(
            "cannot verify: downstream arm starts before the F arm end on a linear coordinate. "
            "Check circular origin / arm orientation; do not invent a junction."
        )
    dropped = backbone[f_end:d_at]
    product = backbone[:f_end] + insert + backbone[d_at:]
    if not circular:
        return product
    # already a full circle representation with origin unchanged
    _ = dropped
    return product


def first_in_frame_stop(aa: str) -> int | None:
    i = aa.find("*")
    return None if i < 0 else i


def run_check(
    *,
    backbone: str,
    insert: str,
    f_arm: str,
    r_arm: str,
    translate_from: str,
    expect_motif: str | None,
    keep_native_stop: bool,
    circular: bool,
    title: str,
    orf_start: int | None = None,
) -> int:
    product = rebuild(backbone, f_arm, r_arm, insert, circular=circular)
    print(f"== {title} ==")
    print(f"backbone {len(backbone)} bp; insert {len(insert)} bp; product {len(product)} bp")
    print(f"F arm {f_arm} ({len(f_arm)} nt)")
    print(f"R arm {r_arm} ({len(r_arm)} nt); plus-strand downstream arm {rc(r_arm)}")

    query = translate_from.upper()
    if orf_start is not None:
        if orf_start < 0 or orf_start >= len(product):
            print(f"FAIL: orf_start {orf_start} is outside the product (0..{len(product)-1}).")
            return 1
        if product[orf_start:orf_start + len(query)] != query:
            print(
                f"FAIL: product does not have {query!r} at orf_start {orf_start} "
                f"(found {product[orf_start:orf_start + len(query)]!r}). Cannot verify frame."
            )
            return 1
        start = orf_start
        print(f"using declared orf_start {orf_start}")
    else:
        start = product.find(query)
        if start < 0:
            print("FAIL: translation start motif not found on the product. Cannot verify frame.")
            return 1
        if product.find(query, start + 1) >= 0:
            f_at = product.find(f_arm.upper())
            f_end = f_at + len(f_arm)
            insert_end = f_end + len(insert)
            starts = []
            at = -1
            while True:
                at = product.find(query, at + 1)
                if at < 0:
                    break
                starts.append(at)
            in_insert = [s for s in starts if f_end <= s < insert_end]
            if query == "ATG" and len(in_insert) == 1:
                start = in_insert[0]
                print(f"note: start motif occurs {len(starts)} times; using the unique ATG inside the insert (0-based {start})")
            else:
                print(
                    f"FAIL: translation start motif is not unique ({len(starts)} hits at {starts}). "
                    "Pass --orf-start (0-based on the product) rather than guessing."
                )
                return 1

    window = product[start:]
    if len(window) < 3:
        print("FAIL: not enough sequence after the start motif to translate.")
        return 1
    window = window[: len(window) - len(window) % 3]
    aa = translate(window)
    stop_at = first_in_frame_stop(aa)
    shown = aa if stop_at is None else aa[: stop_at + 1]
    print(f"translation start at 0-based {start}")
    print(f"protein (through first stop): {shown}")

    if keep_native_stop:
        if stop_at is None:
            print("FAIL: expected the insert to provide a stop codon, but none was found in frame.")
            return 1
        if stop_at == 0:
            print("FAIL: in-frame stop immediately at the start motif.")
            return 1
        print(f"OK: first in-frame stop at amino-acid position {stop_at + 1} (1-based).")
    else:
        # C-terminal fusion: no stop until after the tag
        if stop_at is None:
            print("FAIL: no in-frame stop found after the fusion (tag + linker should end with a stop).")
            return 1
        if expect_motif and expect_motif.upper() not in aa.upper().replace("*", ""):
            # also allow motif with stop stripped
            if expect_motif.upper() not in aa.upper():
                print(
                    f"FAIL: expected tag motif {expect_motif!r} was not found in the translation. "
                    "Frame may be wrong. Stop; do not invent a correction."
                )
                return 1
        # motif should appear before the first stop
        prefix = aa[: stop_at + 1]
        if expect_motif and expect_motif.upper() not in prefix.upper():
            print(
                f"FAIL: tag motif {expect_motif!r} is not before the first stop. "
                "Possible frameshift or premature stop."
            )
            return 1
        if expect_motif:
            print(f"OK: tag motif {expect_motif} occurs before the first stop (aa {stop_at + 1}, 1-based).")
        else:
            print(f"OK: open reading frame extends to amino-acid position {stop_at + 1} (1-based).")
    return 0


def load_scheme(path: Path, scheme_id: str) -> dict:
    import json
    text = path.read_text(encoding="utf-8")
    if path.suffix.lower() == ".json":
        data = json.loads(text)
    else:
        # Minimal loader for the example YAML (no PyYAML dependency).
        try:
            import yaml  # type: ignore
            data = yaml.safe_load(text) or {}
        except ImportError:
            json_path = path.with_suffix(".json")
            if json_path.exists():
                data = json.loads(json_path.read_text(encoding="utf-8"))
            else:
                sys.exit(
                    f"cannot parse {path.name}: install PyYAML, or use the sibling .json scheme file."
                )
    schemes = data.get("schemes") or data
    if scheme_id not in schemes:
        sys.exit(f"scheme {scheme_id!r} not in {path}. Known: {sorted(schemes)}")
    return schemes[scheme_id]


def resolve(path_value: str, base: Path) -> Path:
    p = Path(path_value)
    return p if p.is_absolute() else (base / p).resolve()


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--scheme", help="YAML file with named IF schemes")
    parser.add_argument("--id", help="scheme id inside --scheme")
    parser.add_argument("--backbone")
    parser.add_argument("--insert")
    parser.add_argument("--f-arm")
    parser.add_argument("--r-arm")
    parser.add_argument("--translate-from", default="ATG", help="DNA motif at which to start translation")
    parser.add_argument("--expect-motif", help="amino-acid motif that must appear in-frame (e.g. DYKD)")
    parser.add_argument("--keep-native-stop", action="store_true")
    parser.add_argument("--linear", action="store_true")
    parser.add_argument("--orf-start", type=int, help="0-based start of the ORF on the rebuilt product")
    args = parser.parse_args(argv)

    title = "manual"
    if args.scheme:
        if not args.id:
            sys.exit("--id is required with --scheme")
        scheme = load_scheme(Path(args.scheme), args.id)
        base = Path(args.scheme).resolve().parent
        backbone_path = resolve(scheme["backbone"], base)
        insert_path = resolve(scheme["insert"], base)
        f_arm = scheme["f_arm"]
        r_arm = scheme["r_arm"]
        translate_from = scheme.get("translate_from", args.translate_from)
        expect_motif = scheme.get("expect_motif", args.expect_motif)
        keep_native_stop = scheme.get("keep_native_stop", args.keep_native_stop)
        circular = not scheme.get("linear", args.linear)
        orf_start = scheme.get("orf_start", args.orf_start)
        title = args.id
    else:
        missing = [n for n in ("backbone", "insert", "f_arm", "r_arm") if not getattr(args, n.replace("-", "_"))]
        if missing:
            sys.exit("provide --scheme/--id or all of --backbone --insert --f-arm --r-arm")
        backbone_path = Path(args.backbone)
        insert_path = Path(args.insert)
        f_arm, r_arm = args.f_arm, args.r_arm
        translate_from = args.translate_from
        expect_motif = args.expect_motif
        keep_native_stop = args.keep_native_stop
        circular = not args.linear
        orf_start = args.orf_start

    _, backbone = read_fasta(backbone_path)
    _, insert = read_fasta(insert_path)
    for label, arm in (("F", f_arm), ("R", r_arm)):
        if any(ch not in "ACGTNacgtn" for ch in arm):
            sys.exit(f"non-ACGTN bases in {label} arm")
    return run_check(
        backbone=backbone,
        insert=insert,
        f_arm=f_arm,
        r_arm=r_arm,
        translate_from=translate_from,
        expect_motif=expect_motif,
        keep_native_stop=keep_native_stop,
        circular=circular,
        title=title,
        orf_start=orf_start,
    )


if __name__ == "__main__":
    sys.exit(main())

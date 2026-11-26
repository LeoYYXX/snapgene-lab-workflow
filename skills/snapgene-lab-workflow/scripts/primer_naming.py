#!/usr/bin/env python3
"""Generic primer naming helper.

Unsynthesized primers use a functional name:
    {function}-{gene}-{position}
    e.g. UTR-EXA1-F
         IF-EVFLAG-EXA1-R-TAA
         IF-EVGUS-proEXA1-F

After a vendor order is actually sent, prefix a vendor serial:
    Pr0001-UTR-EXA1-F

The primer table keeps a stable internal id (default pattern PID-000001).
Ids and vendor serials are assigned once and never reused.

Usage, from a work root that contains primers/primers.csv:
    python3 primer_naming.py apply
    python3 primer_naming.py mark-synth PID-000003 PID-000004
    python3 primer_naming.py --root /path/to/work apply

mark-synth must run only after the user confirms the order was sent.
A draft spreadsheet is not a sent order.

When names change, matching primer features in plasmids/**/*.dna are updated
via sgffp (preserve=True) and sgffp_primer_write.ensure_primer_display_defaults.
Never raw-byte-replace inside a SnapGene package. Each modified .dna is copied
to plasmids/_archive/ first.
"""
from __future__ import annotations

import argparse
import csv
import glob
import re
import shutil
import sys
import time
from pathlib import Path

FUNCS = {"UTR", "CDS", "IF", "QPCR", "SEQ", "MUTANT", "CRISPR", "USER", "ID", "PROM"}
VENDOR_PREFIX_RE = re.compile(r"^Pr(\d{4})-")
ID_RE = re.compile(r"^PID-(\d+)$", re.I)


def work_root(cli_root: str | None) -> Path:
    return Path(cli_root).resolve() if cli_root else Path.cwd()


def csv_path(root: Path) -> Path:
    return root / "primers" / "primers.csv"


def base_name(name: str) -> str:
    """Strip a vendor serial if present; keep the functional name including ±TAA."""
    name = VENDOR_PREFIX_RE.sub("", name)
    parts = name.split("-")
    if parts and parts[0].upper() in FUNCS:
        return name
    # Legacy order {gene}-{function}-{F|R} → {function}-{gene}-{F|R}
    if (
        len(parts) >= 3
        and parts[1].upper() in FUNCS
        and parts[2].upper().split("+")[0].split("-")[0] in {"F", "R"}
        or (len(parts) >= 3 and parts[1].upper() in FUNCS and parts[2].upper() in {"F", "R"})
    ):
        rest = "-".join(parts[2:])
        return f"{parts[1].upper() if parts[1].upper() in FUNCS else parts[1]}-{parts[0]}-{rest}"
    return name


def load(path: Path):
    if not path.exists():
        sys.exit(f"primer table not found: {path}")
    with path.open(encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    if not rows:
        sys.exit(f"primer table is empty: {path}")
    fields = list(rows[0].keys())
    for key in ("legacy_name", "synth_status", "vendor_id"):
        if key not in fields:
            fields.append(key)
    for row in rows:
        row.setdefault("legacy_name", "")
        row["legacy_name"] = row["legacy_name"] or row.get("name", "")
        if not row.get("synth_status"):
            row["synth_status"] = "ordered" if row.get("status") == "ordered" else "unsynthesized"
        row.setdefault("vendor_id", "")
    return rows, fields


def next_vendor_serial(rows) -> int:
    nums = []
    for row in rows:
        vid = row.get("vendor_id") or ""
        m = re.match(r"^Pr(\d{4})$", vid)
        if m:
            nums.append(int(m.group(1)))
    return max(nums, default=0) + 1


def is_frozen(row, frozen_before: int) -> bool:
    """Optional: freeze names of internal ids below a cutoff (default: freeze none)."""
    if frozen_before <= 0:
        return False
    m = ID_RE.match(row.get("id") or "")
    if not m:
        return False
    return int(m.group(1)) < frozen_before


def rename(rows, frozen_before: int):
    changes = []
    for row in rows:
        if is_frozen(row, frozen_before):
            continue
        b = base_name(row["name"])
        new = f"{row['vendor_id']}-{b}" if row.get("vendor_id") else b
        if new != row["name"]:
            changes.append((row["name"], new))
            row["name"] = new
    return changes


def dna_paths(root: Path):
    patterns = [
        "plasmids/genes/*.dna",
        "plasmids/empty/*.dna",
        "plasmids/constructs/*.dna",
        "plasmids/constructs/*/*.dna",
    ]
    out = []
    for pat in patterns:
        out.extend(glob.glob(str(root / pat)))
    return out


def sync_dna(root: Path, changes):
    """Rename primer features via sgffp. Packet lengths stay valid. Never raw-byte-replace."""
    if not changes:
        return
    try:
        from sgffp import SgffReader, SgffWriter

        scripts = str(Path(__file__).resolve().parent)
        if scripts not in sys.path:
            sys.path.insert(0, scripts)
        from sgffp_primer_write import ensure_primer_display_defaults
    except ImportError as exc:
        print("sgffp not available; skip .dna sync:", exc)
        return

    mapping = dict(changes)
    ts = time.strftime("%Y%m%d_%H%M%S")
    archive = root / "plasmids" / "_archive"
    archive.mkdir(parents=True, exist_ok=True)

    for path in dna_paths(root):
        sgff = SgffReader.from_file(path)
        _ = sgff.primers.items
        hit = []
        for primer in sgff.primers.items:
            if primer.name in mapping:
                hit.append((primer.name, mapping[primer.name]))
                primer.name = mapping[primer.name]
        if not hit:
            continue
        backup = archive / f"{Path(path).stem}_{ts}.dna"
        shutil.copy2(path, backup)
        ensure_primer_display_defaults(sgff)
        SgffWriter(path).write(sgff, preserve=True)
        print("dna synced", path, [a for a, _ in hit], "backup", backup)


def save(path: Path, rows, fields):
    bak_dir = path.parent
    bak = bak_dir / f"_bak_primers_{time.strftime('%Y%m%d_%H%M%S')}.csv"
    if path.exists():
        shutil.copy2(path, bak)
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fields)
        writer.writeheader()
        writer.writerows(rows)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--root", help="work root that contains primers/primers.csv (default: cwd)")
    parser.add_argument(
        "--frozen-before-id",
        type=int,
        default=0,
        help="do not rename internal ids below this number (0 = freeze none)",
    )
    parser.add_argument("cmd", nargs="?", default="apply", choices=["apply", "mark-synth"])
    parser.add_argument("ids", nargs="*", help="internal ids for mark-synth")
    args = parser.parse_args(argv)

    root = work_root(args.root)
    path = csv_path(root)
    rows, fields = load(path)

    if args.cmd == "apply":
        pending = [
            row
            for row in rows
            if not is_frozen(row, args.frozen_before_id)
            and row.get("synth_status") == "ordered"
            and not row.get("vendor_id")
        ]
        for row in pending:
            row["vendor_id"] = f"Pr{next_vendor_serial(rows):04d}"
    elif args.cmd == "mark-synth":
        if not args.ids:
            sys.exit("mark-synth requires one or more internal ids, e.g. PID-000003")
        wanted = set(args.ids)
        for row in rows:
            if row.get("id") in wanted and not is_frozen(row, args.frozen_before_id) and not row.get("vendor_id"):
                row["synth_status"] = "ordered"
                row["status"] = "ordered"
                row["vendor_id"] = f"Pr{next_vendor_serial(rows):04d}"
        missing = wanted - {row.get("id") for row in rows}
        if missing:
            sys.exit(f"unknown primer id(s): {sorted(missing)}")

    changes = rename(rows, args.frozen_before_id)
    save(path, rows, fields)
    sync_dna(root, changes)
    if changes:
        for old, new in changes:
            print(f"{old} -> {new}")
    else:
        print("no name changes")


if __name__ == "__main__":
    main()

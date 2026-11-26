# Fictional example project

All sequences, names, and contacts in this folder are **synthetic**. They are not from a real organism, plasmid library, or person.

本目录全部为**虚构**示例，不是真实物种、质粒库或个人数据。

## Layout

```
examples/
  primers/primers.csv
  orders/templates/          empty vendor sheet (CSV + xlsx)
  orders/drafts/             filled IF-EVFLAG-EXA1 pair (OPC, 2 OD, 1 tube; delivery/notes empty)
  sequences/                 FASTA for every molecule
  annotations/               feature table + In-Fusion schemes
  plasmids/genes/            DNA001 EXA1, DNA002 EXA1L
  plasmids/empty/            EV-FLAG, EV-GUS
  plasmids/constructs/       pOE001_EV-FLAG_EXA1, pPRO001_pro-EXA1-GUS
  seq/runs/RUN000_example/   folder layout only (no real traces)
```

The `.dna` files were built from scratch with `sgffp.SgffObject.new`, not exported from any lab library. Regenerate everything with:

```bash
python3 skills/snapgene-lab-workflow/scripts/build_examples.py
```

(run from the repository root; `build_examples.py` never reads a live lab tree).

## Identifiers

| kind | example | notes |
|------|---------|--------|
| Organism | Example plant | fictional |
| Genes | EXA1 (`EXApl.01G000100`), EXA1L (`EXApl.03G000200`) | paralog pair |
| Researcher | Jane Researcher | `researcher@example.org`, phone `00000000000`, Example University |
| Primer ids | `PID-000001` … | not a real inventory |
| Vendor serials | `Pr0001`, `Pr0002` already “sent” on the EXA1L UTR pair (illustrates naming after a real send) | |
| Vendor | Example Vendor | `orders@example-vendor.com` |

See `project.md` for a short walk-through of the two designs.

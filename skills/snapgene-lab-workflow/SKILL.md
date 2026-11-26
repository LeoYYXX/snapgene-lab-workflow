---
name: snapgene-lab-workflow
description: "Use when doing SnapGene-based molecular cloning — gene .dna workfiles, primer design and naming, In-Fusion reading-frame checks, promoter–reporter fusions, vendor primer-order drafts, or Sanger write-back. Never send vendor email until the user explicitly confirms. 在用 SnapGene 做基因工作文件、引物设计与命名、In-Fusion 读框、启动子-报告基因、引物订购草稿或测序回写时使用；未经用户明确确认不得代发供应商邮件。"
---

# SnapGene lab workflow

A reusable workflow for AI agents that read and write SnapGene files, design cloning primers, and prepare vendor orders. This skill is generic. Lab-specific inventories, contacts, and sequences live in the user's work root, not in this recipe.

Companion files (relative to this skill directory):

- `scripts/sgffp_primer_write.py` — required display defaults when writing primers
- `scripts/primer_naming.py` — stable internal ids + vendor serials
- `scripts/infusion_frame_check.py` — rebuild a product from backbone arms and check tag frame
- `config/` — templates (copy and edit; do not commit real contacts)
- `../../examples/` — fictional Example plant project (EXA1 / EXA1L)

## Hard rules / 硬规则

These are the point of the skill. Follow them even if a user request is informal.

1. **Prefer `sgffp` to read and write** SnapGene `.dna` / `.rna` / `.prot`. **Never raw-byte-replace** inside a SnapGene package. When rewriting primer names or features, use `preserve=True` and apply display defaults (`ensure_primer_display_defaults`).
2. **Backup** a `.dna` file to `plasmids/_archive/` before modifying it.
3. **Do not invent biology.** If a reading frame, junction, or homology arm cannot be verified against the actual backbone sequence, stop and say so.
4. **Orders: draft, then wait.** Fill a vendor sheet and draft the email. **Never send vendor email** (and never call a send API) until the user explicitly confirms. A draft spreadsheet is not a sent order.
5. **Constructs are in-silico design files**, distinct from wet-lab clones or glycerol stocks. Missing wet-lab tables are **unknown / unrecorded**, not “not done”.
6. Internal file and primer **IDs are assigned once and never reused.**

## Directory convention / 目录约定

Use these generic names at the lab work root (the example project mirrors them under `examples/`):

```
plasmids/genes/          gene genomic workfiles
plasmids/empty/          empty backbones
plasmids/constructs/     in-silico constructs
plasmids/_archive/       timestamped backups taken before a write
primers/                 primer table (e.g. primers.csv)
orders/templates/        unused vendor templates (do not overwrite)
orders/drafts/           filled sheets waiting for confirmation
orders/placed/           copies kept after an order was actually sent
seq/inbox/               incoming .ab1 zips
seq/runs/                one folder per sequencing run
logs/
config/
docs/
```

## SnapGene I/O / 读写规则

- Install `sgffp>=0.22.1`. Create new files with `SgffObject.new(sequence, topology=...)`. Write with `SgffWriter(path).write(sgff, preserve=True)` when updating an existing file.
- New files may be written with `preserve=False`.
- After changing a history tree, call the sgffp history sync (for example `history._sync_tree()`) before writing, so document name and history stay consistent.
- **Primer display defaults** (write these on the Primers wrapper or primers are hidden on the map until the user loosens “minimum continuous match” to 8):

```yaml
HybridizationParams:
  minContinuousMatchLen: "8"
  allowMismatch: "1"
  minMeltingTemperature: "40"
  showAdditionalFivePrimeMatches: "1"
  minimumFivePrimeAnnealing: "15"
```

Each `BindingSite` also needs `Component.hybridizedRange` equal to its 1-based inclusive `location`, plus `Component.bases` = annealed bases. Use `scripts/sgffp_primer_write.py`.

## Gene workfiles / 基因工作文件

Naming template (configure in `config/naming.yaml`):

```
DNA{NNN}_{symbol}_{accession}.dna
```

Example (fictional): `DNA001_EXA1_EXApl.01G000100.dna`. Numbers start at 001 in a new lab, increment per file, and are never reused.

Sequence: usually the genomic locus, plus-strand 5′→3′ to the right (reverse-complement a minus-strand locus before writing). Annotate **one isoform** (the longest CDS, unless the user names another), feature by feature:

`upstream? → 5′UTR → CDS1 → intron1 → … → 3′UTR → downstream?`

Colors (write on the feature **and** each segment; otherwise SnapGene falls back to type defaults):

| feature | color | 说明 |
|---------|-------|------|
| 5′UTR | `#008000` | green 绿 |
| all CDS | `#ff9900` | orange 橙 |
| all intron, upstream, downstream | `#808080` | gray 灰 |
| 3′UTR | `#3366ff` | blue 蓝 |

Description (Info / Notes) is HTML, fields separated by `<br/>`, `key: value` with a space after the colon, only fields that exist, in the configured order. Content must belong to **this** gene. Suggested fields: symbol, accession, organism, function, annotation_versions, isoform (plus a longest-CDS length line), region, created (include the lab timezone label).

Accession in Notes = the chosen gene id. Type is typically `Natural` for genomic workfiles and `Synthetic` for designed constructs.

If the lab maintains **two annotation versions**, cross-check the longest CDS before designing primers. Do not assume two identifiers that merely look similar are the same gene. If they match, keep one file and note `annotation_versions: identical`. If they diverge, keep two files (primary + a suffix for the other version), design from the agreed primary assembly, and still use the other sequence as an alignment check. This skill does not ship any species-specific mapping table.

These genomic colors do **not** automatically apply to recombinant constructs that have no UTR/CDS/intron partition.

## Primer naming / 引物命名

Keep a **stable internal id** (template `PID-000001`; a lab may choose another prefix in `naming.yaml`). Never reuse an id.

- **Unsynthesized:** a functional name such as `{function}-{gene}-{position}`.
  Examples: `UTR-EXA1-F`, `IF-EVFLAG-EXA1-R-TAA`, `IF-EVGUS-proEXA1-F`.
- **After a vendor order is actually sent:** prefix a vendor serial, e.g. `Pr0001-UTR-EXA1L-F`. Serials increment in send order and never change.
- Known function tokens: `UTR`, `CDS`, `IF`, `QPCR`, `SEQ`, `MUTANT`, `CRISPR`, `USER`, `ID`, `PROM`.
- Optional: a lab may freeze names of ids below a cutoff (`primer_frozen_before_id`; default 0 = freeze none). Do not invent a cutoff from another lab's history.

`scripts/primer_naming.py apply` normalizes names. `mark-synth PID-…` assigns the next serial **only after the user confirms the order was sent**, then syncs matching primer features in `.dna` files (backup first, `preserve=True`, display defaults on).

## UTR–cDNA cloning / UTR 锚定克隆

Default path when the goal is a full-length CDS clone from mixed or tissue cDNA:

1. Open the genomic workfile and splice a logic cDNA: `5′UTR + all CDS segments + 3′UTR`. Wet-lab template is **cDNA**, not gDNA (genomic PCR would include introns).
2. **Hard success criterion:** the amplicon contains the **complete longest CDS** (ATG through the stop). UTR flanks are optional and may be 0 nt.
3. Prefer F in the 5′UTR and R in the 3′UTR. Primers **may cross ATG or Stop** as long as ATG and Stop remain inside the product. Never truncate the CDS.
4. Tm / GC windows: start at 57–63 °C (opt 60, ΔTm ≤ 3), GC 40–60%, core 18–28 nt (opt 22). If an In-Fusion arm is added, **Tm is computed on the core only**.
5. Specificity: search a transcript / cDNA library first, then review a genomic library for gDNA-contamination products. NCBI spot-checks are off by default unless the user asks.

## Paralogs / 旁系同源

If mixed cDNA would co-amplify close paralogs (high CDS identity, weak UTR identity is common in polyploid or recently duplicated plants):

- **Require a UTR-specific amplicon**, clone, and **sequence to confirm copy identity** before designing a CDS In-Fusion pair that sits in the conserved CDS.
- Put distinguishing mismatches at the primer **3′ end** (last 1–3 nt, ideally the last base). Strong mismatches (A:A, A:G, G:A, C:C) outrank wobble pairs (G:T / T:G). Two differences, or both primers specific, are more reliable than one internal SNP.
- Default preference for cloning: at least one primer (often a UTR primer) is copy-specific. qPCR should require both ends to be specific.
- If a strictly specific pair cannot be designed, a **non-default downgrade** may allow possible paralog amplicons, but only when the report flags `degraded_specificity` and lists possible products, expected lengths, and a follow-up way to tell copies apart (Sanger, diagnostic digest, clone genotyping). Do not silently treat a relaxed pair as target-specific.
- The fictional pair **EXA1 / EXA1L** in `examples/` illustrates the pattern: CDS identity is high (~93%), UTR identity is low (~13%). UTR primers are copy-specific; a primer pair placed only in the conserved CDS would not be trusted on mixed cDNA.

## In-Fusion arms and reading frame / 同源臂与读框

- Homology arms are **rebuilt from the actual backbone sequence** around the chosen window (typically 15–20 bp at the 5′ end of each insert primer). F arm = plus-strand sequence that remains upstream of the insert; R arm = the **reverse complement** of the plus-strand sequence that remains downstream. If the arms include a complete restriction site, In-Fusion **rebuilds that site**.
- **Do not** compute phase from sticky-end cut offsets (a “cut +1” compensation). That model can invent a spurious extra base and a premature stop (for example a rebuilt SpeI becoming `TAG`).
- Re-check the frame through any **C-terminal tag**. Simulate the product (or run `scripts/infusion_frame_check.py`) and translate from the insert ATG through the tag. If the junction cannot be verified, stop.
- Compensation nucleotides, when truly required, sit **between** the homology arm and the gene-specific core, and must be justified by the rebuilt sequence, not by a sticky-end table.
- Name the reverse primer with **`+TAA` or `-TAA`** so the stop-codon decision is visible:
  - **`-TAA`**: a downstream C-terminal fusion tag is in frame (FLAG, HA, Myc, fluorescent protein, TurboID, C-His, …). **Remove** the native stop from the insert / reverse core.
  - **`+TAA`**: no C-terminal tag, or the protein must stop at the GOI (untagged OE, N-terminal fusions, etc.). **Keep or add** a stop (TAA / TAG / TGA).
- The forward primer does **not** take a ±TAA suffix.
- **Exception — transcriptional promoter–reporter fusions** (e.g. promoter–GUS / promoter–REP): do **not** add ±TAA. The reporter starts from **its own ATG**. The reverse core usually ends **before** the gene's ATG so a second start codon is not placed just upstream of the reporter. This is not a protein fusion; do not reuse a CDS-tag arm set on a reporter backbone.

## Empty backbones: enzyme choice / 空载选酶

Goal: a unique, directional insertion window that does not break resistance markers, origins, or T-DNA borders (unless replacing that marker is the point of the construct).

Prefer a **double digest** (two different enzymes, ideally one cut each, or a second cut only inside a discarded stuffer). Avoid single-cut + In-Fusion unless the arms already fix direction and the user accepts the risk. Do not cut up a C-terminal tag (the 3′ enzyme should stay upstream of the tag ATG) or an N-terminal tag (the 5′ enzyme should stay downstream of the tag). USER backbones use their native chemistry (often PacI + Nt.BbvCI); do not mix USER tails and In-Fusion arms as if they were one scheme.

Record, per backbone, after verifying against the `.dna` sequence: chemistry, enzyme pair, F/R arms taken from the backbone, compensation nt, and whether the reverse primer is +TAA, -TAA, or neither. **Do not copy another lab's empty-vector inventory into this skill.**

## Promoter–reporter designs / 启动子-报告基因

Typical path (see `examples/` `EV-GUS` + `pPRO001_pro-EXA1-GUS`):

1. Extend the gene workfile with an upstream promoter region (often ~1–3 kb, species-dependent). Keep the same DNA number; do not assign a new one. Color the new upstream gray.
2. The PCR insert is often an inner slice: the F core need not start at base 1 of the added flank. Record workfile coordinates so “the file has 3 kb” is not confused with “the construct has 3 kb”.
3. Rebuild arms from that reporter backbone. Check that the reporter ATG is still present and that the insert did not introduce a second ATG immediately upstream. NcoI-style sites (`CCATGG`) often contribute the reporter ATG; confirm the start coordinate rather than assuming “ATG is at offset +1”.
4. Internal cut sites in a promoter do not block In-Fusion, but they can break a diagnostic digest — pick identification enzymes accordingly.
5. If the 5′UTR is included, uORFs are a biological feature of that UTR, not a frameshift into the reporter. Still note them.

## Construct naming / 构建命名

Constructs are design files. A workable template:

```
{class}{seq:03d}_{backbone}_{insert}.dna
```

The example project uses `pOE001_EV-FLAG_EXA1` (CDS fusion, `-TAA`) and `pPRO001_pro-EXA1-GUS` (transcriptional). Labs may use other class prefixes. Assign the next number in that class and never reuse it. Document Name should match the file stem.

When reporting status, separate: in-silico design (this file) vs cloned isolate vs glycerol stock. If the wet-lab table is missing, write **unknown / unrecorded**.

## Vendor orders / 订购（安全规则）

1. Only primers whose status is “to order” (or the lab's equivalent) go onto a draft sheet.
2. Copy `orders/templates/` to a new file in `orders/drafts/` (timestamped). Never overwrite the template.
3. Default fill (unless the user says otherwise): **purity OPC**, **2 OD / tube**, **1 tube**. Leave **delivery-form** and **notes** columns **empty**.
4. Sequences: uppercase A/C/G/T, no spaces. Primer names match the primer table.
5. Draft the vendor email (subject `{date}-{group}-{purpose}`) and show the user the **recipients, subject, body, and attachment path**.
6. **Stop.** Wait for an explicit confirmation. If there is no confirmation, do not send.
7. Only then may the same attachment and text go to the vendor address the **user** has configured. After a real send, copy the sheet to `orders/placed/` and run `primer_naming.py mark-synth` for those ids.
8. Do not invent a vendor mailbox. The addresses in this repo (`orders@example-vendor.com`, `researcher@example.org`) are fictional examples.

## Sequencing write-back / 测序回写

1. Receive a zip of `.ab1` traces into `seq/inbox/`.
2. Unzip into a new `seq/runs/<run_id>/` folder (keep the zip or a checksum note).
3. Align traces to the **construct** (the in-silico design), not only to a gene fragment.
4. Write a run report, then write confirmed sequence / disagreements into the construct `.dna` (backup first; sgffp only). Prefer Notes + a report first if the alignment is messy.
5. This skill does not include real traces. Do not invent base calls.

## Enzyme / In-Fusion checklist before you call a design done

- Arms found uniquely on the backbone?
- Product length = backbone − dropped stuffer + insert?
- Translation from the intended ATG reaches the tag (or the intended stop) with no surprise stop?
- Reverse primer name carries `+TAA`, `-TAA`, or neither, and that choice matches the biology?
- Any internal sites that would mislead a diagnostic digest noted?
- `.dna` backup + `preserve=True` + primer display defaults applied?

If any check fails, report the failure and stop. Do not “fix” a frame by guessing a base.

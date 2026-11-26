# Example project walk-through

Timezone used in the example Notes: Asia/Shanghai (UTC+8).

## 1. Two gene workfiles

`DNA001_EXA1_EXApl.01G000100` and `DNA002_EXA1L_EXApl.03G000200` are short fictional genomic molecules (333 bp) with the same structure: 60 bp upstream, 30 bp 5′UTR, CDS1, a tiny intron, CDS2 (in-frame XhoI `CTCGAG`), 3′UTR, downstream.

Colors follow the skill: 5′UTR green `#008000`, CDS orange `#ff9900`, intron and flanks gray `#808080`, 3′UTR blue `#3366ff`.

The two CDS sequences are about 93% identical; the UTRs are about 13% identical. A primer pair that sits only in the conserved CDS would likely co-amplify both copies from mixed cDNA. The skill therefore requires a **UTR-specific** amplicon, clone, and sequence **before** a CDS In-Fusion on this kind of gene.

`UTR-EXA1-F/R` (unsynthesized functional names) and `Pr0001-UTR-EXA1L-F/R` (already given a vendor serial, as if that order had been sent) are the illustration.

Logic cDNA = 5′UTR + CDS + 3′UTR (213 bp including stop). The wet-lab template would be cDNA, not gDNA.

## 2. C-terminal FLAG fusion (`-TAA`)

Empty backbone `EV-FLAG` (401 bp, circular, fictional) carries a synthetic promoter, an XhoI–SpeI window, and a public 3xFLAG-like tag (amino-acid motif DYKD…).

The insert is EXA1 CDS **without** the native TAA (150 bp). In-Fusion arms are taken from the backbone itself (F ends with `CTCGAG`, R is the reverse complement of the SpeI + linker window). Compensation is 0 / 0. The reverse primer is named `IF-EVFLAG-EXA1-R-TAA`.

Rebuilt translation (insert ATG through the tag):

```
MKLAIGAPFAKIEARLNGAIPAAYKEIANLLEDFIEKAPGNLAIDGKIEL-TSGG-DYKDHDGDYKDHDIDYKDDDDK*
```

The in-frame SpeI + linker contributes `TSGG`. A sticky-end “+1” model must not be used to invent an extra base here.

The construct file `pOE001_EV-FLAG_EXA1` is an **in-silico design**. Its Notes record the wet-lab clone as **unknown / unrecorded**.

A filled vendor draft for this primer pair is in `orders/drafts/` (purity OPC, 2 OD, 1 tube, delivery form and notes left empty, customer Jane Researcher). It is a draft only.

## 3. Promoter–reporter transcriptional fusion (no ±TAA)

Empty backbone `EV-GUS` (406 bp, circular, fictional) has an EcoRI–NcoI window. The reporter CDS is a short synthetic ORF, **not** uidA. The reporter ATG is the `ATG` inside NcoI `CCATGG` (the A is at offset +2 of the site, not +1).

The insert is EXA1 `upstream[10:60] + 5′UTR` (80 bp), ending **before** the gene ATG. The F core is allowed to start inward of the 60 bp flank. Reverse primer: `IF-EVGUS-proEXA1-R` — **no** `+TAA` or `-TAA`. The reporter uses its own start.

Translation from the declared reporter ATG:

```
MAGSSYRAPFAKIEARLNGAIPAAYKE*
```

## 4. Commands

From the repository root:

```bash
python3 skills/snapgene-lab-workflow/scripts/infusion_frame_check.py \
  --scheme examples/annotations/if_schemes.json --id EVFLAG-EXA1

python3 skills/snapgene-lab-workflow/scripts/infusion_frame_check.py \
  --scheme examples/annotations/if_schemes.json --id EVGUS-proEXA1
```

Both should print `OK`.

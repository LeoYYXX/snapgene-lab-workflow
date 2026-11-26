# Sanitization check

This repository was assembled as a generic skill with **fictional** example data only. It was not copied from a live plasmid library, primer table, or order-mail config.

## What was searched

The working tree, `.dna` and `.xlsx` bytes, and `git log -p` were searched for:

- real personal names, personal or institutional email domains, and personal phone numbers from the source lab
- real gene symbols, accession prefixes, empty-backbone codes, construct accession codes, and primer-table id prefixes from that lab
- real vendor mailbox addresses and customer blocks
- cloud folder or spreadsheet identifiers, and machine or hosting labels that would identify the source lab

The exact strings are **not** repeated in this file, so that a history-wide search of the public repository stays clean.

**Result: CLEAN — no matches in the working tree, in packaged binaries, or in commit history** (other than this generic description).

## What *is* in the tree (all fictional)

- Emails: `researcher@example.org`, `orders@example-vendor.com`, `modified@example-vendor.com`, `sequencing@example-vendor.com`
- Phone: `00000000000`
- Person / place: Jane Researcher, Example University, Example Lab, Example Vendor
- Genes: `EXA1` / `EXA1L` with accessions `EXApl.01G000100` / `EXApl.03G000200`
- Primer ids: prefix `PID-` (not a real inventory prefix)
- Constructs: `pOE001_EV-FLAG_EXA1`, `pPRO001_pro-EXA1-GUS`
- Backbones: `EV-FLAG`, `EV-GUS`

Example `.dna` files are 1–3 KB packages created with `sgffp.SgffObject.new` from the short FASTA in `examples/sequences/`. They were not exported from a live plasmids tree.

## Intentionally not included

- Any real primer sequence table or vendor serial list
- Any real gene workfile, empty-backbone inventory, or construct
- Real vendor mailboxes and customer blocks
- Cloud folder / spreadsheet identifiers
- Species-specific genome mapping tables and historical rename maps
- Lab progress / website / kanban plans

If you fork this skill for a real lab, keep private contacts in a gitignored `*.local.yaml` and never commit the live inventory to a public repository.

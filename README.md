# snapgene-lab-workflow

A reusable **SnapGene / molecular cloning skill** for coding agents (Cursor, Grok Bot, and similar). It condenses how to keep gene workfiles, name primers, rebuild In-Fusion junctions and reading frames, draft vendor primer orders, and write Sanger results back into `.dna` files.

This repository is **generic**. It ships **fictional** example data only (Example plant, genes EXA1 / EXA1L, backbones EV-FLAG and EV-GUS). It does **not** contain a real lab inventory, personal contacts, or real plasmid sequences.

一个可供 Cursor、Grok Bot 等编程代理复用的 **SnapGene / 分子克隆技能**：如何维护基因工作文件、命名引物、按真实骨架序列重建 In-Fusion 接合与读框、起草引物订购表，以及把 Sanger 结果写回 `.dna`。

本仓库是**通用模板**，只含**虚构**示例数据（Example plant，基因 EXA1 / EXA1L，空载 EV-FLAG、EV-GUS），**不含**真实实验室库存、个人联系方式或真实质粒序列。

## Install / 安装

Clone this repository, then copy the skill folder into the agent’s skills directory (name may vary by product):

```bash
git clone <this-repository-url>
cp -r snapgene-lab-workflow/skills/snapgene-lab-workflow /path/to/your/skills/
```

For Cursor, user skills typically live under a `~/.cursor/skills/` (or project `.cursor/skills/`) folder; after copying, the skill should appear when you type `/` or `@`. For Grok Bot, place `SKILL.md` (and the companion `scripts/` and `config/` directories) in the bot’s skills store, or point the agent at this clone and ask it to open `skills/snapgene-lab-workflow/SKILL.md`.

Python helpers (optional but recommended):

```bash
pip install -r requirements.txt
```

`sgffp` is required to read or write SnapGene packages. Biopython is optional (the frame-check script has a built-in translator). `openpyxl` is optional and only used to emit example Excel order sheets.

把本仓库中的 `skills/snapgene-lab-workflow/` 复制到代理的 skills 目录即可。读写 `.dna` 请安装 `sgffp`。

## What’s in the skill / 技能里有什么

The recipe is `skills/snapgene-lab-workflow/SKILL.md`. Sections:

1. Hard rules (sgffp only, backup before write, no invented biology, never auto-send vendor mail)
2. Directory convention
3. SnapGene I/O and primer display defaults
4. Gene workfile naming and feature colors
5. Primer naming (stable internal id; vendor serial only after a real send)
6. UTR–cDNA cloning
7. Paralog / copy-specific design
8. In-Fusion arms, reading frame, and `+TAA` / `-TAA`
9. Promoter–reporter transcriptional fusions (no ±TAA)
10. Vendor order drafts
11. Sequencing write-back

Helpers:

| file | role |
|------|------|
| `skills/snapgene-lab-workflow/scripts/sgffp_primer_write.py` | HybridizationParams + binding-site component when writing primers |
| `skills/snapgene-lab-workflow/scripts/primer_naming.py` | apply / mark-synth naming |
| `skills/snapgene-lab-workflow/scripts/infusion_frame_check.py` | rebuild a product from backbone arms and translate the junction |
| `skills/snapgene-lab-workflow/scripts/build_examples.py` | regenerate the fictional example project |
| `skills/snapgene-lab-workflow/config/` | templates (copy; do not commit real contacts) |

## Example project / 虚构示例

See `examples/README.md` and `examples/project.md`.

- Organism: Example plant
- Genes: **EXA1** and its paralog **EXA1L** (high CDS identity, divergent UTRs)
- Empty backbones: **EV-FLAG** (C-terminal FLAG-like tag) and **EV-GUS** (transcriptional reporter; the reporter CDS is synthetic, not uidA)
- Constructs: `pOE001_EV-FLAG_EXA1` (CDS fusion, reverse primer `-TAA`) and `pPRO001_pro-EXA1-GUS` (promoter-reporter, no ±TAA)
- Order draft customer: Jane Researcher, `researcher@example.org`, phone `00000000000`, Example University

Verify the two In-Fusion frames:

```bash
python3 skills/snapgene-lab-workflow/scripts/infusion_frame_check.py \
  --scheme examples/annotations/if_schemes.json --id EVFLAG-EXA1
python3 skills/snapgene-lab-workflow/scripts/infusion_frame_check.py \
  --scheme examples/annotations/if_schemes.json --id EVGUS-proEXA1
```

## Safety / 安全

Agents following this skill **must not** send a vendor email until the user explicitly confirms the draft. That rule is independent of which vendor or mailbox the lab uses.

Do not add real primer tables, real `.dna` libraries, phone numbers, or personal emails to a public fork. Put private contacts in a gitignored `order_email.local.yaml`.

使用本技能的代理：**在用户明确确认之前，不得向供应商发送邮件**。公开 fork 中不要加入真实引物表、真实质粒库或个人联系方式。

## License

MIT. See [LICENSE](LICENSE).

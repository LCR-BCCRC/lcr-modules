# Container Support Status

Tracks progress on adding Apptainer/Singularity `container:` directives to all modules so pipelines can run with `--use-singularity`.

**Status definitions:**
- ✅ **Complete** — every rule has a non-`None` container directive
- ⚠️ **Partial** — some rules have `container: None` (typically R scripts or custom tools without a suitable public image)
- ❌ **None** — all rules have `container: None`; module is conda-only

Last updated: 2026-05-19 (hmftools 1.0/1.1: all compute containers added; snpeff via custom lcr-scripts/snpeff:1.0; wget still null)

---

## Future Work: Pre-flight container check

**Idea:** Add an `onstart:` block to the demo (and eventually production) Snakefiles that aborts early when `--use-singularity` is active but one or more rules in the workflow have `container: None`.

Currently, rules with `container: None` silently fall back to the host environment under `--use-singularity`, so failures only surface when the tool is missing — not at startup.

```python
onstart:
    if workflow.use_singularity:
        uncontainerized = [
            r.name for r in workflow.rules
            if r.name not in workflow.localrules
            and r.container_img is None
        ]
        if uncontainerized:
            sys.exit(
                "ERROR: --use-singularity is active but the following rules have no "
                "container directive:\n  " + "\n  ".join(sorted(uncontainerized))
            )
```

**Open questions before implementing:**
- Hard error (`sys.exit`) vs. warning (`logger.warning`) — or make it configurable via a top-level config flag.
- Where to place it: demo Snakefile only, or a shared include (e.g., `lcr-snakefiles`) so production pipelines also benefit.
- Partial modules (e.g., `hmftools`, `sequenza`) intentionally leave some R-based rules as `container: None`. A known-exceptions list or downgrade to warning may be needed for those.

---

## ✅ Complete (57 module versions)

| Module | Version |
|--------|---------|
| bam2fastq | 1.0, 1.1, 1.2 |
| bwa_mem | 1.0, 1.1 |
| cnvkit | 1.0 |
| controlfreec | 1.0, 1.1, 1.2 |
| cutadapt | 1.0 |
| cutesv | 1.0 |
| freebayes | 1.0 |
| gridss | 1.0, 1.1, 2.0 |
| ichorcna | 1.0, 1.1 |
| liftover | 1.0, 1.1, 1.2, 2.0 |
| lofreq | 1.0, 1.1 |
| manta | 2.0, 2.1, 2.2, 2.3 |
| mixcr | 1.1, 1.2 |
| modkit | 1.0 |
| mutect2 | 1.0, 2.0 |
| picard_qc | 1.0 |
| qc | 1.0 |
| sage | 1.0 |
| sequenza | 1.4 |
| slms_3 | 1.0 |
| sniffles | 1.0 |
| star | 1.0, 1.1, 1.2, 1.3, 1.4 |
| starfish | 1.0, 2.0 |
| strelka | 1.0, 1.1 |
| stringtie | 1.0 |
| utils | 2.0, 2.1 |
| varscan | 1.0, 1.1 |
| vcf2maf | 1.0, 1.1, 1.2 |
| whatshap | 1.0, 2.0 |

---

## ⚠️ Partial (21 module versions)

Rules with `container: None` are typically R-based rules or scripts without a suitable public image.

| Module | Version | None rules / Total | Notes |
|--------|---------|-------------------|-------|
| battenberg | 1.1 | 2/6 | |
| battenberg | 1.2 | 1/6 | wget (reference download only); battenberg compute container added via lcr-scripts/battenberg:1.0 |
| clairs | 1.0 | 2/5 | |
| clairs_to | 1.0 | 2/5 |
| controlfreec | 1.3 | 4/15 |
| hmftools | 1.0 | 1/16 | wget (reference downloads only) |
| hmftools | 1.1 | 1/14 | wget (reference downloads only) |
| hotmaps | 1.0 | 13/16 |
| hotmaps | 1.1 | 13/16 |
| manta | 1.0 | 1/4 |
| pathseq | 1.0 | 1/7 |
| sage | 1.1 | 1/4 |
| salmon | 1.0 | 1/2 |
| salmon | 1.1 | 1/4 |
| salmon | 1.2 | 1/2 |
| sequenza | 1.0, 1.1, 1.2, 1.3 | 2/3 each |
| svar_master | 1.0 | 2/3 |
| vcf2maf | 1.3 | 2/6 |

---

## ❌ None (34 module versions)

These modules have no container support and will not run under `--use-singularity`.

| Module | Versions |
|--------|---------|
| battenberg | 1.0 |
| clair3 | 1.0 |
| cnv_master | 1.0 |
| dlbclass | 1.0 |
| dnds | 1.0, 1.1, 1.2 |
| ecotyper | 1.0 |
| ega_download | 1.0 |
| fishhook | 1.0, 1.1, 1.2 |
| gistic2 | 1.0, 1.1, 1.2 |
| igcaller | 1.0 |
| lymphgen | 1.0, 1.1, 2.0, 2.1 |
| mutsig | 1.0, 1.1, 1.2 |
| nanomethphase | 1.0, 1.1 |
| oncodriveclustl | 1.0, 1.1 |
| oncodrivefml | 1.0, 1.1 |
| rainstorm | 1.0, 1.1, 1.2 |
| sigprofiler | 1.0 |
| spechla | 1.0 |

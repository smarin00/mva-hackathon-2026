# LLM usage declaration

Required by the organisers' 28 Aug 2026 update ("record the provider, the plan
or tier, and the relevant data-handling setting").

- **Provider and plan:** Anthropic, Claude via claude.ai, consumer Pro plan
  (Claude Opus 5.5).
- **Data-handling setting:** model training disabled ("Help improve Claude"
  off); no output ratings or feedback given, since feedback falls outside the
  training opt-out.
- **Scope:** pipeline design, code and troubleshooting (Track 1 and Track 2),
  HPO term mapping, interpretation of named candidate variants, literature
  search for mechanism and drug candidates; drafting of the gene knowledge
  file `config/genes/bub1b.yaml` (domain boundaries, calibration alleles,
  pathway panel, cancer regimens), checked by the author against the cited
  papers; drafting of the two reports, the methods forms and the pitch video
  script and slides. All analysis code was run by the author on the analysis
  server; results and literature claims were reviewed by the author.
- **What was shared:** code, tool logs, named candidate variants (two BUB1B
  alleles plus Track 1 candidate rows), per-chromosome aneuploidy summaries,
  and drug/protein tables derived from public data. No genome-scale variant
  tables, reads or alignments were shared; the assistant never had access to
  the server.
- **Basis:** Sage Privacy & Compliance Officer guidance (Space Discussion #2):
  the provider acts as a processor.
- **Deletion:** chats containing variant rows are deleted with the data by
  23 Nov 2026 (`docs/DELETION.md`).
- **Pitch video narration:** synthetic voice from Kokoro v1.0 (open-source
  text-to-speech model, voice `af_heart`), run locally on the narration text
  only; no third-party service and no patient data involved.

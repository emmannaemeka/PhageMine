# Start refinement and reviewer exports

Alternative short SD motifs are attenuated to 25% above 58% genomic GC. Caller-supported starts retain their existing motif scores, while canonical 5/6-base motifs remain unchanged. Sequence-confirmed 1/4-base translational coupling supplies independent initiation evidence. No reference coordinates are used.

`--phrogs-start-db` accepts a DIAMOND database containing PHROGs protein members. Candidate and subject N-termini must both align within the first two residues, with at least 40% identity, 80% query and subject coverage, and E <=1e-5. General PHROGs family/domain matches alone do not establish a start boundary. Swiss-Prot anchoring remains available.

`--strict-precision` excludes unsupported single-caller CDSs shorter than 150 nt. Default mode retains them. Structural confidence is recorded separately from product confidence in GFF3, QC and `cds_confidence.json`; excluded calls remain in the review queue.

Metadata can be JSON or YAML via `--metadata`. CLI overrides include `--organism`, `--host`, `--isolation-source`, `--country`, `--topology` and `--submitter-template`. Topology defaults to linear. Source qualifiers are included in the feature table and FASTA, and the GenBank source feature is synchronized. Missing source/submitter fields are reported; none are invented. Split-homology alerts propose `ribosomal_slippage` for review without asserting recoding or constructing a join from insufficient evidence.

Overlap screening includes same-strand +1/+2 frames, excludes RNA intervals, and retains RBS or spanin-context candidates. For SD-free candidates, a first-pass search of the longest peptide at each shared stop preserves a protein-family evidence exception. Weak preliminary family hits only qualify candidates for confirmation; promotion still requires the original stringent E-value and coverage thresholds. The earlier +1 policy is evaluated before newly eligible +2 candidates. Searches run sequentially in batches of at most 128 queries with at most two threads per overlap adapter; this bounds concurrent work, not total database memory.

Every completed run writes `final_annotation_package/` with final sequence/annotation files, QC reports, and four-sheet `<prefix>_manual_verification.xlsx`. Annotation strings are written as Excel text, not formulas. The benchmark uses the default sensitivity mode; precision-mode exclusions are tested independently.

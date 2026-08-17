# External databases and executables

PhageMine Core does not bundle Pfam, VOGDB, PHROGs, Swiss-Prot, or PMFDB.
These independently installed, versioned resources are supplied only when
evidence-enhanced annotation or external PMF validation is requested.

Optional components include HMMER (`hmmscan`) with Pfam/VOGDB HMM files,
PHROGs profiles and annotations with MMseqs2, DIAMOND with Swiss-Prot and
matching metadata, and MMseqs2 with a PMFDB target database/index.

Record paths, versions, and checksums in run provenance. Missing optional
resources are recorded as unavailable, not as biological absence. PMFDB is an
independently versioned external knowledgebase.

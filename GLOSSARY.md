# Glossary

Shared terms across hgvs-weaver-data docs and code.

## Reference data

- **Assembly** — one version of the reference genome, GRCh38 or GRCh37. Coordinates mean nothing without the assembly
  they are on, so every alignment and every store names its assembly.
- **Annotation release** — one publication of a publisher's gene and transcript models against an assembly, named by the
  publisher (`RS_2024_08` for RefSeq).
- **Transcript record** — the publisher's own sequence of a transcript (an `NM_` or `NR_` accession with its version).
  It usually matches the genome spliced at the transcript's exons, but not always.
- **Alignment** (or **placement**) — where a transcript sits on one sequence of an assembly: its exons, strand, and the
  per-exon CIGAR that says how the record and the genome differ. Published by NCBI; never derived here. A transcript can
  have several — on X and on Y, on a chromosome and on a patch, or on each alternate locus that carries its gene.
- **Alternate locus** — a scaffold holding another haplotype of a region too variable for one sequence (the MHC, the KIR
  cluster). A gene the primary assembly's haplotype lacks is placed only on alternate loci. Alternate loci are `NT_`
  accessions, as are the unlocalized and unplaced scaffolds of the primary assembly.
- **Patch** — an `NW_` scaffold released between assembly versions: a fix patch corrects the primary assembly (and is
  folded into the next major release), a novel patch adds alternate sequence.
- **MANE** — the NCBI/EMBL-EBI project that picks one representative transcript per protein-coding gene (MANE Select)
  and matches it across RefSeq and Ensembl. It ranks transcripts; it does not place them.
- **HGNC** — the committee that assigns human gene symbols; the source of each gene's current, previous and alias
  symbols.
- **refget digest** — a sequence's identity computed from its residues (sha512t24u, written `SQ.…`), so two copies of
  the same sequence share an identifier whatever they are called.

## Project

- **Gene bundle** (or **bundle**) — the unit of storage: one gene with every transcript, protein, alignment and sequence
  that belongs to it (`GeneBundle` in `proto/`).
- **Shard** — one annotation release's bundles in one bagz file, immutable once written and named for its contents.
- **Store** — an index over an ordered list of shards: the key tables (symbol, accession, protein, digest), the interval
  tables (bundles by genome position), and the **manifest** that names the shards, written last.
- **Genome** — a derived copy of an assembly: the **blocks file**, every sequence cut into fixed-size compressed blocks,
  and the **catalogue** of names, lengths and refget digests, written last.
- **Format version** — the number every manifest and catalogue carries; a reader refuses one it does not know.
- **Provider** — the object weaver reads through: weaver's `DataProvider` over one assembly's store and genome.

## Testing

- **Test double** — anything standing in for a production dependency in a test. The kinds below differ in fidelity and
  in what a test can do with them
  ([`docs/style/writing-tests.md`](docs/style/writing-tests.md#test-doubles-in-order-of-fidelity)).
- **Fake** — a working, lightweight implementation of an interface, unsuitable for production and held to the real
  implementation's contract.
- **Stub** — returns canned values to put the code under test in a state; makes no claim about being called.
- **Mock** — carries expectations about how it is called and fails the test when they are not met; the only double a
  test verifies interactions on.
- **Change detector** — a test that pins the code's current shape rather than an invariant, so it fails on every
  legitimate change and catches no defect.
- **Flaky test** — one that passes and fails on the same code; a defect to root-cause, never a reason to rerun.

# GenBank and UniProt records

**Related:** [`../PRODUCT.md`](../PRODUCT.md) — the principles this design applies; [`store.md`](store.md) — bundles,
shards stacked behind one index, and the key tables these records are found by; [`placements.md`](placements.md) — what
a placement is and why it is always the publisher's; [`retired-versions.md`](retired-versions.md) — the other shard
stacked under a release, and the precedent for records left out and counted; [`../../GLOSSARY.md`](../../GLOSSARY.md) —
bundle, shard, store, refget digest.

## Overview

Some variant names cite a sequence that no RefSeq or Ensembl release holds: a GenBank mRNA such as BRCA1's `U14680.1`,
cited by thousands of ClinVar submissions, its protein `AAA73985.1`, or a UniProt protein such as `P38398`. weaver can
parse every one of them, but the store has nothing to answer with.

The decision is to **serve these records for what they state, and never to place them on the genome.**

- **What a record states is served.** A GenBank mRNA carries its sequence, its CDS and its protein. A UniProt entry
  carries its proteins. That is enough to check a name's stated bases or residues, to normalise it on the record, and to
  go from `c.` to `p.` within a GenBank record. It is not enough to project to the genome or to another transcript.
- **Projection is left to the caller.** No publisher places these records on a current assembly. A placement this
  package computed would be a choice of its own presented as data, which the product rules out. The same residues under
  another publisher's name are visible for free, because sequences are stored by refget digest; what a caller does with
  that is its policy.
- **Each publisher is a shard of its own, with no placements.** It stacks under any store's index, GRCh37 or GRCh38,
  like the historical set does, and the inputs given to the build decide which records it holds.
- **These records are found by accession, never by gene.** A gene's transcripts and proteins, as weaver enumerates them
  for a name that cites only a gene, stay the reference publishers' records.
- **Versions are stored as published, and a bare accession is not resolved by the provider.** A UniProt name almost
  never carries a version. The provider lists the versions it holds; weaver reads the name on each
  ([hgvs-weaver#145](https://github.com/populationgenomics/hgvs-weaver/issues/145)), and choosing among them is the
  caller's policy.

Both need hgvs-weaver 2.0, which reads a transcript with no placement on its own sequence and lists the versions of an
unversioned accession.

## Background

### Who cites these records

HGVS recommends RefSeq, Ensembl and LRG as reference sequences. GenBank accessions are not forbidden, only not
recommended, and UniProt is not mentioned. In practice both appear:

- **GenBank mRNAs** were the reference standards for a few genes before RefSeq existed. BRCA1's `U14680.1` and BRCA2's
  `U43746.1` carry the Breast Cancer Information Core's numbering, which ClinVar still holds about 3,700 submissions in.
  MECP2's `AJ132917.1` and TREX1's `AF319569.1` follow. NCBI itself lists nine such records as reference standards for
  their genes. Across ClinVar, names on GenBank accessions are about 0.07% of submitted expressions, concentrated on
  these few records ([Appendix](#who-cites-genbank-records)).
- **UniProt accessions** are how structural biology and protein-function papers name residues, and how UniProt's own
  variant catalogue writes them (`P38398:p.Cys61Gly`).

Other tools mostly refuse these names. VariantValidator and LOVD accept only RefSeq, Ensembl and LRG. Mutalyzer
describes a GenBank name on the record itself and does not project it. Ensembl VEP accepts a UniProt accession only as a
way to find a gene, then tries each of the gene's transcripts until one has the stated residue, which is the guess this
package exists not to make.

### What a GenBank record states, and what it does not

A GenBank mRNA record is a sequence, a CDS location in the record's own coordinates, and on the CDS a protein accession
with its translation. That determines `c.` to `p.` exactly, as a RefSeq record's CDS does
([`placements.md`](placements.md)). Two caveats, both from the INSDC feature table:

- A CDS marked with an exception has a translation that intentionally differs from the conceptual one.
- An mRNA record need not have a CDS feature at all. That means the CDS was not annotated, not that the RNA is
  non-coding.

What a record never states is where it sits on a genome. NCBI aligned GenBank mRNAs to GRCh38 when building its
annotation, and published those alignments for two releases up to 2022, then stopped. It published none for GRCh37.
UCSC's alignments of GenBank mRNAs are its own BLAT runs, last updated in 2020, on unversioned accessions. RefSeq's
records state which spans of which GenBank records they were built from, but that covers only the records RefSeq cites:
BRCA2's `U43746.1` is not among the sources of `NM_000059.4`.

### What UniProt identifies, and how stably

A UniProt entry has an accession (`P38398`) and one sequence, its canonical isoform, with a sequence version that
increments when the sequence changes. Accession and sequence version together name one sequence, retrievable from
UniProt's archive. More than half of human Swiss-Prot entries have changed sequence at least once, and a change can
shift the numbering: haemoglobin beta's first sequence version lacked the initiator methionine, so its sickle-cell
residue was Glu6 there and is Glu7 now.

Each isoform, the canonical included, has an ID of its own (`P38398-2`), but no version, and an isoform ID has been
reassigned to a different sequence: `P38398-2` named an 1,863-residue sequence until 2013 and a 63-residue one since.
Citations of either almost never carry a version.

UniProt's canonical isoform is the MANE Select protein for about 91% of genes, and UniProt's entries cite their MANE
Select transcript. Where a UniProt cross-reference names the isoform, the sequences were identical in every case
measured; where it names only the entry, it is a correspondence, not an identity.

## Non-goals

- **No placement of a GenBank record on the genome**, by any route: an alignment this package runs, NCBI's frozen 2022
  alignments, or one inherited through a RefSeq transcript. Each is weighed under
  [Alternatives considered](#alternatives-considered).
- **No choice of version for a bare accession.** The provider lists versions; weaver reads the name on each, and the
  caller chooses.
- **No mapping between publishers beyond sequence identity.** That `AAA73985.1` has the same residues as `NP_009225.1`
  is a fact the store holds; that a `U14680.1` position "is" an `NM_007294.4` position is a mapping, and mapping is
  projection.
- **No downloading.** A build reads local files; which GenBank records, which UniProt release, and whether TrEMBL is
  included are decided by the files the caller provides.
- **GenBank chromosome accessions** (`CM000663.2` for chromosome 1) are not these records. They are other names for
  sequences the genome already holds, and belong to the genome's catalogue.

## Design

### A record is served for what it states

The provider answers three things about these records: their residues, the CDS of a GenBank transcript, and the protein
a GenBank transcript encodes. weaver then validates a name against the record's own sequence, normalises it there, and
predicts a GenBank transcript's protein consequence from its own CDS. Everything that needs a genome position is refused
with an error that names the missing placement: `c.` to `g.`, the transcripts over a region, an intronic offset such as
`c.4096+1`. The provider gives the record without a placement whatever genome the caller asks for, so the error weaver
raises says the placement is missing rather than that the lookup failed.

A shift in normalisation runs the length of the record. A placed transcript's shift stops at an exon junction, so that
the change stays one change on the genome; a GenBank mRNA has no genome to disagree with, and states no exons.

Refusing is the honest answer because the alternative is a placement nobody published. Where a GenBank mRNA differs from
the genome by an insertion or deletion inside a repeat, where the gap goes is a choice, and the choice moves every
position downstream of it. [`placements.md`](placements.md) explains why the store takes that choice from the publisher,
and GenBank's publisher made none it still publishes.

### Each publisher is a shard, and the build's inputs choose its records

GenBank and UniProt each get a builder that cuts one shard, as RefSeq and Ensembl do. The shard stacks behind a store's
index with the releases already there, so one provider serves all of them.

A shard of these records has no placements, so it is on no assembly. The index already checks assembly per placement,
not per shard, so the same GenBank or UniProt shard can be listed in a GRCh37 store and a GRCh38 store alike.

Which records a shard holds is decided by the files the build is given, and the shard records which files those were:

- **GenBank:** any GenBank flat file, plus NCBI's table of accessions by gene. A caller wanting the nine reference
  standards fetches nine records; one wanting every human mRNA with a protein fetches those. The builder treats them
  alike.
- **UniProt:** one release's reviewed human entries with their isoforms, and, if the caller also gives it, the release's
  unreviewed entries, which add about 127 thousand entries to the 20 thousand reviewed.

### What each publisher's record becomes

| Record                        | Stored as                                    | Notes                                               |
| ----------------------------- | -------------------------------------------- | --------------------------------------------------- |
| GenBank mRNA with a CDS       | a transcript: sequence, CDS, no placement    | the CDS as the record states it, open ends included |
| its CDS's protein             | a protein, linked to that transcript         | the record's translation                            |
| GenBank mRNA without a CDS    | left out, counted                            | no CDS stated is not "non-coding"                   |
| GenBank CDS with an exception | left out, counted                            | its translation is not the CDS's                    |
| UniProt canonical isoform     | a protein, versioned by its sequence version | also keyed under its isoform ID; no transcript      |
| UniProt other isoform         | a protein, under its bare ID                 | one sequence per shard, as the release states it    |

Records left out are counted in the build's report with the reason, as the historical set's undecipherable records are
([`retired-versions.md`](retired-versions.md)). The provider then answers "not in the reference data" for them: a
missing answer, not a wrong one.

A UniProt protein has no transcript. Its cross-references to RefSeq and Ensembl transcripts are not stored as links:
where they mean identity, the shared digest already says so; where they do not, a stored link would claim more than
UniProt does.

### A record belongs to the genes its publisher names

A bundle is one gene ([`store.md`](store.md)), so every record needs a gene. It takes the one its publisher states:
NCBI's gene for a GenBank accession, UniProt's gene cross-references for an entry. The gene's identity then comes from
HGNC, as for every other publisher.

Two cases depart from one record per gene:

- **A record naming several genes goes into each of their bundles.** UniProt has one entry for a protein several genes
  encode identically. The histone H4 entry names more than a dozen. The key tables already map one key to every bundle
  holding it.
- **A record naming no gene is left out, counted.** Guessing a gene from a symbol in free text would be the publisher's
  statement made up.

### Found by accession, never by gene

A name that cites only a gene, `BRCA1:c.68_69del`, has weaver try each of the gene's transcripts. Those are the gene's
reference transcripts: RefSeq's and Ensembl's. A GenBank mRNA joining them would add a candidate weaver cannot place,
and a historical sequence nobody meant. A UniProt protein joining the gene's proteins would add a second name for a
sequence the list already holds.

So these records are not found by gene. The index does not key their bundles under gene symbols. They are found by
transcript accession, protein accession and digest only. The gene each belongs to is still recorded, so a caller holding
`U14680.1` can ask whose it is.

### Identity across publishers is the digest

Every sequence is stored once per bundle and keyed by its refget digest. When `AAA73985.1` has exactly the residues of
`NP_009225.1`, the two share a digest, and a lookup by that digest finds both. That is the whole crosswalk, and it
claims nothing beyond identical residues.

Where a digest-to-accession answer has to name one accession, the provider names a reference publisher's first: RefSeq
or Ensembl, then GenBank, then UniProt. Any of them is correct, but the order makes the answer stable as shards are
added.

### Versions are stored as published; the provider does not choose one

A GenBank record is versioned like a RefSeq one, and is keyed the same way. A UniProt canonical is keyed by its
accession and sequence version. No publisher writes that pair as `P38398.3`, but it names one sequence, as UniProt's
archive does. The provider's answers therefore carry the version, which UniProt's own spelling does not.

Every UniProt isoform, the canonical included, has an isoform ID and no version, so it is keyed by its bare ID,
`P38398-2`, with the sequence the shard's release gives it. The canonical is keyed twice: under its sequence version,
and under the isoform ID the release marks as the one it displays. That ID is not always `-1`. In release 2026_03 it is
another number in 481 of the 10,684 human entries with isoforms, and 116 of them assign no `-1` at all. An ID the
release does not assign is not keyed: `Q96AP0-1` is not in the reference data, although UniProt's own web service
answers it with the canonical. A shard cut from one release holds one sequence per isoform ID; where shards from two
releases are stacked, the later one answers, as for every record. The store does not number the sequences an isoform ID
has named over time. A `P38398-2.2` would be a name no publisher issued, which a caller could then write and nobody else
resolve.

A name with no version is not resolved by the provider. It answers a versioned name and refuses a bare one, unless the
record has no version of its own. It also lists what it holds under a bare accession: the versions, newest first, or for
an isoform the bare ID alone. weaver reads a name on each
([hgvs-weaver#145](https://github.com/populationgenomics/hgvs-weaver/issues/145)), and returns every reading consistent
with its sequence, with the store's status and tags for each and whether the name stated anything that was checked. The
caller chooses among them: MANE Select first, the newest version, or only an answer one sequence supports. The reason
the choice lives with the caller and not here: a choice made inside the provider would be invisible to everything
downstream. And "the newest version whose residue matches" is sometimes right only by chance, after an insertion
upstream shifted the numbering.

These records carry no status or tags from their publisher. A GenBank or UniProt record is current in the release the
shard was cut from, and has no MANE or Select tag.

### The provider asks the store what kind an accession is

The provider decides whether an accession is a genomic sequence, a transcript or a protein from a fixed list of
prefixes. That list cannot name GenBank or UniProt accessions, which follow no prefix scheme a reader could rely on. The
provider instead asks its key tables: an accession held as a transcript is a transcript, one held as a protein is a
protein, a sequence of the genome is genomic.

### A store holding these shards is format version 2

A reader of the current format decides an accession's kind by its prefix. From a store holding `U14680.1` it would
answer that the record is not in the reference data, which is an answer, and false. So a store with GenBank or UniProt
shards is stamped format version 2, which the readers deployed so far refuse at open, and a store with only RefSeq and
Ensembl shards stays at version 1, which they still read. Each shard records when it is cut whether it needs version 2,
so the index builder takes the stamp from its shards rather than from a list of publishers. The range of versions a
reader opens, and how it narrows, is the store's policy ([`store.md`](store.md#a-store-is-files-and-it-opens-once)).

### Consequences

- **The provider needs hgvs-weaver 2.0.** A protein is served through a sequence lookup, which any weaver reads; a
  transcript with no placement and the list of an accession's versions are 2.0's.
- **The BIC names are answered on the record, as written.** ClinVar's BRCA1/2 names on `U14680.1` and `U43746.1` are
  written `n.`, which HGVS reserves for non-coding records, but BIC numbering counts from the record's first base, which
  is what `n.` on the record means: `U14680.1`'s CDS starts at base 120, and BIC's `185delAG` deletes the record's bases
  185 and 186, an AG. weaver reads an `n.` position on a coding record, so these names validate and normalise as
  written. Reading them as `c.`, or projecting them to `NM_007294.4`, is still the caller's step.
- **A canonical that moves leaves its old name behind.** In release 2026_01 UniProt made `P08195-2` the canonical of
  P08195 in place of `P08195-1`, and the sequence version went from 3 to 4. Each isoform ID kept its sequence; the
  accession's did not. A citation of `P08195:p.…` written in 2025 meant the sequence now keyed `P08195-1`, but `P08195`
  lists only `P08195.4`, and nothing in one release connects the two names. Connecting them needs UniProt's archive,
  which a build does not read.
- **Equivalence on these records is refused**, even between two spellings on one record, because weaver compares on the
  genome. Comparing on the record alone is weaver's to add
  ([hgvs-weaver#103](https://github.com/populationgenomics/hgvs-weaver/issues/103)).
- **A store's key and interval tables are unchanged in shape.** The new records add keys, and no intervals.

## Alternatives considered

### Place a GenBank record on the genome ourselves

Splice-align the record to the assembly, project its `c.` to `g.`, then to every RefSeq and Ensembl transcript through
their published placements.

The information is there, and most positions leave the aligner no choice. Rejected because the rest do. An mRNA matches
its processed pseudogenes without introns, and paralogs nearly as well. A 1990s cDNA carries polymorphisms and
sequencing errors: `U14680.1` has four mismatches and a gap against GRCh38. And where it differs from the genome by an
insertion or deletion in a repeat, the aligner chooses where the gap goes. Each of those is a placement this package
would choose and then present as data.

### Inherit a placement through a RefSeq transcript

Align the GenBank record to a RefSeq transcript of the same gene, end to end with no gap, and carry over that
transcript's published exons. With no gap there is no indel to place, and the gene, not a genome search, picks the
candidates.

This is the strongest of the placement routes, and the records it covers are the ones cited. It was set aside, not
refuted. It is still a mapping this package composes, and the product's line puts mapping with the caller. A caller
wanting it now has what it needs: both sequences, both CDSs, and the RefSeq transcript's placement. If it is wanted here
later, NCBI's 2022 alignments and the Locus Reference Genomic project's per-exon differences are independent oracles to
test it against, built by neither this package nor its inputs.

### Use NCBI's 2022 alignments, or RefSeq's source spans

NCBI's last published alignments of GenBank mRNAs are a publisher's placement, but frozen, on GRCh38 only, and not
maintained. A store built from them would answer GRCh38 and refuse GRCh37 for the same record, for no reason a caller
could see in the record. RefSeq's statements of which GenBank spans each record was built from are current, but cover
only the records RefSeq cites, and say only that spans correspond, not that their bases are identical.

### Choose a version inside the provider

Resolve a bare accession to the only version held, or to the newest one consistent with the name.

Rejected because the choice would be invisible. weaver would receive one sequence and report its answer as though the
name had cited that version. The newest-consistent rule is also wrong in exactly the case it cannot see, a numbering
shift between versions. weaver lists every consistent reading instead, and the caller's choice among them is visible in
its code.

### Number an isoform's sequences

Give each sequence an isoform ID has named a version of the store's own, from UniProt's sequence archive or from the
releases that held it, so that a reassigned ID like `P38398-2` keeps both sequences and weaver can report that they
disagree. Rejected because the version would be this package's invention, issued under UniProt's accession: a caller
would see `P38398-2.2`, write it, and find no other resolver that knows it. One sequence per isoform per shard is what
the release states.

### Key the canonical as `-1`

UniProt's web service answers `P38398-1` with BRCA1's canonical, and answers `-1` with the canonical even for an entry
that assigns no `-1`. Rejected because the release says otherwise for one entry in twenty: `P08195-1` is a non-canonical
isoform with a sequence of its own, so keying the canonical under it would serve the wrong residues for a name that
cites the isoform.

### Leave the format version alone, or raise it for every store

Leaving it alone lets a current reader answer "not in the reference data" for records a store holds. Raising it for
every store refuses RefSeq-only stores that current readers read correctly, and makes every deployed reader upgrade for
a change most stores do not carry.

### Find these records by gene too

List a gene's GenBank mRNAs and UniProt proteins with its RefSeq and Ensembl ones. Rejected for the reasons under
[Found by accession, never by gene](#found-by-accession-never-by-gene): a name citing only a gene means the gene's
reference sequences, and every extra candidate is one more way for weaver to answer with a sequence nobody meant.

### Separate stores

Build GenBank and UniProt records into stores of their own, beside the RefSeq store. Rejected because weaver reads
through one provider, and a caller would have to write the provider that merges them. Shards stacked behind one index
already merge them, by construction.

## Appendix

### Who cites GenBank records

ClinVar's submitted HGVS expressions, by reference, measured over its tab-delimited release in November 2026. Of 5.4
million expressions, 3,882 cite a GenBank accession:

| Reference    | Gene               | Expressions | Form                |
| ------------ | ------------------ | ----------- | ------------------- |
| `U43746.1`   | BRCA2              | 1,950       | `n.`, BIC numbering |
| `U14680.1`   | BRCA1              | 1,713       | `n.`, BIC numbering |
| `AJ132917.1` | MECP2              | 92          | `c.`                |
| `CM0006xx.1` | GRCh37 chromosomes | ~90         | `g.`                |
| `AF319569.1` | TREX1              | 14          |                     |
| `X95694.1`   | TFAP2B             | 8           | marked previous     |
| five others  |                    | 1 each      |                     |

NCBI's RefSeqGene project lists nine GenBank mRNAs as reference standards for their genes: those of BRCA1, BRCA2, GALT,
MECP2, PAH, RB1, SCN1A, SRD5A2 and TREX1.

# hgvs-weaver-data — product north star

## What it is

[hgvs-weaver](https://github.com/populationgenomics/hgvs-weaver) is an HGVS engine: it parses a variant name, checks it
against the reference it cites, and projects it between transcript, protein and genome. What it computes is only as
correct as what it is told about each transcript. hgvs-weaver-data is what tells it:

- **a builder** (`weaver-data-build`) that turns a publisher's annotation release and assembly into files on local disk,
  and
- **a reader** that serves those files, from local disk or `gs://`, as weaver's `DataProvider`.

```
publisher's release ──build──▶ store + genome ──read──▶ DataProvider ──▶ weaver
   (GFF3, FASTA,                 (bagz files)
    alignment BAMs, …)
```

## Why these principles

1. **The publisher's placement, not a derived one.** A RefSeq transcript's sequence is not always the genome's — about
   one record in forty differs — and where they differ by an insertion or deletion, where the difference sits is a
   choice. NCBI publishes its choice as an alignment of every transcript to the assembly. The builder takes that
   alignment exon by exon and never splices the genome to stand in for a record, so the positions weaver computes are
   the ones ClinVar's names and VariantValidator's projections follow. Every placement NCBI publishes is one, on a
   chromosome, an alternate locus or a patch ([`design/placements.md`](design/placements.md)).
1. **The files at rest are a contract.** A store is built once and read for years, from buckets this repo does not own.
   Its format is defined by the protos under `proto/`, evolves additively, and carries a format version the reader
   checks before it reads anything else. A shard of bundles is immutable at its path and named for its contents, so an
   index can never silently point at different bytes.
1. **An answer or an error, never a guess.** "Not in the reference data" is a valid answer. A substitute is not: the
   provider does not fall back to another version of a transcript, or to the genome in place of a record — choosing a
   stand-in is the caller's policy, stated where it is applied. A file that disagrees with its own manifest raises.

## Shape

- **Store**: one record per gene — every transcript, its protein, its alignments and the sequences they cite — in bagz
  shards, one per release, with key and interval tables read into memory at open. A lookup is an in-memory search and
  one ranged read.
- **Genome**: the assembly's sequences cut into compressed blocks, with a catalogue of names, lengths and refget
  digests. A slice, or a transcript's exons together, is one batched read.
- **Provider**: weaver's `DataProvider`, `Refget` and `TranscriptSearch` over one assembly's store and genome.
- **Language**: Python, with the at-rest format in protobuf.

## Relationship to weaver

At runtime the dependency runs one way: this package depends on hgvs-weaver, and weaver never imports it. Weaver's own
test suite may use this package at test time to check itself against real reference data; that cycle is confined to
tests and is acceptable. If changes to the protocol between the two keep needing to land in both repos at once, they
belong in one repo.

## Non-goals

- Not a downloader or mirror. Every builder input is a local file and every output a local directory; fetching the
  publisher's files and uploading the result are the caller's. The one fetch the builder makes itself is each RefSeq
  version's status from Entrez, which no published file states; it is a command of its own and writes a local table.
- Not an HGVS engine. Parsing, normalisation and projection are weaver's.
- Not a resolution policy. Which version to try when a cited one is absent, or how to read a legacy numbering, is
  decided by the application using weaver (Themis, for one).
- No liftover between assemblies; a provider serves one assembly, and a caller working across builds holds one per
  build.

## Scope

RefSeq and Ensembl annotation releases on GRCh38 and GRCh37, with MANE and HGNC as the crosswalks that name genes and
rank transcripts.

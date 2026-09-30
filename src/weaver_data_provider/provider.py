"""The weaver data provider over a store and its genome: what an HGVS engine asks for, answered locally.

weaver's `VariantMapper` takes a `DataProvider`: transcript models with per-exon alignments, sequence
slices of transcripts, proteins and the genome, symbol and accession crosswalks, and refget
identities. One provider serves one assembly; a resolver working across builds holds one per build.
Every answer comes from the store's in-memory tables and one bundle read, or one indexed read of the
genome.
"""

from __future__ import annotations

import weaver

from weaver_data_provider import genome as genome_mod
from weaver_data_provider import store as store_mod
from weaver_data_provider.v1 import bundle_pb2

_GENOMIC = ('NC_', 'NT_', 'NW_')
_TRANSCRIPT = ('NM_', 'NR_', 'XM_', 'XR_', 'ENST')
_PROTEIN = ('NP_', 'XP_', 'YP_', 'ENSP')


def _kind_name(kind: str | weaver.IdentifierType) -> str:
    return str(kind).lower()


def _sequence(bundle: bundle_pb2.GeneBundle, digest: str) -> bytes:
    for sequence in bundle.sequences:
        if sequence.digest == digest:
            return sequence.residues
    raise weaver.DataProviderError(f'{bundle.gene.symbol}: bundle cites sequence {digest} it does not carry')


class BundleProvider:
    """A `weaver.DataProvider` (and `Refget`, `TranscriptSearch`) over one assembly's store and genome."""

    def __init__(self, store: store_mod.BundleStore, genome: genome_mod.Genome) -> None:
        """One assembly's provider.

        Raises:
            ValueError: If the store and the genome are on different assemblies.
        """
        if store.assembly != genome.assembly:
            name = bundle_pb2.Assembly.Name
            raise ValueError(f'store is on {name(store.assembly)}, genome on {name(genome.assembly)}')
        self._store = store
        self._genome = genome

    # ---- what weaver reads ------------------------------------------------------------------------

    def _transcript(self, accession: str) -> tuple[bundle_pb2.GeneBundle, bundle_pb2.Transcript]:
        if '.' not in accession:
            raise weaver.DataProviderError(
                f'{accession}: a transcript model is per version; give the versioned accession'
            )
        base, _, version = accession.partition('.')
        # Records run across shards oldest release first, so the newest holding this version is last.
        for bundle in reversed(self._store.by_accession(accession)):
            for transcript in bundle.transcripts:
                if transcript.accession == base and str(transcript.version) == version:
                    return bundle, transcript
        raise weaver.DataProviderError(f'{accession}: not in the reference data')

    def get_transcript(self, transcript_ac: str, reference_ac: str | None) -> weaver.TranscriptData:
        bundle, transcript = self._transcript(transcript_ac)
        alignments = [a for a in transcript.alignments if reference_ac is None or a.chromosome == reference_ac]
        if not alignments:
            where = f' on {reference_ac}' if reference_ac else ''
            raise weaver.DataProviderError(f'{transcript_ac}: no alignment{where}')
        alignment = alignments[0]
        strand = 1 if alignment.strand == bundle_pb2.STRAND_PLUS else -1
        coding = transcript.HasField('cds')
        return {
            'ac': transcript_ac,
            'gene': bundle.gene.symbol,
            'cds_start_index': transcript.cds.start_index if coding else None,
            'cds_end_index': transcript.cds.end_index_inclusive if coding else None,
            'strand': strand,
            'reference_accession': alignment.chromosome,
            'exons': [
                {
                    'transcript_start': exon.transcript_start,
                    'transcript_end': exon.transcript_end,
                    'reference_start': exon.genome_start,
                    'reference_end': exon.genome_end_inclusive,
                    'alt_strand': strand,
                    'cigar': exon.cigar,
                }
                for exon in alignment.exons
            ],
        }

    def _residues(self, accession: str) -> bytes:
        if accession.startswith(_PROTEIN):
            for bundle in reversed(self._store.by_protein(accession)):
                for protein in bundle.proteins:
                    if f'{protein.accession}.{protein.version}' == accession:
                        return _sequence(bundle, protein.sequence_digest)
            raise weaver.DataProviderError(f'{accession}: not in the reference data')
        bundle, transcript = self._transcript(accession)
        return _sequence(bundle, transcript.sequence_digest)

    def get_seq(self, ac: str, start: int, end: int | None, kind: str | weaver.IdentifierType) -> str:
        """Residues over the 0-based half-open range, clipped to the sequence at both ends."""
        if ac in self._genome:
            return self._genome.fetch(ac, start, end)
        if 'genomic' in _kind_name(kind):
            raise weaver.DataProviderError(f'{ac}: not a sequence of this genome')
        residues = self._residues(ac)
        start = max(start, 0)
        return residues[start:].decode('ascii') if end is None or end < 0 else residues[start:end].decode('ascii')

    def get_symbol_accessions(
        self, symbol: str, source_kind: str, target_kind: str
    ) -> list[tuple[weaver.IdentifierType, str]]:
        """`gene -> c`: a gene's transcripts; `c -> p` or `gene -> p`: proteins; anything else is empty."""
        if source_kind == 'gene' and target_kind == 'c':
            transcripts = (f'{t.accession}.{t.version}' for b in self._store.by_symbol(symbol) for t in b.transcripts)
            return [(weaver.IdentifierType.TranscriptAccession, ac) for ac in dict.fromkeys(transcripts)]
        if source_kind == 'c' and target_kind == 'p':
            try:
                _, transcript = self._transcript(symbol)
            except weaver.DataProviderError:
                return []
            if not transcript.protein_accession:
                return []
            return [
                (weaver.IdentifierType.ProteinAccession, f'{transcript.protein_accession}.{transcript.protein_version}')
            ]
        if source_kind == 'gene' and target_kind == 'p':
            proteins = (f'{p.accession}.{p.version}' for b in self._store.by_symbol(symbol) for p in b.proteins)
            return [(weaver.IdentifierType.ProteinAccession, ac) for ac in dict.fromkeys(proteins)]
        return []

    def get_identifier_type(self, identifier: str) -> weaver.IdentifierType:
        if identifier.startswith(_GENOMIC) or identifier in self._genome:
            return weaver.IdentifierType.GenomicAccession
        if identifier.startswith(_TRANSCRIPT):
            return weaver.IdentifierType.TranscriptAccession
        if identifier.startswith(_PROTEIN):
            return weaver.IdentifierType.ProteinAccession
        if self._store.by_symbol(identifier):
            return weaver.IdentifierType.GeneSymbol
        return weaver.IdentifierType.Unknown

    def get_transcripts_for_region(self, chrom: str, start: int, end: int) -> list[str]:
        """Versioned accessions of transcripts placed over the 0-based closed interval, introns included."""
        found: dict[str, None] = {}
        for bundle in self._store.overlapping(chrom, start, end):
            for transcript in bundle.transcripts:
                if any(
                    a.chromosome == chrom
                    and min(e.genome_start for e in a.exons) <= end
                    and max(e.genome_end_inclusive for e in a.exons) >= start
                    for a in transcript.alignments
                ):
                    found[f'{transcript.accession}.{transcript.version}'] = None
        return list(found)

    # ---- refget ------------------------------------------------------------------------------------

    def get_refget_accession(self, ac: str) -> str | None:
        if ac in self._genome:
            return self._genome.digest(ac)
        try:
            if ac.startswith(_PROTEIN):
                for bundle in self._store.by_protein(ac):
                    for protein in bundle.proteins:
                        if f'{protein.accession}.{protein.version}' == ac:
                            return protein.sequence_digest
                return None
            _, transcript = self._transcript(ac)
        except weaver.DataProviderError:
            return None
        return transcript.sequence_digest

    def get_accession_for_refget(self, refget: str) -> str | None:
        """The accession carrying a digest: a chromosome first, then a transcript, then a protein.

        A `ga4gh:` prefix is accepted; a string that is not an `SQ.` digest names nothing, so answers None.
        """
        refget = refget.removeprefix('ga4gh:')
        if not refget.isascii() or not refget.startswith('SQ.'):
            return None
        chromosome = self._genome.name_for_digest(refget)
        if chromosome is not None:
            return chromosome
        try:
            holding = self._store.by_digest(refget)
        except ValueError:  # not a well-formed digest
            return None
        for bundle in holding:
            for transcript in bundle.transcripts:
                if transcript.sequence_digest == refget:
                    return f'{transcript.accession}.{transcript.version}'
            for protein in bundle.proteins:
                if protein.sequence_digest == refget:
                    return f'{protein.accession}.{protein.version}'
        return None

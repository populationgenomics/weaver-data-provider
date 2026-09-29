r"""`weaver-data-build`: cut a release into a shard, index shards into a store, cut a genome into blocks.

    weaver-data-build refseq --assembly GRCh38 --release RS_2024_08 \
        --annotation genomic.gff.gz --transcripts rna.fna.gz --proteins protein.faa.gz \
        --alignments knownrefseq_alns.bam --alignments modelrefseq_alns.bam \
        --hgnc hgnc_complete_set.txt --mane MANE.summary.txt.gz --records rna.gbff.gz --shards shards/
    weaver-data-build index --assembly GRCh38 --out store/ \
        shards/RS_2023_10-1a2b3c4d.bagz shards/RS_2024_08-5e6f7a8b.bagz
    weaver-data-build genome --assembly GRCh38 --fasta genomic.fna.gz --out genome/

Every input is a local file and every output a local directory; fetching and uploading are the caller's.
"""

from __future__ import annotations

import argparse
import pathlib
import sys
import warnings

from weaver_data_provider import build
from weaver_data_provider.build import genome as genome_build
from weaver_data_provider.build import refseq
from weaver_data_provider.build import store as store_build
from weaver_data_provider.v1 import bundle_pb2

_ASSEMBLIES = {'GRCh38': bundle_pb2.ASSEMBLY_GRCH38, 'GRCh37': bundle_pb2.ASSEMBLY_GRCH37}


def _refseq(args: argparse.Namespace) -> None:
    release = refseq.Release(
        assembly=_ASSEMBLIES[args.assembly],
        release=args.release,
        annotation=args.annotation,
        transcripts=args.transcripts,
        proteins=args.proteins,
        alignments=tuple(args.alignments),
        hgnc=args.hgnc,
        mane=args.mane,
        records=args.records,
    )
    path = store_build.write_shard(
        refseq.bundles(release), args.shards, release=args.release, inputs=release.inputs(), prefix=args.prefix
    )
    record = store_build.shard_record(path)
    print(f'{path}  {record.records} bundles', file=sys.stderr)
    print(path)


def _index(args: argparse.Namespace) -> None:
    manifest = store_build.write_index(args.out, args.shards, assembly=_ASSEMBLIES[args.assembly])
    total = sum(shard.records for shard in manifest.shards)
    print(f'{args.out}: {len(manifest.shards)} shards, {total} bundles', file=sys.stderr)


def _genome(args: argparse.Namespace) -> None:
    catalogue = genome_build.build_genome(args.fasta, args.out, assembly=_ASSEMBLIES[args.assembly])
    print(f'{args.out}: {len(catalogue.sequences)} sequences', file=sys.stderr)


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(
        prog='weaver-data-build', description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    commands = parser.add_subparsers(dest='command', required=True)

    cut = commands.add_parser('refseq', help='cut one RefSeq annotation release into a shard')
    cut.add_argument('--assembly', choices=sorted(_ASSEMBLIES), required=True)
    cut.add_argument('--release', required=True, help="the publisher's release identifier, RS_2024_08")
    cut.add_argument('--annotation', type=pathlib.Path, required=True, help='GFF3, gzipped')
    cut.add_argument('--transcripts', type=pathlib.Path, required=True, help='RNA FASTA, gzipped')
    cut.add_argument('--proteins', type=pathlib.Path, required=True, help='protein FASTA, gzipped')
    cut.add_argument(
        '--alignments',
        type=pathlib.Path,
        action='append',
        required=True,
        help='an alignment BAM with its .bai; repeat for model alignments',
    )
    cut.add_argument('--hgnc', type=pathlib.Path, required=True, help="HGNC's complete set, TSV")
    cut.add_argument('--mane', type=pathlib.Path, required=True, help='the MANE summary, gzipped TSV')
    cut.add_argument(
        '--records', type=pathlib.Path, required=True, help="the transcripts' GenBank records (rna.gbff), gzipped"
    )
    cut.add_argument('--shards', type=pathlib.Path, required=True, help='the directory the shard is written into')
    cut.add_argument(
        '--prefix', default='', help='prepended to the shard name, for a layout that orders shards by name'
    )
    cut.set_defaults(run=_refseq)

    index = commands.add_parser('index', help='write the index and manifest over shards, in the order given')
    index.add_argument('--assembly', choices=sorted(_ASSEMBLIES), required=True)
    index.add_argument('--out', type=pathlib.Path, required=True)
    index.add_argument('shards', type=pathlib.Path, nargs='+')
    index.set_defaults(run=_index)

    genome = commands.add_parser('genome', help='cut a genome FASTA into blocks and catalogue its sequences')
    genome.add_argument('--assembly', choices=sorted(_ASSEMBLIES), required=True)
    genome.add_argument('--fasta', type=pathlib.Path, required=True)
    genome.add_argument('--out', type=pathlib.Path, required=True)
    genome.set_defaults(run=_genome)

    args = parser.parse_args(argv)
    failure: build.BuildError | None = None
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always', build.BuildWarning)
        try:
            args.run(args)
        except build.BuildError as error:
            failure = error
    # Reported only once recording has stopped: a warning shown while recording is recorded again.
    _report(caught)
    if failure is not None:
        print(f'FAILED: {failure}', file=sys.stderr)
        raise SystemExit(1) from failure


def _report(caught: list[warnings.WarningMessage]) -> None:
    """Each build warning as one line on stderr, then their count; any other warning as Python shows it."""
    ours = 0
    for w in caught:
        if issubclass(w.category, build.BuildWarning):
            ours += 1
            print(f'warning: {w.message}', file=sys.stderr)
        else:
            warnings.showwarning(w.message, w.category, w.filename, w.lineno)
    if ours:
        print(f'{ours} build warnings', file=sys.stderr)


if __name__ == '__main__':
    main()

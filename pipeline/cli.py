"""One command-line interface for inputs, computation and manuscript outputs."""
import argparse
from pathlib import Path
import subprocess
import sys
from . import __version__


def parser():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--version', action='version', version=__version__)
    sub = p.add_subparsers(dest='action', required=True)
    v = sub.add_parser('verify', help='Verify the upload package, or a completed scientific run')
    v.add_argument('--work', type=Path)
    inp = sub.add_parser('prepare-inputs', help='Reuse verified caches; optionally download pinned sources')
    inp.add_argument('--output', type=Path, required=True)
    inp.add_argument('--cache-root', type=Path, action='append', default=[])
    inp.add_argument('--sources', type=Path, help='Manifest containing real URLs for the exact prepared inputs')
    inp.add_argument('--allow-downloads', action='store_true')
    inp.add_argument('--archive', type=Path, help='Import a pack-inputs tar archive without network access')
    check = sub.add_parser('check-inputs', help='Verify the exact scientific inputs before fitting')
    check.add_argument('--inputs', type=Path, required=True)
    check.add_argument('--bulk-cache', type=Path)
    pack = sub.add_parser('pack-inputs', help='Create a verified local data archive for deposition or transfer')
    pack.add_argument('--inputs', type=Path, required=True)
    pack.add_argument('--output', type=Path, required=True)
    run = sub.add_parser('run', help='Run the manuscript experiments and export their new results')
    from .run import STAGES
    run.add_argument('--stage', choices=('all', *STAGES), default='all')
    run.add_argument('--inputs', type=Path, required=True)
    run.add_argument('--work', type=Path, required=True)
    run.add_argument('--bulk-cache', type=Path)
    run.add_argument('--workers', type=int, default=36)
    run.add_argument('--rscript', default='Rscript')
    run.add_argument('--r-library', type=Path, required=True)
    run.add_argument('--resume', action='store_true')
    run.add_argument('--dry-run', action='store_true')
    run.add_argument('--output', type=Path, help='Publication output; defaults to WORK/final')
    run.add_argument('--formats', nargs='+', choices=('pdf','png','eps','tiff'), default=['pdf'])
    run.add_argument('--workbook', action='store_true', help='Also generate the optional Excel supplement')
    for name in ('export', 'figures', 'tables'):
        out = sub.add_parser(name, help='Export explicitly selected computed or published data')
        source = out.add_mutually_exclusive_group(required=True)
        source.add_argument('--work', type=Path, help='A complete pipeline run; every stage receipt is checked')
        source.add_argument('--published', action='store_true', help='Replay the archived manuscript data, without refitting')
        out.add_argument('--output', type=Path, required=True)
        if name != 'tables':
            out.add_argument('--formats', nargs='+', choices=('pdf','png','eps','tiff'), default=['pdf'])
        out.add_argument('--workbook', action='store_true')
    state = sub.add_parser('status', help='Read run progress without changing the run')
    state.add_argument('--work', type=Path, required=True)
    package = sub.add_parser('package', help='Build a checksum-verified GitHub/source ZIP')
    package.add_argument('--output', type=Path, required=True)
    return p


def main(argv=None):
    args = parser().parse_args(argv)
    from .package import verify
    if args.action == 'verify':
        verify()
        if args.work:
            from .run import validate_run
            print('Verified scientific run:', validate_run(args.work))
    elif args.action == 'prepare-inputs':
        from .inputs import prepare, unpack_inputs
        if args.archive:
            unpack_inputs(args.archive, args.output)
        prepare(args.output, args.cache_root, args.sources, args.allow_downloads)
    elif args.action == 'check-inputs':
        from .inputs import verify_inputs
        verify_inputs(args.inputs, args.bulk_cache)
    elif args.action == 'pack-inputs':
        from .inputs import pack_inputs
        print(pack_inputs(args.inputs, args.output))
    elif args.action == 'run':
        if not args.dry_run:
            verify()
        from .run import execute
        execute(args)
    elif args.action == 'status':
        from .run import status
        status(args.work)
    elif args.action == 'package':
        from .package import archive
        archive(args.output)
    else:
        verify()
        from .publication import export_published, export_run
        formats = [] if args.action == 'tables' else args.formats
        if args.published:
            export_published(args.output, formats, args.workbook)
        else:
            export_run(args.work, args.output, formats, args.workbook)
    return 0


def entry():
    try:
        return main()
    except (ValueError, OSError, subprocess.CalledProcessError) as error:
        print('STOP:', error, file=sys.stderr)
        return 1
    except KeyboardInterrupt:
        print('Interrupted. Resume with the same --inputs, --work and --resume.', file=sys.stderr)
        return 130

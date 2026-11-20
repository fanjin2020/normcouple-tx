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
    raw = sub.add_parser('rebuild-inputs', help='Reconstruct exact banks from pinned author raw-count objects')
    raw.add_argument('--source-dir', type=Path, required=True)
    raw.add_argument('--output', type=Path, required=True)
    raw.add_argument('--allow-downloads', action='store_true')
    raw.add_argument('--rscript', default='Rscript')
    raw.add_argument('--r-library', type=Path, action='append', default=[])
    check = sub.add_parser('check-inputs', help='Verify the exact scientific inputs before fitting')
    check.add_argument('--inputs', type=Path, required=True)
    check.add_argument('--bulk-cache', type=Path)
    pack = sub.add_parser('pack-inputs', help='Create a verified local data archive for deposition or transfer')
    pack.add_argument('--inputs', type=Path, required=True)
    pack.add_argument('--output', type=Path, required=True)
    run = sub.add_parser('run', help='Run the manuscript experiments and export their new results')
    from .run import STAGES
    run.add_argument('--stage', choices=('all','estimators', *STAGES), default='all')
    run.add_argument('--profile', choices=('full','acceptance'), default='full',
                     help='Acceptance uses real inputs with two repetitions; its outputs are not manuscript evidence')
    run.add_argument('--inputs', type=Path, required=True)
    run.add_argument('--work', type=Path, required=True)
    run.add_argument('--bulk-cache', type=Path)
    run.add_argument('--external-bulk-cache', type=Path)
    run.add_argument('--allow-downloads', action='store_true', help='Download only the pinned external clinical cohort if missing')
    run.add_argument('--workers', type=int, default=36)
    run.add_argument('--rscript', default='Rscript')
    run.add_argument('--r-library', type=Path, required=True)
    run.add_argument('--music-library', type=Path, help='Additional existing R library containing the pinned MuSiC implementation')
    run.add_argument('--resume', action='store_true')
    run.add_argument('--dry-run', action='store_true')
    run.add_argument('--output', type=Path, help='Publication output; defaults to WORK/final')
    run.add_argument('--formats', nargs='+', choices=('pdf','png','eps','tiff','svg'), default=['pdf'])
    run.add_argument('--workbook', action='store_true', help='Also generate the optional Excel supplement')
    cross = sub.add_parser('cross-method', help='Fixed-response dtangle, official MuSiC and linear SVR comparison')
    cross.add_argument('--inputs', type=Path, required=True)
    cross.add_argument('--work', type=Path, required=True)
    cross.add_argument('--rscript', default='Rscript')
    cross.add_argument('--r-library', type=Path, action='append', required=True)
    cross.add_argument('--workers', type=int, default=4)
    cross.add_argument('--repetitions', type=int, default=200)
    cross.add_argument('--batch-size', type=int, default=2)
    cross.add_argument('--smoke', action='store_true', help='First donor per setting, at most two replicates')
    cross.add_argument('--resume', action='store_true')
    clinical = sub.add_parser('clinical', help='Fixed-outcome bulk analysis with MuSiC gene policies and linear SVR')
    clinical.add_argument('--inputs', type=Path, required=True)
    clinical.add_argument('--bulk-cache', type=Path, required=True)
    clinical.add_argument('--external-bulk-cache', type=Path, required=True)
    clinical.add_argument('--work', type=Path, required=True)
    clinical.add_argument('--rscript', default='Rscript')
    clinical.add_argument('--r-library', type=Path, action='append', required=True)
    clinical.add_argument('--allow-downloads', action='store_true')
    clinical.add_argument('--resume', action='store_true')
    clinical.add_argument('--estimates-from',type=Path,
                          help='Reuse sealed estimator predictions only after exact raw-count, R-code and plan checks; refit all downstream models')
    audit = sub.add_parser('sensitivity', help='Audit fixed outcomes across supplied composition estimates')
    audit.add_argument('--outcomes', type=Path, required=True)
    audit.add_argument('--design', type=Path, required=True)
    audit.add_argument('--fractions', type=Path, required=True)
    audit.add_argument('--baseline', required=True)
    audit.add_argument('--exposure', default='exposure')
    audit.add_argument('--covariates', nargs='*', default=[])
    audit.add_argument('--output', type=Path, required=True)
    evidence=sub.add_parser('evidence',help='Run all six targeted evidence checks from the sealed full estimator run')
    evidence.add_argument('--root',type=Path,default=Path.home()/'git/repository')
    evidence.add_argument('--clinical',type=Path,help='Completed clinical_final directory (read only)')
    evidence.add_argument('--cross',type=Path,help='Completed cross_full directory (read only)')
    evidence.add_argument('--imaging-cache',type=Path,help='Cache of the three checksum-pinned upstream files')
    evidence.add_argument('--output',type=Path,help='New evidence directory; verified checkpoints resume automatically')
    evidence.add_argument('--allow-downloads',action='store_true',help='Only the three public imaging study files')
    evidence.add_argument('--rscript')
    evidence.add_argument('--r-library',type=Path,action='append',default=[])
    example=sub.add_parser('example',help='Run the small fixed-response software demonstration only')
    example.add_argument('--output',type=Path,required=True)
    for name in ('export', 'figures', 'tables'):
        out = sub.add_parser(name, help='Export explicitly selected computed or published data')
        source = out.add_mutually_exclusive_group(required=True)
        source.add_argument('--work', type=Path, help='A complete pipeline run; every stage receipt is checked')
        source.add_argument('--published', action='store_true', help='Replay the archived manuscript data, without refitting')
        out.add_argument('--evidence',type=Path,help='Regenerate Tables S20 and S21 from verified clinical probe and matched imaging records')
        out.add_argument('--output', type=Path, required=True)
        if name != 'tables':
            out.add_argument('--formats', nargs='+', choices=('pdf','png','eps','tiff','svg'), default=['pdf'])
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
    elif args.action == 'rebuild-inputs':
        from .rebuild import execute
        execute(args)
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
    elif args.action == 'cross-method':
        from .cross_method import execute
        execute(args)
    elif args.action == 'clinical':
        from .clinical_estimators import execute
        execute(args)
    elif args.action == 'sensitivity':
        from .sensitivity import execute
        execute(args)
    elif args.action == 'evidence':
        from .evidence import execute
        execute(args)
    elif args.action == 'example':
        from .evidence import demo
        from .io import separate
        print(demo(separate(args.output)))
    elif args.action == 'package':
        from .package import archive
        archive(args.output)
    else:
        verify()
        from .publication import export_published, export_run
        formats = [] if args.action == 'tables' else args.formats
        if args.published:
            export_published(args.output, formats, args.workbook, evidence=args.evidence)
        else:
            export_run(args.work, args.output, formats, args.workbook, evidence=args.evidence)
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

"""Dependency-ordered execution of the experiments used in the final manuscript."""
import importlib.metadata
import os
from pathlib import Path
import shutil
import signal
import subprocess
import sys
import time
from . import __version__
from .inputs import DEFINITIONS, stage_bulk_cache, verify_inputs
from .io import ROOT, atomic_json, digest, inside, lock, read_json, receipt, separate, sha, verify_receipt, writer_active
from .runtime import materialize

STAGES = ('normalization', 'mechanism', 'marker_grid', 'random_panels',
          'coupling', 'common_support', 'marker_associations', 'cross_method', 'clinical_estimators')
DEPENDENCIES = {'coupling': ('random_panels',), 'common_support': ('coupling',)}
ENTRIES = {'normalization': 'normalization/run_revision.py', 'mechanism': 'mechanism/run_extension.py',
           'marker_grid': 'marker_grid/run.py', 'random_panels': 'random_panels/run.py',
           'coupling': 'common_support/baseline/run_coupling.py',
           'common_support': 'common_support/finish.py', 'marker_associations': 'marker_associations/run.py'}
REQUIRED = {
    'normalization': ('tables/source_summary.csv', 'protocol.json', 'RUN_STATUS.json'),
    'mechanism': ('final_results/summary.csv', 'final_results/donor_summary.csv', 'protocol.json', 'RUN_STATUS.json'),
    'marker_grid': ('final_results/grid_summary.csv', 'final_results/grid_paired_vs_TMM.csv',
                    'final_results/power_mechanism_gene_audit.csv', 'protocol.json', 'RUN_STATUS.json'),
    'random_panels': ('protocol.json', 'selected_markers.json', 'RUN_STATUS.json'),
    'coupling': ('final/analysis_scope.json', 'final/method_summary.csv', 'protocol.json', 'STATE.json',
                 'clinical/complete.json', 'clinical/LIBD/jobs.json', 'clinical/SEAAD/jobs.json'),
    'common_support': ('protocol.json', 'STATE.json', 'common_support/method_summary.csv',
        'common_support/donor_subset_results.csv', 'common_support/dose_retention.csv',
        'common_support/budget_by_class.csv.gz', 'common_support/complete.json',
        'matching/marker_matching_balance.csv', 'corrections/complete.json',
        'corrections/clinical_disease_associations.csv.gz', 'corrections/clinical_method_stability.csv',
        'corrections/clinical_design_diagnostics.csv', 'corrections/diagnostic_qualification_counts.csv',
        'corrections/diagnostic_qualified_performance.csv', 'clinical_multistep/complete.json',
        'clinical_multistep/finite_perturbations.csv', 'clinical_multistep/derivative_steps.csv',
        'clinical_multistep/linear_approximation.csv'),
    'marker_associations': ('protocol.json', 'STATUS.json', 'final_results/marker_gene_associations.csv.gz'),
    'cross_method': ('protocol.json', 'complete.json', 'final/receipt.json', 'final/summary.csv',
                     'final/metrics.csv.gz', 'final/contrasts.csv.gz','final/paired_precision.csv'),
    'clinical_estimators': ('protocol.json','complete.json','final/method_summary.csv',
                            'final/paired_sensitivity.csv','final/validation.json')}
STATUS_FILES = {s: 'RUN_STATUS.json' for s in STAGES}
STATUS_FILES.update(coupling='STATE.json', common_support='STATE.json', marker_associations='STATUS.json', cross_method='complete.json')
STATUS_FILES['clinical_estimators'] = 'complete.json'
SUCCESS = {'normalization': {'COMPLETE'}, 'mechanism': {'COMPLETE_REQUESTED_STAGES'},
           'marker_grid': {'COMPLETE_REQUESTED_STAGES'}, 'random_panels': {'COMPLETE'},
           'coupling': {'CONFIRMATION_COMPLETE'}, 'common_support': {'FOLLOWUP_COMPLETE'},
           'marker_associations': {'COMPLETE', 'COMPLETE_WITH_DECLARED_GAPS'}}


def selected_stages(stage):
    wanted = set(STAGES if stage == 'all' else
                 ('cross_method','clinical_estimators') if stage == 'estimators' else [stage])
    def include(name):
        for dependency in DEPENDENCIES.get(name, ()):
            wanted.add(dependency); include(dependency)
    for name in list(wanted):
        include(name)
    return [s for s in STAGES if s in wanted]


def commands(args, runtime):
    work = Path(args.work).expanduser().resolve()
    inputs = Path(args.inputs).expanduser().resolve()
    cache = work / '.input-cache/bulk'
    common = ['--project-root', str(work), '--workers', str(args.workers), '--rscript', str(args.rscript),
              '--r-library', str(Path(args.r_library).expanduser().resolve())]
    banks = ['--libd-dir', str(inputs / 'LIBD'), '--seaad-dir', str(inputs / 'SEAAD'),
             '--gene-map', str(inputs / 'gene_identity_map.csv'), '--external-cache', str(inputs / 'external')]
    plan = []
    smoke = getattr(args,'profile','full') == 'acceptance'
    for stage in selected_stages(args.stage):
        if stage == 'clinical_estimators':
            external = getattr(args,'external_bulk_cache',None) or inputs/'external_bulk'
            cmd = [sys.executable,str(ROOT/'reproduce.py'),'clinical','--inputs',str(inputs),
                   '--bulk-cache',str(cache),'--external-bulk-cache',str(external),
                   '--work',str(work/stage),'--rscript',str(args.rscript),
                   '--r-library',str(Path(args.r_library).expanduser().resolve())]
            if getattr(args,'music_library',None):
                cmd += ['--r-library',str(Path(args.music_library).expanduser().resolve())]
            if getattr(args,'allow_downloads',False):
                cmd += ['--allow-downloads']
            if getattr(args,'resume',False) and (work/stage/'protocol.json').exists():
                cmd += ['--resume']
            plan.append((stage,cmd))
            continue
        if stage == 'cross_method':
            cmd = [sys.executable, str(ROOT / 'reproduce.py'), 'cross-method', '--inputs', str(inputs),
                   '--work', str(work / stage), '--workers', str(min(args.workers, 8)),
                   '--repetitions', '200', '--batch-size', '4', '--rscript', str(args.rscript),
                   '--r-library', str(Path(args.r_library).expanduser().resolve())]
            if getattr(args, 'music_library', None):
                cmd += ['--r-library', str(Path(args.music_library).expanduser().resolve())]
            if smoke:
                cmd += ['--smoke']
            if getattr(args, 'resume', False) and (work / stage / 'protocol.json').exists():
                cmd += ['--resume']
            plan.append((stage, cmd))
            continue
        cmd = [sys.executable, str(Path(runtime) / ENTRIES[stage]), *common]
        if stage == 'marker_associations':
            cmd += ['--counts', str(inputs / 'external/GSE174367_snRNA-seq_filtered_feature_bc_matrix.h5'),
                    '--metadata', str(inputs / 'external/GSE174367_snRNA-seq_cell_meta.csv.gz'), '--offline']
        else:
            cmd += banks
        cmd += ['--output', str(work / stage)]
        if stage == 'mechanism':
            cmd += ['--only', 'simulations', '--no-input-export', '--offline']
        elif stage == 'marker_grid':
            cmd += ['--only', 'grid', '--offline']
        elif stage == 'random_panels':
            cmd += ['--offline']
        elif stage == 'coupling':
            cmd += ['--phase', 'confirm', '--mode', 'all', '--random-results', str(work / 'random_panels'),
                    '--cache', str(cache), '--offline', '--python', sys.executable]
        elif stage == 'common_support':
            cmd += ['--previous-results', str(work / 'coupling'), '--random-results', str(work / 'random_panels'),
                    '--cache', str(cache), '--python', sys.executable]
        if smoke and stage != 'marker_associations':
            cmd += ['--smoke']
            if stage == 'normalization':
                cmd += ['--repetitions','2']
        plan.append((stage, cmd))
    return plan


def framework_identity():
    paths = [ROOT / 'reproduce.py', ROOT / 'kernels/layout.json']
    paths += [p for p in (ROOT / 'pipeline').iterdir() if p.suffix in ('.py', '.R')]
    paths += list((ROOT / 'publication').glob('*.py'))
    paths += [ROOT / 'publication/table_schema.json', ROOT / 'publication/figure_captions.json']
    paths += sorted((ROOT/'publication/input_definitions').glob('*'))
    return digest({p.relative_to(ROOT).as_posix(): sha(p) for p in sorted(paths)})


def python_environment():
    expected = read_json(DEFINITIONS / 'computational_specification.json')['software']
    problems = []
    for name, version in expected.items():
        try:
            actual = importlib.metadata.version(name)
            __import__(name)
            if actual != version:
                problems.append(name + ': expected ' + version + ', found ' + actual)
        except (ImportError, OSError, importlib.metadata.PackageNotFoundError) as error:
            problems.append(name + ': ' + str(error))
    if problems:
        raise ValueError('Select the archived analysis Python (requirements-archived.txt); '
                         'no packages were installed.\n' + '\n'.join(problems))
    return expected


def child_environment(library):
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE='1')
    for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS',
                'BLIS_NUM_THREADS', 'RCPP_PARALLEL_NUM_THREADS'):
        env[key] = '1'
    env['R_LIBS'] = str(library) + os.pathsep + env.get('R_LIBS', '')
    env['R_LIBS_USER'] = str(library)
    return env


def run_child(command, log, env):
    log.parent.mkdir(parents=True, exist_ok=True)
    with log.open('a', encoding='utf-8') as stream:
        proc = subprocess.Popen(command, cwd=Path(command[1]).parent, env=env,
                                stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                text=True, errors='replace', start_new_session=os.name != 'nt')
        try:
            for line in proc.stdout:
                print(line, end='', flush=True); stream.write(line); stream.flush()
            if proc.wait():
                raise subprocess.CalledProcessError(proc.returncode, command)
        except BaseException:
            if proc.poll() is None:
                if os.name != 'nt':
                    os.killpg(proc.pid, signal.SIGTERM)
                else:
                    proc.terminate()
                try:
                    proc.wait(timeout=30)
                except subprocess.TimeoutExpired:
                    if os.name != 'nt':
                        os.killpg(proc.pid, signal.SIGKILL)
                    else:
                        proc.kill()
                    proc.wait()
            raise


def stage_files(work, stage, software_acceptance=False):
    folder = Path(work) / stage
    state = read_json(folder / STATUS_FILES[stage])
    if stage == 'cross_method':
        protocol = read_json(folder / 'protocol.json')
        units, reps, batch = len(protocol['units']), protocol['repetitions'], protocol['batch_size']
        if (reps < (2 if software_acceptance else 200) or state.get('datasets') != units*reps*2 or
                state.get('batches') != units*((reps+batch-1)//batch)):
            raise ValueError('Cross-method manuscript scope is incomplete')
        verify_receipt(folder / 'final', read_json(folder / 'final/receipt.json')['files'])
    elif stage == 'clinical_estimators':
        if any(state.get(k)!=v for k,v in dict(cohorts=2,references=2,methods=6,probes=3).items()):
            raise ValueError('Clinical estimator scope is incomplete')
        verify_receipt(folder,state['files'])
    elif state.get('status') not in SUCCESS[stage] | ({'SOFTWARE_SMOKE_ONLY'} if software_acceptance else set()):
        raise ValueError('Stage did not finish its manuscript scope: ' + stage + ': ' + str(state.get('status')))
    protocol = read_json(folder / 'protocol.json')
    if protocol.get('settings', {}).get('software_fixture') or protocol.get('software_fixture'):
        raise ValueError('A software fixture cannot be exported as research evidence')
    smoke = bool(protocol.get('software_smoke') or protocol.get('settings',{}).get('software_smoke'))
    if smoke and not software_acceptance:
        raise ValueError('Software acceptance results cannot be exported as research evidence')
    names = list(REQUIRED[stage])
    # Later stages consume sealed controls and clinical arrays, not just summary CSVs.
    if stage in ('random_panels', 'coupling'):
        for p in folder.rglob('*'):
            rel = p.relative_to(folder).as_posix()
            if p.is_file() and ('/batches/' in '/' + rel or rel.startswith('clinical/') or
                                p.name in ('complete.json', 'runtime.json', 'files_sha256.json')):
                names.append(rel)
    return receipt(folder, names)


def completed(work, stage, identity, strict=False):
    marker = Path(work) / '.pipeline/stages' / (stage + '.json')
    if not marker.is_file():
        return False
    record = read_json(marker)
    if record.get('run_signature') != identity:
        raise ValueError('Stage receipt belongs to another run: ' + stage)
    try:
        verify_receipt(Path(work) / stage, record['files'])
    except (ValueError, OSError):
        if strict:
            raise
        return False
    return True


def validate_run(work, allow_acceptance=False):
    work = Path(work).resolve()
    cfg = read_json(work / 'pipeline.json')
    identity = cfg.pop('signature')
    if digest(cfg) != identity:
        raise ValueError('Pipeline configuration checksum differs')
    if cfg.get('framework_source') != framework_identity():
        raise ValueError('Run was made with different framework code; export with that same code package')
    if cfg.get('profile') == 'acceptance' and not allow_acceptance:
        raise ValueError('Software acceptance is not a full scientific run')
    for stage in STAGES:
        if not completed(work, stage, identity, strict=True):
            raise ValueError('Required experiment is incomplete: ' + stage)
    return identity


def execute(args):
    work = separate(args.work, args.inputs)
    bulk_source = Path(args.bulk_cache).expanduser().resolve() if args.bulk_cache else Path(args.inputs).resolve() / 'bulk'
    separate(work, bulk_source)
    separate(args.output or work / 'final', args.inputs, bulk_source,
             *(work / s for s in STAGES), work / '.runtime', work / '.pipeline', work / '.input-cache')
    library = Path(args.r_library).expanduser().resolve()
    if work == library or work.is_relative_to(library) or library.is_relative_to(work):
        raise ValueError('Keep the R library outside the run directory')
    if args.workers < 1 or args.workers > 72:
        raise ValueError('Use 1 to 72 workers')
    if args.dry_run:
        for stage, cmd in commands(args, work / '.runtime'):
            print(stage, subprocess.list2cmdline(cmd))
        return
    if args.workbook:
        try:
            __import__('openpyxl')
        except ImportError as error:
            raise ValueError('--workbook requires optional requirements-export.txt; CSV/figures do not') from error
    if os.name == 'nt':
        raise ValueError('Run full scientific fits on Linux. Published-data export and software tests support Windows.')
    python = python_environment()
    inputs = verify_inputs(args.inputs, args.bulk_cache)
    external = read_json(DEFINITIONS/'external_bulk_sources.json')
    if 'clinical_estimators' in selected_stages(args.stage):
        from .cohorts import external_sources
        external = external_sources(args.external_bulk_cache or Path(args.inputs)/'external_bulk',args.allow_downloads)
    if work.exists() and any(work.iterdir()) and not args.resume:
        raise ValueError('Run directory is not empty. Use --resume for this same experiment.')
    work.mkdir(parents=True, exist_ok=True)
    with lock(work / '.pipeline.lock'):
        runtime = materialize(work / '.runtime')
        stage_bulk_cache(work / '.input-cache/bulk', bulk_source, inputs)
        env = child_environment(library)
        rscript = shutil.which(str(args.rscript))
        if not rscript:
            raise ValueError('Rscript was not found; no packages were installed')
        args.rscript = rscript
        if {'cross_method','clinical_estimators'} & set(selected_stages(args.stage)):
            import json
            from .cross_method import MUSIC_COMMIT
            libs = [str(library)] + ([str(Path(args.music_library).resolve())] if args.music_library else [])
            expression = '.libPaths(c(' + ','.join(json.dumps(x) for x in libs) + ',.libPaths()));stopifnot(identical(packageDescription("MuSiC")$RemoteSha,' + json.dumps(MUSIC_COMMIT) + '))'
            subprocess.run([rscript, '--vanilla', '-e', expression], env=env, check=True)
        r_record = work / '.pipeline/R_versions.json'
        r_record.parent.mkdir(parents=True, exist_ok=True)
        subprocess.run([rscript, '--vanilla', str(runtime / 'normalization/r_bridge.R'), '--check',
                        str(library), str(r_record)], env=env, check=True)
        r_versions = read_json(r_record)
        expected_r = read_json(DEFINITIONS / 'computational_specification.json')['R']
        if any(r_versions.get(k) != v for k, v in expected_r.items()):
            raise ValueError('R/package versions differ from the manuscript protocol: ' + str(r_versions))
        cfg = dict(framework=__version__, framework_source=framework_identity(), kernel_layout=sha(ROOT / 'kernels/layout.json'),
                   input_hashes=inputs, python=python, R=r_versions,profile=args.profile,
                   external_bulk=external)
        identity = digest(cfg)
        stored = dict(cfg, signature=identity)
        if (work / 'pipeline.json').exists() and read_json(work / 'pipeline.json') != stored:
            raise ValueError('Code, inputs or software changed. Use a new run directory.')
        atomic_json(work / 'pipeline.json', stored)
        try:
            for stage, cmd in commands(args, runtime):
                if completed(work, stage, identity):
                    print('VERIFIED COMPLETE', stage, flush=True)
                    continue
                atomic_json(work / 'pipeline_status.json', dict(status='RUNNING', stage=stage, pid=os.getpid(),
                            signature=identity, updated_unix=time.time()))
                print('RUN', stage, flush=True)
                run_child(cmd, work / '.pipeline/logs' / (stage + '.log'), env)
                files = stage_files(work, stage,args.profile=='acceptance')
                atomic_json(work / '.pipeline/stages' / (stage + '.json'),
                            dict(run_signature=identity, files=files, completed_unix=time.time()))
            if args.stage == 'all':
                validate_run(work,allow_acceptance=args.profile=='acceptance')
                from .publication import export_run
                export_run(work, args.output or work / 'final', args.formats, args.workbook,
                           software_acceptance=args.profile=='acceptance')
                state = 'SOFTWARE_ACCEPTANCE_COMPLETE' if args.profile=='acceptance' else 'COMPLETE'
            elif args.stage == 'estimators':
                from .publication import export_estimators
                export_estimators(work,args.output or work/'final',args.formats,args.workbook)
                state = 'SELECTED_STAGES_COMPLETE'
            else:
                state = 'SELECTED_STAGES_COMPLETE'
            atomic_json(work / 'pipeline_status.json', dict(status=state, signature=identity, updated_unix=time.time()))
        except BaseException as error:
            atomic_json(work / 'pipeline_status.json', dict(status='INTERRUPTED_OR_FAILED', signature=identity,
                        error=type(error).__name__ + ': ' + str(error), updated_unix=time.time()))
            raise


def status(work):
    work = Path(work).expanduser().resolve()
    if not (work / 'pipeline.json').exists():
        raise ValueError('No pipeline configuration found: ' + str(work))
    print('RUN', work)
    state = read_json(work / 'pipeline_status.json') if (work / 'pipeline_status.json').exists() else {}
    print('LAST RECORDED STATE', state.get('status', 'NOT_STARTED'))
    active = writer_active(work / '.pipeline.lock')
    print('WRITER', 'ACTIVE' if active else 'INACTIVE')
    for stage in STAGES:
        sealed = (work / '.pipeline/stages' / (stage + '.json')).is_file()
        item = work / stage / STATUS_FILES[stage]
        detail = read_json(item) if item.is_file() else {}
        label = 'RECEIPT_RECORDED' if sealed else detail.get('status', 'PENDING')
        counts = ' '.join(str(k) + '=' + str(detail[k]) for k in
                         ('completed_batches', 'planned_batches', 'total_batches', 'completed', 'total') if k in detail)
        print(stage, label, counts)
    if not active and state.get('status') not in ('COMPLETE', 'SELECTED_STAGES_COMPLETE','SOFTWARE_ACCEPTANCE_COMPLETE'):
        print('No active writer. Resume with the same work/input arguments; do not delete lock files.')
    print('Stage records include reused batches. Export verifies their file checksums.')

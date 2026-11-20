import hashlib
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
import io
from pathlib import Path
import socket
import sys
import tarfile
import tempfile
import threading
import unittest
from unittest.mock import patch

from pipeline.cli import parser
from pipeline.download import download
from pipeline.inputs import prepare, source_manifest, stage_bulk_cache, unpack_inputs
from pipeline.io import ROOT, atomic_json, digest, inside, lock, read_json, separate, sha, verify_receipt, writer_active
from pipeline.publication import AUXILIARY, collect, load_frame, table_frames
from pipeline.run import STAGES, commands, completed, framework_identity, selected_stages, stage_files, validate_run
from pipeline.runtime import materialize


class PlanTests(unittest.TestCase):
    def args(self, root):
        return parser().parse_args(['run', '--inputs', str(root/'inputs'), '--work', str(root/'work'),
                                   '--r-library', str(root/'R'), '--workers', '36'])

    def test_all_includes_parent_and_common_support(self):
        self.assertEqual(selected_stages('all'), list(STAGES))
        self.assertEqual(selected_stages('common_support'), ['random_panels','coupling','common_support'])
        self.assertNotIn('equal_budget', STAGES)

    def test_parent_is_created_by_same_run(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            plans = dict(commands(self.args(root), root/'runtime'))
            cmd = plans['common_support']
            self.assertEqual(Path(cmd[cmd.index('--previous-results')+1]), root/'work/coupling')
            self.assertEqual(Path(cmd[cmd.index('--random-results')+1]), root/'work/random_panels')
            own_paths = {sys.executable, str(ROOT/'reproduce.py')}
            self.assertTrue(all('/home/research/' not in token for cmd in plans.values()
                                for token in cmd if token not in own_paths))

    def test_export_requires_explicit_data_origin(self):
        with self.assertRaises(SystemExit):
            parser().parse_args(['export','--output','somewhere'])
        self.assertTrue(parser().parse_args(['export','--published','--output','somewhere']).published)

    def test_incomplete_run_never_falls_back_to_published_data(self):
        with tempfile.TemporaryDirectory() as tmp:
            body = {'framework':'fixture', 'framework_source':framework_identity(), 'input_hashes':{}}
            atomic_json(Path(tmp)/'pipeline.json', dict(body, signature=digest(body)))
            with self.assertRaisesRegex(ValueError, 'incomplete'):
                validate_run(tmp)

    def test_changed_framework_cannot_export_an_old_run(self):
        with tempfile.TemporaryDirectory() as tmp:
            body = {'framework_source':'changed'}
            atomic_json(Path(tmp)/'pipeline.json', dict(body, signature=digest(body)))
            with self.assertRaisesRegex(ValueError, 'different framework'):
                validate_run(tmp)

    def test_missing_stage_source_is_an_error(self):
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaisesRegex(ValueError, 'Missing analysis output'):
                collect({'normalization': Path(tmp)})

    def test_success_message_without_files_does_not_complete_stage(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            atomic_json(root/'normalization/RUN_STATUS.json', {'status':'COMPLETE'})
            atomic_json(root/'normalization/protocol.json', {'settings':{}})
            with self.assertRaises(ValueError):
                stage_files(root, 'normalization')

    def test_fixture_is_not_research_completion(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            atomic_json(root/'normalization/RUN_STATUS.json', {'status':'COMPLETE'})
            atomic_json(root/'normalization/protocol.json', {'settings':{'software_fixture':True}})
            with self.assertRaisesRegex(ValueError, 'fixture'):
                stage_files(root, 'normalization')

    def test_resume_rechecks_result_hash(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            p = root/'normalization/result.csv'; p.parent.mkdir(); p.write_text('valid')
            atomic_json(root/'.pipeline/stages/normalization.json',
                        {'run_signature':'test','files':{'result.csv':sha(p)}})
            self.assertTrue(completed(root, 'normalization', 'test'))
            p.write_text('changed')
            self.assertFalse(completed(root, 'normalization', 'test'))
            with self.assertRaises(ValueError):
                completed(root, 'normalization', 'test', strict=True)


class StorageTests(unittest.TestCase):
    def test_lock_status_does_not_confuse_stale_file_with_writer(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)/'run.lock'
            self.assertFalse(writer_active(path))
            self.assertFalse(path.exists())
            with lock(path):
                self.assertTrue(writer_active(path))
            self.assertFalse(writer_active(path))

    def test_kernel_files_have_single_canonical_copy(self):
        layout = read_json(ROOT/'kernels/layout.json')['files']
        self.assertGreater(len(layout), len({v['source'] for v in layout.values()}))
        by_sha = {}
        for item in layout.values():
            by_sha.setdefault(item['sha256'], set()).add(item['source'])
        self.assertTrue(all(len(v)==1 for v in by_sha.values()))

    def test_materialized_layout_preserves_frozen_source_hashes(self):
        with tempfile.TemporaryDirectory() as tmp:
            runtime = materialize(Path(tmp)/'runtime')
            for p in runtime.rglob('*.json'):
                if p.name in ('source_manifest.json','source_sha256.json','bundle_manifest.json','bundle_sha256.json'):
                    verify_receipt(p.parent, read_json(p))
            chosen = runtime/'normalization/revision_core.py'
            chosen.write_text('tampered')
            with self.assertRaisesRegex(ValueError, 'modified'):
                materialize(runtime)

    def test_no_path_traversal(self):
        with tempfile.TemporaryDirectory() as tmp:
            for name in ('../escape', '/absolute', 'a/../../escape', 'bad\\name'):
                with self.assertRaises(ValueError):
                    inside(tmp, name)
        with self.assertRaises(ValueError):
            separate(ROOT/'tests/results')

    def test_table_and_figure_share_one_physical_source(self):
        catalog = read_json(ROOT/'publication/data/catalog.json')['views']
        self.assertEqual(catalog['normalization.csv']['file'], catalog['TableS2.csv']['file'])
        self.assertFalse((ROOT/'publication/data/normalization.csv').exists())
        self.assertEqual({p.name for p in (ROOT/'publication/data').glob('*.csv*')},
                         {f'TableS{i}.csv' for i in range(1,22)} | set(AUXILIARY))

    def test_null_outcome_is_not_parsed_as_missing(self):
        f = load_frame(ROOT/'publication/data', 'common_support_summary.csv')
        self.assertTrue(f.outcome.eq('null').any())
        self.assertFalse(f.outcome.isna().any())

    def test_source_manifest_cannot_change_frozen_input(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)/'sources.json'
            atomic_json(path, {'files': {'LIBD/cells.csv': {'sha256':'0'*64,'url':'https://example.org/file'}}})
            with self.assertRaisesRegex(ValueError, 'changes a fixed'):
                source_manifest(path)

    def test_cached_import_needs_no_network_or_old_receipt(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); cache = root/'cache'; cache.mkdir()
            p = cache/'counts.mtx'; p.write_bytes(b'pinned input')
            spec = {'LIBD/counts.mtx': {'sha256':sha(p),'url':None}}
            with patch('pipeline.inputs.specification', return_value=spec), patch('pipeline.inputs.download') as network:
                prepare(root/'inputs', [cache])
                prepare(root/'inputs', [cache])
                network.assert_not_called()
                self.assertEqual((root/'inputs/LIBD/counts.mtx').read_bytes(), p.read_bytes())

    def test_interrupted_cache_copy_can_resume(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); cache = root/'cache'; cache.mkdir()
            source = cache/'counts.mtx'; source.write_bytes(b'complete verified matrix')
            partial = root/'inputs/LIBD/counts.mtx.importing'; partial.parent.mkdir(parents=True)
            partial.write_bytes(b'complete')
            spec = {'LIBD/counts.mtx':{'sha256':sha(source),'url':None}}
            with patch('pipeline.inputs.specification', return_value=spec):
                prepare(root/'inputs', [cache])
            self.assertEqual((partial.parent/'counts.mtx').read_bytes(), source.read_bytes())

    def test_bulk_cache_adapter_is_accepted_by_frozen_kernel(self):
        import importlib.util
        module_spec = importlib.util.spec_from_file_location('bulk_common_test', ROOT/'kernels/coupling/common.py')
        legacy = importlib.util.module_from_spec(module_spec); module_spec.loader.exec_module(legacy)
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); source = root/'old'; source.mkdir()
            p = source/'bulk.tar'; p.write_bytes(b'verified data')
            url = 'https://example.org/bulk.tar'
            spec = {'bulk/bulk.tar':{'sha256':sha(p), 'url':url}}
            with patch('pipeline.inputs.specification', return_value=spec):
                cache = stage_bulk_cache(root/'work/cache', source, {'bulk/bulk.tar':sha(p)})
            result = legacy.acquire(url, cache, 'bulk.tar', offline=True, expected=sha(p))
            self.assertEqual(result.read_bytes(), p.read_bytes())
            self.assertEqual([x.name for x in source.iterdir()], ['bulk.tar'])

    def test_unsafe_input_archive_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); archive = root/'bad.tar'
            with tarfile.open(archive, 'w') as tar:
                item = tarfile.TarInfo('../escape'); item.size=1
                tar.addfile(item, io.BytesIO(b'x'))
            with self.assertRaisesRegex(ValueError, 'Unexpected'):
                unpack_inputs(archive, root/'inputs')
            self.assertFalse((root/'escape').exists())


class DownloadTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        self.payload = b'verified bytes\n' * 4096
        self.expected = hashlib.sha256(self.payload).hexdigest()
        self.ranges = []
        self.ignore_range = False
        self.drop_first = False
        outer = self
        class Handler(BaseHTTPRequestHandler):
            def log_message(self, *args):
                pass
            def do_GET(self):
                value = self.headers.get('Range')
                outer.ranges.append(value)
                start = int(value.split('=')[1].split('-')[0]) if value and not outer.ignore_range else 0
                self.send_response(206 if start else 200)
                if start:
                    self.send_header('Content-Range', f'bytes {start}-{len(outer.payload)-1}/{len(outer.payload)}')
                body = outer.payload[start:]
                self.send_header('Content-Length', str(len(body))); self.end_headers()
                if outer.drop_first and len(outer.ranges)==1:
                    self.wfile.write(body[:2048]); self.wfile.flush()
                    self.connection.shutdown(socket.SHUT_RDWR); self.connection.close()
                else:
                    self.wfile.write(body)
        self.server = ThreadingHTTPServer(('127.0.0.1', 0), Handler)
        self.thread = threading.Thread(target=self.server.serve_forever, daemon=True)
        self.thread.start()
        self.url = 'http://127.0.0.1:%d/data' % self.server.server_port

    def tearDown(self):
        self.server.shutdown(); self.server.server_close(); self.thread.join(); self.tmp.cleanup()

    def test_existing_valid_file_reused_without_receipt(self):
        p = self.root/'data'; p.write_bytes(self.payload)
        download(self.url, p, self.expected)
        self.assertEqual(self.ranges, [])

    def test_partial_download_resumes(self):
        p = self.root/'data'
        p.with_name('data.partial').write_bytes(self.payload[:4096])
        download(self.url, p, self.expected)
        self.assertEqual(self.ranges[0], 'bytes=4096-')
        self.assertEqual(p.read_bytes(), self.payload)

    def test_server_ignoring_range_restarts_safely(self):
        self.ignore_range = True
        p = self.root/'data'; p.with_name('data.partial').write_bytes(self.payload[:4096])
        download(self.url, p, self.expected)
        self.assertEqual(p.read_bytes(), self.payload)

    def test_connection_loss_retries_from_saved_bytes(self):
        self.drop_first = True
        p = self.root/'data'
        with patch('pipeline.download.time.sleep'):
            download(self.url, p, self.expected)
        self.assertGreater(len(self.ranges), 1)
        self.assertIsNotNone(self.ranges[1])
        self.assertEqual(p.read_bytes(), self.payload)

    def test_wrong_checksum_never_becomes_final_file(self):
        p = self.root/'data'
        with self.assertRaisesRegex(ValueError, 'checksum differs'):
            download(self.url, p, '0'*64)
        self.assertFalse(p.exists())
        self.assertTrue(p.with_name('data.partial').exists())


if __name__ == '__main__':
    unittest.main()

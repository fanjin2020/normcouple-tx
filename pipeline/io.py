"""Checked paths, atomic records and run ownership."""
from contextlib import contextmanager
import hashlib
import json
import os
from pathlib import Path
import tempfile

ROOT = Path(__file__).resolve().parents[1]


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()).hexdigest()


def read_json(path):
    return json.loads(Path(path).read_text(encoding='utf-8-sig'))


def atomic_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, name = tempfile.mkstemp(prefix='.' + path.name + '.', dir=path.parent)
    try:
        with os.fdopen(fd, 'w', encoding='utf-8', newline='\n') as stream:
            json.dump(value, stream, indent=2, ensure_ascii=False, allow_nan=False)
            stream.write('\n')
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(name, path)
    finally:
        if os.path.exists(name):
            os.unlink(name)


def inside(root, relative):
    if not isinstance(relative, str) or '\\' in relative or not relative:
        raise ValueError('Invalid relative path: ' + str(relative))
    root = Path(root).resolve()
    path = (root / relative).resolve()
    if Path(relative).is_absolute() or path == root or not path.is_relative_to(root):
        raise ValueError('Path leaves its declared root: ' + relative)
    return path


def separate(output, *protected):
    out = Path(output).expanduser().resolve()
    for value in (ROOT, *protected):
        p = Path(value).expanduser().resolve()
        if out == p or out.is_relative_to(p) or p.is_relative_to(out):
            raise ValueError('Output overlaps a source/input tree: ' + str(p))
    return out


@contextmanager
def lock(path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('a+b') as stream:
        try:
            if os.name == 'nt':
                import msvcrt
                stream.seek(0)
                if not stream.read(1):
                    stream.write(b'0'); stream.flush()
                stream.seek(0)
                msvcrt.locking(stream.fileno(), msvcrt.LK_NBLCK, 1)
            else:
                import fcntl
                fcntl.flock(stream.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except OSError as error:
            raise ValueError('Another writer is using ' + str(path)) from error
        try:
            yield
        finally:
            if os.name == 'nt':
                stream.seek(0)
                msvcrt.locking(stream.fileno(), msvcrt.LK_UNLCK, 1)
            else:
                fcntl.flock(stream.fileno(), fcntl.LOCK_UN)


def receipt(root, paths):
    root = Path(root).resolve()
    result = {}
    for relative in sorted(set(paths)):
        p = inside(root, relative)
        if not p.is_file():
            raise ValueError('Missing required output: ' + str(p))
        result[relative] = sha(p)
    return result


def writer_active(path):
    path = Path(path)
    if not path.is_file():
        return False
    with path.open('rb') as stream:
        if os.name == 'nt':
            import msvcrt
            try:
                msvcrt.locking(stream.fileno(), msvcrt.LK_NBLCK, 1)
            except OSError:
                return True
            msvcrt.locking(stream.fileno(), msvcrt.LK_UNLCK, 1)
        else:
            import fcntl
            try:
                fcntl.flock(stream.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
            except BlockingIOError:
                return True
            fcntl.flock(stream.fileno(), fcntl.LOCK_UN)
    return False


def verify_receipt(root, files):
    if not files:
        raise ValueError('An empty output receipt cannot establish completion')
    for name, expected in files.items():
        p = inside(root, name)
        if not p.is_file() or sha(p) != expected:
            raise ValueError('Missing or changed output: ' + str(p))

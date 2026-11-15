"""Small, atomic, path-checked helpers for the independent continuation."""
import hashlib
import json
import os
from pathlib import Path
import shutil


def require(ok, message):
    if not ok:
        raise ValueError(message)


def read_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8-sig"))


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_name(path.name + ".partial")
    temp.write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False), encoding="utf-8")
    os.replace(temp, path)


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def inside(root, name):
    root = Path(root).resolve()
    require(not Path(name).is_absolute(), "Absolute path in portable manifest: " + str(name))
    path = (root / name).resolve()
    require(path != root and path.is_relative_to(root), "Path escapes manifest directory: " + str(name))
    return path


def copy_verified(src, dst, expected=None):
    src, dst = Path(src), Path(dst)
    expected = expected or sha(src)
    if dst.exists():
        require(dst.is_file() and sha(dst) == expected, "Refusing to overwrite different file: " + str(dst))
        return
    dst.parent.mkdir(parents=True, exist_ok=True)
    temp = dst.with_name(dst.name + ".partial")
    shutil.copyfile(src, temp)
    require(sha(temp) == expected, "Copy checksum mismatch: " + str(src))
    os.replace(temp, dst)


def resolve_results(path):
    path = Path(path).expanduser().resolve()
    direct = [path, path / "final_results"]
    for candidate in direct:
        if all((candidate / n).is_file() for n in ("summary.csv", "protocol.json", "gene_summary.csv")):
            return candidate
    candidates = []
    if path.is_dir():
        for child in path.iterdir():
            if child.is_dir() and (child / "summary.csv").is_file() and (child / "protocol.json").is_file():
                candidates.append(child)
    require(len(candidates) == 1, "No unique completed result folder below " + str(path))
    return candidates[0]

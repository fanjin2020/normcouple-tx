"""I/O shared by the two separately identified Paper 3 additions."""
import hashlib
import json
import os
from pathlib import Path
import time
import urllib.request

HERE = Path(__file__).resolve().parent


def require(ok, message):
    if not ok:
        raise ValueError(message)


def sha(path):
    with Path(path).open("rb") as f:
        return hashlib.file_digest(f, "sha256").hexdigest()


def read_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8-sig"))


def save_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".partial")
    tmp.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")
    os.replace(tmp, path)


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()).hexdigest()


def acquire(url, cache, filename, offline=False, expected=None, max_bytes=150_000_000):
    cache = Path(cache)
    cache.mkdir(parents=True, exist_ok=True)
    dest, receipt = cache / filename, cache / (filename + ".download.json")
    if dest.exists():
        require(receipt.exists(), "Unverified cache file; move it aside or use a new cache: " + str(dest))
        meta = read_json(receipt)
        require(meta["url"] == url and sha(dest) == meta["sha256"], "Cached source identity changed: " + filename)
        require(not expected or meta["sha256"] == expected, "Cache does not match pinned source: " + filename)
        print("REUSE DOWNLOAD", filename, flush=True)
        return dest
    require(not offline, "Missing public file in offline mode: " + filename)
    tmp = dest.with_name(filename + ".partial")
    for attempt in range(4):
        try:
            print("DOWNLOAD", url, flush=True)
            request = urllib.request.Request(url, headers={"User-Agent": "Paper3-reproducible-research/1.0"})
            with urllib.request.urlopen(request, timeout=90) as source, tmp.open("wb") as target:
                declared = int(source.headers.get("Content-Length", 0))
                require(declared <= max_bytes, "Download exceeds size limit")
                size = 0
                for block in iter(lambda: source.read(1024*1024), b""):
                    target.write(block)
                    size += len(block)
                    require(size <= max_bytes, "Download exceeds size limit")
                require(not declared or size == declared, "Truncated source download")
                headers = {k: source.headers.get(k) for k in ("ETag", "Last-Modified", "Content-Length")}
            checksum = sha(tmp)
            require(not expected or checksum == expected, "Public source changed; review instead of silently updating: " + filename)
            os.replace(tmp, dest)
            save_json(receipt, dict(url=url, sha256=checksum, bytes=size, retrieved_unix=time.time(), headers=headers))
            return dest
        except (OSError, ValueError):
            if attempt == 3:
                raise
            time.sleep(2**attempt)


def protocol(path, body):
    value = dict(body, signature=digest(body))
    if Path(path).exists():
        require(read_json(path) == value, "Inputs, code, environment or settings changed; use a new output directory.")
    save_json(path, value)
    return value


def file_manifest(folder):
    folder = Path(folder)
    save_json(folder / "files_sha256.json", {p.relative_to(folder).as_posix(): sha(p)
        for p in sorted(folder.rglob("*")) if p.is_file() and p.name != "files_sha256.json"})

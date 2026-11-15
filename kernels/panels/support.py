"""Standard-library protocol, checkpoint and progress utilities."""
import hashlib
import json
import math
import os
from pathlib import Path
import socket
import statistics
import time

PACKAGE = "paper3_random_marker"
VERSION = "1.0.0"


def read_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def atomic_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_name(path.name + ".%s.partial" % os.getpid())
    with temp.open("w", encoding="utf-8", newline="\n") as stream:
        json.dump(value, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temp, path)


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"),
                                    allow_nan=False).encode()).hexdigest()


def file_sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def tasks_for(protocol):
    s = protocol["settings"]
    for start in range(0, s["repetitions"], s["batch_size"]):
        for block in range(math.ceil(len(s["subset_ids"]) / s["subset_batch_size"])):
            for unit in protocol["units"]:
                yield unit, start, block


def task_folder(output, task):
    unit, start, block = task
    return Path(output) / "units" / unit["unit_id"] / ("s%03d_r%05d" % (block, start))


def task_work(protocol, task):
    _, start, block = task
    s = protocol["settings"]
    count = min(s["subset_batch_size"], len(s["subset_ids"]) - block * s["subset_batch_size"])
    return min(s["batch_size"], s["repetitions"] - start) * count


def checkpoint(folder, signature, verify=True):
    folder = Path(folder)
    receipt = folder / "complete.json"
    if not receipt.exists():
        return None
    value = read_json(receipt)
    if value["signature"] != signature or not value.get("files"):
        raise ValueError("Checkpoint protocol mismatch: " + str(folder))
    for name, meta in value["files"].items():
        path = (folder / name).resolve()
        if not path.is_relative_to(folder.resolve()) or not path.is_file():
            raise ValueError("Missing/unsafe checkpoint file: " + str(path))
        if path.stat().st_size != meta["bytes"] or (verify and file_sha(path) != meta["sha256"]):
            raise ValueError("Damaged checkpoint (not silently reused): " + str(path))
    return value


def seal(folder, names, signature, **metadata):
    folder = Path(folder)
    files = {name: {"bytes": (folder / name).stat().st_size,
                    "sha256": file_sha(folder / name)} for name in names}
    atomic_json(folder / "complete.json", dict(signature=signature, files=files, **metadata))


def state(output, status, **extra):
    atomic_json(Path(output) / "RUN_STATUS.json", dict(status=status, pid=os.getpid(),
                hostname=socket.gethostname(), updated_unix=time.time(), **extra))


def lock_is_active(output):
    path = Path(output) / "run.lock"
    if os.name == "nt" or not path.exists():
        return None
    import fcntl
    with path.open("rb") as stream:
        try:
            fcntl.flock(stream, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            return True
        else:
            fcntl.flock(stream, fcntl.LOCK_UN)
            return False


def progress(output, verify=False):
    output = Path(output)
    if not (output / "protocol.json").exists():
        print("No protocol yet: preparation has not finished, or this output path is wrong.")
        if (output / "RUN_STATUS.json").exists():
            print(read_json(output / "RUN_STATUS.json"))
        return 0
    cfg = read_json(output / "protocol.json")
    current = read_json(output / "RUN_STATUS.json") if (output / "RUN_STATUS.json").exists() else {}
    total = done = work_done = work_total = 0
    finished, damaged, units = [], [], {}
    for task in tasks_for(cfg):
        total += 1
        work_total += task_work(cfg, task)
        unit = task[0]["unit_id"]
        units.setdefault(unit, [0, 0])[1] += 1
        try:
            result = checkpoint(task_folder(output, task), cfg["signature"], verify)
            if result:
                done += 1
                units[unit][0] += 1
                work_done += result["work_items"]
                finished.append(result)
        except (OSError, ValueError, KeyError) as exc:
            damaged.append(str(exc))
    active = lock_is_active(output)
    print("Recorded status:", current.get("status", "UNKNOWN"))
    print("Writer lock:", "ACTIVE" if active else "INACTIVE" if active is False else "not checked on this OS")
    print("%s batches: %d/%d (%.1f%%)" % ("Hash-verified" if verify else "Receipt/size-checked", done, total, 100 * done / total))
    print("Completed repetition-subset pairs: %d/%d (%.1f%%)" % (work_done, work_total, 100 * work_done / work_total))
    print("Complete units: %d/%d; damaged receipts: %d" % (sum(a == b for a, b in units.values()), len(units), len(damaged)))
    print("Workers in last launch:", current.get("workers", "unknown"))
    if finished and work_done < work_total and active is not False:
        recent = sorted(finished, key=lambda x: x["finished_unix"])[-72:]
        seconds = statistics.median(x["elapsed_seconds"] / x["work_items"] for x in recent)
        workers = max(1, min(current.get("workers", 1), total - done))
        estimate = (work_total - work_done) * seconds / workers / 3600
        print("Rough simulation-only remaining time: %.2f h at %d workers." % (estimate, workers))
        print("Based on completed work, not a guarantee; excludes preparation and summaries.")
    elif work_done < work_total and active is False:
        print("No active writer. Re-run the same command to resume; do not delete run.lock.")
    for message in damaged[:5]:
        print(message)
    return 2 if damaged else 0

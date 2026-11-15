"""Pinned downloads with verified cache reuse and HTTP range recovery."""
import os
import http.client
from pathlib import Path
import re
import time
import urllib.error
import urllib.parse
import urllib.request
from .io import atomic_json, lock, read_json, sha


def download(url, destination, expected, attempts=3, timeout=180):
    parsed = urllib.parse.urlsplit(url)
    if parsed.scheme != 'https' and not (parsed.scheme == 'http' and parsed.hostname in ('127.0.0.1', 'localhost')):
        raise ValueError('Use HTTPS source URLs')
    if parsed.username or parsed.password or not re.fullmatch('[a-f0-9]{64}', expected):
        raise ValueError('A public URL and a pinned SHA256 are required')
    dest = Path(destination)
    dest.parent.mkdir(parents=True, exist_ok=True)
    with lock(dest.with_name('.' + dest.name + '.lock')):
        if dest.exists():
            if sha(dest) != expected:
                raise ValueError('Existing file has a different checksum: ' + str(dest))
            print('REUSE', dest.name, flush=True)
            return dest
        part = dest.with_name(dest.name + '.partial')
        metadata = dest.with_name(dest.name + '.partial.json')
        identity = {'url': url, 'sha256': expected}
        if metadata.exists() and read_json(metadata) != identity:
            raise ValueError('Partial file belongs to another source: ' + str(part))
        atomic_json(metadata, identity)
        if part.exists() and sha(part) == expected:
            os.replace(part, dest); metadata.unlink()
            return dest
        for attempt in range(attempts):
            offset = part.stat().st_size if part.exists() else 0
            headers = {'Accept-Encoding': 'identity', 'User-Agent': 'paper3-reproduction/2.0'}
            if offset:
                headers['Range'] = 'bytes=%d-' % offset
            try:
                request = urllib.request.Request(url, headers=headers)
                with urllib.request.urlopen(request, timeout=timeout) as response:
                    status = response.status
                    if status == 206:
                        match = re.fullmatch(r'bytes (\d+)-(\d+)/(\d+|\*)', response.headers.get('Content-Range', ''))
                        if not match or int(match[1]) != offset:
                            raise ValueError('Server returned an inconsistent byte range')
                        mode = 'ab'
                    elif status == 200:
                        mode = 'wb'  # The server ignored Range: safely restart this pinned object.
                    else:
                        raise OSError('Unexpected HTTP status ' + str(status))
                    print('DOWNLOAD', dest.name, 'offset', offset if mode == 'ab' else 0, flush=True)
                    received = 0
                    with part.open(mode) as stream:
                        reader = getattr(response, 'read1', response.read)
                        for block in iter(lambda: reader(8 * 1024 * 1024), b''):
                            stream.write(block); received += len(block)
                        stream.flush(); os.fsync(stream.fileno())
                    length = response.headers.get('Content-Length')
                    if length is not None and received != int(length):
                        raise OSError('Incomplete response body')
                if sha(part) != expected:
                    raise ValueError('Downloaded checksum differs; partial file retained: ' + str(part))
                os.replace(part, dest)
                metadata.unlink()
                return dest
            except (OSError, urllib.error.URLError, http.client.HTTPException) as error:
                if isinstance(error, urllib.error.HTTPError) and error.code == 416 and part.exists() and sha(part) == expected:
                    os.replace(part, dest); metadata.unlink()
                    return dest
                if attempt + 1 == attempts:
                    raise OSError('Download interrupted; rerun to resume ' + str(part)) from error
                time.sleep(min(2 ** attempt, 8))
    raise RuntimeError('Download did not complete')

"""Verified compact result export, available without third-party Python modules."""
import os
from pathlib import Path
import zipfile
from common import require,read_json,sha,save_json,file_manifest


def export_review(output):
    output=Path(output).resolve()
    final=output/'final'
    require((final/'STATUS.json').is_file(),'No finished compact result set; resume the original run first')
    status=read_json(final/'STATUS.json')
    require(status['status'] in ('SOFTWARE_FIXTURE_ONLY','PILOT_COMPLETE','CONFIRMATION_COMPLETE',
                                 'SIMULATION_ONLY_COMPLETE','CLINICAL_ONLY_COMPLETE'),
            'Cannot export incomplete results')
    manifest=read_json(final/'files_sha256.json')
    for name,expected in manifest.items():
        path=(final/name).resolve()
        require(path.is_relative_to(final) and path.is_file() and sha(path)==expected,'Result identity mismatch: '+name)
    cfg=read_json(output/'protocol.json')
    destination=output.parent/('p3r_'+cfg['signature'][:10]+'.zip')
    temporary=destination.with_name(destination.name+'.partial')
    names=[*manifest,'files_sha256.json']
    require(sum((final/name).stat().st_size for name in names)<250_000_000,'Compact output exceeded 250 MB; inspect before export')
    with zipfile.ZipFile(temporary,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as archive:
        for name in sorted(names):
            archive.write(final/name,arcname='p3_review/'+name)
    os.replace(temporary,destination)
    save_json(destination.with_suffix('.sha256.json'),dict(file=destination.name,sha256=sha(destination),bytes=destination.stat().st_size))
    print('REVIEW ZIP',destination,flush=True)
    print('FULL BATCHES AND RAW MATRICES ARE NOT IN THIS ZIP',flush=True)
    return 0

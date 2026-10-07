"""Import downloaded AlphaFold Server ZIPs without scoring against truth.

Run from the repository root, for example:
  python agent2/baselines/server_submission_2026-10-06/import_server_outputs.py \
    --archive T1147=/path/to/T1147.zip --job-url T1147=https://alphafoldserver.com/fold/ID \
    --archive T1133=/path/to/T1133.zip --job-url T1133=https://alphafoldserver.com/fold/ID \
    --archive T1183=/path/to/T1183.zip --job-url T1183=https://alphafoldserver.com/fold/ID

Alternatively --inputs-json file.json accepts a target-keyed object whose values
contain archive, job_url, and optional job_id. Output is always a fresh directory.
All original archives and all archive members, including terms, MSAs and templates,
are preserved. Primary selection uses only exported ranking_score, then lowest
sample index for rounded ties; this script never reads experimental coordinates.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
from io import StringIO
import json
import math
from pathlib import Path, PurePosixPath
import re
import shutil
import stat
from urllib.parse import urlparse
import zipfile

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
TARGETS = ('T1147', 'T1133', 'T1183')
SCHEMA_SOURCE = 'https://github.com/google-deepmind/alphafold/blob/main/server/README.md'


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def relative(path):
    return Path(path).resolve().relative_to(ROOT).as_posix()


def save_json(path, value):
    with Path(path).open('x', encoding='utf-8') as handle:
        json.dump(value, handle, indent=2, allow_nan=False)
        handle.write('\n')


def fasta_sequence(path):
    lines = Path(path).read_text(encoding='utf-8').splitlines()
    require(sum(line.startswith('>') for line in lines) == 1, 'Expected one FASTA record')
    sequence = ''.join(line.strip() for line in lines if line.strip() and not line.startswith('>'))
    require(sequence and not set(sequence) - set('ACDEFGHIKLMNPQRSTVWY'), 'Noncanonical FASTA')
    return sequence


def keyed(values):
    result = {}
    for value in values or []:
        target, separator, payload = value.partition('=')
        require(separator and target in TARGETS and target not in result, 'Use unique TARGET=value arguments')
        result[target] = payload
    return result


def validate_request(archive, target, sequence):
    names = [name for name in archive.namelist() if name.lower().endswith('_job_request.json')]
    require(len(names) == 1, f'{target}: expected one exported job request')
    payload = json.loads(archive.read(names[0]))
    jobs = payload if isinstance(payload, list) else [payload]
    require(len(jobs) == 1 and isinstance(jobs[0], dict), f'{target}: expected one job in archive')
    job = jobs[0]
    expected_name = f'P2_{target}_AF3server_20261006'
    require(str(job.get('name', '')).casefold() == expected_name.casefold(), f'{target}: unexpected job name')
    require(job.get('dialect') == 'alphafoldserver', f'{target}: archive is not a server request')
    require(type(job.get('version')) is int and job['version'] >= 1, f'{target}: unknown exported schema version')
    # Actual server archives export version 3; preserve it rather than rewriting to prepared version 1.
    require(job.get('modelSeeds') == ['1'], f'{target}: actual modelSeeds must be exactly ["1"]')
    entities = job.get('sequences', [])
    require(len(entities) == 1 and set(entities[0]) == {'proteinChain'}, f'{target}: expected one protein entity')
    protein = entities[0]['proteinChain']
    require(protein.get('sequence') == sequence and type(protein.get('count')) is int and protein['count'] == 1,
            f'{target}: downloaded request sequence/count differs from exact target FASTA')
    for key in ('unpairedMsa', 'templates', 'glycans', 'modifications'):
        require(not protein.get(key), f'{target}: unexpected supplied {key} in exported request')
    settings = {key: protein.get(key) for key in ('useStructureTemplate', 'maxTemplateDate')}
    differences = []
    if job['version'] != 1:
        differences.append({'field': 'exported_schema_version', 'prepared': 1, 'actual': job['version']})
    for key, intended in [('useStructureTemplate', True), ('maxTemplateDate', '2021-09-30')]:
        if key not in protein:
            differences.append({'field': key, 'prepared': intended, 'actual': 'omitted in exported request'})
        elif protein[key] != intended:
            differences.append({'field': key, 'prepared': intended, 'actual': protein[key]})
    settings['maxTemplateDate_record'] = (protein['maxTemplateDate'] if 'maxTemplateDate' in protein
                                        else 'server default; not serialized in downloaded request')
    settings['effective_max_template_date_verified'] = 'maxTemplateDate' in protein
    settings['documented_useStructureTemplate_default'] = True
    settings['documented_maximum_allowed_template_date'] = '2025-02-03'
    settings['documented_default_template_cutoff'] = None
    settings['documentation_note'] = ('The cited JSON README documents the maximum allowed date, not the default '
                                       'maxTemplateDate. An omitted effective cutoff cannot be inferred from that maximum.')
    return names[0], job, settings, differences


def archive_plan(path, target, sequence):
    with zipfile.ZipFile(path) as archive:
        require(archive.testzip() is None, f'{target}: ZIP CRC validation failed')
        seen = set()
        for info in archive.infolist():
            name = info.filename.replace('\\', '/')
            parts = PurePosixPath(name)
            require(not parts.is_absolute() and '..' not in parts.parts and not any(':' in p for p in parts.parts),
                    f'{target}: unsafe ZIP member {name}')
            require(not stat.S_ISLNK(info.external_attr >> 16), f'{target}: ZIP symlink is not allowed')
            require(name.casefold() not in seen, f'{target}: duplicate ZIP member {name}')
            seen.add(name.casefold())
        request_name, request, settings, differences = validate_request(archive, target, sequence)
        terms = [name for name in archive.namelist() if PurePosixPath(name).name.lower() == 'terms_of_use.md']
        require(terms, f'{target}: terms_of_use.md missing from server archive')
        samples = []
        for name in archive.namelist():
            match = re.fullmatch(r'(.*)summary_confidences_(\d+)\.json', name)
            if not match:
                continue
            model = f'{match.group(1)}model_{match.group(2)}.cif'
            full_data = f'{match.group(1)}full_data_{match.group(2)}.json'
            require(model in archive.namelist() and full_data in archive.namelist(), f'{target}: incomplete sample')
            confidence = json.loads(archive.read(name))
            score = confidence.get('ranking_score')
            require(type(score) in (int, float) and math.isfinite(score), f'{target}: invalid ranking_score')
            samples.append({'sample_index': int(match.group(2)), 'ranking_score': score, 'model_file': model,
                            'confidence_file': name, 'full_data_file': full_data})
        require(samples and len({s['sample_index'] for s in samples}) == len(samples), f'{target}: no unique ranked samples')
        samples.sort(key=lambda sample: (-sample['ranking_score'], sample['sample_index']))
        primary = samples[0]
        # Validate the selected model's sequence/chain after selection; no reference structure is opened.
        from Bio.PDB import MMCIFParser
        from Bio.SeqUtils import seq1
        structure = MMCIFParser(QUIET=True).get_structure(target, StringIO(archive.read(primary['model_file']).decode('utf-8')))
        require(len(structure) == 1, f'{target}: selected CIF must contain one coordinate model')
        chains = [(chain, [r for r in chain if 'CA' in r]) for chain in structure[0]]
        chains = [(chain, residues) for chain, residues in chains if residues]
        require(len(chains) == 1, f'{target}: selected CIF must have one protein chain')
        chain, residues = chains[0]
        require(''.join(seq1(r.resname) for r in residues) == sequence, f'{target}: selected CIF sequence differs')
        numbered = all(r.id[1] == i and not r.id[2].strip() for i, r in enumerate(residues, 1))
        return {'request_file': request_name, 'request': request, 'actual_settings': settings,
                'differences_from_prepared': differences, 'samples_ranked': samples, 'primary': primary,
                'chain': chain.id, 'mapping_mode': 'author_number_equals_target' if numbered else 'sequence',
                'terms_files': terms, 'member_count': len(archive.infolist())}


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--archive', action='append', help='TARGET=/absolute/path/to/download.zip')
    parser.add_argument('--job-url', action='append', help='TARGET=https://alphafoldserver.com/fold/JOB_ID')
    parser.add_argument('--job-id', action='append', help='TARGET=JOB_ID; otherwise derived from job URL')
    parser.add_argument('--inputs-json', type=Path, help='Target-keyed archive/job_url/optional job_id object')
    parser.add_argument('--config', type=Path, default=ROOT / 'agent2/config.json')
    parser.add_argument('--out', type=Path, help='New directory inside this submission folder')
    args = parser.parse_args()
    inputs = json.loads(args.inputs_json.read_text()) if args.inputs_json else {}
    for field, values in [('archive', args.archive), ('job_url', args.job_url), ('job_id', args.job_id)]:
        for target, value in keyed(values).items():
            inputs.setdefault(target, {})[field] = value
    require(set(inputs) == set(TARGETS), 'Supply exactly T1147, T1133 and T1183')
    config_path = args.config.resolve()
    require(config_path.is_relative_to(ROOT), 'Configuration must remain inside the repository')
    config = json.loads(config_path.read_text())
    entries = {entry['target']: entry for entry in config['targets']}
    require(set(entries) == set(TARGETS), 'Base configuration target set differs')
    plans = {}
    for target in TARGETS:
        item = inputs[target]
        archive_path = Path(item['archive']).resolve()
        require(archive_path.is_file() and zipfile.is_zipfile(archive_path), f'{target}: ZIP missing/invalid')
        url = urlparse(item['job_url'])
        require(url.scheme == 'https' and url.hostname == 'alphafoldserver.com' and url.path.startswith('/fold/'),
                f'{target}: use the actual AlphaFold Server job URL')
        job_id = item.get('job_id') or url.path.rstrip('/').split('/')[-1]
        require(job_id and job_id == url.path.rstrip('/').split('/')[-1], f'{target}: job ID/URL disagree')
        fasta = ROOT / 'targets' / target / (target + '.fasta')
        prepared = HERE / (target + '.fasta')
        require(fasta.read_bytes() == prepared.read_bytes(), f'{target}: prepared FASTA changed')
        plan = archive_plan(archive_path, target, fasta_sequence(fasta))
        plans[target] = {**plan, 'archive_source': str(archive_path), 'archive_sha256': sha(archive_path),
                         'job_url': item['job_url'], 'job_id': job_id, 'source_fasta': relative(fasta),
                         'source_fasta_sha256': sha(fasta), 'sequence_length': len(fasta_sequence(fasta))}
    output = (args.out or HERE / ('imported_' + datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ'))).resolve()
    require(output.is_relative_to(HERE) and output != HERE, 'Output must be a NEW child of the submission folder')
    output.mkdir(parents=True, exist_ok=False)
    # Every archive passes validation and primary selection before any output directory is created.
    records = {}
    for target in TARGETS:
        plan = plans[target]
        directory = output / target
        directory.mkdir()
        archive_copy = directory / Path(plan['archive_source']).name
        shutil.copyfile(plan['archive_source'], archive_copy)
        require(sha(archive_copy) == plan['archive_sha256'], f'{target}: preserved ZIP hash mismatch')
        extracted = directory / 'outputs'
        extracted.mkdir()
        with zipfile.ZipFile(archive_copy) as archive:
            archive.extractall(extracted)
        fasta_copy = directory / 'submitted.fasta'
        shutil.copyfile(ROOT / plan['source_fasta'], fasta_copy)
        plan['preserved_archive'] = relative(archive_copy)
        plan['output_files_sha256'] = {p.relative_to(extracted).as_posix(): sha(p) for p in sorted(extracted.rglob('*')) if p.is_file()}
        plan['selection_rule'] = 'Highest exported server ranking_score; lowest sample index for rounded ties; fixed before truth scoring'
        plan['status'] = 'VALIDATED_IMPORTED_NOT_REFERENCE_SCORED'
        entries[target]['af3_server'] = {
            'path': relative(extracted / plan['primary']['model_file']), 'submitted_fasta': relative(fasta_copy),
            'chain': plan['chain'], 'model_index': 0, 'mapping_mode': plan['mapping_mode'],
            'job_id': plan['job_id'], 'source_url': plan['job_url'],
            'import_provenance': relative(directory / 'provenance.json'),
            'ranking_score': plan['primary']['ranking_score'], 'sample_index': plan['primary']['sample_index'],
            'actual_settings': plan['actual_settings'], 'differences_from_prepared': plan['differences_from_prepared']
        }
        save_json(directory / 'provenance.json', plan)
        records[target] = plan
    save_json(output / 'config_with_server.json', config)
    save_json(output / 'import_provenance.json', {
        'created_utc': datetime.now(timezone.utc).isoformat(), 'status': 'VALIDATED_IMPORTED_NOT_REFERENCE_SCORED',
        'base_config': relative(config_path), 'base_config_sha256': sha(config_path),
        'import_helper_sha256': sha(Path(__file__)), 'documentation': SCHEMA_SOURCE,
        'selection_before_truth': True, 'targets': records
    })
    notice = ('# AlphaFold Server output notice\n\nAll downloaded ZIPs and all extracted files are preserved unchanged. '
              'Each target outputs/ folder contains the server terms_of_use.md; read and retain those terms when reusing outputs.\n\n'
              'The import adds provenance and configuration. Later C-alpha extraction, superposition, figures and metrics '
              'are derived analysis written separately; they do not replace the original server outputs.\n\n'
              'Parameter omissions/differences appear in provenance.json and config_with_server.json. An omitted '
              'maxTemplateDate is a server default not serialized in its request; its effective cutoff is not asserted.\n')
    (output / 'OUTPUT_TERMS_NOTICE.md').write_text(notice, encoding='utf-8')
    print(json.dumps({'output': str(output), 'config_with_server': str(output / 'config_with_server.json'),
                      'primary_samples': {t: records[t]['primary'] for t in TARGETS}}, indent=2))


if __name__ == '__main__':
    main()

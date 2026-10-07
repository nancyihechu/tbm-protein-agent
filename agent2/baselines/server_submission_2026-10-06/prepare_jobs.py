"""Prepare exact, locally validated AlphaFold Server inputs; submit nothing."""
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

DEST = Path(__file__).resolve().parent
ROOT = DEST.parents[2]
TARGETS = [('T1147', 103), ('T1133', 585), ('T1183', 200)]

def digest(data):
    return hashlib.sha256(data).hexdigest()

def write_new(path, data):
    with path.open('xb') as handle:
        handle.write(data)

jobs = []
entries = []
for target, expected_length in TARGETS:
    source = ROOT / 'targets' / target / (target + '.fasta')
    raw = source.read_bytes()
    lines = raw.decode('utf-8').splitlines()
    assert sum(line.startswith('>') for line in lines) == 1
    sequence = ''.join(line.strip() for line in lines if line.strip() and not line.startswith('>'))
    assert len(sequence) == expected_length
    assert not (set(sequence) - set('ACDEFGHIKLMNPQRSTVWY'))
    job = {
        'name': 'P2_' + target + '_AF3server_20261006',
        'modelSeeds': ['1'],
        'sequences': [{'proteinChain': {
            'sequence': sequence,
            'count': 1,
            'useStructureTemplate': True,
            'maxTemplateDate': '2021-09-30'
        }}],
        'dialect': 'alphafoldserver',
        'version': 1
    }
    assert set(job) == {'name', 'modelSeeds', 'sequences', 'dialect', 'version'}
    assert len(job['sequences']) == 1 and set(job['sequences'][0]) == {'proteinChain'}
    chain = job['sequences'][0]['proteinChain']
    assert set(chain) == {'sequence', 'count', 'useStructureTemplate', 'maxTemplateDate'}
    assert chain['sequence'] == sequence and chain['count'] == 1
    assert isinstance(chain['useStructureTemplate'], bool) and chain['useStructureTemplate']
    assert job['modelSeeds'] == ['1'] and 0 <= int(job['modelSeeds'][0]) < 2**32
    fasta_copy = DEST / (target + '.fasta')
    job_file = DEST / (target + '_job.json')
    json_bytes = (json.dumps([job], indent=2) + '\n').encode('utf-8')
    write_new(fasta_copy, raw)
    write_new(job_file, json_bytes)
    assert fasta_copy.read_bytes() == raw
    assert json.loads(job_file.read_text()) == [job]
    entries.append({
        'target': target,
        'job_name': job['name'],
        'status': 'PREPARED_NOT_SUBMITTED',
        'source_fasta': str(source.relative_to(ROOT)),
        'submitted_fasta': str(fasta_copy.relative_to(ROOT)),
        'source_fasta_sha256': digest(raw),
        'submitted_fasta_sha256': digest(fasta_copy.read_bytes()),
        'sequence_sha256': digest(sequence.encode('ascii')),
        'sequence_length': len(sequence),
        'job_json': str(job_file.relative_to(ROOT)),
        'job_json_sha256': digest(json_bytes),
        'server_job_id': None,
        'server_job_url': None,
        'submitted_utc': None
    })
    jobs.append(job)
batch_bytes = (json.dumps(jobs, indent=2) + '\n').encode('utf-8')
batch_file = DEST / 'jobs.json'
write_new(batch_file, batch_bytes)
assert len(json.loads(batch_file.read_text())) == 3
manifest = {
    'schema_version': 1,
    'created_utc': datetime.now(timezone.utc).isoformat(),
    'status': 'PREPARED_NOT_SUBMITTED',
    'schema_source': 'https://github.com/google-deepmind/alphafold/blob/main/server/README.md',
    'validation': 'Checked official documented job fields/types, uint32-string seed, list-of-jobs envelope, single protein chain, canonical sequence alphabet, declared lengths, byte-identical FASTA copies, and JSON round trips.',
    'batch_json': str(batch_file.relative_to(ROOT)),
    'batch_json_sha256': digest(batch_bytes),
    'settings': {'dialect': 'alphafoldserver', 'version': 1, 'modelSeeds': ['1'], 'count': 1, 'useStructureTemplate': True, 'maxTemplateDate': '2021-09-30', 'custom_msa': False, 'custom_templates': False},
    'primary_model_selection_rule': 'Choose the server-ranked primary/top model using server ranking evidence before opening or computing its experimental-reference accuracy scores. Preserve all downloaded samples and the ranking evidence. Do not select by truth scores.',
    'comparison_note': 'Seed 1 and template cutoff match Nancy\'s declared local AF3 setup. Server implementation/database differences may remain; equal seed/cutoff does not guarantee identical predictions.',
    'jobs': entries
}
write_new(DEST / 'manifest.json', (json.dumps(manifest, indent=2) + '\n').encode('utf-8'))
print(json.dumps({'batch_upload': str(batch_file), 'individual_uploads': [str(DEST / (t + '_job.json')) for t, _ in TARGETS], 'validated_lengths': {t: n for t, n in TARGETS}, 'status': manifest['status']}, indent=2))

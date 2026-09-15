import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import hpo.pipeline.trial as pipeline
from hpo.pipeline.trial import normalize_pairs
from hpo.pipeline.dssr import parse_pairs
from hpo.pipeline.scoring import score_pairs, aggregate

FIXTURE = json.loads(
    (Path(__file__).parent / 'fixtures' / 'dssr_example.json').read_text()
)
SEQ = 'AGUAGUAGUAGUAGA'


class Tests(unittest.TestCase):
    def test_dssr_json_parser(self):
        predicted = parse_pairs(FIXTURE, SEQ)
        self.assertEqual(predicted, [(0,14),(1,13),(2,12),(3,11),(4,10),(5,9)])
        score = score_pairs({(0,14),(1,13),(7,12)}, predicted)
        self.assertEqual((score['tp'], score['fp'], score['fn']), (2,4,1))
        self.assertAlmostEqual(score['loss'], 5/9)

    def test_empty_duplicate_invalid(self):
        self.assertEqual(score_pairs([],[])['loss'], 0)
        self.assertEqual(score_pairs([(0,1)],[])['loss'], 1)
        self.assertEqual(normalize_pairs([[0,1],[1,0]],2), {(0,1)})
        for pairs in ([[0,0]], [[-1,2]], [[0,4]], [[0.0,1]], [[True,1]]):
            with self.assertRaises(ValueError): normalize_pairs(pairs,3)
        with self.assertRaises(ValueError): aggregate([])

    def test_mapping_rejects_missing_and_mismatch(self):
        with self.assertRaises(ValueError): parse_pairs(FIXTURE, 'A' * len(SEQ))
        broken = dict(FIXTURE, nts=FIXTURE['nts'][:-1])
        with self.assertRaises(ValueError): parse_pairs(broken, SEQ)

    def test_parse_accounting(self):
        accounting = '\n'.join([
            '7_0|COMPLETED|0:0|',
            '7_1|FAILED+|1:0|',
            '',
            'malformed',
            '7_2||0:0|',
        ])
        self.assertEqual(
            pipeline.parse_accounting(accounting),
            {
                '7_0': ('COMPLETED', '0:0'),
                '7_1': ('FAILED', '1:0'),
            },
        )

    def test_array_accounting(self):
        cfg = dict(wait_timeout_seconds=10,poll_seconds=0)
        with patch.object(pipeline, 'command', side_effect=['', '7_0|COMPLETED|0:0|\n7_1|COMPLETED|0:0|']):
            pipeline.wait_job('7',2,cfg)
        with patch.object(pipeline, 'command', return_value='7_0|FAILED|1:0|'):
            with self.assertRaises(RuntimeError): pipeline.wait_job('7',1,cfg)

    def test_resume_and_failure_no_partial_loss(self):
        # Real orchestration + real DSSR parsing; external binaries substituted.
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            cfg = dict(fixed={}, dssr='fake', dssr_timeout_seconds=10)
            rows = [dict(id='one', sequence=SEQ,pairs=[[0,14],[1,13],[7,12]])]
            pipeline.write_json(directory/'trial_metadata.json',dict(fingerprint='abc',parameters={}))
            pipeline.write_json(directory/'af3_job_id.json',dict(job_id='7'))
            model = directory/'rna_00000/af3/model.cif'
            model.parent.mkdir(parents=True)
            model.write_text('mock model')
            pipeline.write_json(model.parent/'af3_success.json',dict(model=str(model)))
            with patch.object(pipeline, 'wait_job') as wait, patch.object(pipeline.dssr, 'run_dssr', return_value=parse_pairs(FIXTURE, SEQ)):
                loss = pipeline.run_pipeline(cfg,rows,directory,'abc')
                self.assertAlmostEqual(loss,5/9)
                self.assertAlmostEqual(pipeline.run_pipeline(cfg,rows,directory,'abc'),loss)
                self.assertEqual(wait.call_count,1)
            with self.assertRaises(ValueError): pipeline.run_pipeline(cfg,rows,directory,'abc',{'N':3})
            (directory/'trial_result.json').unlink()
            with patch.object(pipeline,'wait_job',side_effect=RuntimeError('failed')):
                with self.assertRaises(RuntimeError): pipeline.run_pipeline(cfg,rows,directory,'abc')
            self.assertFalse((directory/'trial_result.json').exists())

    def test_fresh_submission_and_ambiguous_guard(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            cfg = dict(fixed={}, shs_python=sys.executable,
                shs_seed=42, af3_seed=1, model_dir='~/af3/models',
                generator_timeout_seconds=5, module='bio/alphafold/3.0.1',
                gpu=dict(concurrency=2,partition='gpu-single',cpus=8,memory='20G',gres='gpu:A100:1',time='00:15:00'),
                dssr='fake',dssr_timeout_seconds=5)
            rows=[dict(id='one', sequence=SEQ,pairs=[[0,14]])]
            def generate(args, **kwargs):
                Path(args[-1]).write_text('{}')
            def submit(args):
                self.assertIn('--array=0-0%2',args)
                task=json.loads((directory/'af3_tasks.json').read_text())[0]
                model=Path(task['output'])/'rna_00000/rna_00000_model.cif'
                model.parent.mkdir(parents=True)
                model.write_text('mock model')
                pipeline.write_json(Path(task['output'])/'af3_success.json',dict(model=str(model)))
                return '17'
            with patch.object(pipeline.subprocess, 'run', side_effect=generate), patch.object(pipeline, 'command', side_effect=submit) as submit_mock, patch.object(pipeline, 'wait_job'), patch.object(pipeline.dssr, 'run_dssr', return_value={(0,14)}):
                self.assertEqual(pipeline.run_pipeline(cfg,rows,directory,'abc'),0)
                self.assertEqual(submit_mock.call_count,1)
            (directory/'af3_job_id.json').unlink()
            (directory/'trial_result.json').unlink()
            with patch.object(pipeline.subprocess,'run',side_effect=generate), patch.object(pipeline,'command') as submit_mock:
                with self.assertRaisesRegex(RuntimeError,'submission-intent'):
                    pipeline.run_pipeline(cfg,rows,directory,'abc')
                submit_mock.assert_not_called()

    def test_folder_inputs(self):
        with tempfile.TemporaryDirectory() as tmp:
            folder = Path(tmp)
            for i in range(2):
                (folder/f'{i}.csv').write_text(f'id,sequence,pairs\nrna{i},ACGU,"[[0,3]]"\n')
            self.assertEqual(len(pipeline.create_dataset(folder)),2)
            (folder/'duplicate.csv').write_text('id,sequence,pairs\nrna0,ACGU,[]\n')
            with self.assertRaises(ValueError): pipeline.create_dataset(folder)


if __name__ == '__main__':
    unittest.main()

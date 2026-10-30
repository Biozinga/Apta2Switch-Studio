"""Guide/chemistry/indexing checks; native preparation has a separate runtime."""
import json
import tempfile
import time
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch
import numpy as np
from aptaswitch_core.oxdna_guided import guide_pairs, seed_strand, prepare_guided_system, _backbone, _configuration
from aptaswitch_core.scientific_process import CalculationCancelled
from aptaswitch_core.oxdna_guided_worker import Stage, validation_from_bonds, separate_monomers


class GuideTests(unittest.TestCase):
    def test_preparation_checks_actual_strand_separation_not_only_target_register(self):
        from tests.test_oxdna_ffs_worker import toehold_request
        request = toehold_request()
        request.update(switch_sequence='AC', trigger_sequence='AT', switch_material='rna', trigger_material='dna',
                       switch_guide='..', trigger_guide='..',
                       guides={'switch': {'pairs': []}, 'trigger': {'pairs': []}})
        particles = [SimpleNamespace(pos=x) for x in (0., 8., 12., 9.)]
        stage = Stage.__new__(Stage)
        stage.hb = SimpleNamespace(hb_list=lambda: [])
        from tests.test_oxdna_kinetics_worker import frame
        stage.frame = frame
        stage.info = SimpleNamespace(interaction=SimpleNamespace(begin_energy_computation=lambda: None),
            particles=lambda: particles,
            box=SimpleNamespace(sqr_min_image_distance=lambda p, q: (q.pos-p.pos)**2))
        physical = stage.observe(request)
        self.assertEqual(physical['basin_distance'], 1.)
        self.assertFalse(physical['ffs_eligible'])
        for key in ('distance_scope', 'distance_pairs', 'separation_pairs'):
            request['plan'].pop(key)
        historical = stage.observe(request)
        self.assertEqual(historical['basin_distance'], 4.)
        self.assertTrue(historical['ffs_eligible'])

    def test_wobble_is_rna_only(self):
        self.assertEqual(guide_pairs('GAAAU','(...)','rna'),[(0,4)])
        with self.assertRaises(ValueError):
            guide_pairs('GAAAT','(...)','dna')

    def test_reject_incomplete_wrong_or_interstrand_guides(self):
        for guide in ('((..)', '....+', '[[.]]', '(...)('):
            with self.subTest(guide=guide), self.assertRaises(ValueError):
                guide_pairs('GAAAC',guide,'rna')
        with self.assertRaises(ValueError):
            guide_pairs('AAAAT','(...)','rna')
        with self.assertRaises(ValueError):
            guide_pairs('GAC','(.)','dna')

    def test_hairpin_seed_preserves_pairs_and_orthonormal_axes(self):
        for material,sequence in [('rna','GGGGAAAACCCC'),('dna','GGGGAAAACCCC')]:
            with self.subTest(material=material):
                rows=seed_strand(sequence,'((((....))))',material)
                self.assertEqual(rows.shape,(12,15))
                self.assertTrue(np.isfinite(rows).all())
                np.testing.assert_allclose(np.linalg.norm(rows[:,3:6],axis=1),1)
                np.testing.assert_allclose(np.linalg.norm(rows[:,6:9],axis=1),1)
                np.testing.assert_allclose(np.sum(rows[:,3:6]*rows[:,6:9],axis=1),0,atol=1e-14)
                for i,j in guide_pairs(sequence,'((((....))))',material):
                    self.assertAlmostEqual(np.linalg.norm(rows[i,:3]-rows[j,:3]),1.2)
                    np.testing.assert_allclose(rows[i,3:9],-rows[j,3:9])
                # 5-prime to 3-prime is against a3, as the native topology expects.
                self.assertLess(np.dot(rows[1,:3]-rows[0,:3],rows[0,6:9]),0)

    def test_linear_guide_has_valid_backbone_and_no_paired_seed(self):
        for material,sequence in [('rna','GAUCAUCG'),('dna','GATCATCG')]:
            rows=seed_strand(sequence,'.'*len(sequence),material)
            backs=_backbone(rows,material)
            np.testing.assert_allclose(np.linalg.norm(np.diff(backs,axis=0),axis=1),.76)

    def test_coaxial_bulge_multistem_is_finite(self):
        guide='(((...(((....)))...)))..(((....)))'
        seq=['A']*len(guide)
        stack=[]
        for i,c in enumerate(guide):
            if c=='(':stack.append(i)
            elif c==')':j=stack.pop();seq[j]='G';seq[i]='C'
        rows=seed_strand(''.join(seq),guide,'rna')
        self.assertTrue(np.isfinite(rows).all())

    def test_multiloop_is_finite(self):
        guide='((..((....))..((....))..))'
        seq=['A']*len(guide);stack=[]
        for i,c in enumerate(guide):
            if c=='(':stack.append(i)
            elif c==')':j=stack.pop();seq[j]='G';seq[i]='C'
        self.assertTrue(np.isfinite(seed_strand(''.join(seq),guide,'rna')).all())

    def test_cancel_before_any_native_process_retains_diagnostics(self):
        with tempfile.TemporaryDirectory() as folder:
            with self.assertRaisesRegex(RuntimeError,'cancelled'):
                prepare_guided_system({},Path(folder),cancel_check=lambda:True)
            report=json.loads((Path(folder)/'preparation.json').read_text())
            self.assertEqual(report['status'],'cancelled')
            self.assertNotIn('configuration_file',report)

    def test_diagnostic_report_preserves_exact_recorded_guide_provenance(self):
        for model, source, provenance in (
            ('experimental-Sugimoto95-Vienna04-loopmean', 'native_hybrid_mfe_isolated_strands',
             {'engine': 'native-hybrid', 'parameter_digest': 'fixture-only-digest'}),
            ('external-archive', 'archived_mfe_isolated_strands', {}),
        ):
            request = {'guide_model': model, 'guide_source': source, 'guide_thermodynamics': provenance}
            with self.subTest(model=model), tempfile.TemporaryDirectory() as folder:
                with self.assertRaisesRegex(RuntimeError, 'cancelled'):
                    prepare_guided_system(request, Path(folder), cancel_check=lambda: True)
                report = json.loads((Path(folder) / 'preparation.json').read_text())
                self.assertEqual(report['method'], 'sequence_and_recorded_2d')
                self.assertEqual(report['guide_model'], model)
                self.assertEqual(report['guide_source'], source)
                self.assertEqual(report['guide_thermodynamics'], provenance)

    def test_expired_budget_never_starts_native_process(self):
        with tempfile.TemporaryDirectory() as folder:
            with patch('aptaswitch_core.oxdna_guided.run_stage') as native:
                with self.assertRaises(TimeoutError):
                    prepare_guided_system({},Path(folder),deadline=time.monotonic()-1)
                native.assert_not_called()
            report=json.loads((Path(folder)/'preparation.json').read_text())
            self.assertEqual(report['status'],'failed')
            self.assertNotIn('configuration_file',report)

    def test_cancel_wins_over_a_completed_worker_race(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder)
            parameter=root/'params.txt';parameter.write_text('test')
            runtime={'python':'unused',**{'seq_dep_file_'+chem:str(parameter) for chem in ('DNA','RNA','NA')}}
            request={'switch_sequence':'GGGGAAAACCCC','switch_material':'rna','switch_guide':'((((....))))',
                     'trigger_sequence':'AAAT','trigger_material':'dna','trigger_guide':'....',
                     'temperature_c':25,'salt_m':.5,'seed':1,'runtime':runtime}
            output=root/'job'
            def completed_then_cancelled(*args, **kwargs):
                (output/'guided_result.json').write_text(json.dumps({'status':'completed','configuration_file':'candidate.dat'}))
                raise CalculationCancelled('cancelled just after native completion')
            with patch('aptaswitch_core.oxdna_guided.run_stage',side_effect=completed_then_cancelled):
                with self.assertRaisesRegex(RuntimeError,'cancelled'):
                    prepare_guided_system(request,output)
            report=json.loads((output/'preparation.json').read_text())
            self.assertEqual(report['status'],'cancelled')
            self.assertFalse(report['ffs_eligible'])
            self.assertNotIn('configuration_file',report)

    def test_missing_guide_is_inspectable_but_not_ffs_eligible(self):
        request={'switch_sequence':'GGGGAAAACCCC','guides':{'switch':{'pairs':[[0,11],[1,10],[2,9],[3,8]]},'trigger':{'pairs':[]}},
                 'plan':{'stem_pairs':[[1,10],[2,9],[3,8]],'basin_min_stem':3,'basin_distance':4}}
        observation=validation_from_bonds(request,{(1,10),(2,9),(3,8)},6.)
        self.assertTrue(observation['off_basin_eligible'])
        self.assertFalse(observation['guide_accepted'])
        self.assertFalse(observation['ffs_eligible'])
        self.assertEqual(observation['guides']['switch']['missing_pairs'],[[0,11]])

    def test_full_guide_with_extra_pair_is_explicitly_not_exact(self):
        request={'switch_sequence':'GGGGAAAACCCC','guides':{'switch':{'pairs':[[0,11],[1,10],[2,9]]},'trigger':{'pairs':[]}},
                 'plan':{'stem_pairs':[[1,10],[2,9]],'basin_min_stem':2,'basin_distance':4}}
        observation=validation_from_bonds(request,{(0,11),(1,10),(2,9),(3,8)},6.)
        self.assertFalse(observation['guide_accepted'])
        self.assertFalse(observation['ffs_eligible'])
        self.assertFalse(observation['guide_exact_match'])
        self.assertEqual(observation['unexpected_pairs'],[[3,8]])
        self.assertFalse(validation_from_bonds(request,{(0,11),(1,10),(2,9)},3.)['ffs_eligible'])

    def test_final_separation_preserves_internal_geometry_and_has_periodic_gap(self):
        first=seed_strand('GGGGAAAACCCC','((((....))))','rna')
        second=seed_strand('AAAT','....','dna')
        original=np.vstack([first,second])
        original[:,:3]+=20
        original[len(first):,0]+=1  # Deliberately close after diffusion.
        with tempfile.TemporaryDirectory() as folder:
            source=Path(folder)/'source.dat';target=Path(folder)/'placed.dat'
            _configuration(source,original,40)
            placement=separate_monomers(source,target,len(first),6.)
            placed=np.array([[float(value) for value in line.split()] for line in target.read_text().splitlines()[3:]])
        np.testing.assert_allclose(placed[:,3:],original[:,3:],atol=1e-14)
        for a,b in ((0,len(first)),(len(first),len(original))):
            np.testing.assert_allclose(np.diff(placed[a:b,:3],axis=0),np.diff(original[a:b,:3],axis=0),atol=1e-12)
        self.assertAlmostEqual(placed[len(first):,0].min()-placed[:len(first),0].max(),6.)
        self.assertGreaterEqual(placement['minimum_periodic_x_gap'],6.)

    def test_nonempty_output_is_not_overwritten(self):
        with tempfile.TemporaryDirectory() as folder:
            marker=Path(folder)/'keep.dat';marker.write_text('original')
            with self.assertRaisesRegex(ValueError,'empty'):
                prepare_guided_system({},Path(folder))
            self.assertEqual(marker.read_text(),'original')

if __name__=='__main__':unittest.main()

"""Shared native oxNA frame, validation and checkpoint primitives."""
from pathlib import Path
from copy import deepcopy
import os
import tempfile
from types import SimpleNamespace
import unittest

from aptaswitch_core.oxdna_ffs_worker import FrameBuffer, OxpyAdapter, _same_checkpoint, distance_pair_groups, validate_request


def request(**sampling):
    return {"prepared": {"topology_file": "unused.top", "configuration_file": "initial.dat"},
            "sampling": {"temperature_c": 37, "salt_m": .15, "frame_interval": 1,
                         "max_movie_frames": 20, **sampling},
            "plan": {"switch_length":2, "nucleotide_count":4,
                     "binding_pairs":[[0,3],[1,2]], "stem_pairs":[[0,1]], "opening_pairs":[[0,2]],
                     "interfaces":[{"name":"opening", "min_binding":1,"max_stem":0,"max_opening":0}]}}


def toehold_request():
    data = request()
    data['plan'].update(distance_scope='trigger_to_toehold', toehold_pairs=[[0, 3]],
        distance_pairs=[[0, 2], [0, 3]], separation_pairs=[[0, 2], [0, 3], [1, 2], [1, 3]])
    return data


def terminal_toehold_request():
    data = toehold_request()
    data['plan'].update(distance_scope='trigger_3prime_to_toehold_start', distance_pairs=[[0, 3]])
    return data


class NativePrimitiveTests(unittest.TestCase):
    def test_terminal_distance_rejects_wrong_end_and_other_toehold_bases(self):
        valid = validate_request(terminal_toehold_request())
        self.assertEqual(distance_pair_groups(valid['plan']), {
            'distance': ((0, 3),), 'separation': ((0, 2), (0, 3), (1, 2), (1, 3)),
        })
        for mutation in (
            {'distance_pairs': [[0, 2]]},
            {'distance_pairs': [[0, 2], [0, 3]]},
            {'distance_pairs': [[1, 3]]},
            {'toehold_pairs': [[1, 2]], 'distance_pairs': [[1, 2]]},
            {'separation_pairs': [[0, 3]]},
        ):
            with self.subTest(mutation=mutation), self.assertRaises(ValueError):
                data = terminal_toehold_request()
                data['plan'].update(mutation)
                validate_request(data)

    def test_distance_scope_validates_region_register_and_all_contact_pairs(self):
        valid = validate_request(toehold_request())
        self.assertEqual(distance_pair_groups(valid['plan']), {
            'distance': ((0, 2), (0, 3)), 'separation': ((0, 2), (0, 3), (1, 2), (1, 3)),
        })
        for mutation in (
            {'toehold_pairs': []}, {'toehold_pairs': [[0, 2]]},
            {'distance_scope': 'centres_of_strands'}, {'distance_pairs': [[0, 3]]},
            {'distance_pairs': [[0, 2], [0, 4]]}, {'separation_pairs': [[0, 3], [1, 2]]},
        ):
            with self.subTest(mutation=mutation), self.assertRaises(ValueError):
                data = toehold_request()
                data['plan'].update(mutation)
                validate_request(data)
        data = toehold_request()
        data.update(switch_sequence='AC', trigger_sequence='GT')
        validate_request(data)
        data['trigger_sequence'] = 'TT'
        with self.assertRaisesRegex(ValueError, 'sequence register'):
            validate_request(data)

    def test_legacy_distance_scope_and_identity_are_not_silently_redefined(self):
        old = validate_request(request())
        self.assertEqual(distance_pair_groups(old['plan']), {'distance': ((0, 3), (1, 2))})
        for field in ('distance_scope', 'distance_pairs', 'separation_pairs'):
            self.assertNotIn(field, old['plan'])
        data = deepcopy(old)
        data['plan']['distance_pairs'] = [[0, 3]]
        with self.assertRaises(ValueError):
            validate_request(data)

    def test_observed_distances_include_off_register_contacts_and_periodic_images(self):
        def observe(positions, plan, box_size):
            particles = [SimpleNamespace(pos=[x, 0., 0.]) for x in positions]
            def sqr_distance(p, q):
                delta = q.pos[0] - p.pos[0]
                return (delta - box_size * round(delta / box_size)) ** 2
            adapter = OxpyAdapter.__new__(OxpyAdapter)
            adapter.box_policy = None
            adapter.info = SimpleNamespace(interaction=SimpleNamespace(begin_energy_computation=lambda: None),
                box=SimpleNamespace(sqr_min_image_distance=sqr_distance))
            adapter.hb_threshold = -.1
            adapter.hb_observable = SimpleNamespace(hb_list=lambda: [])
            adapter.switch_length = 2
            adapter.pair_indices = {key: set(tuple(p) for p in plan[key]) for key in ('binding_pairs', 'stem_pairs', 'opening_pairs')}
            adapter.distance_pairs = {name: [(particles[i], particles[j]) for i, j in pairs]
                                      for name, pairs in distance_pair_groups(plan).items()}
            adapter.manager = SimpleNamespace(current_step=100)
            return adapter.observe()
        new = observe([0., 5., 19., 12.], toehold_request()['plan'], 20.)
        old = observe([0., 5., 19., 12.], request()['plan'], 20.)
        self.assertEqual((new.distance, new.separation, new.binding), (1., 1., 0))
        self.assertEqual((old.distance, old.separation), (6., None))
        body_contact = observe([0., 8., 9., 13.], toehold_request()['plan'], 40.)
        self.assertEqual((body_contact.distance, body_contact.separation), (9., 1.))
        terminal = terminal_toehold_request()['plan']
        # Bringing only 5′ close cannot cross even the first distance interface.
        wrong_end = observe([0., 5., 1., 12.], terminal, 40.)
        self.assertEqual((wrong_end.distance, wrong_end.separation), (12., 1.))
        right_end = observe([0., 5., 12., 1.], terminal, 40.)
        self.assertEqual((right_end.distance, right_end.separation), (1., 1.))
        wrapped = observe([0., 5., 12., 19.], terminal, 20.)
        self.assertEqual((wrapped.distance, wrapped.separation), (1., 1.))

    def test_terminal_coordinate_changes_checkpoint_identity(self):
        from aptaswitch_core.oxdna_ffs_worker import request_identity
        self.assertNotEqual(request_identity(toehold_request()),
                            request_identity(terminal_toehold_request()))

    def test_checkpoint_comparison_allows_native_roundoff_but_not_changed_phase_space(self):
        first = "t = 100\nb = 20 20 20\nE = 1 2 3\n5 5 5 1 0 0 0 0 1 .1 .2 .3 .4 .5 .6\n"
        rounded = first.replace("E = 1 2 3", "E = 1.0000000000001 2 3").replace(".1 .2", ".10000000000001 .2")
        self.assertTrue(_same_checkpoint(first, rounded))
        for altered in (first.replace("t = 100", "t = 101"),
                        first.replace("b = 20", "b = 30"),
                        first.replace("5 5 5", "5.1 5 5"),
                        first.replace(".1 .2", ".2 .2"),
                        first.replace(".4 .5 .6", ".4 .5 .7"),
                        first.replace(".1 .2", "nan .2")):
            with self.subTest(altered=altered):
                self.assertFalse(_same_checkpoint(first, altered))


    def test_frame_buffer_has_bounded_storage_and_keeps_actual_endpoints(self):
        buffer = FrameBuffer(1, 6)
        for i in range(100):
            buffer.add(i, lambda i=i: f"actual {i}")
        buffer.add(100, lambda: "actual 100", force=True)
        self.assertLessEqual(len(buffer.frames), 6)
        self.assertEqual(buffer.frames[0], (0, "actual 0"))
        self.assertEqual(buffer.frames[-1], (100, "actual 100"))


    def test_validation_normalizes_pairs_and_rejects_invalid_md_and_censoring(self):
        data=request()
        data['plan']['binding_pairs'].append([3,0])
        self.assertEqual(len(validate_request(data)['plan']['binding_pairs']),2)
        for sampling in ({'salt_m':0},{'dt':-1},{'budget_seconds':1},
                         {'max_steps_per_trial':100},{'diversity':{'version':2}},
                         {'temperature_c':float('nan')}):
            with self.subTest(sampling=sampling),self.assertRaises(ValueError):
                validate_request(request(**sampling))
        data=request()
        data['plan']['binding_pairs'].append([0,100])
        with self.assertRaises(ValueError):
            validate_request(data)

    def test_request_validation_does_not_inject_unused_sampling_controls(self):
        sampling=validate_request(request())['sampling']
        for key in ('budget_seconds','max_steps_per_trial','flux_crossings','trials_per_interface','diversity'):
            self.assertNotIn(key,sampling)


@unittest.skipUnless(os.environ.get("APTA_OXDNA_SOURCE"), "requires external oxpy and NA parameter files")
class NativeCheckpointSmokeTests(unittest.TestCase):
    def test_mixed_NA_checkpoint_preserves_phase_space_and_step_counter(self):
        from aptaswitch_core.oxdna_box import estimate_ffs_box
        source = Path(os.environ["APTA_OXDNA_SOURCE"])
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            top, conf = root / "pair.top", root / "pair.dat"
            top.write_text("2 2 5->3\nG type=DNA\nC type=RNA\n")
            conf.write_text("t = 0\nb = 40 40 40\nE = 0 0 0\n"
                            "5 5 5 1 0 0 0 0 1 0 0 0 0 0 0\n"
                            "15 5 5 -1 0 0 0 0 -1 0 0 0 0 0 0\n")
            settings = request()
            settings["prepared"] = {"topology_file": str(top), "configuration_file": str(conf)}
            settings["plan"]["binding_pairs"] = [[0, 1]]
            settings["plan"]["switch_length"] = 1
            settings["plan"]["nucleotide_count"] = 2
            settings["plan"]["stem_pairs"] = []
            settings["plan"]["opening_pairs"] = []
            settings['box_policy'] = estimate_ffs_box(1, 1, 'dna', 'rna', 37, .15)
            settings['box_size'] = 40.
            settings["runtime"] = {
                "seq_dep_file_DNA": str(source / "oxDNA2_sequence_dependent_parameters.txt"),
                "seq_dep_file_RNA": str(source / "rna_sequence_dependent_parameters.txt"),
                "seq_dep_file_NA": str(source / "NA_sequence_dependent_parameters.txt"),
            }
            settings = validate_request(settings)
            adapter = OxpyAdapter(settings, conf, 31, True)
            try:
                self.assertEqual(adapter.observe().binding, 0)
                self.assertEqual(adapter.step().step, 1)
                self.assertAlmostEqual(adapter.periodic_clearance['minimum_clearance'], 40.)
                before = adapter.frame()
                checkpoint = root / "checkpoint.dat"
                checkpoint.write_text(before)
            finally:
                adapter.close()
            adapter = OxpyAdapter(settings, checkpoint, 32, False)
            try:
                after = adapter.frame()
                self.assertEqual(adapter.observe().step, 1)
                for original, restored in zip(before.splitlines()[3:], after.splitlines()[3:]):
                    values = list(map(float, original.split()))
                    self.assertEqual(len(values), 15)  # positions, orientation, velocity, angular momentum
                    for expected, actual in zip(values, map(float, restored.split())):
                        self.assertAlmostEqual(expected, actual, places=10)
                self.assertEqual(adapter.step().step, 2)
                self.assertGreater(adapter.periodic_clearance['minimum_clearance'],
                                   adapter.periodic_clearance['required_clearance'])
            finally:
                adapter.close()



if __name__ == "__main__":
    unittest.main()

"""Association help preserves recorded physical inputs in either language."""
from dataclasses import replace
import unittest

from aptaswitch_core.localization import language_context, tr
from aptaswitch_core.models import ScoringWeights, ThermoConditions
from aptaswitch_studio.scoring_help import (
    ASSOCIATION_WITH_APTAMER_LABEL, association_with_aptamer_help,
)


class AssociationHelpTests(unittest.TestCase):
    def test_new_defaults_disable_aptamer_association_and_increase_linearity(self):
        weights = ScoringWeights()
        self.assertEqual(weights.leak, 0.)
        self.assertEqual(weights.trigger_linearity, 3.)

    def test_physical_definition_and_recorded_weight_in_both_languages(self):
        conditions = replace(ThermoConditions(), trigger_concentration_m=2e-6,
                             switch_concentration_m=7e-6, validation_max_complex_size=5)
        for language in ('fr', 'en'):
            with self.subTest(language=language), language_context(language):
                text = association_with_aptamer_help(conditions, weight=.2)
                self.assertIn('2 µM', text)
                self.assertIn('7 µM', text)
                self.assertIn('0.285714:0.285714:1', text)
                self.assertIn('5 brins' if language == 'fr' else '5 strands', text)
                self.assertIn('exactement 1 trigger et 1 switch' if language == 'fr'
                              else 'exactly 1 trigger and 1 switch', text)
                self.assertIn('Poids enregistré : 0.2' if language == 'fr'
                              else 'Recorded weight: 0.2', text)
                self.assertIn('pas une fuite de traduction' if language == 'fr'
                              else 'neither OFF translation leakage', text)
                self.assertNotIn('Poids 0 :' if language == 'fr' else 'Weight 0:', text)
                label = tr(ASSOCIATION_WITH_APTAMER_LABEL)
                self.assertIn('aptamère' if language == 'fr' else 'aptamer', label)

    def test_unknown_historical_weight_is_not_replaced_by_the_current_default(self):
        with language_context('fr'):
            text = association_with_aptamer_help(ThermoConditions())
            self.assertNotIn('Poids 0', text)
            self.assertNotIn('Poids enregistré', text)
            self.assertIn('n’intervient pas dans le score',
                          association_with_aptamer_help(ThermoConditions(), weight=0.))

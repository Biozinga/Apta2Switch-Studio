"""Full reporters and short design prefixes must remain distinct and portable."""

import hashlib
import unittest
from dataclasses import asdict
from unittest.mock import patch

from aptaswitch_core.models import REPORTERS, ReporterContext
from aptaswitch_core import reporter_library as library
from aptaswitch_core.reporter_sequences import (
    SFGFP_CDS_DNA, SFGFP_CDS_SHA256, MNEONGREEN_CDS_DNA, MNEONGREEN_CDS_SHA256,
    MCHERRY_CDS_DNA, MCHERRY_CDS_SHA256, MCHERRY_PROTEIN_SHA256,
    MRFP1_CDS_DNA, MRFP1_CDS_SHA256, MRFP1_PROTEIN_SHA256,
    LACZ_CDS_DNA, LACZ_CDS_SHA256, LACZ_CDS_RANGE, LACZ_CDS_STRAND, LACZ_PROTEIN_SHA256,
)


class CompleteReporterTests(unittest.TestCase):
    def test_builtin_full_sequences_match_versioned_genbank_cds_and_existing_prefixes(self):
        cases = (
            ("sfgfp", SFGFP_CDS_DNA, SFGFP_CDS_SHA256, 717, "HQ873313.1", "CGUAAAGGCGAGGAGCUGUUC"),
            ("mng", MNEONGREEN_CDS_DNA, MNEONGREEN_CDS_SHA256, 711, "KC295282.1", "GUGAGCAAGGGCGAGGAGGAU"),
            ("mcherry", MCHERRY_CDS_DNA, MCHERRY_CDS_SHA256, 711, "AY678264.1", "GUGAGCAAGGGCGAGGAGGAU"),
            ("mrfp1", MRFP1_CDS_DNA, MRFP1_CDS_SHA256, 678, "AF506027.1", "GCCUCCUCCGAGGACGUC"),
            ("lacz", LACZ_CDS_DNA, LACZ_CDS_SHA256, 3075, "NC_000913.3", "ACCAUGAUUACGGAUUCA"),
        )
        for key, dna, digest, length, accession, historical in cases:
            with self.subTest(key=key):
                self.assertEqual(len(dna), length)
                self.assertEqual(hashlib.sha256(dna.encode()).hexdigest(), digest)
                self.assertEqual(dna[:3], "ATG")
                self.assertIn(dna[-3:], {"TAA", "TGA", "TAG"})
                reporter = REPORTERS[key]
                reporter.validate()
                self.assertEqual(reporter.sequence_after_start_rna, dna[3:33].replace("T", "U"))
                self.assertEqual(reporter.full_sequence_after_start_rna, dna[3:].replace("T", "U"))
                self.assertTrue(reporter.sequence_after_start_rna.startswith(historical))
                self.assertEqual(reporter.source_accession, accession)
                self.assertEqual(reporter.source_url, "https://www.ncbi.nlm.nih.gov/nuccore/" + accession)

    def test_new_reporters_translate_exactly_to_the_annotated_source_proteins(self):
        # Standard genetic code; these CDS annotations use translation table 11,
        # whose internal codons and ATG initiation have the same translation.
        codons = (a+b+c for a in "TCAG" for b in "TCAG" for c in "TCAG")
        code = dict(zip(codons, "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"))
        for dna, digest, length in ((MCHERRY_CDS_DNA, MCHERRY_PROTEIN_SHA256, 236),
                                    (MRFP1_CDS_DNA, MRFP1_PROTEIN_SHA256, 225),
                                    (LACZ_CDS_DNA, LACZ_PROTEIN_SHA256, 1024)):
            protein = "".join(code[dna[i:i+3]] for i in range(0, len(dna), 3))
            self.assertEqual(protein[0], "M")
            self.assertEqual(protein[-1], "*")
            self.assertNotIn("*", protein[:-1])
            self.assertEqual(len(protein)-1, length)
            self.assertEqual(hashlib.sha256(protein[:-1].encode()).hexdigest(), digest)
        self.assertEqual(LACZ_CDS_RANGE, (363231, 366305))
        self.assertEqual(LACZ_CDS_STRAND, -1)

    def test_full_custom_sequence_removes_only_initial_aug_and_preserves_terminal_stop(self):
        full = "atg" + "acg" * 10 + "atg" + "cct" * 10 + "taa"
        reporter = library.reporter_from_full_sequence("test", "Rapporteur", full)
        self.assertEqual(reporter.sequence_after_start_rna, "ACG" * 10)
        self.assertEqual(reporter.full_sequence_after_start_rna, full[3:].upper().replace("T", "U"))
        self.assertIn("AUG", reporter.full_sequence_after_start_rna)
        self.assertTrue(reporter.full_sequence_after_start_rna.endswith("UAA"))
        self.assertEqual(ReporterContext(**asdict(reporter)), reporter)

    def test_no_terminal_stop_is_accepted_for_explicit_full_fusion_sequence(self):
        reporter = library.reporter_from_full_sequence("tag", "Fusion", "ATG" + "ACG" * 12)
        self.assertEqual(reporter.full_sequence_after_start_rna, "ACG" * 12)

    def test_invalid_frames_internal_stops_and_missing_start_are_rejected(self):
        for sequence in ("ACG" * 12, "ATG" + "ACG" * 9, "ATG" + "ACG" * 10 + "A", "ATG" + "ACG" * 10 + "TAA" + "ACG"):
            with self.subTest(sequence=sequence), self.assertRaises(ValueError):
                library.reporter_from_full_sequence("bad", "Bad", sequence)


    def test_full_sequence_roundtrip_and_prefix_consistency(self):
        settings = {}

        def save(new):
            settings.clear()
            settings.update(new)

        with patch.object(library, "load_settings", side_effect=lambda: dict(settings)), patch.object(library, "save_settings", side_effect=save):
            original = library.add_custom_reporter("Full", "ATG" + "ACG" * 15 + "TAA")
            self.assertEqual(library.load_custom_reporters()[original.key], original)
            self.assertEqual(settings["reporters"][0]["full_sequence_after_start_rna"], original.full_sequence_after_start_rna)
            settings["reporters"][0]["full_sequence_after_start_rna"] = "GGG" * 15
            self.assertNotIn(original.key, library.load_custom_reporters())

    def test_custom_reporter_named_like_new_builtin_is_saved_separately(self):
        settings = {}
        def save(new):
            settings.clear()
            settings.update(new)
        with patch.object(library, "load_settings", side_effect=lambda: dict(settings)), patch.object(library, "save_settings", side_effect=save):
            custom = library.add_custom_reporter("mCherry", "ATG" + "GCT" * 15 + "TAG")
            self.assertEqual(custom.key, "mcherry-2")
            reporters = library.all_reporters()
            self.assertEqual(len(reporters), 6)
            self.assertEqual(reporters["mcherry"], REPORTERS["mcherry"])
            self.assertEqual(reporters[custom.key].full_sequence_after_start_rna, "GCU" * 15 + "UAG")
            library.remove_custom_reporter(custom.key)
            self.assertEqual(library.all_reporters(), REPORTERS)


if __name__ == "__main__":
    unittest.main()

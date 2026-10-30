import importlib.util
from pathlib import Path
import unittest
ROOT = Path(__file__).resolve().parents[1]

class RuntimeNoticeTests(unittest.TestCase):
    def test_binary_notice_retains_complete_biopython_license(self):
        spec = importlib.util.spec_from_file_location("prepare_thermodynamics_runtime", ROOT / "scripts/prepare_thermodynamics_runtime.py")
        thermodynamics = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(thermodynamics)
        notice = thermodynamics.third_party_notices()
        self.assertIn("The Biopython Contributors", notice)
        self.assertIn("Redistributions in binary form must reproduce", notice)
        self.assertIn("POSSIBILITY OF SUCH DAMAGE.", notice)


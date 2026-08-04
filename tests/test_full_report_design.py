import base64
import re
import unittest
from pathlib import Path


REPORT_HTML = (
    Path(__file__).resolve().parents[1]
    / "design"
    / "report"
    / "full-report"
    / "index.html"
)


class FullReportDesignTest(unittest.TestCase):
    def test_html_is_self_contained(self) -> None:
        html = REPORT_HTML.read_text(encoding="utf-8")

        self.assertNotIn("../../../assets/", html)
        embedded_svgs = re.findall(
            r"data:image/svg\+xml;base64,([A-Za-z0-9+/=]+)", html
        )
        self.assertEqual(len(embedded_svgs), 3)

        for encoded_svg in embedded_svgs:
            svg = base64.b64decode(encoded_svg, validate=True)
            self.assertIn(b"<svg", svg[:1000])


if __name__ == "__main__":
    unittest.main()

"""Create application icons from the central app logo in the assets package."""

from __future__ import annotations

import base64
from io import BytesIO
from pathlib import Path

from PIL import Image
import resvg_py


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "src" / "aptaswitch_studio" / "assets" / "logo.png"
APP_ICON = SOURCE.with_name("app_icon.png")
WINDOWS_ICON = ROOT / "packaging" / "icons" / "AppIcon.ico"
LINUX_ICON = ROOT / "packaging" / "icons" / "aptaswitch-studio.png"
MACOS_ICON = ROOT / "packaging" / "icons" / "AppIcon.icns"
PNG_ICON = ROOT / "packaging" / "icons" / "AppIcon.png"


def icon_svg(logo_png: bytes) -> str:
    """Put the unmodified transparent logo on a softly lit macOS-style tile."""

    logo = base64.b64encode(logo_png).decode("ascii")
    # Keep 64 px around the tile so its soft shadow remains inside the canvas.
    # The logo has its own inset and cannot reach the rounded corners.
    return f'''<svg xmlns="http://www.w3.org/2000/svg"
        xmlns:xlink="http://www.w3.org/1999/xlink"
        width="1024" height="1024" viewBox="0 0 1024 1024">
      <defs>
        <path id="tile" d="M300 64 H724 C880 64 960 144 960 300
          V724 C960 880 880 960 724 960 H300
          C144 960 64 880 64 724 V300 C64 144 144 64 300 64 Z"/>
        <linearGradient id="surface" x1="0.12" y1="0" x2="0.88" y2="1">
          <stop offset="0" stop-color="#fffefa"/>
          <stop offset="0.42" stop-color="#f9fcff"/>
          <stop offset="1" stop-color="#e2ecf6"/>
        </linearGradient>
        <radialGradient id="light" cx="0.18" cy="0.02" r="0.9">
          <stop offset="0" stop-color="#ffffff" stop-opacity="0.8"/>
          <stop offset="1" stop-color="#ffffff" stop-opacity="0"/>
        </radialGradient>
        <linearGradient id="rim" x1="0" y1="0" x2="0.6" y2="1">
          <stop offset="0" stop-color="#ffffff" stop-opacity="0.95"/>
          <stop offset="0.55" stop-color="#ffffff" stop-opacity="0.55"/>
          <stop offset="1" stop-color="#abc1d7" stop-opacity="0.5"/>
        </linearGradient>
        <filter id="shadow" x="0" y="0" width="1024" height="1024"
            filterUnits="userSpaceOnUse" color-interpolation-filters="sRGB">
          <feGaussianBlur in="SourceAlpha" stdDeviation="11"/>
          <feOffset dy="10" result="blurredAlpha"/>
          <feFlood flood-color="#1a3048" flood-opacity="0.19"/>
          <feComposite in2="blurredAlpha" operator="in"/>
        </filter>
      </defs>
      <use xlink:href="#tile" fill="#000000" filter="url(#shadow)"/>
      <use xlink:href="#tile" fill="url(#surface)"/>
      <use xlink:href="#tile" fill="url(#light)"/>
      <use xlink:href="#tile" fill="none" stroke="url(#rim)" stroke-width="2"/>
      <image x="132" y="132" width="760" height="760"
          preserveAspectRatio="xMidYMid meet"
          xlink:href="data:image/png;base64,{logo}"/>
    </svg>'''


def main() -> int:
    if not SOURCE.is_file():
        raise SystemExit(f"Icône source introuvable : {SOURCE}")

    # Rasterize the vector tile at Retina resolution, then downsample once.
    # No opaque page background: corners and the outer shadow keep their alpha.
    png = resvg_py.svg_to_bytes(svg_string=icon_svg(SOURCE.read_bytes()), zoom=2.0)
    with Image.open(BytesIO(png)) as rendered:
        image = rendered.convert("RGBA").resize((1024, 1024), Image.Resampling.LANCZOS)

    WINDOWS_ICON.parent.mkdir(parents=True, exist_ok=True)
    image.save(APP_ICON, format="PNG", optimize=True)
    image.save(PNG_ICON, format="PNG", optimize=True)
    image.save(MACOS_ICON, format="ICNS")
    image.resize((512, 512), Image.Resampling.LANCZOS).save(
        LINUX_ICON,
        format="PNG",
        optimize=True,
    )
    image.save(
        WINDOWS_ICON,
        format="ICO",
        sizes=[(16, 16), (24, 24), (32, 32), (48, 48), (64, 64), (128, 128), (256, 256)],
    )

    print(f"Logo source : {SOURCE}")
    print(f"Icône d'application : {APP_ICON}")
    print(f"Icône Windows : {WINDOWS_ICON}")
    print(f"Icône macOS : {MACOS_ICON}")
    print(f"Icône Linux : {LINUX_ICON}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

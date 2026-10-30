"""Maintainer-only offline browser verification of an oxDNA film or 3D view."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('html', type=Path)
    parser.add_argument('--browser-executable', required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--preparation', action='store_true', help='Check a static prepared configuration')
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    from playwright.sync_api import sync_playwright

    errors, requests, checks = [], [], []
    video = None
    with sync_playwright() as playwright:
        browser = playwright.chromium.launch(executable_path=args.browser_executable, headless=True,
            args=['--use-angle=swiftshader', '--enable-unsafe-swiftshader'])
        context = browser.new_context(offline=True, viewport={'width': 980, 'height': 740}, accept_downloads=True)
        page = context.new_page()
        page.on('pageerror', lambda error: errors.append(str(error)))
        page.on('request', lambda request: requests.append(request.url))
        page.set_content(args.html.read_text(), wait_until='load')
        page.wait_for_function("document.querySelector('#info').textContent.length > 0")
        first = page.locator('#view').screenshot(path=str(args.output / 'first-frame.png'))
        atoms = page.evaluate('viewer.getModel().selectedAtoms({}).length')
        frames = int(page.locator('#slider').get_attribute('max')) + 1
        assert atoms > 0, atoms
        if args.preparation:
            assert frames == 1, frames
            assert not page.locator('#playback').is_visible(), 'Film controls on a static configuration'
            assert page.locator('#fit').is_visible(), 'Missing reset-view control'
            page.evaluate('viewer.rotate(40, "y"); viewer.render()')
            rotated = page.locator('#view').screenshot(path=str(args.output / 'rotated.png'))
            assert hashlib.sha256(first).digest() != hashlib.sha256(rotated).digest(), '3D rotation did not change the view'
            checks.extend(['static prepared configuration', '3D rotation', 'no film controls'])
        else:
            assert frames > 1, frames
            page.locator('#slider').focus()
            page.locator('#slider').press('End')
            page.wait_for_function('Number(document.querySelector("#slider").value) === Number(document.querySelector("#slider").max)')
            last = page.locator('#view').screenshot(path=str(args.output / 'last-frame.png'))
            assert hashlib.sha256(first).digest() != hashlib.sha256(last).digest(), 'The trajectory did not move'
            page.locator('#play').click()
            page.wait_for_function('Number(document.querySelector("#slider").value) < Number(document.querySelector("#slider").max)')
            page.locator('#play').click()
            page.locator('#speed').select_option('24')
            with page.expect_download(timeout=60000) as event:
                page.locator('#record').click()
            download = event.value
            video = args.output / download.suggested_filename
            download.save_as(video)
            assert video.stat().st_size > 1000, 'Empty video'
            checks.extend(['actual trajectory playback', 'frame slider', 'video export'])
        page.set_viewport_size({'width': 380, 'height': 740})
        page.screenshot(path=str(args.output / 'compact.png'), full_page=True)
        assert page.evaluate('document.documentElement.scrollWidth <= innerWidth'), 'Horizontal overflow'
        assert not errors, errors
        assert not requests, requests
        browser.close()
    report = {'ok': True, 'offline': True, 'frames': frames, 'coarse_grain_sites': atoms,
              'video': str(video) if video else None, 'video_bytes': video.stat().st_size if video else 0,
              'javascript_errors': errors, 'network_requests': requests,
              'checks': checks + ['compact viewport', 'no network requests']}
    (args.output / 'report.json').write_text(json.dumps(report, indent=2))
    print(json.dumps(report))


if __name__ == '__main__':
    main()

# Maintainer-only native Windows launcher build. Scientific code stays in WSL2.
param([Parameter(Mandatory=$true)][string]$BackendArchive)
$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest
Set-Location (Split-Path -Parent $PSScriptRoot)
if ($env:OS -ne 'Windows_NT' -or [System.Runtime.InteropServices.RuntimeInformation]::OSArchitecture -ne 'X64') {
    throw 'Build this launcher on Windows x64.'
}
$backend = (Resolve-Path $BackendArchive).Path
$python = if ($env:PYTHON_BIN) { $env:PYTHON_BIN } else { 'python' }
function Invoke-Python {
    & $script:python @args
    if ($LASTEXITCODE -ne 0) { throw "Python command failed (exit $LASTEXITCODE)." }
}
# Do not install the scientific Python dependencies into the Windows launcher.
Invoke-Python -m pip install 'pyinstaller==6.21.0' 'pillow==12.3.0'
Invoke-Python scripts/check_bundled_app.py --archive $backend
$buildRoot = Join-Path ([System.IO.Path]::GetTempPath()) ('aptaswitch-windows-' + [guid]::NewGuid())
New-Item -ItemType Directory -Path $buildRoot | Out-Null
try {
    $icon = Join-Path $buildRoot 'AppIcon.ico'
    Invoke-Python -c 'from PIL import Image; import sys; Image.open(sys.argv[1]).save(sys.argv[2], sizes=[(16,16),(32,32),(48,48),(64,64),(128,128),(256,256)])' 'src/aptaswitch_studio/assets/app_icon.png' $icon
    Invoke-Python -m PyInstaller --noconfirm --clean --onedir --windowed --name Apta2Switch-Studio --icon $icon --distpath (Join-Path $buildRoot 'dist') --workpath (Join-Path $buildRoot 'build') --specpath $buildRoot packaging/windows/launcher.py
    $package = Join-Path $buildRoot 'dist/Apta2Switch-Studio'
    $runtime = Join-Path $package 'runtime'
    New-Item -ItemType Directory -Path $runtime | Out-Null
    Copy-Item $backend (Join-Path $runtime 'Apta2Switch-Studio-linux-x64.tar.gz')
    Copy-Item packaging/windows/install_backend.py $runtime
    Copy-Item LICENSE $package
    Copy-Item packaging/windows/START-HERE.txt (Join-Path $package 'START-HERE.txt')
    Invoke-Python scripts/windows_package.py manifest $package
    $check = Join-Path $buildRoot 'launcher-check.json'
    $process = Start-Process -FilePath (Join-Path $package 'Apta2Switch-Studio.exe') -ArgumentList @('--check-package', ('"' + $check + '"')) -PassThru
    if (-not $process.WaitForExit(180000)) { $process.Kill(); throw 'Windows launcher self-check timed out.' }
    if ($process.ExitCode -ne 0 -or -not (Test-Path $check)) { throw 'The frozen Windows launcher did not complete its self-check.' }
    if (-not (Get-Content -Raw $check | ConvertFrom-Json).ok) { throw (Get-Content -Raw $check) }
    New-Item -ItemType Directory -Force dist | Out-Null
    # Python ZIP avoids Compress-Archive's omissions of hidden runtime files.
    Invoke-Python scripts/windows_package.py archive $package 'dist/Apta2Switch-Studio-windows-x64.zip'
    Write-Host 'Archive: dist/Apta2Switch-Studio-windows-x64.zip'
} finally {
    Remove-Item -LiteralPath $buildRoot -Recurse -Force
}

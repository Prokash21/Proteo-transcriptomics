param(
  [string]$RscriptPath = "C:\Program Files\R\R-4.5.1\bin\Rscript.exe",
  [switch]$InstallPackages
)

$ErrorActionPreference = "Stop"
$AnalysisDir = Resolve-Path (Join-Path $PSScriptRoot "..")

if (-not (Test-Path $RscriptPath)) {
  throw "Rscript not found at: $RscriptPath"
}

if ($InstallPackages) {
  $installLog = Join-Path $AnalysisDir "logs\install_clusterProfiler_packages.log"
  $installErr = Join-Path $AnalysisDir "logs\install_clusterProfiler_packages.err.log"
  $installProc = Start-Process -FilePath $RscriptPath -ArgumentList @((Join-Path $PSScriptRoot "install_clusterProfiler_packages.R"), $AnalysisDir) -Wait -PassThru -NoNewWindow -RedirectStandardOutput $installLog -RedirectStandardError $installErr
  Get-Content -Path $installLog -ErrorAction SilentlyContinue
  Get-Content -Path $installErr -ErrorAction SilentlyContinue
  if ($installProc.ExitCode -ne 0) {
    throw "Package installation failed. See $installLog"
  }
}

$runLog = Join-Path $AnalysisDir "logs\run_clusterProfiler_turquoise_module.log"
$runErr = Join-Path $AnalysisDir "logs\run_clusterProfiler_turquoise_module.err.log"
$runProc = Start-Process -FilePath $RscriptPath -ArgumentList @((Join-Path $PSScriptRoot "run_clusterProfiler_turquoise_module.R"), $AnalysisDir) -Wait -PassThru -NoNewWindow -RedirectStandardOutput $runLog -RedirectStandardError $runErr
Get-Content -Path $runLog -ErrorAction SilentlyContinue
Get-Content -Path $runErr -ErrorAction SilentlyContinue
if ($runProc.ExitCode -ne 0) {
  throw "clusterProfiler analysis failed. See $runLog"
}

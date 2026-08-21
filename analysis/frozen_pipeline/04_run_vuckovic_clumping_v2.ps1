[CmdletBinding()]
param()

Set-StrictMode -Version Latest
$ErrorActionPreference = 'Stop'

$ProjectRoot = (Get-Location).Path
$ProtocolPath = Join-Path $ProjectRoot 'docs\protocol\analysis_plan_v1.1.md'
$PlinkExe = Join-Path $ProjectRoot 'tools\plink2\plink2.exe'
$ReferencePrefix = Join-Path $ProjectRoot 'resources\ld\1kg_v3\EUR'
$IncludedClumpInputPath = Join-Path $ProjectRoot 'data_derived\clumping_inputs\vuckovic_hb_clump_input_apoe_included_v1.tsv'
$ExcludedClumpInputPath = Join-Path $ProjectRoot 'data_derived\clumping_inputs\vuckovic_hb_clump_input_apoe_excluded_v1.tsv'
$IncludedPrefix = Join-Path $ProjectRoot 'data_derived\clumped\vuckovic_hb_apoe_included_v2'
$ExcludedPrefix = Join-Path $ProjectRoot 'data_derived\clumped\vuckovic_hb_apoe_excluded_v2'
$ValidatorPath = Join-Path $ProjectRoot 'R\04_validate_vuckovic_clumping_v2.R'
$PipelineLog = Join-Path $ProjectRoot 'results\logs\vuckovic_hb_clumping_pipeline_v2.log'

$ExpectedSha256 = @{
    Included = 'b9e1c00f45213dc574c34ff6c11ca548876445a36b653b22a38205cbcb9a5b50'
    Excluded = 'ccb7c3adb6338f444bd01f9731f83f69c858ffb84ceda9349423b813c5280858'
    Plink = '247491bfca7512e070dc99d6565e9fc56f3a52ad5afc01286016271d34c4992f'
}

$V2Targets = @(
    "$IncludedPrefix.clumps", "$IncludedPrefix.log", "$IncludedPrefix.clumps.missing_id", "$IncludedPrefix.clumps.missing_allele",
    "$ExcludedPrefix.clumps", "$ExcludedPrefix.log", "$ExcludedPrefix.clumps.missing_id", "$ExcludedPrefix.clumps.missing_allele",
    (Join-Path $ProjectRoot 'data_derived\instruments\vuckovic_hb_instruments_apoe_included_v2.parquet'),
    (Join-Path $ProjectRoot 'data_derived\instruments\vuckovic_hb_instruments_apoe_included_v2.tsv'),
    (Join-Path $ProjectRoot 'data_derived\instruments\vuckovic_hb_instruments_apoe_excluded_v2.parquet'),
    (Join-Path $ProjectRoot 'data_derived\instruments\vuckovic_hb_instruments_apoe_excluded_v2.tsv'),
    (Join-Path $ProjectRoot 'results\qc\vuckovic_hb_clumping_v2.json'),
    (Join-Path $ProjectRoot 'results\qc\vuckovic_hb_clumping_comparison_v2.csv'),
    (Join-Path $ProjectRoot 'results\qc\vuckovic_hb_clumped_by_chr_v2.csv'),
    $PipelineLog
)

function Assert-LeafFile {
    param([Parameter(Mandatory)][string]$Path, [Parameter(Mandatory)][string]$Label)
    if (-not (Test-Path -LiteralPath $Path -PathType Leaf)) { throw "Missing $($Label): $Path" }
}
function Assert-Sha256 {
    param([Parameter(Mandatory)][string]$Path, [Parameter(Mandatory)][string]$Expected, [Parameter(Mandatory)][string]$Label)
    $Actual = (Get-FileHash -LiteralPath $Path -Algorithm SHA256).Hash.ToLowerInvariant()
    if ($Actual -ne $Expected) { throw "SHA-256 mismatch for $($Label). Expected $Expected; observed $Actual" }
}
function Write-PipelineLog {
    param([Parameter(Mandatory)][string]$Message)
    Add-Content -LiteralPath $PipelineLog -Value "$(Get-Date -Format o) $Message" -Encoding utf8
}
function Invoke-PlinkClump {
    param(
        [Parameter(Mandatory)][string]$AnalysisName,
        [Parameter(Mandatory)][string]$ClumpInputPath,
        [Parameter(Mandatory)][string]$OutputPrefix
    )
    if ($null -eq $ClumpInputPath) { throw "$($AnalysisName) ClumpInputPath is NULL" }
    if ([string]::IsNullOrWhiteSpace($ClumpInputPath)) { throw "$($AnalysisName) ClumpInputPath is empty" }
    Assert-LeafFile -Path $ClumpInputPath -Label "$($AnalysisName) clumping input"
    if ((Get-Item -LiteralPath $ClumpInputPath).Length -le 0) { throw "$($AnalysisName) clumping input has zero bytes" }
    foreach ($OutputFile in @("$OutputPrefix.clumps", "$OutputPrefix.log", "$OutputPrefix.clumps.missing_id", "$OutputPrefix.clumps.missing_allele")) {
        if (Test-Path -LiteralPath $OutputFile) { throw "Refusing to overwrite $($AnalysisName) output: $OutputFile" }
        if (Test-Path -LiteralPath "$OutputFile.partial") { throw "Residual partial $($AnalysisName) output: $OutputFile.partial" }
    }
    $PlinkArguments = @(
        '--bfile', $ReferencePrefix,
        '--clump', $ClumpInputPath,
        '--clump-id-field', 'SNP',
        '--clump-p-field', 'LOG10P',
        '--clump-a1-field', 'effect_allele',
        '--clump-force-a1',
        '--clump-log10',
        '--clump-log10-p1', '7.301029995663981',
        '--clump-log10-p2', '0',
        '--clump-r2', '0.001',
        '--clump-kb', '10000',
        '--clump-unphased',
        '--threads', '8',
        '--memory', '8000',
        '--out', $OutputPrefix
    )
    $ClumpPosition = [Array]::IndexOf([string[]]$PlinkArguments, '--clump')
    if ($ClumpPosition -lt 0 -or $ClumpPosition -ge ($PlinkArguments.Count - 1) -or $PlinkArguments[$ClumpPosition + 1] -ne $ClumpInputPath) {
        throw "$($AnalysisName) parameter array does not place --clump immediately before ClumpInputPath"
    }
    Write-PipelineLog "$($AnalysisName) parameter array begins."
    for ($Index = 0; $Index -lt $PlinkArguments.Count; $Index++) { Write-PipelineLog "$($AnalysisName) arg[$Index]=$($PlinkArguments[$Index])" }
    Write-PipelineLog "$($AnalysisName) parameter array ends."
    $Started = Get-Date
    Write-PipelineLog "$($AnalysisName) started."
    & $PlinkExe @PlinkArguments
    $PlinkExitCode = $LASTEXITCODE
    $Finished = Get-Date
    Write-PipelineLog "$($AnalysisName) finished. exit_code=$PlinkExitCode duration=$((New-TimeSpan -Start $Started -End $Finished).ToString())"
    if ($PlinkExitCode -ne 0) { throw "$($AnalysisName) PLINK exit code was $PlinkExitCode; stopping pipeline." }
}

try {
    if ((Split-Path -Leaf $ProjectRoot) -ne 'hb_delirium_bidir_mr') { throw "Run from E:\Research\hb_delirium_bidir_mr; current directory is $ProjectRoot" }
    foreach ($Pair in @(
        @{Path=$ProtocolPath;Label='current protocol'}, @{Path=$PlinkExe;Label='PLINK executable'},
        @{Path="$ReferencePrefix.bed";Label='EUR.bed'}, @{Path="$ReferencePrefix.bim";Label='EUR.bim'}, @{Path="$ReferencePrefix.fam";Label='EUR.fam'},
        @{Path=$IncludedClumpInputPath;Label='APOE-included input'}, @{Path=$ExcludedClumpInputPath;Label='APOE-excluded input'},
        @{Path=$ValidatorPath;Label='V2 R validator'}
    )) { Assert-LeafFile -Path $Pair.Path -Label $Pair.Label }
    foreach ($Target in $V2Targets) {
        if (Test-Path -LiteralPath $Target) { throw "Refusing to overwrite V2 target: $Target" }
        if (Test-Path -LiteralPath "$Target.partial") { throw "Residual partial V2 target: $Target.partial" }
    }
    Assert-Sha256 -Path $IncludedClumpInputPath -Expected $ExpectedSha256.Included -Label 'APOE-included input'
    Assert-Sha256 -Path $ExcludedClumpInputPath -Expected $ExpectedSha256.Excluded -Label 'APOE-excluded input'
    Assert-Sha256 -Path $PlinkExe -Expected $ExpectedSha256.Plink -Label 'PLINK executable'
    New-Item -ItemType Directory -Force -Path (Split-Path -Parent $IncludedPrefix), (Split-Path -Parent $PipelineLog) | Out-Null
    New-Item -ItemType File -Path $PipelineLog -ErrorAction Stop | Out-Null
    $PlinkVersion = (& $PlinkExe --version 2>&1 | Out-String).Trim()
    Write-PipelineLog 'Protocol=docs/protocol/analysis_plan_v1.1.md'
    Write-PipelineLog "PLINK=$PlinkVersion"
    Write-PipelineLog 'Frozen parameters: r2=0.001 kb=10000 log10p1=7.301029995663981 log10p2=0 threads=8 memory_mb=8000'
    Write-PipelineLog 'No --clump-allow-overlap; no proxy, liftOver, FinnGen, harmonisation, or MR.'
    Invoke-PlinkClump -AnalysisName 'APOE included' -ClumpInputPath $IncludedClumpInputPath -OutputPrefix $IncludedPrefix
    Invoke-PlinkClump -AnalysisName 'APOE excluded' -ClumpInputPath $ExcludedClumpInputPath -OutputPrefix $ExcludedPrefix
    $RscriptExe = (Get-Command Rscript.exe -ErrorAction SilentlyContinue).Source
    if (-not $RscriptExe) { $RscriptExe = 'C:\Program Files\R\R-4.6.1\bin\Rscript.exe' }
    Assert-LeafFile -Path $RscriptExe -Label 'Rscript executable'
    & $RscriptExe $ValidatorPath -Root $ProjectRoot -PlinkPath $PlinkExe -PlinkVersion $PlinkVersion
    $ValidatorExitCode = $LASTEXITCODE
    Write-PipelineLog "V2 validator finished. exit_code=$ValidatorExitCode"
    if ($ValidatorExitCode -ne 0) { throw "V2 validator exit code was $ValidatorExitCode" }
    Write-PipelineLog 'V2 pipeline completed successfully.'
    exit 0
}
catch {
    if (Test-Path -LiteralPath $PipelineLog) { Write-PipelineLog "FAILED: $($_.Exception.Message)" }
    Write-Error $_
    exit 1
}


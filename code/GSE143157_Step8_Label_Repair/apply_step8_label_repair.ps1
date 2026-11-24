param(
    [string]$Target = "D:\Paper\MSA_Research\02_Epigenetics\GSE143157\02_Preprocessing\Step8_BiasAware_Pathway_Enrichment",
    [switch]$NoOpenReport
)

$ErrorActionPreference = "Stop"
$Payload = $PSScriptRoot

$ReplaceFiles = @(
    "gometh_GO_5mC.csv",
    "gometh_KEGG_5mC.csv",
    "gometh_GO_5hmC.csv",
    "gometh_KEGG_5hmC.csv",
    "top_terms_all.csv",
    "significant_terms_FDR05.csv",
    "top_enrichment_5mC.png",
    "top_enrichment_5hmC.png",
    "STEP8_REPORT.html",
    "RUN_CONFIGURATION.json",
    "step8_summary.json",
    "SUCCESS.txt"
)

$NewFiles = @(
    "label_repair_audit.csv",
    "LABEL_REPAIR_APPLIED.txt"
)

$ExpectedOriginal = @{
    "gometh_GO_5mC.csv" = "348CE85DDF0BDBAC2D25B3F304F8962AB3D13621EE7069B3085B1C0452D22F0E"
    "gometh_KEGG_5mC.csv" = "86CF47BC1CEAD3BEE258527E9422E58E0C81D015FA34AF4592CED66E41CC59E5"
    "gometh_GO_5hmC.csv" = "04A944EAE779A008240D3692DE8A35CEF56EA9D21EECB319BC131FEDFF2F4636"
    "gometh_KEGG_5hmC.csv" = "5C3C767A607CAA9052D6C2E358F341D899FFFBB8C93D781DDEC1FE281A3A3148"
    "top_terms_all.csv" = "BAAAFCE8392CCD4A74F24AE2CB5B93E42AF3CD5F219A7A36A8B56BB284F84B1D"
    "significant_terms_FDR05.csv" = "DB2AD219A456A830C8D4FBD05B3AAD6F95446C8FF283A598B3A470C016E31A8C"
    "top_enrichment_5mC.png" = "E9A5E14588ABDBFBA7E06EE5293C37B8FBCD2997BFCF826519B9A2092AC051B0"
    "top_enrichment_5hmC.png" = "5D8616ACAEE3D2140BCB796CEE68BDB088B8E0565743EC492599CA2721713038"
    "STEP8_REPORT.html" = "7985D00F9BF2B159D6B1F5D4FF922B0AEDBD48AA57A4A6FE9660F0456DE7E54D"
    "RUN_CONFIGURATION.json" = "11F4662749BB61723A766CD820492EF4FD76901A8941BF2BF7E4085E7669A16E"
    "step8_summary.json" = "BB2A4B924E571957D14D302F340B9B34A2B6A567D28B0F74B4663043467FE596"
    "SUCCESS.txt" = "F30B71F6CC85EAF0981E6FCC94072A147732353C1F88618ABE57B40755B42E14"
}

$ExpectedPatch = @{
    "gometh_GO_5mC.csv" = "153ADB2E81ED01445E96361122174A4978E03380D8CF31882BC3E86C3F2E206E"
    "gometh_KEGG_5mC.csv" = "EF4F5E39DEB922E3E23D0D7E0B3177CC3F3CB7A5E52863836CD7E0E55EEDADD7"
    "gometh_GO_5hmC.csv" = "CAE7DE71E5F01DD20F8299A5A676D2808CCF1F2606F2CC79E371DEDA5D961E9C"
    "gometh_KEGG_5hmC.csv" = "968D245086D8EE41B56136EE4844609D1D91F65D99D4069FE5F7A87C9143BAF8"
    "top_terms_all.csv" = "6D32AF944F15DE42B3DF3B9F0A1DE8874CB099F19D121036E824DD38DB23F8EA"
    "significant_terms_FDR05.csv" = "3280CD7990EF92BE8716D7147EA3281B8C9BCEDB96B9147FBB081D0777DCE089"
    "top_enrichment_5mC.png" = "14E94FFE68F06F160584A401078FA52A6E2DC8C3A7A3426838B74F285EC7A04A"
    "top_enrichment_5hmC.png" = "3C7D5A1CB1FD4291D61A52670D95F4DE7E7C8B0A4CCF3A0551A9E94973C99BA0"
    "STEP8_REPORT.html" = "596573370EF110FBB869EF5EFB7C86F1A6BB6C79C3C09546E2C0DB3DC7D3A6F6"
    "RUN_CONFIGURATION.json" = "C1F799C74AF1E8C0D309FABA7F205C267ABE4BDCCA5CB114E669075652D28A6A"
    "step8_summary.json" = "7B9DAF3B82F37DF9D069AC2E292FD713397CDA3590BCEF5FF0E982910DEB92EC"
    "SUCCESS.txt" = "7D399E028C7B40F53F55D46BB3FA22844A833A78779A9DD8EF7F006C07449F95"
    "label_repair_audit.csv" = "EF0375D997B0B5594A65471419D9995D9786E779564C89CB5806D5975E2CC8CE"
    "LABEL_REPAIR_APPLIED.txt" = "847942AF6B4F3C485BFAF49839D0004DE923D3ECDF7CFA199D2B9116E51E38FC"
}

function Get-Sha256([string]$Path) {
    return (Get-FileHash -LiteralPath $Path -Algorithm SHA256).Hash.ToUpperInvariant()
}

Write-Host ""
Write-Host "============================================================"
Write-Host "GSE143157 STEP 8 - VERIFIED LABEL-ONLY REPAIR"
Write-Host "============================================================"
Write-Host "GOmeth will NOT be rerun. Term IDs, P-values, FDR, and genes stay unchanged."
Write-Host ""

if (-not (Test-Path -LiteralPath $Target -PathType Container)) {
    throw "Step 8 output folder was not found: $Target"
}

$allFiles = $ReplaceFiles + $NewFiles
foreach ($name in $allFiles) {
    $source = Join-Path $Payload $name
    if (-not (Test-Path -LiteralPath $source -PathType Leaf)) {
        throw "Patch payload is incomplete: $name"
    }
    if ((Get-Sha256 $source) -ne $ExpectedPatch[$name]) {
        throw "Patch payload hash failed: $name"
    }
}

$alreadyApplied = Test-Path -LiteralPath (Join-Path $Target "LABEL_REPAIR_APPLIED.txt") -PathType Leaf
if ($alreadyApplied) {
    $allCurrentMatch = $true
    foreach ($name in $allFiles) {
        $current = Join-Path $Target $name
        if (-not (Test-Path -LiteralPath $current -PathType Leaf) -or (Get-Sha256 $current) -ne $ExpectedPatch[$name]) {
            $allCurrentMatch = $false
            break
        }
    }
    if ($allCurrentMatch) {
        Write-Host "Step 8 label repair was already applied and all hashes match."
        if (-not $NoOpenReport) {
            Start-Process (Join-Path $Target "STEP8_REPORT.html")
        }
        exit 0
    }
    throw "A repair marker exists, but the repaired-file hashes do not match. No files were changed."
}

foreach ($name in $ReplaceFiles) {
    $current = Join-Path $Target $name
    if (-not (Test-Path -LiteralPath $current -PathType Leaf)) {
        throw "Original Step 8 file is missing: $name"
    }
    if ((Get-Sha256 $current) -ne $ExpectedOriginal[$name]) {
        throw "Original hash differs for $name. This patch refuses to overwrite a changed Step 8 run."
    }
}

$stamp = Get-Date -Format "yyyyMMdd_HHmmss"
$staging = Join-Path $Target ("_label_repair_staging_" + $stamp)
$backup = Join-Path $Target ("_backup_before_label_repair_" + $stamp)

New-Item -ItemType Directory -Path $staging -ErrorAction Stop | Out-Null
foreach ($name in $allFiles) {
    Copy-Item -LiteralPath (Join-Path $Payload $name) -Destination (Join-Path $staging $name) -ErrorAction Stop
}

New-Item -ItemType Directory -Path $backup -ErrorAction Stop | Out-Null
foreach ($name in $ReplaceFiles) {
    Copy-Item -LiteralPath (Join-Path $Target $name) -Destination (Join-Path $backup $name) -ErrorAction Stop
}

foreach ($name in $ReplaceFiles) {
    Move-Item -LiteralPath (Join-Path $staging $name) -Destination (Join-Path $Target $name) -Force -ErrorAction Stop
}
foreach ($name in $NewFiles) {
    Move-Item -LiteralPath (Join-Path $staging $name) -Destination (Join-Path $Target $name) -ErrorAction Stop
}
Remove-Item -LiteralPath $staging -ErrorAction Stop

foreach ($name in $allFiles) {
    $published = Join-Path $Target $name
    if ((Get-Sha256 $published) -ne $ExpectedPatch[$name]) {
        throw "Published-file verification failed: $name. Originals are preserved in $backup"
    }
}

$localMarker = Join-Path $Target "LABEL_REPAIR_LOCAL_APPLY.txt"
@(
    "STEP 8 VERIFIED LABEL REPAIR APPLIED",
    "Local apply time: $(Get-Date -Format 'yyyy-MM-dd HH:mm:ss zzz')",
    "Backup folder: $backup",
    "Core GOmeth statistics changed: FALSE"
) | Set-Content -LiteralPath $localMarker -Encoding UTF8

Write-Host ""
Write-Host "STEP 8 LABEL REPAIR FINISHED SUCCESSFULLY" -ForegroundColor Green
Write-Host "Backup: $backup"
Write-Host "All 14 corrected files passed SHA-256 verification."
Write-Host ""
if (-not $NoOpenReport) {
    Start-Process (Join-Path $Target "STEP8_REPORT.html")
}

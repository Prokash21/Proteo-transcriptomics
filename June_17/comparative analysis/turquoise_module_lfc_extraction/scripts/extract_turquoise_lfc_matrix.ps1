param(
  [string]$AnalysisDir = ".."
)

$ErrorActionPreference = "Stop"

function Resolve-FromScriptDir([string]$PathValue) {
  if ([System.IO.Path]::IsPathRooted($PathValue)) { return $PathValue }
  return Join-Path $PSScriptRoot $PathValue
}

function Split-GeneField([string]$Value) {
  @($Value -split ";" | ForEach-Object { $_.Trim() } | Where-Object { $_ })
}

function To-NullableDouble([string]$Value) {
  $parsed = 0.0
  if ([double]::TryParse($Value, [System.Globalization.NumberStyles]::Float, [System.Globalization.CultureInfo]::InvariantCulture, [ref]$parsed)) {
    return $parsed
  }
  return $null
}

function Get-PadjScore($Row) {
  $padj = To-NullableDouble ([string]$Row.padj)
  if ($null -ne $padj) { return $padj }
  $pvalue = To-NullableDouble ([string]$Row.pvalue)
  if ($null -ne $pvalue) { return $pvalue }
  return [double]::PositiveInfinity
}

$analysisRoot = Resolve-FromScriptDir $AnalysisDir
$inputDir = Join-Path $analysisRoot "inputs"
$lfcDir = Join-Path $inputDir "mapped_resLFC"
$resultsDir = Join-Path $analysisRoot "results"
New-Item -ItemType Directory -Force -Path $resultsDir | Out-Null

$moduleFile = Join-Path $inputDir "turquoise_module_genes.csv"
$depMappedFile = Join-Path $inputDir "dep_without_outliers_02_gene_names_mapped.csv"

$moduleRows = Import-Csv $moduleFile
$moduleGenes = New-Object System.Collections.Generic.List[string]
$moduleLookup = [System.Collections.Generic.Dictionary[string,string]]::new([System.StringComparer]::OrdinalIgnoreCase)
foreach ($row in $moduleRows) {
  foreach ($gene in (Split-GeneField ([string]$row.gene))) {
    if (-not $moduleLookup.ContainsKey($gene)) {
      $moduleLookup.Add($gene, $gene)
      $moduleGenes.Add($gene)
    }
  }
}

$sourceFiles = @(
  [pscustomobject]@{
    Column = "GSE_40_colon_cancer_vs_normal_log2FC"
    File = "mapped_GSE_40_02_resLFC_Treatment_colon_cancer_vs_normal_mygene.csv"
  },
  [pscustomobject]@{
    Column = "GSE_45_tumor_vs_normal_log2FC"
    File = "mapped_GSE_45_OR_02_resLFC_condition_tumor_vs_normal_mygene.csv"
  },
  [pscustomobject]@{
    Column = "GSE_60_tumor_vs_normal_log2FC"
    File = "mapped_GSE_60_resLFC_raw_wolTreatment_tumor_vs_normal_mygene.csv"
  },
  [pscustomobject]@{
    Column = "GSE_79_tumor_vs_normal_log2FC"
    File = "mapped_GSE_79_02_resLFC_Treatment_tumor_vs_normal_mygene.csv"
  },
  [pscustomobject]@{
    Column = "TCGA_COAD_Tumor_vs_Normal_log2FC"
    File = "mapped_tcga_coad_resLFC_Treatment_Tumor_vs_Normal_mygene.csv"
  },
  [pscustomobject]@{
    Column = "TCGA_READ_Tumor_vs_Normal_log2FC"
    File = "mapped_tcga_read_resLFC_SampleType_Tumor_vs_Normal_mygene.csv"
  }
)
$depColumn = "DEP_without_outliers_02_logFC"
$lfcColumns = @($sourceFiles.Column) + $depColumn

$matrixRows = [ordered]@{}
foreach ($gene in $moduleGenes) {
  $rowData = [ordered]@{ gene = $gene }
  foreach ($column in $lfcColumns) { $rowData[$column] = "" }
  $matrixRows[$gene] = [pscustomobject]$rowData
}

$selectedScore = @{}
$matchedDetails = New-Object System.Collections.Generic.List[object]

foreach ($source in $sourceFiles) {
  $path = Join-Path $lfcDir $source.File
  if (-not (Test-Path $path)) { throw "Missing LFC file: $path" }

  foreach ($row in (Import-Csv $path)) {
    $geneId = ([string]$row.gene_id).Trim()
    if (-not $geneId) { continue }
    if (-not $moduleLookup.ContainsKey($geneId)) { continue }

    $canonicalGene = $moduleLookup[$geneId]
    $lfc = ([string]$row.log2FoldChange).Trim()
    if (-not $lfc) { continue }

    $scoreKey = "$canonicalGene|$($source.Column)"
    $score = Get-PadjScore $row
    if (-not $selectedScore.ContainsKey($scoreKey) -or $score -lt $selectedScore[$scoreKey]) {
      $selectedScore[$scoreKey] = $score
      $matrixRows[$canonicalGene].$($source.Column) = $lfc
    }

    $matchedDetails.Add([pscustomobject]@{
      source_column = $source.Column
      source_file = $source.File
      module_gene = $canonicalGene
      matched_gene_id = $geneId
      lfc = $lfc
      baseMean = $row.baseMean
      lfcSE = $row.lfcSE
      pvalue = $row.pvalue
      padj = $row.padj
      gene_name = $row.gene_name
      entrez_id = $row.entrez_id
      ensembl_gene_id = $row.ensembl_gene_id
    })
  }
}

if (-not (Test-Path $depMappedFile)) { throw "Missing DEP mapped file: $depMappedFile" }
foreach ($row in (Import-Csv $depMappedFile)) {
  foreach ($symbol in (Split-GeneField ([string]$row.resolved_symbol))) {
    if (-not $moduleLookup.ContainsKey($symbol)) { continue }

    $canonicalGene = $moduleLookup[$symbol]
    $lfc = ([string]$row.logFC).Trim()
    if (-not $lfc) { continue }

    $scoreKey = "$canonicalGene|$depColumn"
    $scoreValue = To-NullableDouble $lfc
    $score = if ($null -ne $scoreValue) { -[Math]::Abs($scoreValue) } else { 0 }
    if (-not $selectedScore.ContainsKey($scoreKey) -or $score -lt $selectedScore[$scoreKey]) {
      $selectedScore[$scoreKey] = $score
      $matrixRows[$canonicalGene].$depColumn = $lfc
    }

    $matchedDetails.Add([pscustomobject]@{
      source_column = $depColumn
      source_file = "dep_without_outliers_02_gene_names_mapped.csv"
      module_gene = $canonicalGene
      matched_gene_id = $symbol
      lfc = $lfc
      baseMean = ""
      lfcSE = ""
      pvalue = ""
      padj = ""
      gene_name = $row.resolved_gene_name
      entrez_id = $row.entrezgene
      ensembl_gene_id = ""
    })
  }
}

$matrix = @($moduleGenes | ForEach-Object { $matrixRows[$_] })
$matrixPath = Join-Path $resultsDir "turquoise_module_lfc_matrix.csv"
$matrix | Export-Csv -Path $matrixPath -NoTypeInformation

$tsvPath = Join-Path $resultsDir "turquoise_module_lfc_matrix.tsv"
$header = @("gene") + $lfcColumns
Set-Content -Path $tsvPath -Value ($header -join "`t") -Encoding UTF8
foreach ($row in $matrix) {
  $values = @($row.gene) + @($lfcColumns | ForEach-Object { $row.$_ })
  Add-Content -Path $tsvPath -Value ($values -join "`t") -Encoding UTF8
}

$matchedDetails | Export-Csv -Path (Join-Path $resultsDir "matched_lfc_long_details.csv") -NoTypeInformation

$coverage = New-Object System.Collections.Generic.List[object]
foreach ($column in $lfcColumns) {
  $matchedCount = @($matrix | Where-Object { [string]$_.($column) }).Count
  $coverage.Add([pscustomobject]@{
    lfc_column = $column
    matched_genes = $matchedCount
    missing_genes = $moduleGenes.Count - $matchedCount
    total_module_genes = $moduleGenes.Count
    percent_matched = [Math]::Round(($matchedCount / [double]$moduleGenes.Count) * 100, 2)
  })
}
$coverage | Export-Csv -Path (Join-Path $resultsDir "lfc_coverage_summary.csv") -NoTypeInformation

$missingRows = foreach ($gene in $moduleGenes) {
  $missingColumns = @($lfcColumns | Where-Object { -not [string]$matrixRows[$gene].$_ })
  if ($missingColumns.Count -gt 0) {
    [pscustomobject]@{
      gene = $gene
      missing_lfc_columns = ($missingColumns -join ";")
      missing_count = $missingColumns.Count
    }
  }
}
$missingRows | Export-Csv -Path (Join-Path $resultsDir "missing_lfc_by_gene.csv") -NoTypeInformation

Write-Output "Module genes: $($moduleGenes.Count)"
Write-Output "LFC columns: $($lfcColumns.Count)"
Write-Output "Matrix: $matrixPath"
$coverage | Format-Table -AutoSize

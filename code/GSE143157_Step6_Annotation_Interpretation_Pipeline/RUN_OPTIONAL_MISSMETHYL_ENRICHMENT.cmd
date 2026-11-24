@echo off
setlocal
title OPTIONAL GSE143157 bias-aware enrichment
set "PIPELINE_DIR=%~dp0"
set "RESULTS=D:\Paper\MSA_Research\02_Epigenetics\GSE143157\02_Preprocessing\Step6_Annotation_Interpretation_Primary77"
cd /d "%PIPELINE_DIR%"

echo.
echo OPTIONAL: bias-aware GO/KEGG enrichment with R and missMethyl
echo This is exploratory and may require 30-90 minutes for first package installation.
echo.
if not exist "%RESULTS%\SUCCESS.txt" (
  echo ERROR: Complete the main Step 6 first.
  pause
  exit /b 2
)
where Rscript >nul 2>nul
if errorlevel 1 (
  echo Rscript was not found. Do not install anything yet.
  echo Send the main Step 6 outputs to ChatGPT first for review.
  pause
  exit /b 2
)
if exist "%RESULTS%\MISSMETHYL_SUCCESS.txt" (
  echo Optional enrichment already completed. Nothing was overwritten.
  pause
  exit /b 0
)
Rscript run_bias_aware_enrichment.R "%RESULTS%"
if errorlevel 1 (
  echo ERROR: Optional enrichment failed. Main Step 6 results are still valid.
  pause
  exit /b 2
)
echo Optional bias-aware enrichment finished.
pause
exit /b 0

@echo off
setlocal EnableExtensions
title GSE143157 Step 8 Bias-aware Pathway Enrichment

set "PIPELINE_DIR=%~dp0"
set "STEP6=D:\Paper\MSA_Research\02_Epigenetics\GSE143157\02_Preprocessing\Step6_Annotation_Interpretation_Primary77"
set "RESULTS=D:\Paper\MSA_Research\02_Epigenetics\GSE143157\02_Preprocessing\Step8_BiasAware_Pathway_Enrichment"

cd /d "%PIPELINE_DIR%"

echo.
echo ============================================================
echo GSE143157 STEP 8 - BIAS-AWARE PATHWAY ENRICHMENT
echo ============================================================
echo Pipeline version: 1.0.4
echo Method: missMethyl::gometh for the Illumina EPIC array
echo Separate tested-probe backgrounds: 5mC and 5hmC
echo Interpretation: exploratory; FDR below 0.05 is the significance rule
echo.

if not exist "%STEP6%\SUCCESS.txt" (
  echo ERROR: Step 6 SUCCESS.txt was not found:
  echo %STEP6%
  pause
  exit /b 2
)

if exist "%RESULTS%\SUCCESS.txt" (
  echo A successful Step 8 result already exists. Nothing was overwritten.
  echo %RESULTS%
  start "" "%RESULTS%\STEP8_REPORT.html"
  pause
  exit /b 0
)

set "RSCRIPT="
where Rscript.exe >nul 2>nul
if not errorlevel 1 set "RSCRIPT=Rscript.exe"

if not defined RSCRIPT if exist "%ProgramFiles%\R" (
  for /f "delims=" %%R in ('dir /b /ad /o-n "%ProgramFiles%\R\R-*" 2^>nul') do if not defined RSCRIPT (
    if exist "%ProgramFiles%\R\%%R\bin\Rscript.exe" set "RSCRIPT=%ProgramFiles%\R\%%R\bin\Rscript.exe"
  )
)

if not defined RSCRIPT (
  echo ERROR: R was not found on this computer.
  echo.
  echo Install the current Windows version of R from:
  echo https://cran.r-project.org/bin/windows/base/
  echo Keep the default installation folder, then run this CMD file again.
  echo.
  pause
  exit /b 2
)

set "PYTHON_LAUNCHER="
where py >nul 2>nul
if not errorlevel 1 set "PYTHON_LAUNCHER=py -3"
if not defined PYTHON_LAUNCHER (
  where python >nul 2>nul
  if not errorlevel 1 set "PYTHON_LAUNCHER=python"
)
if not defined PYTHON_LAUNCHER (
  echo ERROR: Python 3 was not found. It is needed for the final integrity check.
  pause
  exit /b 2
)

echo [1/3] R found:
"%RSCRIPT%" --version
echo.
echo [2/3] Installing any missing R/Bioconductor packages and running GOmeth...
echo Packages will be installed in your writable AppData folder; Administrator access is not needed.
echo Source annotation packages will be installed sequentially, so Rtools is not required.
echo KEGGREST will be installed and verified before its dependent annotation packages.
echo The first run can take 20-60 minutes because annotation packages are large.
echo Keep the internet connected and do not close this window.
"%RSCRIPT%" run_gse143157_step8.R "%STEP6%" "%RESULTS%"
if errorlevel 1 (
  echo.
  echo ERROR: Step 8 stopped. Read the last lines above or send a screenshot here.
  echo Partial files are in: %RESULTS%
  pause
  exit /b 2
)

echo [3/3] Verifying method locks, inputs, P-values, FDR values, and result files...
%PYTHON_LAUNCHER% check_step8_outputs.py "%STEP6%" "%RESULTS%"
if errorlevel 1 (
  echo ERROR: Step 8 integrity check failed. Send a screenshot here.
  pause
  exit /b 2
)

echo.
echo ============================================================
echo STEP 8 FINISHED SUCCESSFULLY
echo ============================================================
echo Report: %RESULTS%\STEP8_REPORT.html
echo.
start "" "%RESULTS%\STEP8_REPORT.html"
pause
exit /b 0

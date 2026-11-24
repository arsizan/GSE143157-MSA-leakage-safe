@echo off
setlocal EnableExtensions
title GSE143157 Step 9 Final Audit and Paper Package

set "PIPELINE_DIR=%~dp0"
set "BASE=D:\Paper\MSA_Research\02_Epigenetics\GSE143157\02_Preprocessing"
set "STEP1=%BASE%\Step1_results_full78"
set "STEP2=%BASE%\Step2_PCA_Outliers_full78"
set "STEP2B=%BASE%\Step2B_Cohorts"
set "STEP3=%BASE%\Step3_Confounders_CV"
set "STEP4=%BASE%\Step4_Covariate_Benchmark"
set "STEP5=%BASE%\Step5_Molecular_Benchmark_Primary77"
set "STEP6=%BASE%\Step6_Annotation_Interpretation_Primary77"
set "STEP7=%BASE%\Step7_Molecular_Sensitivity"
set "STEP8=%BASE%\Step8_BiasAware_Pathway_Enrichment"
set "RESULTS=%BASE%\Step9_Final_Audit_Paper_Package"

cd /d "%PIPELINE_DIR%"

echo.
echo ============================================================
echo GSE143157 STEP 9 - FINAL AUDIT AND PAPER PACKAGE
echo ============================================================
echo This step reads Steps 1-8, verifies the locked decisions,
echo and creates the final source-of-truth tables, figures, hashes,
echo exact paper numbers, and claim guardrails.
echo.

for %%D in ("%STEP1%" "%STEP2%" "%STEP2B%" "%STEP3%" "%STEP4%" "%STEP5%" "%STEP6%" "%STEP7%" "%STEP8%") do (
  if not exist "%%~D\SUCCESS.txt" (
    echo ERROR: Required SUCCESS.txt was not found in:
    echo %%~D
    echo Finish that earlier step before running Step 9.
    pause
    exit /b 2
  )
)

if exist "%RESULTS%\SUCCESS.txt" (
  echo A successful Step 9 freeze already exists. Nothing was overwritten.
  echo %RESULTS%
  start "" "%RESULTS%\STEP9_REPORT.html"
  pause
  exit /b 0
)

set "PYTHON_LAUNCHER="
where py >nul 2>nul
if not errorlevel 1 set "PYTHON_LAUNCHER=py -3"
if not defined PYTHON_LAUNCHER (
  where python >nul 2>nul
  if not errorlevel 1 set "PYTHON_LAUNCHER=python"
)
if not defined PYTHON_LAUNCHER (
  echo ERROR: Python 3 was not found.
  pause
  exit /b 2
)

if not exist ".venv\Scripts\python.exe" (
  echo [1/4] Creating a private Python environment...
  %PYTHON_LAUNCHER% -m venv .venv
  if errorlevel 1 (
    echo ERROR: Python could not create the private environment.
    pause
    exit /b 2
  )
) else (
  echo [1/4] Existing private Python environment found.
)

echo [2/4] Installing the small reporting dependencies...
".venv\Scripts\python.exe" -m pip install --timeout 120 --retries 10 -r requirements_step9.txt
if errorlevel 1 (
  echo ERROR: Package installation failed. Check the internet connection.
  pause
  exit /b 2
)

echo [3/4] Auditing Steps 1-8 and building the final freeze...
".venv\Scripts\python.exe" run_gse143157_step9.py ^
  --step1-dir "%STEP1%" ^
  --step2-dir "%STEP2%" ^
  --step2b-dir "%STEP2B%" ^
  --step3-dir "%STEP3%" ^
  --step4-dir "%STEP4%" ^
  --step5-dir "%STEP5%" ^
  --step6-dir "%STEP6%" ^
  --step7-dir "%STEP7%" ^
  --step8-dir "%STEP8%" ^
  --output "%RESULTS%"
if errorlevel 1 (
  echo ERROR: Step 9 build stopped. Send a screenshot of the error.
  pause
  exit /b 2
)

echo [4/4] Recalculating all frozen hashes and checking final tables...
".venv\Scripts\python.exe" check_step9_outputs.py "%RESULTS%"
if errorlevel 1 (
  echo ERROR: Step 9 integrity check failed. Do not use the package for writing yet.
  pause
  exit /b 2
)

echo.
echo ============================================================
echo STEP 9 FINISHED SUCCESSFULLY
echo ============================================================
echo Final report: %RESULTS%\STEP9_REPORT.html
echo.
start "" "%RESULTS%\STEP9_REPORT.html"
pause
exit /b 0

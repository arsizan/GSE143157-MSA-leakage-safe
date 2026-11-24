@echo off
setlocal
title GSE143157 Step 2B

set "PIPELINE_DIR=%~dp0"
set "STEP1=D:\Paper\MSA_Research\02_Epigenetics\GSE143157\02_Preprocessing\Step1_results_full78"
set "STEP2=D:\Paper\MSA_Research\02_Epigenetics\GSE143157\02_Preprocessing\Step2_PCA_Outliers_full78"
set "RESULTS=D:\Paper\MSA_Research\02_Epigenetics\GSE143157\02_Preprocessing\Step2B_Cohorts"

cd /d "%PIPELINE_DIR%"

echo.
echo ============================================================
echo GSE143157 STEP 2B - FREEZE ANALYSIS COHORTS
echo ============================================================
echo Primary: Clean-77 ^(MSA 30 removed from both modalities^)
echo Sensitivity: Full-78
echo Optional stress-test: 76 ^(MSA 30 and Control 11 removed^)
echo.

if not exist "%STEP1%\SUCCESS.txt" (
  echo ERROR: Step 1 SUCCESS.txt was not found.
  pause
  exit /b 2
)
if not exist "%STEP2%\SUCCESS.txt" (
  echo ERROR: Step 2 SUCCESS.txt was not found.
  pause
  exit /b 2
)
if exist "%RESULTS%\SUCCESS.txt" (
  echo A successful Step 2B result already exists. Nothing was overwritten.
  echo %RESULTS%
  pause
  exit /b 0
)
if exist "%RESULTS%" (
  echo ERROR: The result folder exists but is incomplete or non-empty.
  echo Rename it before trying again: %RESULTS%
  pause
  exit /b 2
)

if not exist ".venv\Scripts\python.exe" (
  echo [1/4] Creating a private Python environment...
  python -m venv .venv
  if errorlevel 1 (
    echo ERROR: Python could not create the environment.
    echo Run python --version in Command Prompt and send a screenshot.
    pause
    exit /b 2
  )
) else (
  echo [1/4] Existing private Python environment found.
)

echo [2/4] Installing required packages...
".venv\Scripts\python.exe" -m pip install -r requirements_step2b.txt
if errorlevel 1 (
  echo ERROR: Package installation failed. Check the internet connection.
  pause
  exit /b 2
)

echo [3/4] Building the memory-safe Clean-77 HDF5 and cohort manifests...
echo The large source file is read in small row chunks. Do not close this window.
".venv\Scripts\python.exe" build_gse143157_step2b.py --step1-dir "%STEP1%" --step2-dir "%STEP2%" --output "%RESULTS%"
if errorlevel 1 (
  echo ERROR: Step 2B stopped. Take a screenshot of the full error and send it here.
  pause
  exit /b 2
)

echo [4/4] Checking donor order, matrix shapes, and selected values...
".venv\Scripts\python.exe" check_step2b_outputs.py "%RESULTS%"
if errorlevel 1 (
  echo ERROR: Step 2B integrity check failed. Send a screenshot here.
  pause
  exit /b 2
)

echo.
echo ============================================================
echo STEP 2B FINISHED SUCCESSFULLY
echo ============================================================
echo Primary cohort: 77 donors ^(37 Control, 40 MSA^)
echo Report: %RESULTS%\STEP2B_REPORT.html
echo.
start "" "%RESULTS%\STEP2B_REPORT.html"
pause
exit /b 0

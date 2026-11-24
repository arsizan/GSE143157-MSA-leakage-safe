@echo off
setlocal
title GSE143157 Step 4

set "PIPELINE_DIR=%~dp0"
set "STEP2B=D:\Paper\MSA_Research\02_Epigenetics\GSE143157\02_Preprocessing\Step2B_Cohorts"
set "STEP3=D:\Paper\MSA_Research\02_Epigenetics\GSE143157\02_Preprocessing\Step3_Confounders_CV"
set "RESULTS=D:\Paper\MSA_Research\02_Epigenetics\GSE143157\02_Preprocessing\Step4_Covariate_Benchmark"

cd /d "%PIPELINE_DIR%"

set "PYTHON_LAUNCHER="
where py >nul 2>nul
if not errorlevel 1 set "PYTHON_LAUNCHER=py -3"
if not defined PYTHON_LAUNCHER (
  where python >nul 2>nul
  if not errorlevel 1 set "PYTHON_LAUNCHER=python"
)

echo.
echo ============================================================
echo GSE143157 STEP 4 - COVARIATE-ONLY BENCHMARK
echo ============================================================
echo Models: age only, non-age covariates, all covariates
echo Frozen nested CV: 5 outer folds x 10 repeats, 4 inner folds
echo.

if not exist "%STEP2B%\SUCCESS.txt" (
  echo ERROR: Step 2B SUCCESS.txt was not found.
  pause
  exit /b 2
)
if not exist "%STEP3%\SUCCESS.txt" (
  echo ERROR: Step 3 SUCCESS.txt was not found.
  pause
  exit /b 2
)
if exist "%RESULTS%\SUCCESS.txt" (
  echo A successful Step 4 result already exists. Nothing was overwritten.
  echo %RESULTS%
  pause
  exit /b 0
)
if exist "%RESULTS%" (
  echo ERROR: The Step 4 result folder exists but is incomplete or non-empty.
  echo Rename that folder before trying again: %RESULTS%
  pause
  exit /b 2
)

if not exist ".venv\Scripts\python.exe" (
  echo [1/4] Creating a private Python environment...
  if not defined PYTHON_LAUNCHER (
    echo ERROR: Python 3 was not found on this computer.
    echo Install Python 3, enable Add Python to PATH, then run this file again.
    pause
    exit /b 2
  )
  %PYTHON_LAUNCHER% -m venv .venv
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
".venv\Scripts\python.exe" -m pip install -r requirements_step4.txt
if errorlevel 1 (
  echo ERROR: Package installation failed. Check the internet connection.
  pause
  exit /b 2
)

echo [3/4] Running leakage-safe covariate benchmarks...
echo This fits many small inner-CV models and can take several minutes.
".venv\Scripts\python.exe" run_gse143157_step4.py --step2b-dir "%STEP2B%" --step3-dir "%STEP3%" --output "%RESULTS%"
if errorlevel 1 (
  echo ERROR: Step 4 stopped. Take a screenshot of the full error and send it here.
  pause
  exit /b 2
)

echo [4/4] Checking every prediction against the frozen Step 3 folds...
".venv\Scripts\python.exe" check_step4_outputs.py "%RESULTS%" --step3-dir "%STEP3%"
if errorlevel 1 (
  echo ERROR: Step 4 integrity check failed. Send a screenshot here.
  pause
  exit /b 2
)

echo.
echo ============================================================
echo STEP 4 FINISHED SUCCESSFULLY
echo ============================================================
echo Report: %RESULTS%\STEP4_REPORT.html
echo.
start "" "%RESULTS%\STEP4_REPORT.html"
pause
exit /b 0

@echo off
setlocal
title GSE143157 Step 3

set "PIPELINE_DIR=%~dp0"
set "STEP2B=D:\Paper\MSA_Research\02_Epigenetics\GSE143157\02_Preprocessing\Step2B_Cohorts"
set "RESULTS=D:\Paper\MSA_Research\02_Epigenetics\GSE143157\02_Preprocessing\Step3_Confounders_CV"

cd /d "%PIPELINE_DIR%"

echo.
echo ============================================================
echo GSE143157 STEP 3 - CONFOUNDERS AND NESTED CV
echo ============================================================
echo Input: %STEP2B%
echo Output: %RESULTS%
echo.

if not exist "%STEP2B%\SUCCESS.txt" (
  echo ERROR: Step 2B SUCCESS.txt was not found.
  pause
  exit /b 2
)
if not exist "%STEP2B%\GSE143157_primary_clean77.h5" (
  echo ERROR: The primary Clean-77 HDF5 was not found.
  pause
  exit /b 2
)
if exist "%RESULTS%\SUCCESS.txt" (
  echo A successful Step 3 result already exists. Nothing was overwritten.
  echo %RESULTS%
  pause
  exit /b 0
)
if exist "%RESULTS%" (
  echo ERROR: The Step 3 result folder exists but is incomplete or non-empty.
  echo Rename that folder before trying again: %RESULTS%
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
".venv\Scripts\python.exe" -m pip install -r requirements_step3.txt
if errorlevel 1 (
  echo ERROR: Package installation failed. Check the internet connection.
  pause
  exit /b 2
)

echo [3/4] Auditing covariates and creating fixed nested-CV folds...
echo No CpG values are used to create these folds.
".venv\Scripts\python.exe" run_gse143157_step3.py --step2b-dir "%STEP2B%" --output "%RESULTS%"
if errorlevel 1 (
  echo ERROR: Step 3 stopped. Take a screenshot of the full error and send it here.
  pause
  exit /b 2
)

echo [4/4] Checking every outer and inner fold...
".venv\Scripts\python.exe" check_step3_outputs.py "%RESULTS%"
if errorlevel 1 (
  echo ERROR: Step 3 integrity check failed. Send a screenshot here.
  pause
  exit /b 2
)

echo.
echo ============================================================
echo STEP 3 FINISHED SUCCESSFULLY
echo ============================================================
echo Report: %RESULTS%\STEP3_REPORT.html
echo.
start "" "%RESULTS%\STEP3_REPORT.html"
pause
exit /b 0

@echo off
setlocal
title GSE143157 Step 2

set "PIPELINE_DIR=%~dp0"
set "STEP1=D:\Paper\MSA_Research\02_Epigenetics\GSE143157\02_Preprocessing\Step1_results_full78"
set "RESULTS=D:\Paper\MSA_Research\02_Epigenetics\GSE143157\02_Preprocessing\Step2_PCA_Outliers_full78"

cd /d "%PIPELINE_DIR%"

echo.
echo ============================================================
echo GSE143157 STEP 2 - PCA AND OUTLIER REVIEW
echo ============================================================
echo Step 1: %STEP1%
echo Results: %RESULTS%
echo.

if not exist "%STEP1%\SUCCESS.txt" (
  echo ERROR: Step 1 SUCCESS.txt was not found.
  echo Expected: %STEP1%\SUCCESS.txt
  pause
  exit /b 2
)

if not exist "%STEP1%\GSE143157_step1_matrices.h5" (
  echo ERROR: The Step 1 HDF5 matrix was not found.
  pause
  exit /b 2
)

if exist "%RESULTS%\SUCCESS.txt" (
  echo A successful Step 2 result already exists.
  echo Nothing was overwritten: %RESULTS%
  pause
  exit /b 0
)

if exist "%RESULTS%" (
  echo ERROR: The Step 2 result folder already exists but is incomplete or non-empty.
  echo Rename that folder before trying again: %RESULTS%
  pause
  exit /b 2
)

if not exist ".venv\Scripts\python.exe" (
  echo [1/4] Creating a private Python environment...
  python -m venv .venv
  if errorlevel 1 (
    echo ERROR: Python could not create the environment.
    echo Run python --version and send a screenshot.
    pause
    exit /b 2
  )
) else (
  echo [1/4] Existing private Python environment found.
)

echo [2/4] Installing the required packages...
".venv\Scripts\python.exe" -m pip install -r requirements_step2.txt
if errorlevel 1 (
  echo ERROR: Package installation failed. Check the internet connection.
  pause
  exit /b 2
)

echo [3/4] Running PCA, correlations, covariate plots, and outlier evidence...
echo No donor will be deleted automatically.
".venv\Scripts\python.exe" run_gse143157_step2.py --step1-dir "%STEP1%" --output "%RESULTS%"
if errorlevel 1 (
  echo ERROR: Step 2 stopped. Take a screenshot of the full error and send it here.
  pause
  exit /b 2
)

echo [4/4] Checking the Step 2 outputs...
".venv\Scripts\python.exe" check_step2_outputs.py "%RESULTS%"
if errorlevel 1 (
  echo ERROR: Step 2 integrity check failed. Send a screenshot here.
  pause
  exit /b 2
)

echo.
echo ============================================================
echo STEP 2 FINISHED SUCCESSFULLY
echo ============================================================
echo Open this report in your browser:
echo %RESULTS%\STEP2_REPORT.html
echo.
pause
exit /b 0

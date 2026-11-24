@echo off
setlocal
title GSE143157 Step 1

set "PIPELINE_DIR=%~dp0"
set "MATRIX=D:\Paper\MSA_Research\GSE143157_BMIQ_normalized_beta.txt\GSE143157_BMIQ_normalized_beta.txt"
set "SOFT=D:\Paper\MSA_Research\GSE143157_BMIQ_normalized_beta.txt\GSE143157_family.soft.gz"
set "RESULTS=D:\Paper\MSA_Research\02_Epigenetics\GSE143157\02_Preprocessing\Step1_results_full78"

cd /d "%PIPELINE_DIR%"

echo.
echo ============================================================
echo GSE143157 STEP 1
echo ============================================================
echo Matrix:  %MATRIX%
echo Metadata: %SOFT%
echo Results: %RESULTS%
echo.

if not exist "%MATRIX%" (
  echo ERROR: The normalized matrix was not found.
  echo Check this path: %MATRIX%
  pause
  exit /b 2
)

if not exist "%SOFT%" (
  echo ERROR: GSE143157_family.soft.gz was not found.
  echo Check this path: %SOFT%
  pause
  exit /b 2
)

if exist "%RESULTS%\SUCCESS.txt" (
  echo A successful Step 1 result already exists.
  echo Nothing was overwritten: %RESULTS%
  pause
  exit /b 0
)

if exist "%RESULTS%" (
  echo ERROR: The result folder already exists but is incomplete or non-empty.
  echo Rename that folder before trying again: %RESULTS%
  pause
  exit /b 2
)

if not exist ".venv\Scripts\python.exe" (
  echo [1/4] Creating a private Python environment...
  python -m venv .venv
  if errorlevel 1 (
    echo ERROR: Python could not create the environment.
    echo Run: python --version
    pause
    exit /b 2
  )
) else (
  echo [1/4] Existing private Python environment found.
)

echo [2/4] Installing the required packages...
".venv\Scripts\python.exe" -m pip install -r requirements_step1.txt
if errorlevel 1 (
  echo ERROR: Package installation failed. Check the internet connection.
  pause
  exit /b 2
)

echo [3/4] Building 5mC and 5hmC matrices. Please keep the computer awake.
".venv\Scripts\python.exe" build_gse143157_step1.py --input "%MATRIX%" --metadata "%SOFT%" --output "%RESULTS%"
if errorlevel 1 (
  echo ERROR: Step 1 stopped. Take a screenshot of the error and send it here.
  pause
  exit /b 2
)

echo [4/4] Checking the output files...
".venv\Scripts\python.exe" check_step1_outputs.py "%RESULTS%"
if errorlevel 1 (
  echo ERROR: Output integrity check failed. Send a screenshot here.
  pause
  exit /b 2
)

echo.
echo ============================================================
echo STEP 1 FINISHED SUCCESSFULLY
echo ============================================================
echo Results are here:
echo %RESULTS%
echo.
pause
exit /b 0

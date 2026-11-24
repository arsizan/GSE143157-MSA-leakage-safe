@echo off
setlocal
title GSE143157 Step 6 Annotation and Interpretation

set "PIPELINE_DIR=%~dp0"
set "STEP2B=D:\Paper\MSA_Research\02_Epigenetics\GSE143157\02_Preprocessing\Step2B_Cohorts"
set "STEP5=D:\Paper\MSA_Research\02_Epigenetics\GSE143157\02_Preprocessing\Step5_Molecular_Benchmark_Primary77"
set "MANIFEST=D:\Paper\MSA_Research\GSE143157_RAW\GPL21145_MethylationEPIC_15073387_v-1-0.csv.gz"
set "RESULTS=D:\Paper\MSA_Research\02_Epigenetics\GSE143157\02_Preprocessing\Step6_Annotation_Interpretation_Primary77"

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
echo GSE143157 STEP 6 - CpG ANNOTATION AND INTERPRETATION AUDIT
echo ============================================================
echo Input: Step 5 Primary Clean-77 stable features
echo Annotation: local Illumina MethylationEPIC v1.0 manifest
echo Expected time: usually 5-20 minutes
echo The large manifest is streamed and is not loaded fully into RAM.
echo.

if not exist "%STEP2B%\SUCCESS.txt" (
  echo ERROR: Step 2B SUCCESS.txt was not found.
  pause
  exit /b 2
)
if not exist "%STEP5%\SUCCESS.txt" (
  echo ERROR: Step 5 SUCCESS.txt was not found.
  pause
  exit /b 2
)
if not exist "%MANIFEST%" (
  echo ERROR: The GPL21145 EPIC manifest was not found:
  echo %MANIFEST%
  pause
  exit /b 2
)
if exist "%RESULTS%\SUCCESS.txt" (
  echo A successful Step 6 result already exists. Nothing was overwritten.
  echo %RESULTS%
  start "" "%RESULTS%\STEP6_REPORT.html"
  pause
  exit /b 0
)

if not exist ".venv\Scripts\python.exe" (
  echo [1/4] Creating a private Python environment...
  if not defined PYTHON_LAUNCHER (
    echo ERROR: Python 3 was not found. Install Python 3 and enable Add Python to PATH.
    pause
    exit /b 2
  )
  %PYTHON_LAUNCHER% -m venv .venv
  if errorlevel 1 (
    echo ERROR: Python could not create the private environment.
    pause
    exit /b 2
  )
) else (
  echo [1/4] Existing private Python environment found.
)

echo [2/4] Installing required scientific packages...
".venv\Scripts\python.exe" -m pip install --timeout 120 --retries 10 -r requirements_step6.txt
if errorlevel 1 (
  echo ERROR: Package installation failed. Check the internet connection.
  pause
  exit /b 2
)

echo [3/4] Streaming the manifest and building annotated evidence tables...
echo Do not close this window.
".venv\Scripts\python.exe" run_gse143157_step6.py --step5-dir "%STEP5%" --step2b-dir "%STEP2B%" --manifest "%MANIFEST%" --published-signals "%PIPELINE_DIR%published_reference_signals.csv" --output "%RESULTS%"
if errorlevel 1 (
  echo ERROR: Step 6 stopped. Send a screenshot of the error.
  pause
  exit /b 2
)

echo [4/4] Verifying candidates, annotations, backgrounds, and outputs...
".venv\Scripts\python.exe" check_step6_outputs.py "%RESULTS%" --step5-dir "%STEP5%"
if errorlevel 1 (
  echo ERROR: Step 6 integrity check failed. Send a screenshot here.
  pause
  exit /b 2
)

echo.
echo ============================================================
echo STEP 6 FINISHED SUCCESSFULLY
echo ============================================================
echo Report: %RESULTS%\STEP6_REPORT.html
echo.
start "" "%RESULTS%\STEP6_REPORT.html"
pause
exit /b 0

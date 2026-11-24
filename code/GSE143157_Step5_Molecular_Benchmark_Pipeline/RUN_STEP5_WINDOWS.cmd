@echo off
setlocal
title GSE143157 Step 5 Molecular Benchmark

set "PIPELINE_DIR=%~dp0"
set "STEP2B=D:\Paper\MSA_Research\02_Epigenetics\GSE143157\02_Preprocessing\Step2B_Cohorts"
set "STEP3=D:\Paper\MSA_Research\02_Epigenetics\GSE143157\02_Preprocessing\Step3_Confounders_CV"
set "STEP4=D:\Paper\MSA_Research\02_Epigenetics\GSE143157\02_Preprocessing\Step4_Covariate_Benchmark"
set "RESULTS=D:\Paper\MSA_Research\02_Epigenetics\GSE143157\02_Preprocessing\Step5_Molecular_Benchmark_Primary77"
set "OMP_NUM_THREADS=1"
set "MKL_NUM_THREADS=1"
set "OPENBLAS_NUM_THREADS=1"

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
echo GSE143157 STEP 5 - PRIMARY MOLECULAR BENCHMARK
echo ============================================================
echo Cohort: Primary Clean-77
echo Models: 5mC, 5hmC, combined x 3 adjustment views
echo Frozen nested CV: 5 outer folds x 10 repeats, 4 inner folds
echo Expected time: often 15-60 minutes; slower computers may take longer
echo Checkpoints are saved after every outer fold.
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
if not exist "%STEP4%\SUCCESS.txt" (
  echo ERROR: Step 4 SUCCESS.txt was not found.
  pause
  exit /b 2
)
if exist "%RESULTS%\SUCCESS.txt" (
  echo A successful Step 5 result already exists. Nothing was overwritten.
  echo %RESULTS%
  start "" "%RESULTS%\STEP5_REPORT.html"
  pause
  exit /b 0
)

if not exist ".venv\Scripts\python.exe" (
  echo [1/4] Creating a private Python environment...
  if not defined PYTHON_LAUNCHER (
    echo ERROR: Python 3 was not found.
    echo Install Python 3 and enable Add Python to PATH.
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
".venv\Scripts\python.exe" -m pip install --timeout 120 --retries 10 -r requirements_step5.txt
if errorlevel 1 (
  echo ERROR: Package installation failed. Check the internet connection.
  pause
  exit /b 2
)

echo [3/4] Running leakage-safe molecular models...
echo Do not close this window. If interrupted, run this file again to resume.
".venv\Scripts\python.exe" run_gse143157_step5.py --step2b-dir "%STEP2B%" --step3-dir "%STEP3%" --step4-dir "%STEP4%" --output "%RESULTS%"
if errorlevel 1 (
  echo ERROR: Step 5 stopped. Completed outer-fold checkpoints are preserved.
  echo Run this file again to resume, or send a screenshot of the error.
  pause
  exit /b 2
)

echo [4/4] Verifying donor predictions, frozen folds, tuning, and outputs...
".venv\Scripts\python.exe" check_step5_outputs.py "%RESULTS%" --step3-dir "%STEP3%"
if errorlevel 1 (
  echo ERROR: Step 5 integrity check failed. Send a screenshot here.
  pause
  exit /b 2
)

echo.
echo ============================================================
echo STEP 5 FINISHED SUCCESSFULLY
echo ============================================================
echo Report: %RESULTS%\STEP5_REPORT.html
echo.
start "" "%RESULTS%\STEP5_REPORT.html"
pause
exit /b 0

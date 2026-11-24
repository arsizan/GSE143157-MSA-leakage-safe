@echo off
setlocal
title GSE143157 Step 7 Sensitivity Validation

set "PIPELINE_DIR=%~dp0"
set "STEP1=D:\Paper\MSA_Research\02_Epigenetics\GSE143157\02_Preprocessing\Step1_results_full78"
set "STEP2B=D:\Paper\MSA_Research\02_Epigenetics\GSE143157\02_Preprocessing\Step2B_Cohorts"
set "STEP3=D:\Paper\MSA_Research\02_Epigenetics\GSE143157\02_Preprocessing\Step3_Confounders_CV"
set "STEP4=D:\Paper\MSA_Research\02_Epigenetics\GSE143157\02_Preprocessing\Step4_Covariate_Benchmark"
set "STEP5=D:\Paper\MSA_Research\02_Epigenetics\GSE143157\02_Preprocessing\Step5_Molecular_Benchmark_Primary77"
set "STEP6=D:\Paper\MSA_Research\02_Epigenetics\GSE143157\02_Preprocessing\Step6_Annotation_Interpretation_Primary77"
set "RESULTS=D:\Paper\MSA_Research\02_Epigenetics\GSE143157\02_Preprocessing\Step7_Molecular_Sensitivity"
set "FULL78=%RESULTS%\sensitivity_full78"
set "CLEAN76=%RESULTS%\sensitivity_no_MSA30_Control11_76"
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
echo GSE143157 STEP 7 - PRESPECIFIED SENSITIVITY VALIDATION
echo ============================================================
echo Cohorts: Full-78 and strict Clean-76
echo Locked models: 5mC, 5hmC, combined + covariates
echo Frozen nested CV: 5 outer folds x 10 repeats, 4 inner folds
echo Expected time: often 30-120 minutes; slower computers may take longer
echo A checkpoint is saved after every outer fold.
echo.

for %%D in ("%STEP1%" "%STEP2B%" "%STEP3%" "%STEP4%" "%STEP5%" "%STEP6%") do (
  if not exist "%%~D\SUCCESS.txt" (
    echo ERROR: Required SUCCESS.txt was not found in:
    echo %%~D
    pause
    exit /b 2
  )
)

if exist "%RESULTS%\SUCCESS.txt" (
  echo A successful Step 7 result already exists. Nothing was overwritten.
  echo %RESULTS%
  start "" "%RESULTS%\STEP7_REPORT.html"
  pause
  exit /b 0
)

if not exist ".venv\Scripts\python.exe" (
  echo [1/6] Creating a private Python environment...
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
  echo [1/6] Existing private Python environment found.
)

echo [2/6] Installing required scientific packages...
".venv\Scripts\python.exe" -m pip install --timeout 120 --retries 10 -r requirements_step7.txt
if errorlevel 1 (
  echo ERROR: Package installation failed. Check the internet connection.
  pause
  exit /b 2
)

if exist "%FULL78%\SUCCESS.txt" (
  echo [3/6] Full-78 result already complete; skipping it.
) else (
  echo [3/6] Running Full-78 sensitivity cohort...
  echo Do not close this window. If interrupted, run this file again to resume.
  ".venv\Scripts\python.exe" run_gse143157_step7.py --step1-dir "%STEP1%" --step2b-dir "%STEP2B%" --step3-dir "%STEP3%" --step4-dir "%STEP4%" --cohort sensitivity_full78 --output "%FULL78%"
  if errorlevel 1 (
    echo ERROR: Full-78 run stopped. Finished checkpoints are preserved.
    echo Run this file again to resume, or send a screenshot of the error.
    pause
    exit /b 2
  )
)

if exist "%CLEAN76%\SUCCESS.txt" (
  echo [4/6] Clean-76 result already complete; skipping it.
) else (
  echo [4/6] Running strict Clean-76 sensitivity cohort...
  echo Do not close this window. If interrupted, run this file again to resume.
  ".venv\Scripts\python.exe" run_gse143157_step7.py --step1-dir "%STEP1%" --step2b-dir "%STEP2B%" --step3-dir "%STEP3%" --step4-dir "%STEP4%" --cohort sensitivity_no_MSA30_Control11_76 --output "%CLEAN76%"
  if errorlevel 1 (
    echo ERROR: Clean-76 run stopped. Finished checkpoints are preserved.
    echo Run this file again to resume, or send a screenshot of the error.
    pause
    exit /b 2
  )
)

echo [5/6] Combining performance and candidate-stability evidence...
".venv\Scripts\python.exe" summarize_gse143157_step7.py --step4-dir "%STEP4%" --step5-dir "%STEP5%" --step6-dir "%STEP6%" --full78-dir "%FULL78%" --clean76-dir "%CLEAN76%" --output "%RESULTS%"
if errorlevel 1 (
  echo ERROR: Step 7 summary generation failed. Send a screenshot here.
  pause
  exit /b 2
)

echo [6/6] Verifying donors, predictions, folds, tuning, candidates, and outputs...
".venv\Scripts\python.exe" check_step7_outputs.py "%RESULTS%"
if errorlevel 1 (
  echo ERROR: Step 7 integrity check failed. Send a screenshot here.
  pause
  exit /b 2
)

echo.
echo ============================================================
echo STEP 7 FINISHED SUCCESSFULLY
echo ============================================================
echo Report: %RESULTS%\STEP7_REPORT.html
echo.
start "" "%RESULTS%\STEP7_REPORT.html"
pause
exit /b 0

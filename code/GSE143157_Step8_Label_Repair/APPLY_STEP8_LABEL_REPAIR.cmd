@echo off
setlocal EnableExtensions
title GSE143157 Step 8 Verified Label Repair
cd /d "%~dp0"
powershell.exe -NoProfile -ExecutionPolicy Bypass -File "%~dp0apply_step8_label_repair.ps1"
if errorlevel 1 (
  echo.
  echo ERROR: The repair was not completed. Send a screenshot of this window.
  pause
  exit /b 1
)
echo.
echo You may now rerun or finish Step 9 using the corrected Step 8 labels.
pause
exit /b 0


@echo off
REM ==============================================================
REM K-CHOPORE Pipeline Launcher (Windows)
REM ==============================================================
REM Usage: run_pipeline.bat E:\kchopore_data [threads]
REM Example: run_pipeline.bat E:\ONT_data 12
REM ==============================================================

setlocal EnableDelayedExpansion

set DATA_DIR=%~1
set THREADS=%~2
set IMAGE=k-chopore:latest

if "%DATA_DIR%"=="" (
    echo ERROR: Please specify the path to your data directory.
    echo.
    echo Usage: %~nx0 E:\kchopore_data [threads]
    echo.
    echo Your data directory must contain:
    echo   data\raw\fastq\WT_C_R1.fastq
    echo   data\raw\fastq\WT_C_R2.fastq
    echo   data\reference\genome\TAIR10_chr_all.fas.fasta
    echo   data\reference\annotations\AtRTDv2_QUASI_19April2016.gtf
    echo   data\raw\summaries\*_sequencing_summary_*.txt
    exit /b 1
)

if "%THREADS%"=="" set THREADS=12

echo ========================================
echo  K-CHOPORE Pipeline Launcher
echo ========================================
echo.
echo Data directory: %DATA_DIR%
echo Threads:        %THREADS%
echo Docker image:   %IMAGE%
echo.

REM Check Docker image exists
docker image inspect %IMAGE% >nul 2>&1
if errorlevel 1 (
    echo ERROR: Docker image '%IMAGE%' not found!
    echo Build it first: docker build -t k-chopore:latest .
    exit /b 1
)

REM Check key files
echo Checking required files...
if exist "%DATA_DIR%\data\reference\genome\TAIR10_chr_all.fas.fasta" (
    echo   OK: Reference genome found
) else (
    echo   MISSING: data\reference\genome\TAIR10_chr_all.fas.fasta
)
if exist "%DATA_DIR%\data\reference\annotations\AtRTDv2_QUASI_19April2016.gtf" (
    echo   OK: GTF annotation found
) else (
    echo   MISSING: data\reference\annotations\AtRTDv2_QUASI_19April2016.gtf
)
if exist "%DATA_DIR%\data\raw\fastq\WT_C_R1.fastq" (
    echo   OK: WT_C_R1.fastq found
) else (
    echo   MISSING: data\raw\fastq\WT_C_R1.fastq
)
if exist "%DATA_DIR%\data\raw\fastq\WT_C_R2.fastq" (
    echo   OK: WT_C_R2.fastq found
) else (
    echo   MISSING: data\raw\fastq\WT_C_R2.fastq
)
echo.

REM Create output dirs
if not exist "%DATA_DIR%\results" mkdir "%DATA_DIR%\results"
if not exist "%DATA_DIR%\logs" mkdir "%DATA_DIR%\logs"

echo Starting K-CHOPORE pipeline in Docker...
echo.

docker run -it --rm ^
    -v "%DATA_DIR%":/workspace ^
    -w /workspace ^
    %IMAGE% ^
    snakemake ^
        --snakefile /workspace/Snakefile ^
        --configfile /workspace/config/config.yml ^
        --cores %THREADS% ^
        --latency-wait 30 ^
        --printshellcmds ^
        --reason ^
        --keep-going ^
        --rerun-incomplete

if %ERRORLEVEL% equ 0 (
    echo.
    echo ========================================
    echo  Pipeline completed successfully!
    echo ========================================
    echo Results: %DATA_DIR%\results\
) else (
    echo.
    echo ========================================
    echo  Pipeline finished with errors!
    echo ========================================
    echo Check logs: %DATA_DIR%\logs\
)

endlocal

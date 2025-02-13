@echo off

setlocal

set "FOLDER_PATH=.\BASE"
if not exist "%FOLDER_PATH%" (
    mkdir "%FOLDER_PATH%"
)

set "FOLDER_PATH=.\INI"
if not exist "%FOLDER_PATH%" (
    mkdir "%FOLDER_PATH%"
)

set "FOLDER_PATH=.\M"
if not exist "%FOLDER_PATH%" (
    mkdir "%FOLDER_PATH%"
)

set "FOLDER_PATH=.\OUTPUT"
if not exist "%FOLDER_PATH%" (
    mkdir "%FOLDER_PATH%"
)

set "FOLDER_PATH=.\PML"
if not exist "%FOLDER_PATH%" (
    mkdir "%FOLDER_PATH%"
)

set "FOLDER_PATH=.\Y"
if not exist "%FOLDER_PATH%" (
    mkdir "%FOLDER_PATH%"
)

endlocal
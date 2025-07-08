# Initial and termination time settings
$current_t=3
$target_t=1000
# Timeout and wait duration settings (unit: seconds)
$timeoutSeconds = 279
$waitSeconds = 70
# Whether initialization is completed
$is_ini = 0
#1: initialized, 0: uninitialized, -1: error in previous run, 2: M matrix calculated, -2: regenerate previous two time steps

$is_error = 0
$currentDir = Get-Location  

# Build target path: parent directory of parent directory
$exePath = Join-Path -Path $currentDir -ChildPath "\x64\Debug\Console1.exe"
$outputFilePath = Join-Path -Path $currentDir -ChildPath "output.txt"
$arguments = "$current_t $target_t $is_ini"

$process = Start-Process -FilePath $exePath -ArgumentList $arguments -WorkingDirectory $currentDir `
    -RedirectStandardOutput $outputFilePath -ErrorAction SilentlyContinue -WindowStyle Hidden -PassThru
<# $process = New-Object System.Diagnostics.Process  
$process.StartInfo = $processInfo  
$process.StartInfo.FileName += " > " + $outputFilePath 
$process.Start() | Out-Null   #>

function ReadAvailableOutput([System.IO.StreamReader]$reader,[int]$lastnum) {  
    $num=$lastnum
    while (!$reader.EndOfStream) {  
        $line = $reader.ReadLine()
        if ($line -match '^\s*\d+\s*$') {  
            #Write-Host  $line
            $num = [int]$line.Trim()  
        }elseif ($line.Contains("finish") ) {
            return 1
        } 
        elseif($line.Contains("error")){
            return 0
        }
        elseif($line.Contains("NaN")){
            Write-Host 'NaN'
        }
        else {  
            #Write-Host  $line
        }
    } 
    if ($num -gt $lastnum){
        $lastnum = $num
        Write-Host $lastnum
    }elseif ($num -eq 3 ){
        write-host $line
    }
    
    return $lastnum  
}  
$startTime = Get-Date
$endTime = Get-Date

$last_t = $current_t
# Timer setup 
$encoding = [System.Text.Encoding]::UTF8  
$file = [System.IO.File]::Open($outputFilePath, [System.IO.FileMode]::Open, [System.IO.FileAccess]::Read, [System.IO.FileShare]::ReadWrite)  
$stdout = New-Object System.IO.StreamReader($file, $encoding) 
while($current_t-le$target_t){
    # Read all currently available outputs

    $current_t = ReadAvailableOutput $stdout $current_t 
    if($current_t -eq 1){break}
    elseif($current_t -eq 3){
        $startTime = Get-Date
        $endTime = Get-Date
    }
    
    $timeSpan = $endTime - $startTime
    $timeSpan = [int]$timeSpan.TotalSeconds
    # Check for new outputs 
    if ($current_t -gt $last_t) {  
        # Update last output 
        $last_t = $current_t  
        # Process new outputs here (e.g., print to console) 
        $is_error = 0
        $startTime = Get-Date
    } elseif ($current_t -eq 0 -or $timeSpan-gt $timeoutSeconds) { # Check for timeout condition
         # Process restart logic
        Write-Host "No new output for $timeoutSeconds seconds, restarting process..."  
        $process.Kill()  
        $process.WaitForExit() 
        if($current_t -eq 0){
            Write-Host "error2"
            $current_t=$last_t
            if($is_error -eq 0){$is_error = 1}
            else{$current_t=$current_t+1}}

        # Restart process (Note: Must reconfigure all required StartInfo properties)
        #$process.StartInfo = $processInfo  
        #$process.Start() | Out-Null   
        $is_ini = -1
        $arguments = "$current_t $target_t $is_ini"
        $process = Start-Process -FilePath $exePath -ArgumentList $arguments -WorkingDirectory $currentDir `
            -RedirectStandardOutput $outputFilePath -ErrorAction SilentlyContinue -WindowStyle Hidden -PassThru
        $file = [System.IO.File]::Open($outputFilePath, [System.IO.FileMode]::Open, [System.IO.FileAccess]::Read, [System.IO.FileShare]::ReadWrite)  
        $stdout = New-Object System.IO.StreamReader($file, $encoding) 
        $last_t = $current_t   # Reset last output tracking after restart  
        $startTime = Get-Date
    } else {  
        Write-Host "wait..."
        Start-Sleep -Seconds $waitSeconds
    } 
    $endTime = Get-Date
}
$process.Kill()  
$process.WaitForExit() 
 

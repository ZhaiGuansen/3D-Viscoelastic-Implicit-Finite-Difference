#初始时刻与中止时刻设置
$current_t=3
$target_t=1000
#超时与等待用时设置（单位：秒）
$timeoutSeconds = 279
$waitSeconds = 70
#是否已经初始化
$is_ini = 0
#1表示已经初始化,0表示未初始化,2:进行初始化但M矩阵已计算,-1表示先前运行出现意外中止,-2表示t-1，t-2时刻值均需重新读取
# 获取当前工作目录  
$is_error = 0
$currentDir = Get-Location  

# 构建目标路径：当前目录的上级目录的上级目录  
$exePath = Join-Path -Path $currentDir -ChildPath "\x64\Debug\Console1.exe"
$outputFilePath = Join-Path -Path $currentDir -ChildPath "output.txt"
$arguments = "$current_t $target_t $is_ini"#exe参数
# 创建ProcessStartInfo对象  
<# $processInfo = New-Object System.Diagnostics.ProcessStartInfo  
$processInfo.FileName = $exePath  
$processInfo.WorkingDirectory = $currentDir
$processInfo.Arguments = $arguments  
$processInfo.RedirectStandardOutput = $false
$processInfo.RedirectStandardError = $true  
$processInfo.UseShellExecute = $false  
$processInfo.CreateNoWindow = $true  # 如果不需要看到控制台窗口，可以设置为true #>

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
            # 转换为整数并输出  
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
# 定时器设置  
$encoding = [System.Text.Encoding]::UTF8  
$file = [System.IO.File]::Open($outputFilePath, [System.IO.FileMode]::Open, [System.IO.FileAccess]::Read, [System.IO.FileShare]::ReadWrite)  
$stdout = New-Object System.IO.StreamReader($file, $encoding) 
while($current_t-le$target_t){
    # 读取当前所有可用的输出  

    $current_t = ReadAvailableOutput $stdout $current_t 
    if($current_t -eq 1){break}
    elseif($current_t -eq 3){
        $startTime = Get-Date
        $endTime = Get-Date
    }
    
    $timeSpan = $endTime - $startTime
    $timeSpan = [int]$timeSpan.TotalSeconds
    # 检查是否有新输出  
    if ($current_t -gt $last_t) {  
        # 更新最后输出  
        $last_t = $current_t  
        # 可以在这里处理新输出，比如打印到控制台  
        $is_error = 0
        $startTime = Get-Date
    } elseif ($current_t -eq 0 -or $timeSpan-gt $timeoutSeconds) {  # 检查是否超时
        # 重启进程逻辑  
        Write-Host "No new output for $timeoutSeconds seconds, restarting process..."  
        $process.Kill()  
        $process.WaitForExit() 
        if($current_t -eq 0){
            Write-Host "error2"
            $current_t=$last_t
            if($is_error -eq 0){$is_error = 1}
            else{$current_t=$current_t+1}}

        # 重新启动进程（注意：这里应该再次设置所有必要的StartInfo属性）  
        #$process.StartInfo = $processInfo  
        #$process.Start() | Out-Null   
        $is_ini = -1
        $arguments = "$current_t $target_t $is_ini"
        $process = Start-Process -FilePath $exePath -ArgumentList $arguments -WorkingDirectory $currentDir `
            -RedirectStandardOutput $outputFilePath -ErrorAction SilentlyContinue -WindowStyle Hidden -PassThru
        $file = [System.IO.File]::Open($outputFilePath, [System.IO.FileMode]::Open, [System.IO.FileAccess]::Read, [System.IO.FileShare]::ReadWrite)  
        $stdout = New-Object System.IO.StreamReader($file, $encoding) 
        $last_t = $current_t  # 重置最后输出，因为我们已经重启了进程  
        $startTime = Get-Date
    } else {  
        Write-Host "wait..."
        Start-Sleep -Seconds $waitSeconds
    } 
    $endTime = Get-Date
}
$process.Kill()  
$process.WaitForExit() 
 
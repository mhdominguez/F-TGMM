$root=Join-Path -Resolve -Path $PSScriptRoot .. 
$srcdirs=@(
    "src"
)

function doit() {
    try{
        #write-host "doing the build"
        . (join-path $PSScriptRoot "build.ps1")
    } catch {
        write-host $_.Exception.Message
    }
    write-host "---------------------------"
}

write-host "--- W A T C H ---" -BackgroundColor green -ForegroundColor black

$srcdirs | ForEach-Object {
    $watcher=new-object System.IO.FileSystemWatcher
    $watcher.Path = (Join-Path $root $_)
    $watcher.Filter = '*.*'
    $watcher.IncludeSubDirectories=$true
    $watcher.NotifyFilter="LastAccess,LastWrite,FileName,DirectoryName"
    $watcher.EnableRaisingEvents=$true

    write-host $watcher.Path -ForegroundColor White

    Register-ObjectEvent -InputObject $watcher -EventName "Changed" -Action $function:doit
    Register-ObjectEvent -InputObject $watcher -EventName "Created" -Action {write-host "Deleted"} #$function:doit
    Register-ObjectEvent -InputObject $watcher -EventName "Deleted" -Action {write-host "Deleted"}
    Register-ObjectEvent -InputObject $watcher -EventName "Disposed" -Action {write-host "Disposed"}
    Register-ObjectEvent -InputObject $watcher -EventName "Error" -Action {write-host "Error"}
    Register-ObjectEvent -InputObject $watcher -EventName "Renamed" -Action {write-host "Renamed"}    
}

try {
    Write-Host "Press any key to stop ..." -ForegroundColor Cyan

    $running=$true
    while($running) {
        if([console]::KeyAvailable) {
            $running=$false
        }
    }
} finally {
    # FIXME: this removes ALL registered events.
    #        usually there aren't any just hanging around but there could be
    Get-EventSubscriber | Unregister-Event
}


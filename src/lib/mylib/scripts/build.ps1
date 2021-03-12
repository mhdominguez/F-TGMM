$config="Release"
$build=(join-path  $PSScriptRoot "../build")
$logfile=join-path $build "log.txt"

function msbuild($ver) {
    if( $ver -ge 15 ) {
        $vswhere = "${env:ProgramFiles(x86)}\Microsoft Visual Studio\Installer\vswhere"
        $path = . $vswhere -latest -products * -requires Microsoft.Component.MSBuild -property installationPath
        join-path $path 'MSBuild\15.0\Bin\MSBuild.exe'        
    } else {
        $regKey = "HKLM:\software\Microsoft\MSBuild\ToolsVersions\$ver"
        $regProperty = "MSBuildToolsPath"
        join-path -path (Get-ItemProperty $regKey).$regProperty -childpath "msbuild.exe"
    }
}

if (-Not (Test-Path build)){
    $installdir=join-path (pwd) 'win64'
    mkdir build -ErrorAction SilentlyContinue
    cd build
    cmake -G "Visual Studio 12 2013 Win64" -DCMAKE_INSTALL_PREFIX="${installdir}" ..    
    cd ..
    $conf_ret_code=$LASTEXITCODE
    if($conf_ret_code -ne 0){
        throw "CMake failed with exit code $LASTEXITCODE."
    }
}

write-host "Building ..."
cd build
try {
    (&(msbuild("15.0")) INSTALL.vcxproj  /m /nologo /clp:"Verbosity=minimal;ShowTimestamp" /property:"Configuration=${config}")>>$logfile
    $build_ret_code=$LASTEXITCODE
} finally {
    cd ..
}

function color ($line) {
    if($line -match "error C\d+")       {"red"}
    elseif($line -match "error:")       {"red"}
    elseif($line -match "error LNK\d+") {"red"}
    elseif($line -match "warning C\d+") {"yellow"}
    elseif($line -match "warning:")     {"yellow"}
    else {"white"}
}
@(get-content $logfile) | %{write-host $_ -ForegroundColor (color $_)}

remove-item -Path $logfile -Force

if($build_ret_code -ne 0){
    throw "MSBUILD failed with exit code $LASTEXITCODE."
}
write-host "...Done"

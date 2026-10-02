# Start the SpatialViewer tile server (run from anywhere).
# Password: -Password, or $env:SPATIALVIZ_PASSWORD, or a .password file next to server.py (gitignored).
param([string]$Cache = "D:\SpatialVizCache", [int]$Port = 8760, [string]$Password = "")
$here = Split-Path -Parent $MyInvocation.MyCommand.Path
$args = @("$here\server.py", "--cache", $Cache, "--port", $Port)
if ($Password -ne "") { $args += @("--password", $Password) }
& "$here\.venv\Scripts\python.exe" @args

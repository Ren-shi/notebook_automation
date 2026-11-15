# Drive the GEMINI++ GUI (Gemini.exe from LISE++) for one input file, without a mouse.
#
# Gemini.exe has no silent mode: given a .gemini file on the command line it opens the file in its window and waits.
# This script presses Execute through Windows UI Automation, waits for the results window, saves the text of every
# results tab next to the input (<base>.tab<i>.txt), and closes the program. GEMINI++ itself writes <base>.html and
# <base>.cs4 next to the input when the run finishes. Called by scripts/gemini_reference.py, which starts this script
# on a separate, invisible Windows desktop, so that GEMINI++'s windows never show or take the focus (UI Automation
# works within the desktop the script runs on). Started from a normal console, the windows show for a few seconds.
#
#   powershell -NoProfile -ExecutionPolicy Bypass -File gemini_drive.ps1 -Exe <Gemini.exe> -InputFile <file.gemini>
param(
  [Parameter(Mandatory = $true)][string]$Exe,
  [Parameter(Mandatory = $true)][string]$InputFile,
  [string]$Tab = "",
  [int]$TimeoutSeconds = 3600,
  [int]$StallSeconds = 45,
  [string]$LogFile = ""  # result line or error, for callers that cannot read this script's output
)
$ErrorActionPreference = "Stop"
trap {
  if ($LogFile) { "ERROR: $_" | Out-File -Append -Encoding utf8 $LogFile }
  Write-Error "$_" -ErrorAction Continue
  exit 1
}
Add-Type -AssemblyName UIAutomationClient, UIAutomationTypes
$A = [System.Windows.Automation.AutomationElement]
$TS = [System.Windows.Automation.TreeScope]
$True_ = [System.Windows.Automation.Condition]::TrueCondition

$inputPath = (Resolve-Path $InputFile).Path
$dir = Split-Path $inputPath -Parent
$base = [System.IO.Path]::GetFileNameWithoutExtension($inputPath)

function Get-Windows($procId) {
  $pc = New-Object System.Windows.Automation.PropertyCondition($A::ProcessIdProperty, $procId)
  return $A::RootElement.FindAll($TS::Children, $pc)
}

$p = Start-Process -FilePath $Exe -ArgumentList "`"$inputPath`"" -WorkingDirectory $dir -PassThru
try {
  # wait for the main window
  $deadline = (Get-Date).AddSeconds(120)
  $main = $null
  while (-not $main -and (Get-Date) -lt $deadline) {
    Start-Sleep -Milliseconds 500
    foreach ($w in (Get-Windows $p.Id)) { if ($w.Current.AutomationId -eq "MainWindow") { $main = $w } }
  }
  if (-not $main) { throw "Gemini main window did not appear" }
  $title = $main.Current.Name
  if ($title -ne "Gemini - $base") { throw "Gemini did not open the input file (window title '$title')" }

  if ($Tab) {
    # which calculation runs follows the selected tab of the main window, not the file's mode= key
    $tc = New-Object System.Windows.Automation.AndCondition(
      (New-Object System.Windows.Automation.PropertyCondition($A::NameProperty, $Tab)),
      (New-Object System.Windows.Automation.PropertyCondition($A::ControlTypeProperty, [System.Windows.Automation.ControlType]::TabItem)))
    $ti = $main.FindFirst($TS::Descendants, $tc)
    if (-not $ti) { throw "tab '$Tab' not found" }
    $ti.GetCurrentPattern([System.Windows.Automation.SelectionItemPattern]::Pattern).Select()
    Start-Sleep -Milliseconds 500
  }
  $cond = New-Object System.Windows.Automation.PropertyCondition($A::NameProperty, "Execute")
  $btn = $main.FindFirst($TS::Descendants, $cond)
  if (-not $btn) { throw "Execute button not found" }
  $html = Join-Path $dir "$base.html"
  if (Test-Path $html) { Remove-Item $html }
  $t0 = Get-Date
  try { $btn.GetCurrentPattern([System.Windows.Automation.InvokePattern]::Pattern).Invoke() } catch { }

  # Wait for GEMINI++ to write <base>.html, which it does when the run is done. Do not touch the window through UI
  # Automation meanwhile: querying it while the calculation runs crashed Gemini.exe (access violations in Gemini.exe
  # and Qt6Core.dll) or left it hung behind its "Calculating.." dialog.
  # A run that stops using CPU without writing results has hung (seen for 16O + 144Sm, which later crashes).
  $deadline = (Get-Date).AddSeconds($TimeoutSeconds)
  $lastCpu = -1.0; $stalledSince = Get-Date
  while (-not (Test-Path $html)) {
    Start-Sleep -Seconds 2
    if ($p.HasExited) { throw "Gemini.exe exited during the run (exit code $($p.ExitCode))" }
    if ((Get-Date) -gt $deadline) { throw "no results after $TimeoutSeconds s" }
    $p.Refresh(); $cpu = $p.TotalProcessorTime.TotalSeconds
    if ($cpu -gt $lastCpu + 0.5) { $lastCpu = $cpu; $stalledSince = Get-Date }
    elseif (((Get-Date) - $stalledSince).TotalSeconds -gt $StallSeconds) {
      throw "Gemini.exe hung: no CPU use for $StallSeconds s and no results (CPU $([math]::Round($cpu, 1)) s)"
    }
  }
  $elapsed = ((Get-Date) - $t0).TotalSeconds
  Start-Sleep -Seconds 3
  $res = $null
  $deadline = (Get-Date).AddSeconds(60)
  while (-not $res -and (Get-Date) -lt $deadline) {
    foreach ($w in (Get-Windows $p.Id)) { if ($w.Current.ClassName -eq "Result_Widget") { $res = $w } }
    if (-not $res) { Start-Sleep -Seconds 2 }
  }
  if (-not $res) { throw "results written but no results window" }

  $tabCond = New-Object System.Windows.Automation.PropertyCondition($A::ControlTypeProperty, [System.Windows.Automation.ControlType]::TabItem)
  $tabs = $res.FindAll($TS::Descendants, $tabCond)
  $i = 0
  foreach ($t in $tabs) {
    $t.GetCurrentPattern([System.Windows.Automation.SelectionItemPattern]::Pattern).Select()
    Start-Sleep -Milliseconds 800
    $texts = @()
    foreach ($e in $res.FindAll($TS::Descendants, $True_)) {
      try {
        $tp = $e.GetCurrentPattern([System.Windows.Automation.TextPattern]::Pattern)
        if (-not $e.Current.IsOffscreen) { $texts += $tp.DocumentRange.GetText(-1) }
      } catch { }
    }
    [System.IO.File]::WriteAllText((Join-Path $dir "$base.tab$i.txt"), ($texts -join "`n"), [System.Text.Encoding]::UTF8)
    $i++
  }
  $msg = "OK $base elapsed_s=$([math]::Round($elapsed, 1)) tabs=$i"
  if ($LogFile) { $msg | Out-File -Append -Encoding utf8 $LogFile }
  $msg
}
finally {
  if (-not $p.HasExited) { Stop-Process -Id $p.Id -Force }
}

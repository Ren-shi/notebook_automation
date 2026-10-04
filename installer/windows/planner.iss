; Inno Setup script for the physim experiment planner on Windows (backlog item 44).
;
; Packs the folder that installer/windows/build.py assembled (default build\planner) into one setup program:
;   iscc /DAppVersion=0.1.0 /DBundle=..\..\build\planner installer\windows\planner.iss
;
; It installs for the current user only (no administrator rights, no UAC prompt) into
; %LOCALAPPDATA%\Programs\physim planner, adds a Start-menu shortcut (and optionally one on the desktop) that runs
; "pythonw.exe -m physim app --desktop" from the bundled Python, and registers an uninstaller. Installing a newer
; version over an older one replaces it.

#ifndef AppVersion
  #define AppVersion "0.0.0"
#endif
#ifndef Bundle
  #define Bundle "..\..\build\planner"
#endif
#ifndef OutputDir
  #define OutputDir "..\..\dist"
#endif

[Setup]
AppId={{6C1F4E0B-8B1D-4C55-9C2E-8E3A1F5B7D44}
AppName=physim planner
AppVersion={#AppVersion}
AppVerName=physim planner {#AppVersion}
AppPublisher=physim
AppPublisherURL=https://github.com/Ren-shi/notebook_automation
DefaultDirName={autopf}\physim planner
DefaultGroupName=physim
DisableProgramGroupPage=yes
PrivilegesRequired=lowest
ArchitecturesAllowed=x64compatible
ArchitecturesInstallIn64BitMode=x64compatible
MinVersion=10.0
OutputDir={#OutputDir}
OutputBaseFilename=physim-planner-{#AppVersion}-windows-x64-setup
SetupIconFile={#Bundle}\planner.ico
UninstallDisplayIcon={app}\planner.ico
LicenseFile={#Bundle}\LICENSE.txt
Compression=lzma2/max
SolidCompression=yes
WizardStyle=modern
CloseApplications=yes

[Tasks]
Name: "desktopicon"; Description: "{cm:CreateDesktopIcon}"; GroupDescription: "{cm:AdditionalIcons}"

[Files]
Source: "{#Bundle}\*"; DestDir: "{app}"; Flags: ignoreversion recursesubdirs createallsubdirs

[InstallDelete]
; An older version's packages would otherwise stay next to the new ones.
Type: filesandordirs; Name: "{app}\python"

[UninstallDelete]
; Python writes bytecode caches next to modules it compiles at run time.
Type: filesandordirs; Name: "{app}\python"

[Icons]
Name: "{group}\physim planner"; Filename: "{app}\python\pythonw.exe"; Parameters: "-m physim app --desktop"; \
  WorkingDir: "{userdocs}"; IconFilename: "{app}\planner.ico"; \
  Comment: "Plan a nuclear-physics experiment in your web browser"
Name: "{autodesktop}\physim planner"; Filename: "{app}\python\pythonw.exe"; Parameters: "-m physim app --desktop"; \
  WorkingDir: "{userdocs}"; IconFilename: "{app}\planner.ico"; Tasks: desktopicon

[Run]
Filename: "{app}\python\pythonw.exe"; Parameters: "-m physim app --desktop"; WorkingDir: "{userdocs}"; \
  Description: "Start the physim planner now"; Flags: postinstall nowait skipifsilent

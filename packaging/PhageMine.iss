#define MyAppName "PhageMine"
#ifndef MyAppVersion
  #error MyAppVersion must be supplied from phagemine.__version__
#endif

[Setup]
AppId={{9C8E0934-99D5-48D2-853D-F4B8A56C1339}
AppName={#MyAppName}
AppVersion={#MyAppVersion}
VersionInfoVersion=1.1.0.1
VersionInfoDescription=PhageMine Desktop {#MyAppVersion} installer
DefaultDirName={autopf}\PhageMine
DefaultGroupName=PhageMine
OutputDir=..\installer-dist
OutputBaseFilename=PhageMine-Windows-Setup
Compression=lzma2
SolidCompression=yes
ArchitecturesAllowed=x64compatible
ArchitecturesInstallIn64BitMode=x64compatible
UninstallDisplayIcon={app}\PhageMine.exe
WizardStyle=modern

[Files]
Source: "..\dist\PhageMine\*"; DestDir: "{app}"; Flags: ignoreversion recursesubdirs createallsubdirs

[Icons]
Name: "{group}\PhageMine"; Filename: "{app}\PhageMine.exe"
Name: "{autodesktop}\PhageMine"; Filename: "{app}\PhageMine.exe"; Tasks: desktopicon

[Tasks]
Name: "desktopicon"; Description: "Create a desktop shortcut"; GroupDescription: "Additional icons:"

[Run]
Filename: "{app}\PhageMine.exe"; Description: "Launch PhageMine"; Flags: nowait postinstall skipifsilent

schema = {
	"packages": [{
		'_vsInfo': Optional({
			'ignoreChanges': True
		}),
		'breadcrumbTemplate': Optional({
			'templates': [{
				'projectSortOrder': 2 | 900,
				'selects': [
					'Component.Unreal' |
					'Microsoft.VisualStudio.ComponentGroup.UWP.VC' |
					'Microsoft.VisualStudio.ComponentGroup.UWP.Xamarin' |
					'Component.Cocos'],
				'id': (
					'Microsoft.VisualStudio.Breadcrumb.Graphics.GetCocos2d' |
					'Microsoft.VisualStudio.Breadcrumb.Universal.GetVCSupport' |
					'Microsoft.VisualStudio.Breadcrumb.Graphics.GetUnreal' |
					'Microsoft.VisualStudio.Breadcrumb.Xamarin.UWPSupport'),
				'projectSubTypeSortOrder': 4 | None,
				'projectType': 'CSharp' | 'VC' | 'Game',
				'projectSubType': 'Windows Universal' | None
			}],
			'localizedResources': [{
				'description': str,
				'language': 'fr-fr' | 'en-us' | 'ja-jp' | 'ko-kr' | 'ru-ru' | 'cs-cz' | 'zh-cn' | 'de-de' | 'it-it' | 'tr-tr' | 'pl-pl' | 'zh-tw' | 'pt-br' | 'es-es',
				'projectTypeDisplayName': 'Visual C#' | 'Oyun' | 'Game' | 'Gra' | 'Jogo' | 'Hra' | 'Jeu' | 'ゲーム' | 'Gioco' | 'Visual C++' | 'Игра' | 'Juego' | 'Spiel' | '遊戲' | '游戏' | '게임',
				'projectSubTypeDisplayName': 'Universelles Windows' | 'Windows Evrensel' | 'Universale di Windows' | 'Windows Universel' | 'Windows 通用' | 'Windows ユニバーサル' | 'Windows universel' | 'Windows Universal' | 'Windows 유니버설' | 'Windows universale' | 'Universal do Windows' | 'Universal de Windows' | 'Projekt uniwersalny systemu Windows' | 'Универсальные приложения для Windows' | 'Univerzální aplikace pro Windows' | 'Универсальное приложение для Windows' | 'Platforma uniwersalna systemu Windows' | None,
				'templateId': 'Microsoft.VisualStudio.Breadcrumb.Graphics.GetCocos2d' | 'Microsoft.VisualStudio.Breadcrumb.Universal.GetVCSupport' | 'Microsoft.VisualStudio.Breadcrumb.Graphics.GetUnreal' | 'Microsoft.VisualStudio.Breadcrumb.Xamarin.UWPSupport',
				'title': str
			}]
		}),
		'chip': 'x86' | 'neutral' | 'x64' | 'arm64' | None,
		'defaultInstallDirectory': (
			'[ProgramFilesFolder]Microsoft Visual Studio\\2022\\BuildTools' |
			'[ProgramFiles64]Microsoft Visual Studio\\2022\\Professional' |
			'[ProgramFiles64]Microsoft Visual Studio\\2022\\TeamExplorer' |
			'[ProgramFilesFolder]Microsoft Visual Studio\\2022\\TestController' |
			'[ProgramFilesFolder]Microsoft Visual Studio\\2022\\TestAgent' |
			'[ProgramFiles64]Microsoft Visual Studio\\2022\\Enterprise' |
			'[ProgramFiles64]Microsoft Visual Studio\\2022\\Community' |
			None),
		'defaultProgram': Optional({
			'descriptionPosition': -1005 | -1004 | -1002 | -1000,
			'id': 'VisualStudio.[InstanceId]',
			'name': '' | 'Microsoft Visual Studio 2022',
			'registrationPath': 'Software\\Microsoft\\VisualStudio_[InstanceId]\\Capabilities',
			'descriptionPath': '[InstallDir]\\Common7\\IDE\\devenvdesc.dll'
		}),
		'dependencies': Optional({
			str: '' | '{version}' | '[{version},{version})' | {
				'chip': 'x64' | 'x86' | None,
				'id': str | None,
				'when': Optional(['Microsoft.VisualStudio.Product.Professional' | 'Microsoft.VisualStudio.Product.Enterprise' | 'Microsoft.Visualstudio.Product.WDExpress' | 'Microsoft.VisualStudio.Product.ProfessionalX86' | 'Microsoft.VisualStudio.Product.WDExpress' | 'Microsoft.VisualStudio.Product.TestController' | 'Microsoft.VisualStudio.Product.TestAgent' | 'Microsoft.VisualStudio.Product.EnterpriseX86' | 'Microsoft.VisualStudio.Product.BuildTools' | 'Microsoft.VisualStudio.Product.TeamExplorerX86' | 'Microsoft.VisualStudio.Product.SQL' | 'Microsoft.VisualStudio.Product.CommunityX86' | 'Microsoft.VisualStudio.Product.TeamExplorer' | 'Microsoft.VisualStudio.Product.Community']),
				'behaviors': 'ignoreApplicabilityFailures' | None,
				'language': 'de-DE' | 'en-US' | 'zh-CN' | 'ru-RU' | 'ko-KR' | 'zh-TW' | 'pt-BR' | 'pl-PL' | 'es-ES' | 'ja-JP' | 'fr-FR' | 'tr-TR' | 'it-IT' | 'cs-CZ' | None,
				'type': 'Recommended' | 'Optional' | None,
				'machineArch': 'x86' | 'x64' | 'ARM64' | 'arm64' | None,
				'version': '' | "{version}" | "[{version},{version})"
			}
		}),
		'detectConditions': {
			'expression': """an expression composed of the case insensitive operators "or", "and", "not", parenthetical expressions, and at least these keywords ['UnityRegKey', 'VCRedistUpgradeBuild', 'UwpManagedFeatureRegKeyNative', 'IncredibuildRegKey', 'Win10Cx86PendingRebootIdentityKey', 'Win10Bx86PackageIdentityLatestFx', 'ita', 'Win8amd64PendingRebootIdentityKey', 'UwpManagedFeatureRegKeyWoW', 'CocosDashboardRegKey', 'Win10Aamd64PackageIdentityLatestFx', 'PermanentIdentityKey', 'NetFx4ExtendedASPNET45Selected', 'Win81x86PackageIdentityLatestFx', 'Win8x86IdentityKey', 'Win81amd64PackageIdentityLatestFx', 'Win10Ax86IdentityKey', 'DotNetZDPMSPIdentityKey', 'VCRedistInstalled', 'MIERegKey', 'UWPNativeFeatureRegKeyNative', 'InstalledIdentityKey', 'Win81x86IdentityKey', 'csy', 'PendingRebootIdentityKey', 'chs', 'Win10Bx86IdentityKey', 'trk', 'Win81amd64PendingRebootIdentityKey', 'And', 'ITraceReloggerKey', 'Win8amd64IdentityKey', 'Win10Camd64IdentityKey', 'Win10Carm64IdentityKey', '140.VC.Installed', 'deu', 'VstoRKey', 'UnrealRegKey', 'Win10Aamd64PendingRebootIdentityKey', 'DontWantToInstall', 'IpOverUsbFeatureRegKeyWoW', 'Win10x86IdentityKey', 'DesktopArm64FeatureRegKeyNative', 'Win81amd64IdentityKey', 'FileConditionIdentityKey', 'Win10Bamd64PendingRebootIdentityKey', 'Win10Bamd64IdentityKey', 'Win10Ax86PackageIdentityLatestFx', 'Win10Aamd64IdentityKey', 'Win81x86PendingRebootIdentityKey', 'Win10Cx86IdentityKey', 'kor', 'Desktopx64FeatureRegKeyNative', 'WIFRegKey', 'BundleRegKeyWoW', 'Win10Ax86PendingRebootIdentityKey', 'Win10Camd64PendingRebootIdentityKey', 'DesktopArm64FeatureRegKeyWoW', 'Win8x86PackageIdentityLatestFx', 'jpn', 'UWPNativeFeatureRegKeyWoW', 'Win10amd64IdentityKey', 'esn', 'Win8x86PendingRebootIdentityKey', 'BundleRegKeyNative', 'Win10Bamd64PackageIdentityLatestFx', 'IpOverUsbFeatureRegKeyNative', 'VCRedistUpgradeMinor', 'plk', 'Win10Bx86PendingRebootIdentityKey', 'cht', 'VCRedistUpgradeMajor', 'Win8amd64PackageIdentityLatestFx', 'rus', 'enu', 'ptb', 'DotNetRegKey', 'Win10Carm64PendingRebootIdentityKey', 'CumulativeTargetingPackRegKey', 'Desktopx64FeatureRegKeyWoW', 'fra']""",
			'conditions': Optional([{
				'registryKey': str | None,
				'id': str("one of the keywords present in 'expression'") | None,
				'registryData': "{version}" | "[{version},)" | None,
				'fileVersionRange': '4.8.4375.0' | None,
				'registryType': 'Integer' | None,
				'registryValue': 'Major' | 'Bld' | 'OptionId.DesktopCPPx64' | 'version' | 'Version' | 'BundleVersion' | 'FullVersion' | 'Minor' | 'pv' | 'ThisVersionInstalled' | 'exe' | 'Selection' | 'Install Path' | 'OptionId.IpOverUsb' | 'OptionId.UWPCPP' | 'OptionId.UWPManaged' | 'release' | 'OptionId.DesktopCPPARM64' | 'InstallPath' | 'CurrentState' | None,
				'filePath': str | None,
				'productCode': str,
				'join': 'productCode' | 'or'
			}])
		},
		'extensionDir': str | None,
		'fileAssociations': Optional([{
			'perceivedType': 'text' | 'string' | 'image' | None,
			'defaultProgramRegistrationPath': 'Software\\Microsoft\\VisualStudio_[InstanceId]\\Capabilities' | None,
			'contentType': 'Application/xml' | 'application/xml' | 'string' | 'application/xml-dtd' | 'text/plain' | 'image/bmp' | 'application/javascript' | 'text/xml' | None,
			'extension': str,
			'isIconOnly': True | None,
			'progId': str
		}]),
		'folderMappings': Optional(['$ReferenceAssemblies' | '$VCTargets' | '$Schemas' | '$MSBuild' | '$RemoteDebugger' | '$Licenses' | '$PublicAssemblies']),
		'icon': Optional({
			'mimeType': 'image/svg+xml',
			'base64': str
		}),
		'iconPath': '[InstallDir]\\Common7\\IDE\\devenv.ico' | None,
		'id': str,
		'initInstallParams': Optional({
			'parameters': '-Operation install -InstallationID [InstanceId] -InstallationName [InstallationName] -InstallationVersion [InstallationVersion] -InstallationWorkloads [InstallationWorkloads] -InstallationPackages [InstallationPackages] -InstallationPath """[InstallDir]""" -ComponentId [ComponentId] -ChannelsPath """[ChannelsPath]""" -SetupEngineFilePath """[SetupEngineFilePath]""" -Log """[logfile]"""',
			'fileName': '[InstallDir]\\Common7\\IDE\\VSInitializer.exe'
		}),
		'initRepairParams': Optional({
			'parameters': '-Operation repair -InstallationID [InstanceId] -InstallationName [InstallationName] -InstallationVersion [InstallationVersion] -InstallationWorkloads [InstallationWorkloads] -InstallationPackages [InstallationPackages] -InstallationPath """[InstallDir]""" -ComponentId [ComponentId] -ChannelsPath """[ChannelsPath]""" -SetupEngineFilePath """[SetupEngineFilePath]""" -Log """[logfile]"""',
			'fileName': '[InstallDir]\\Common7\\IDE\\VSInitializer.exe'
		}),
		'initUninstallParams': Optional({
			'parameters': '-Operation uninstall -InstallationID [InstanceId] -InstallationName [InstallationName] -InstallationVersion [InstallationVersion] -InstallationWorkloads [InstallationWorkloads] -InstallationPackages [InstallationPackages] -InstallationPath """[InstallDir]""" -ComponentId [ComponentId] -ChannelsPath """[ChannelsPath]""" -SetupEngineFilePath """[SetupEngineFilePath]""" -Log """[logfile]"""',
			'fileName': '[InstallDir]\\Common7\\IDE\\VSInitializer.exe'
		}),
		'installParams': Optional({
			'parameters': str | None,
			'fileName': 'WinSdkInstaller.exe' | 'MicrosoftEdgeWebView2RuntimeInstallerX64.exe' | 'MicrosoftEdgeWebView2RuntimeInstallerX86.exe' | '[Payload]' | '[InstallDir]\\Common7\\IDE\\TestToolsFinalizer.exe' | 'MicrosoftEdgeWebView2RuntimeInstallerARM64.exe' | '[SystemFolder]\\WindowsPowerShell\\v1.0\\powershell.exe' | '[InstallDir]\\Common7\\IDE\\VSFinalizer.exe'
		}),
		'installSizes': Optional({
			'win10KitsInstallDir': int | None,
			'targetDrive': int | None,
			'sharedDrive': int | None,
			'systemDrive': int | None
		}),
		'isUiGroup': True | None,
		'language': 'de-DE' | 'en-US' | 'neutral' | 'zh-CN' | 'ru-RU' | 'ko-KR' | 'zh-TW' | 'pt-BR' | 'pl-PL' | 'es-ES' | 'ja-JP' | 'fr-FR' | 'tr-TR' | 'it-IT' | 'cs-CZ' | None,
		'launchParams': Optional({
			'fileName': 'Common7\\IDE\\devenv.exe' | 'Common7\\IDE\\TestAgentConfigUI.exe' | 'Common7\\Tools\\LaunchDevCmd.bat' | 'Common7\\IDE\\TestControllerConfigUI.exe'
		}),
		'license': (
			'https://incredibuild.com/eula.html' |
			'https://go.microsoft.com/fwlink/?LinkID=733796' |
			'https://go.microsoft.com/fwlink/?LinkID=2180154' |
			'https://go.microsoft.com/fwlink/?linkid=823951' |
			'https://github.com/aspnet/aspnetcore' |
			'https://go.microsoft.com/fwlink/?LinkId=691992' |
			'https://go.microsoft.com/fwlink/?LinkID=616920' |
			'https://go.microsoft.com/fwlink/?LinkId=708610' |
			'https://go.microsoft.com/fwlink/?linkid=859432' |
			'https://www.microsoft.com/web/webpi/eula/Azure_Data_Lake_Tools_for_VS.htm' |
			'https://github.com/dotnet/core-sdk' |
			'https://github.com/dotnet/core-setup' |
			'https://go.microsoft.com/fwlink/?linkid=617019' |
			'https://go.microsoft.com/fwlink/?LinkID=627227' |
			'https://go.microsoft.com/fwlink/?LinkID=832025&clcid=0x409' |
			'https://go.microsoft.com/fwlink/?LinkID=617016&clcid=0x409' |
			'https://go.microsoft.com/fwlink/?LinkID=691682' |
			None),
		'localizedResources': Optional([{
			'longDescription': str | None,
			'description': str,
			'license': 'https://go.microsoft.com/fwlink/?LinkId=2179811' | 'https://go.microsoft.com/fwlink/?LinkId=2179911' | 'https://go.microsoft.com/fwlink/?LinkId=2180117' | 'https://go.microsoft.com/fwlink/?LinkId=2229259' | None,
			'language': 'ja-jp' | 'ru-ru' | 'pt-BR' | 'zh-cn' | 'es-ES' | 'it-it' | 'ja-JP' | 'zh-tw' | 'en-us' | 'ko-kr' | 'zh-CN' | 'ko-KR' | 'cs-cz' | 'de-de' | 'tr-TR' | 'it-IT' | 'cs-CZ' | 'es-es' | 'de-DE' | 'zh-TW' | 'ru-RU' | 'pl-PL' | 'pl-pl' | 'fr-fr' | 'en-US' | 'tr-tr' | 'fr-FR' | 'pt-br',
			'category': str | None,
			'keywords': Optional([str]),
			'title': str
		}]),
		'logFiles': Optional([{
			'pattern': str | None,
			'directory': '[Windows]\\logs\\CBS' | 'StandaloneSDK\\UnionWinmdWorkingFolder\\logs' | '[ProgramData]\\Microsoft\\EdgeUpdate\\Log' | 'windowssdk' | '[Windows]\\logs\\DISM' | '[Windows]\\logs\\cbs' | '[Windows]\\temp' | None
		}]),
		'machineArch': 'x86' | 'x64' | 'ARM64' | 'arm64' | None,
		'msiProperties': Optional({
			'EXTUI': '1' | None,
			'FEEDBACKOPTIN': '1' | None,
			'INSTALLDIR': '[ProgramFilesOrSharedDrive]\\Epic Games' | '[SharedInstallDir]' | '[SharedInstallDir]\\Entity Framework Tools' | None,
			'DOTNETHOME': '[ProgramFilesX64]dotnet' | '[ProgramFilesX86]dotnet' | None,
			'ALLUSERS': '1' | None,
			'MSIFASTINSTALL': '7' | None,
			'APPLICATIONFOLDER': '[SharedInstallDir]\\SDKs\\Azure' | None,
			'SKIPPENDINGREBOOTCHECK': '1' | None,
			'IACCEPTMSSQLCMDLNUTILSLICENSETERMS': 'YES' | None,
			'INSTALLLEVEL': '2' | None,
			'IACCEPTSQLLOCALDBLICENSETERMS': 'YES' | None,
			'VSVERSION': '17.0.33626.349' | '17.0.33605.316' | None,
			'ALLOWMSIINSTALL': 'True' | None,
			'ALLOWMSIUNINSTALL': 'True' | None,
			'VSINSTALLER': '1' | None,
			'USING_EXUIH': '1' | None,
			'IACCEPTMSODBCSQLLICENSETERMS': 'YES' | None,
			'VS7.71E8CD30_5AB3_4754_907A_EF452D9A12F5': '[CustomInstallPath]' | None,
			'VSEXTUI': '1' | None,
			'VS7.3643236F_FC70_11D3_A536_0090278A1BB8': '[SharedInstallDir]' | '[SharedInstallDrive]\\Program Files (x86)\\Microsoft Visual Studio 14.0\\' | '[CustomInstallPath]' | None,
			'TARGETDIR': '[SharedInstallDir]' | None,
			'ADDLOCAL': 'ALL' | None,
			'ISVSINSTALL': '1' | None,
			'INSTALLEDBY': 'VisualStudio' | None,
			'VSFOLDER': '[SharedInstallDir]' | None,
			'PIDKEY': 'NGKBDRWKQFTT82MTRMPKRM6XM' | None
		}),
		'name': str,
		'nonCriticalProcesses': Optional([
			'IntelliTrace' |
			'VBCSCompiler' |
			'XDesProc' |
			'vcpkgsrv' |
			'mspdbsrv' |
			'VcxprojReader'|
			'vsls-agent' |
			'git' |
			'Microsoft.Alm.Shared.Remoting.RemoteContainer.dll' |
			'Microsoft.ServiceHub.Controller' |
			'vctip' |
			'ScriptedSandbox64.exe' |
			'ScriptedSandbox32.exe']),
		'nuGetPackageId': 'Microsoft.Windows.SDK.BuildTools' | None,
		'nuGetPackageVersion': '10.0.22621.756' | None,
		'outOfSupport': True | None,
		'payloads': Optional([{
			'cache': False | None,
			'fileName': str,
			'url': str,
			'isDynamicEndpoint': True | None,
			'sha256': str,
			'size': int,
			'signer': Optional({
				"$ref": '1' - '15'
			})
		}]),
		'permanent': True | None,
		'productArch': 'x86' | 'neutral' | 'x64' | 'ARM64' | 'Neutral' | 'Arm64' | 'arm64' | None,
		'productCode': str | None,
		'productLanguage': int | None,
		'productVersion': str or None,
		'progIds': Optional([{
			'ddeTopic': 'system' | None,
			'id': str,
			'defaultIconPath': str | None,
			'arguments': '/openuri' | '"%1"' | '/SnapshotDebugger' | '/dde "%1"' | '/openTfsLinkAsIs' | None,
			'dde': True | None,
			'defaultIconPosition': int | None,
			'ddeApplication': 'VsGraphics.17.0' | 'VisualStudio.15.0' | 'VisualStudio.17.0' | 'VisualStudio.17.6' | None,
			'clsid': '{DDD49901-1DD4-4288-A835-F9C413C4EF1F}' | None,
			'appUserModelId': 'Blend.[InstanceId]' | 'VisualStudio.[InstanceId]' | None,
			'path': '[ProgramFiles(x86)]\\Common Files\\Microsoft Shared\\MSEnv\\VSLauncher.exe' | '[InstallDir]\\Common7\\IDE\\VSWebLauncher.exe' | '[InstallDir]\\Common7\\IDE\\blend.exe' | '[InstallDir]\\Common7\\IDE\\devenv.exe' | '[InstallDir]\\Common7\\IDE\\vsga.exe' | None,
			'alwaysShowExtension': True | None,
			'iconHandler': '{9A2B23E4-2A50-48DB-B3C3-F5EA12947CB8}' | None,
			'displayName': str | None
		}]),
		'projectClassifiers': Optional([{
			'matcherId': 'MSBuild.Generic' | None,
			'selects': Optional(['Microsoft.VisualStudio.Component.Windows10SDK.18362' | 'Component.Xamarin' | 'Microsoft.VisualStudio.Component.Windows11SDK.22621' | 'Microsoft.VisualStudio.ComponentGroup.UWP.VC' | 'Microsoft.VisualStudio.Web.Mvc4.ComponentGroup' | 'Microsoft.VisualStudio.Component.VC.ATL' | 'Component.MDD.IOS' | 'Microsoft.VisualStudio.Component.Windows10SDK.19041' | 'Microsoft.VisualStudio.Component.VC.ATLMFC']),
			'id': str,
			'matcherData': Optional([{
				'regExMatchSource': str,
				'type': 'Import' | 'Item' | 'Property',
				'capabilityType': 'CSharpSharedTargetImports' | 'USQLTargetImport' | 'WindowsTargetPlatformVersion' | 'ApplicationType' | 'MicrosoftCSharpReference' | 'MicrosoftAspNetCoreReference' | 'UseOfAtl' | 'WcfConfigValidationEnabled' | 'TargetPlatformVersion' | 'TypeGuid' | 'UseOfMFC' | 'HiveTargetImport' | 'PigTargetImport' | 'SharedCommonDefaultProps' | 'SharedCommonProps' | 'OutputType' | 'TargetPlatformIdentifier' | 'VisualBasicTargetImports' | 'VisualBasicSharingTargetImports' | 'ApplicationTypeRevision' | 'Import',
				'projectPropertyId': 'TargetPlatformVersion' | 'OutputType' | 'Reference' | 'PackageReference' | 'UseOfMFC' | 'TargetPlatformIdentifier' | 'ProjectTypeGuids' | 'WindowsTargetPlatformVersion' | 'ApplicationTypeRevision' | 'UseOfAtl' | 'WcfConfigValidationEnabled' | 'ProjectReference' | 'ApplicationType' | None
			}]),
			'factoryGuid': '{8BC9CEB8-8B4A-11D0-8D11-00A0C91BC942}' | '{F2A71F9B-5D33-465A-A702-920D77279786}' | '{F184B08F-C81C-45F6-A57F-5ABD9991F28F}' | '{416D63FD-0477-49AA-A954-A7C5B95A9B51}' | '{39E2626F-3545-4960-A6E8-258AD8476CE5}' | '888888a0-9f3d-457c-b088-3a5042f75d52' | '{4D4E14FB-86F2-46A5-8BFB-41569A68D9E8}' | '{E24C65DC-7377-472b-9ABA-BC803B73C61A}' | '{81BE5007-4CAB-4670-AEEB-3DA6572B9D12}' | '{12D26815-A329-4BD9-922E-4A35B783B8AA}' | '{778DAE3C-4631-46EA-AA77-85C1314464D9}' | '{C7167F0D-BC9F-4E6E-AFE1-012C56B48DB5}' | '{182E2583-ECAD-465B-BB50-91101D7C24CE}' | '{9A19103F-16F7-4668-BE54-9A1E7A4F7556}' | '{FAE04EC0-301F-11D3-BF4B-00C04F79EFBC}' | '{fae04ec0-301f-11d3-bf4b-00c04f79efbc}' | '{79875B75-C0D1-4DF0-91FB-7D5CF7C909B0}' | '{cc5fd16D-436d-48ad-a40c-5a424c6e3e79}',
			'extension': 'pigproj' | 'wapproj' | 'fsproj' | 'shproj' | 'webproj' | 'ccproj' | 'csproj' | 'androidproj' | 'hiveproj' | 'vcxproj' | 'vbproj' | 'usqlproj' | 'pyproj',
			'appliesTo': 'OutputType+VisualBasicTargetImports' | 'CSharpSharedTargetImports+SharedCommonProps+SharedCommonDefaultProps' | 'USQLTargetImport' | 'ApplicationType+ApplicationTypeRevision+WindowsTargetPlatformVersion' | 'ApplicationType' | 'OutputType+CSharpSharedTargetImports+MicrosoftCSharpReference' | 'TypeGuid+WcfConfigValidationEnabled' | 'MicrosoftAspNetCoreReference' | 'SharedCommonDefaultProps+SharedCommonProps+VisualBasicSharingTargetImports' | 'UseOfAtl' | 'TargetPlatformVersion' | 'TypeGuid' | 'UseOfMFC' | 'HiveTargetImport' | 'PigTargetImport' | 'TypeGuid+TargetPlatformIdentifier+TargetPlatformVersion' | 'Import' | 'ApplicationType+ApplicationTypeRevision' | 'TypeGuid+TargetPlatformIdentifier' | None,
			'priority': 10 | 100 | 101 | 200 | 204 | 300 | 400 | 500 | 600 | 1000 | None
		}]),
		'propertyInitializers': Optional([{
			'value': 'KitsRoot10',
			'key': 'HKEY_LOCAL_MACHINE\\SOFTWARE\\Microsoft\\Windows Kits\\Installed Roots',
			'defaultValue': '[ProgramFilesOrSharedDrive]\\Windows Kits\\10',
			'property': 'Win10KitsInstallDir'
		}]),
		'providerKey': str | None,
		'recommendSelection': True | None,
		'relatedProcessDirectories': Optional([str]),
		'relatedProcessFiles': Optional([str]),
		'relatedServices': 'VSStandardCollectorService150' | None,
		'relativePath': str | None,
		'releaseNotes': 'https://docs.microsoft.com/en-us/visualstudio/releases/2022/release-notes-v17.6#17.6.2' | None,
		'repairParams': Optional({
			'parameters': str | None,
			'fileName': '[Payload]' | '[InstallDir]\\Common7\\IDE\\TestToolsFinalizer.exe' | '[SystemFolder]\\WindowsPowerShell\\v1.0\\powershell.exe' | '[InstallDir]\\Common7\\IDE\\VSFinalizer.exe'
		}),
		'requirements': Optional({
			'functors': Optional({
				'architecture': 'arm64' | 'x64'
			}),
			'supportedOS': "{version}" | "[{version},{version})" | None,
			'conditions': Optional({
				'conditions': Optional([{
					'registryKey': "HKEY_LOCAL_MACHINE\\SYSTEM\\CurrentControlSet\\Control\\MUI\\UILanguages\\{lang}" | False,
					'chip': 'x86' | 'x64' | None,
					'id': 'IISW3SvcInstalledx86' | 'Server-Gui-Mgmt-Enabled' | 'OSLPRegKey' | 'ServerCore' | 'IIScoreWebEngineInstalledx64' | 'IsFunctionDiscoveryPnPProviderPresentWoW' | 'Win10ThresholdOneAndTwoBuildNumber' | 'IIScoreWebEngineInstalledx86' | 'IISW3SvcInstalledx64' | 'Win10ThresholdBuildNumber' | 'IsFunctionDiscoveryPnPProviderPresentNative' | 'Powershell5',
					'registryData': '1049' | 'Server Core' | '2052' | '3082' | '1036' | '[1.0]' | '1040' | '1031' | '1046' | '[10240.0,14393.0)' | '1041' | '1055' | '1' | '1028' | '1045' | '[5.0,)' | '1029' | '1042' | None,
					'registryValue': 'CurrentBuildNumber' | 'Server-Gui-Mgmt' | 'LCID' | 'CoreWebEngine' | 'InstallationType' | '00000000' | 'W3SVC' | 'PowershellVersion',
					'registryType': 'Integer' | None
				}]),
				'expression': """an expression composed of the case insensitive operators "or", "and", "not", parenthetical expressions, and at least these keywords ['ServerCore', 'Server-Gui-Mgmt-Enabled', 'OSLPRegKey', 'Win10ThresholdOneAndTwoBuildNumber', 'IsFunctionDiscoveryPnPProviderPresentNative', 'IsFunctionDiscoveryPnPProviderPresentWoW', 'Powershell5', ''IISCoreWebEngineInstalledx64', 'IISW3SvcInstalledx64', 'Win10ThresholdBuildNumber', 'IISCoreWebEngineInstalledx860', 'IISW3SvcInstalledx86']"""
			})
		}),
		'returnCodes': Optional(['5100' | '1638' | '-1' | '-2147219187' | '50' | '1' | '-2147219970' | '2' | '-2147219416']),
		'shortcuts': Optional([{
			'description': 'Open Visual Studio 2022 Tools Command Prompt for targeting x86 with ARM64-hosted tools' | 'Open Visual Studio 2022 Tools Command Prompt for targeting x64 with x86-hosted tools' | 'Open Visual Studio 2022 Tools Command Prompt' | 'Open Visual Studio 2022 Tools Command Prompt for targeting x86' | 'Open Visual Studio 2022 Tools Command Prompt for targeting x64 with ARM64-hosted tools' | 'Open Visual Studio 2022 Tools Command Prompt for targeting x64' | 'Open Visual Studio 2022 Tools Command Prompt for targeting ARM64' | 'Microsoft Visual Studio 2022' | 'Microsoft Blend for Visual Studio 2022' | 'Test Agent Configuration Tool' | 'Microsoft Visual Studio Debuggable Package Manager PowerShell Session' | 'Open Visual Studio 2022 Tools PowerShell' | 'Load Test Controller Configuration Tool' | 'Open Visual Studio 2022 Tools Command Prompt for targeting x86 with x64-hosted tools' | None,
			'folder': '[startmenu]\\programs\\Visual Studio 2022\\Visual Studio Tools' | '[startmenu]\\programs\\Visual Studio 2022\\Visual Studio Tools\\VC' | '[startmenu]\\programs' | '[startmenu]\\programs\\Visual Studio 2022\\Microsoft Visual Studio SDK\\Tools',
			'arguments': '/k "[installdir]\\VC\\Auxiliary\\Build\\vcvarsall.bat" arm64_x86' | '/rootSuffix Exp' | '/k "[installdir]\\VC\\Auxiliary\\Build\\vcvarsall.bat" arm64_x64' | '-noe -c "&{Import-Module """[installdir]\\Common7\\Tools\\Microsoft.VisualStudio.DevShell.dll"""; Enter-VsDevShell [InstanceId]}"' | '/k "[installdir]\\VC\\Auxiliary\\Build\\vcvarsamd64_x86.bat"' | '/k "[installdir]\\VC\\Auxiliary\\Build\\vcvars32.bat"' | '/k "[installdir]\\VC\\Auxiliary\\Build\\vcvarsx86_amd64.bat"' | '/k "[installdir]\\VC\\Auxiliary\\Build\\vcvars64.bat"' | '/k "[installdir]\\VC\\Auxiliary\\Build\\vcvarsall.bat" arm64' | '/k "[installdir]\\Common7\\Tools\\VsDevCmd.bat"' | '-NoExit -Command "& { Import-Module Appx; Import-Module .\\AppxDebug.dll; Show-AppxDebug}"' | '/C "[installdir]\\VSSDK\\VisualStudioIntegration\\Tools\\Bin\\CreateExpInstance.exe" /Reset /VSInstance=17.0_[InstanceId] /RootSuffix=Exp && PAUSE' | None,
			'targetPath': '[installdir]\\Common7\\IDE\\TestControllerConfigUI.exe' | '[installdir]\\Common7\\IDE\\blend.exe' | '%comspec%' | '[installdir]\\Common7\\IDE\\devenv.exe' | '[installdir]\\Common7\\IDE\\TestAgentConfigUI.exe' | '[systemfolder]WindowsPowerShell\\v1.0\\PowerShell.exe',
			'shellProperties': Optional({
				'System.AppUserModel.ID': 'Blend.[InstanceId]' | 'VisualStudio.[InstanceId]'
			}),
			'workingDirectory': '[installdir]\\Common7\\IDE\\' | '%comspec%' | '[installdir]\\' | '[installdir]\\Common7\\IDE\\Remote Debugger\\Appx\\' | None,
			'displayName': str 
		}]),
		'supportsDownloadThenInstall': True | None,
		'supportsDownloadThenUpdate': True | None,
		'supportsExtensions': bool | None,
		'thirdPartyNotices': 'https://go.microsoft.com/fwlink/?LinkId=661288' | None,
		'type': 'Exe' | 'Vsix' | 'Nupkg' | 'Msu' | 'Component' | 'Zip' | 'Msi' | 'Workload' | 'Group' | 'WindowsFeature' | 'Product',
		'ui': Optional({
			'isAdvertisedPackage': True
		}),
		'uninstallParams': Optional({
			'parameters': str,
			'fileName': 'WinSdkInstaller.exe' | '[Payload]' | '[InstallDir]\\Common7\\IDE\\TestToolsFinalizer.exe' | '[SystemFolder]\\WindowsPowerShell\\v1.0\\powershell.exe' | '[InstallDir]\\Common7\\IDE\\VSFinalizer.exe'
		}),
		'upgradeCode': str | None,
		'urlAssociations': Optional([{
			'protocol': 'vsweb' | 'git-client' | 'vstfs' | 'vssd',
			'progId': 'VisualStudio.git-client.[InstanceId]' | 'VisualStudio.vsweb.[InstanceId]' | 'VisualStudio.vstfs.[InstanceId]' | 'VisualStudio.vssd.[InstanceId]',
			'defaultProgramRegistrationPath': 'Software\\Microsoft\\VisualStudio_[InstanceId]\\Capabilities',
			'displayName': 'TFS Protocol Handler' | 'Git Protocol Handler' | 'Snapshot Debugger Protocol Handler' | 'Web Protocol Handler'
		}]),
		'version': str,
		'visible': True | None,
		'vsixId': str | None
	}]
}

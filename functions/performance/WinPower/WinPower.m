% Shutdown, reboot, sleep, stay awake
% This function controls the power functions of Windows.
%
% Status = WinPower(Command, Force)
% INPUT:
%   Command: String, not case-sensitive. On demand all applications are closed.
%      'poweroff':  Switch power off.
%      'reboot':    Reboot the machine.
%      'logoff':    Logoff the current user.
%      'shutdown':  Shut down the machine to a state, which allows the user to
%                   switch off the power securely by hand.
%      'sleep':     Let the machine fall asleep.
%                   The additional argument 'off' disables the ability to fall
%                   into sleep, 'on' re-enables sleeping again.
%      'hibernate': Write memory to disk and fall into deep sleep.
%      'lock':      Lock the machine, password is required for wake-up.
%      'RebootMatlab': Reboot the machine, restart Matlab when the user is
%                   logged in again, use the optional 2nd input as startup
%                   command. See EXAMPLES. The Force level can be defined as
%                   3rd input. Vista or higher is required.
%      'Monitor':   Set monitor status without stopping the processing.
%                   2nd input: 'off' (default), 'on', 'standby'.
%                   Moving the mouse etc. enables the monitor automatically.
%      'BatteryStatus': Reply the battery related parameters.
%   Force: String, not case-sensitive. Force applications to shut down. This
%      may cause a loss of data, if changed data cannot be saved to disk!
%      For 'sleep', 'hibernate' and 'lock' the Force is ignored.
%      'noforce':   Wait until the applications close voluntarily. Default.
%      'force':     Close waiting applications. Dangerous.
%      'forceifhung': Close waiting and crashed applications. Dangerous.
%
% OUTPUT:
%   Status: For the command 'BatteryStatus' a struct is replied:
%     .ACLine: Double value, 1: connected to AC power, 0: not connected,
%              NaN: unknown status.
%     .BatteryStatus:       String, 'high', 'low', 'critical', 'unknown' or
%                           'no_battery'. BatteryLifePercent is neater.
%     .BatteryCharging:     Logical flag, TRUE if the battery is charging.
%     .BatteryLifePercent:  Remaining life of the battery in percent.
%     .BatteryLifeTime:     Number of remaining seconds. NaN if unknown.
%     .BatteryFullLifeTime: Life time after a full charge. NaN if unknown.
%
% EXAMPLES:
% 1. Try a poweroff, do not kill waiting applications:
%      WinPower('poweroff');
% 2. Force a poweroff (safe open documents before!):
%      WinPower('poweroff', 'forceifhung');
% 3. Do not let the computer fall asleep during a long computation:
%      WinPower('Sleep', 'off');  Long_Calculation();  WinPower('Sleep', 'on');
%    Exiting Matlab restores the ability to sleep automatically.
% 4. Reboot the machine and after the user is logged in Matlab is started with
%    the bench() function:
%      WinPower('RebootMatlab', 'bench(2)');
%    This creates the startup argument for Matlab: -r "bench(2)"
%    The argument string is limited to 250 characters. Use \" to mask double
%    quotes inside the string.
% 5. Get the battery parameters:
%      Status = WinPower('BatteryStatus')
% 6. Switch off the monitor for 5 seconds:
%      WinPower('Monitor', 'off'); pause(5); WinPower('Monitor', 'on');
%
% NOTE: The change of the power status can fail although this function has been
%   executed successfully, when the user or another applications stops the
%   asynchronous process.
%
% WinPower - control power functions of Windows
% COMPILE:
%   The compilation of the C-file is started automatically when WinPower is
%   called the first time. For more details see: WinPower.c->COMPILE section.
%
% Other compilers: The case-insensitive string comparison may vary. Set the
% macro STRCMP_NI accordingly on demand: strnicmp, strncmpi, _strnicmp, ...
%
% Tested: Matlab/64 7.8, 7.13, 8.6, 9.1, Win7/64
%         Compiler: OWC1.8, MSVC2008/2010
% Assumed Compatibility: higher Matlab versions, Vista, Windows8.
% Author: Jan Simon, Heidelberg, (C) 2012-2017 matlab.2010(a)n(MINUS)simon.de

% $JRev: R-m V:012 Sum:LpUfWTizRBsl Date:03-May-2015 20:17:17 $
% $License: BSD (use/copy/change/redistribute on own risk, mention the author) $
% $File: Tools\GLMisc\WinPower.m $
% History:
% 001: 16-Apr-2012 13:04, First version.
% 007: 26-Jul-2012 11:11, Stable, supports OpenWatcom1.8 + Matlab6.5/32.
% 008: 13-Aug-2012 21:54, BatteryStatus.
% 011: 03-May-2015 12:56, Monitor on/off/standby.

% ==============================================================================
% This is just a dummy function to support the HELP command and to start an
% automatic compilation of the C-file. Alternatives:
% - Compile manually by "mex -O WinPower.c",
% - or call "InstallMex WinPower.c" directly,
% - or download the pre-copiled MEX: http:%www.n-simon.de/mex
% - Delete or comment the dummy function below, but keep the help section above!

function WinPower(varargin)  %#ok<VANUS>

warning(['JSimon:', mfilename, ':NoMex'], ...
   'The Mex file is compiled automatically.');

% No unit-test, no special parameters for the compilation:
InstallMex('WinPower.c');

% Now the call could be forwarded to the compiled MEX, but
% varargout = WinPower(varargin{:});

% return;

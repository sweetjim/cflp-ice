// WinPower.c
// Shutdown, reboot, sleep, disable sleep
// This function controls the power functions of Windows.
//
// Status = WinPower(Command, Force)
// INPUT:
//   Command: String, not case-sensitive. On demand all applications are closed.
//      'poweroff':  Switch power off.
//      'reboot':    Reboot the machine.
//      'logoff':    Logoff the current user.
//      'shutdown':  Shut down the machine to a state, which allows the user to
//                   switch off the power securely by hand.
//      'sleep':     Let the machine fall asleep.
//                   The additional argument 'off' disables the ability to fall
//                   into sleep, 'on' re-enables sleeping again.
//      'hibernate': Write memory to disk and fall into deep sleep.
//      'lock':      Lock the machine, password is required for wake-up.
//      'RebootMatlab': Reboot the machine, restart Matlab when the user is
//                   logged in again, use the optional 2nd input as startup
//                   command. See EXAMPLES. The Force level can be defined as
//                   3rd input. Vista or higher is required.
//      'Monitor':   Set monitor status without stopping the processing.
//                   2nd input: 'off' (default), 'on', 'standby'.
//                   Moving the mouse etc. enables the monitor automatically.
//      'Battery':   Reply the battery related parameters.
//   Force: String, not case-sensitive. Force applications to shut down. This
//      may cause a loss of data, if changed data cannot be saved to disk!
//      For 'sleep', 'hibernate' and 'lock' the Force is ignored.
//      'noforce':   Wait until the applications close voluntarily. Default.
//      'force':     Close waiting applications. Dangerous.
//      'forceifhung': Close waiting and crashed applications. Dangerous.
//
// OUTPUT:
//   Status: For the command 'Battery' a struct is replied:
//     .ACLine: Double value, 1: connected to AC power, 0: not connected,
//              NaN: unknown status.
//     .BatteryStatus:       String, 'high', 'low', 'critical', 'unknown' or
//                           'no_battery'. BatteryLifePercent is neater.
//     .BatteryCharging:     Logical flag, TRUE if the battery is charging.
//     .BatteryLifePercent:  Remaining life of the battery in percent.
//     .BatteryLifeTime:     Number of remaining seconds. NaN if unknown.
//     .BatteryFullLifeTime: Life time after a full charge. NaN if unknown.
//
// EXAMPLES:
// 1. Try a poweroff, do not kill waiting applications:
//      WinPower('poweroff');
// 2. Force a poweroff (safe open documents before!):
//      WinPower('poweroff', 'forceifhung');
// 3. Do not let the computer fall asleep during a long computation:
//      WinPower('Sleep', 'off');  Long_Calculation();  WinPower('Sleep', 'on');
//    Exiting Matlab restores the ability to sleep automatically.
// 4. Reboot the machine and after the user is logged in Matlab is started with
//    the bench() function:
//      WinPower('RebootMatlab', 'bench(2)');
//    This creates the startup argument for Matlab: -r "bench(2)"
//    The argument string is limited to 250 characters. Use \" to mask double
//    quotes inside the string.
// 5. Get the battery parameters:
//      Status = WinPower('Battery')
// 6. Switch off the monitor for 5 seconds:
//      WinPower('Monitor', 'off'); pause(5); WinPower('Monitor', 'on');
//
// NOTE: The change of the power status can fail although this function has been
//   executed successfully, when the user or another applications stops the
//   asynchronous process.
//
// COMPILE:
// This function must be compiled before using:
//   mex -O WinPower.c
// Calling the M-function WinPower tries to compile the function automatically.
// MSVC2008/2010 (Express or the corresponding SDKs) can compile the code, but
// BCC5.5 and the LCC shipped with Matlab 32bit fail. Small changes in the code
// allow to compile it with OpenWatcom1.8 for Matlab 6.5. Precompiled MEX files
// can be downloaded: http://www.n-simon.de/mex
//
// Other compilers: The case-insensitive string comparison may vary. Set the
// macro STRCMP_NI accordingly on demand: strnicmp, strncmpi, _strnicmp, ...
//
// Tested: Matlab 6.5, 7.7, 7.8, 7.13, WinXP/32, Win7/64
//         Compiler: OWC1.8, MSVC2008/2010
// Assumed Compatibility: higher Matlab versions, Vista, Windows8.
// Author: Jan Simon, Heidelberg, (C) 2012-2016 matlab.2010(a)n(MINUS)simon.de

/*
% $JRev: R-w V:022 Sum:Z1Xiyw0sgvAP Date:21-Dec-2017 23:53:53 $
% $License: BSD (use/copy/change/redistribute on own risk, mention the author) $
% $File: Tools\Mex\Source\WinPower.c $
% History:
% 001: 16-Apr-2012 13:04, First version.
% 005: 26-Jul-2012 11:11, Stable, supports OpenWatcom1.8 + Matlab6.5/32.
% 008: 13-Aug-2012 21:54, GetPowerStatus.
% 014: 21-Aug-2013 22:57, GetPowerStatus -> BatteryStatus.
% 019: 03-May-2015 12:20, Monitor Off command.
% 021: 18-Jun-2016 02:32, EWX_RESTARTAPPS defined on demand.
%      Thanks to Stuart Rogers.
% 022: 16-Dec-2017 19:43, Monitor commands work under Windows 10.
%      Thanks to Richard Walker.
*/

// Headers: --------------------------------------------------------------------
#if !defined(__WINDOWS__) && !defined(_WIN32) && !defined(_WIN64)
#  error Implemented for Windows only!
#endif

#if defined(__LCC__)
#  error LCC compiler cannot compile WinPower - see WinPower.c->COMPILE
#endif

// To compile with OpenWatcom1.8 under Windows XP:
#if defined(__WATCOMC__)
#  define _WIN32_WINNT 0x0501
#  define WINVER 0x0501
#endif

#pragma comment(lib, "user32.lib")
#pragma comment(lib, "advapi32.lib")
#pragma comment(lib, "PowrProf.lib")

#include <windows.h>
#include <winuser.h>
#include <PowrProf.h>  // Fails for LCC: ULONG causes troubles
#include <string.h>
#include <shellapi.h>
#include "mex.h"

// Compiler dependent settings: ------------------------------------------------
// Case-insensitive string comparison depends on the compiler:
//   strncmpi, strnicmp, _strnicmp, ...
#define STRCMP_NI _strnicmp

// Missing under __WATCOMC__ and MinGW:
#if !defined(EWX_RESTARTAPPS)  // Restart Application: Vista or higher
#  define EWX_RESTARTAPPS 0x00000040L
#endif
#if !defined(EWX_FORCEIFHUNG)
#  define EWX_FORCEIFHUNG 0x00000010L
#endif
        
// Define the default force level:
// 0:               Don't force applications to quit (recommended)
// EWX_FORCE:       Force applications to stop (unsaved data are lost)
// EWX_FORCEIFHUNG: Stop crashed application (unsaved data are lost)
#define DEFAULT_FORCE_LEVEL 0
                  
// Assume 32 bit addressing for Matlab 6.5:
// See MEX option "compatibleArrayDims" for MEX in Matlab >= 7.7.
#ifndef MWSIZE_MAX
#define mwSize  int32_T           // Defined in tmwtypes.h
#define mwIndex int32_T
#define MWSIZE_MAX MAX_int32_T
#endif

// Error messages do not contain the function name in Matlab 6.5! This is not
// necessary in Matlab 7, but it does not bother:
#define ERR_ID    "JSimon:WinPower:"
#define ERR_HEAD  "*** WinPower[mex]: "
#define WARN_HEAD "### WinPower[mex]: "
#define ERROR_2(id,msg)   mexErrMsgIdAndTxt(ERR_ID id, ERR_HEAD msg);
#define ERROR_3(id,msg,p) mexErrMsgIdAndTxt(ERR_ID id, ERR_HEAD msg, p);

#define Command_LEN  15  // Consider this number of characters only
#define ForceStr_LEN 15

// Types: ---------------------------------------------------------------------
// Dynamic pointer to modern function:
typedef HRESULT (WINAPI *REGAPP_FCN)(PCWSTR, DWORD);

// Enumarate different actions:
enum PowerAction{myVOID, myEXIT, mySLEEP, myHIBERNATE};

// GLobals: --------------------------------------------------------------------
// Flag: If TRUE, the sleep mode needs to be enabled:
static BOOL doCleanupSleep = FALSE;

// Prototypes: -----------------------------------------------------------------
void mySetPower(UINT Flag, enum PowerAction Action);
void myAllowSleep(BOOL Enable);
static void myCleanupSleep(void);
BOOL mySetMatlabArg(const mxArray *Arg);
mxArray *myBatteryStatus(void);

// Main function ===============================================================
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  char Command[Command_LEN + 1], ForceStr[ForceStr_LEN + 1];
  UINT Force = DEFAULT_FORCE_LEVEL;
  UINT Flag;
  int  ForceArgIndex = 1;     // 0-based!
  const mxArray *RebootArg;
  BOOL InputIsChar = TRUE, doSetPower = TRUE, canRestart;
  int  nOutput = 0;
  enum PowerAction Action = myEXIT;
  const int MONITOR_ON = -1, MONITOR_OFF = 2, MONITOR_STANDBY  = 1;
  int MonitorAction;
  
  // Check number and type of inputs:
  switch (nrhs) {
     case 3:  InputIsChar &= mxIsChar(prhs[2]);  // fall through!
     case 2:  InputIsChar &= mxIsChar(prhs[1]);  // fall through!
     case 1:  InputIsChar &= mxIsChar(prhs[0]);
              break;
     default: ERROR_2("BadNInput", "1 ro 3 inputs accepted.");
  }
  
  if (!InputIsChar) {
     ERROR_2("BadTypeInput", "Inputs must be a strings.");
  }
  
  // Set flags according to first input:
  mxGetString(prhs[0], Command, Command_LEN);
  if        (STRCMP_NI(Command, "poweroff", 8) == 0) {
     Flag = EWX_POWEROFF;
  } else if (STRCMP_NI(Command, "reboot",   6) == 0) {
     Flag = EWX_REBOOT;
  } else if (STRCMP_NI(Command, "logoff",   6) == 0) {
     Flag = EWX_LOGOFF;
  } else if (STRCMP_NI(Command, "shutdown", 8) == 0) {
     Flag = EWX_SHUTDOWN;
     
  } else if (STRCMP_NI(Command, "rebootmatlab", 12) == 0) {
     // Register startup parameters:
     ForceArgIndex = 2;  // [Force] would be defined as 3rd input:
     if (nrhs >= 2) {
        RebootArg = prhs[1];
     } else {
        RebootArg = NULL;
     }
     
     // Register the arguments for restarting Matlab:
     canRestart = mySetMatlabArg(RebootArg);
     
     // WinXP cannot restart applications with specified arguments:
     if (canRestart) {
        Flag = EWX_RESTARTAPPS;
     } else {
        Flag = EWX_REBOOT;
     }
     
  } else if (STRCMP_NI(Command, "sleep", 5) == 0) {
     if (nrhs == 1) {    // WinPower('sleep')
       Action = mySLEEP;
     } else {            // WinPower('sleep', 'off') or 'on':
       Action = myVOID;
       mxGetString(prhs[1], ForceStr, ForceStr_LEN);
       myAllowSleep((BOOL) (STRCMP_NI(ForceStr, "on", 3) == 0));
     }
     
  } else if (STRCMP_NI(Command, "hibernate", 9) == 0) {
     Action = myHIBERNATE;
     
  } else if (STRCMP_NI(Command, "lock", 4) == 0) {
     LockWorkStation();
     Action = myVOID;
     
  } else if (STRCMP_NI(Command, "battery", 7) == 0) {
     plhs[0] = myBatteryStatus();
     Action  = myVOID;
     nOutput = 1;

  } else if (STRCMP_NI(Command, "monitor", 7) == 0) {
     Action = myVOID;       // Nothing to do for power mode
          
     MonitorAction = MONITOR_OFF;  // Default Monitor action
     if (nrhs == 2) {              // WinPower('monitor', 'off'/'on'):
       mxGetString(prhs[1], ForceStr, ForceStr_LEN);
       if        (STRCMP_NI(ForceStr, "on", 3) == 0) {
          MonitorAction = MONITOR_ON;
       } else if (STRCMP_NI(ForceStr, "standby", 7) == 0) {
          MonitorAction = MONITOR_STANDBY;     // No effect on my desktop?!
       } else if (STRCMP_NI(ForceStr, "off", 4) != 0) {
          ERROR_3("UnknownArgument",
                  "Unknown argument for [Monitor]: [%s].", ForceStr);
       }
     }
     
     // 16-Dec-2017 19:14
     // Richard Walker suggested to use SendNotifyMessage under Windows 10.
     // This worked under Win7:
     // SendMessage(HWND_BROADCAST,
     //             WM_SYSCOMMAND, SC_MONITORPOWER, MonitorAction);
     SendNotifyMessage(HWND_BROADCAST,
                       WM_SYSCOMMAND, SC_MONITORPOWER, MonitorAction);
     
     // Wobble the mouse to wake up the monitor under Windows 10:
     if (MonitorAction == MONITOR_ON) {
        mouse_event(MOUSEEVENTF_MOVE, 0, 1, 0, NULL);
        Sleep(40);
        mouse_event(MOUSEEVENTF_MOVE, 0, -1, 0, NULL);
     }
             
  } else {
     ERROR_3("UnknownCommand", "Unknown Command: [%s].", Command);
  }
  
  // Check number of outputs:
  if (nlhs > nOutput) {
     ERROR_2("BadNOutput", "Bad number of outputs.");
  }
  
  // Call the core to set the power status:
  if (Action != myVOID) {
     // Parse 2nd input argument:
     if (nrhs > ForceArgIndex) {
        mxGetString(prhs[ForceArgIndex], ForceStr, ForceStr_LEN);
        if (STRCMP_NI(ForceStr, "force", 6) == 0) {
           Force = EWX_FORCE;
        } else if (STRCMP_NI(ForceStr, "forceifhung", 12) == 0) {
           Force = EWX_FORCEIFHUNG;
        } else if (STRCMP_NI(ForceStr, "noforce", 8) == 0) {
           Force = 0;
        } else if (!mxIsEmpty(prhs[ForceArgIndex])){
           ERROR_2("BadValueInput2", "2nd input not recognized.");
        }
     }
     
     // Set the power status:
     mySetPower(Flag | Force, Action);
  }
  
  return;
}

// =============================================================================
void mySetPower(UINT Flag, enum PowerAction Action)
{
  // Call ExitWindowsEx or SetSuspendState with the apropriate flags.
  HANDLE hToken;
  TOKEN_PRIVILEGES tkp;
  DWORD  Success;
  BOOL   Exit_Ok = TRUE, Suspend_Ok = TRUE, SetToken_Ok;
  
  // Reason for the action (not used for Sleep and Hibernate):
  // "Other (Planned)": A planned shutdown or restart.
  DWORD Reason = SHTDN_REASON_MAJOR_OTHER | SHTDN_REASON_MINOR_OTHER |
                 SHTDN_REASON_FLAG_PLANNED;
  
  // Get a token for this process:
  if (!OpenProcessToken(GetCurrentProcess(),
      TOKEN_ADJUST_PRIVILEGES | TOKEN_QUERY, &hToken)) {
     ERROR_2("TokenNotOpend", "Cannot open token.");
  }
  
  // Get the LUID for the shutdown privilege:
  LookupPrivilegeValue(NULL, SE_SHUTDOWN_NAME, &tkp.Privileges[0].Luid);
  tkp.PrivilegeCount           = (DWORD) 1;
  tkp.Privileges[0].Attributes = SE_PRIVILEGE_ENABLED;

  // Get the shutdown privilege for this process:
  SetToken_Ok = AdjustTokenPrivileges(hToken, FALSE, &tkp, 0,
                                      (PTOKEN_PRIVILEGES) NULL, 0);
  
  // Set the power status:
  if (SetToken_Ok) {
     switch (Action) {
       case myEXIT:
         Exit_Ok = ExitWindowsEx(Flag, Reason);
         break;
       case mySLEEP:
         Suspend_Ok = SetSuspendState(FALSE, FALSE, FALSE);
         break;
       case myHIBERNATE:
         Suspend_Ok = SetSuspendState(TRUE, FALSE, FALSE);
         break;
       default:
         ERROR_2("BadSwitch", "Bad switch - programming error.");
     }
  }
  Success = GetLastError();
  
  // Remove privilegs before the error handling:
  tkp.Privileges[0].Attributes = 0;
  AdjustTokenPrivileges(hToken, FALSE, &tkp, 0,
                        (PTOKEN_PRIVILEGES) NULL, 0);
  CloseHandle(hToken);  // Not sure, but sounds logical to close the handle
  
  // Handle errors:
  if (!SetToken_Ok) {
     ERROR_3("BadNInput", "Setting token privilegs failed: [%d].", Success);
  } else if (!Exit_Ok) {
     ERROR_3("ExitWindowsExFailed", "ExitWindowsEx failed: [%d].", Success);
  } else if (!Exit_Ok) {
     ERROR_3("SetSuspendFailed", "SetSuspendState failed: [%d].", Success);
  }
  
  return;
}

// =============================================================================
void myAllowSleep(BOOL Enable)
{
  // Toggle the ability to fall into sleep depending on the [flag].
  extern BOOL doCleanupSleep;
  
  if (Enable) {  // Allow sleep mode:
     SetThreadExecutionState(ES_CONTINUOUS);
     doCleanupSleep = FALSE;
     
  } else {       // Disable sleep mode:
     // Allow the screen saver:
     // SetThreadExecutionState(ES_CONTINUOUS | ES_SYSTEM_REQUIRED);
     
     // Disable even the screen saver:
     SetThreadExecutionState(ES_CONTINUOUS | ES_SYSTEM_REQUIRED |
                             ES_AWAYMODE_REQUIRED);
     
     // Care for restoring the sleep state:
     if (!doCleanupSleep) {
        mexAtExit(myCleanupSleep);
        doCleanupSleep = TRUE;
     }
  }
  
  return;
}

// =============================================================================
static void myCleanupSleep(void)
{
  // It would not be friendly to let the sleep mode disabled, when Matlab is
  // shut down unexpectedly.
  // This function is called, when the MEX file is cleared from the memory
  // either by "clear all", "clear mex", "clear WinPower" or when exiting the
  // Matlab session. mexLock() blocks "clear()", but this would have other
  // disadvantages.
  extern BOOL doCleanupSleep;
   
  if (doCleanupSleep) {
     myAllowSleep(TRUE);
     mexPrintf(WARN_HEAD "Sleep mode is restored.\n");
  }
  
  // return;
}

// =============================================================================
BOOL mySetMatlabArg(const mxArray *Arg)
{
  // Set parameters for restarting Matlab. This feature needs Vista or higher.
  wchar_t    ArgString[MAX_PATH];
  mwSize     Len = 0;
  HRESULT    Ok;
  REGAPP_FCN pRegApp;
  
  // Check is function is available - Vista or higher needed:
  pRegApp = (REGAPP_FCN) GetProcAddress(GetModuleHandle(TEXT("kernel32.dll")),
                                        "RegisterApplicationRestart");
  if (pRegApp == NULL) {
     mexPrintf(WARN_HEAD "Restart of Matlab needs Windows Vista or newer.\n");
     return FALSE;
  }
  
  // Get string from 2nd input argument:
  if (Arg != NULL) {
     Len = mxGetNumberOfElements(Arg);
     if (Len > MAX_PATH - 7) {
        ERROR_3("ParamTooLong",
                "Restart parameter longer than %d characters.", MAX_PATH - 7);
     }
  }
  
  // Copy the data array, add double quotes and terminate:
  if (Len > 0) {
     memcpy(ArgString, L"-r \"", 4 * sizeof(mxChar));
     memcpy(ArgString + 4, mxGetData(Arg), Len * sizeof(mxChar));
     ArgString[Len + 4] = L'"';
     ArgString[Len + 5] = L'\0';
  } else {        // Empty string as default:
     ArgString[0] = L'\0';
  }
  
  // Call the API for the actual action:
  // Ok = RegisterApplicationRestart((PCWSTR) ArgString, 0);
  Ok = pRegApp((PCWSTR) ArgString, 0);
  
  // Handle errors:
  if (Ok != S_OK) {
     if (Ok == E_INVALIDARG) {
        ERROR_2("InvalidRegisterArg", "Bad parameter for application restart.");
     } else {
        ERROR_2("UnknownRegisterProblem",
                "RegisterApplicationRestart failed.");
     }
  }
  
  return TRUE;
}

// =============================================================================
mxArray *myBatteryStatus(void)
{
  SYSTEM_POWER_STATUS Status;
  BOOL    Success;
  mxArray *Out;
  double  dValue;
  char    *cValue;
  BYTE    Flag;
  const char *Fields[6] = {"ACLine", "BatteryStatus", "BatteryCharging",
    "BatteryLifePercent", "BatteryLifeTime", "BatteryFullLifeTime"};
  
  // _WIN32_WINNT macro must be defined as >= 0x0400!
  Success = GetSystemPowerStatus(&Status);
  
  // Create the output: --------------------------------------------------------
  Out = mxCreateStructMatrix(1, 1, 6, Fields);

  // ACLine: Connected to AC-power, NaN if unknown:
  switch (Status.ACLineStatus) {
     case 0:   dValue = 0.0;         break;
     case 1:   dValue = 1.0;         break;
     default:  dValue = mxGetNaN();  break;
  }
  mxSetFieldByNumber(Out, 0, 0, mxCreateDoubleScalar(dValue));
  
  // BatteryCharging: Logical flag:
  Flag = Status.BatteryFlag;
  mxSetFieldByNumber(Out, 0, 2,
                     mxCreateLogicalScalar((Flag & 0x08) == 0x08));
  
  // BatteryStatus:
  Flag &= 0xf7;  // ~8, Clear the CHARGING flag
  switch (Flag) {
     case 1:    cValue = "high";       break;
     case 2:    cValue = "low";        break;
     case 3:    cValue = "critical";   break;
     case 128:  cValue = "no_battery"; break;
     case 255:  // Fall through!
     default:   cValue = "unknown";    break;
  }
  mxSetFieldByNumber(Out, 0, 1, mxCreateString(cValue));
    
  // BatteryLifePercent:
  dValue = (double) Status.BatteryLifePercent;
  if (dValue > 100.0) {
     dValue = mxGetNaN();
  }
  mxSetFieldByNumber(Out, 0, 3, mxCreateDoubleScalar(dValue));
  
  // BatteryLifeTime:
  dValue = (double) (int32_T) Status.BatteryLifeTime;
  if (dValue < 0) {
     dValue = mxGetNaN();
  }
  mxSetFieldByNumber(Out, 0, 4, mxCreateDoubleScalar(dValue));
  
  // BatteryFullLifeTime:
  dValue = (double) (int32_T) Status.BatteryFullLifeTime;
  if (dValue < 0) {
     dValue = mxGetNaN();
  }
  mxSetFieldByNumber(Out, 0, 5, mxCreateDoubleScalar(dValue));
  
  return Out;
}

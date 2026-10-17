/* include the COIN-OR-wide system specific configure header */
#include "configall_system.h"

/* this needs to come before the include of config_smi_default.h */
#ifndef SMILIB_EXPORT
#if defined(_WIN32) && defined(DLL_EXPORT)
#define SMILIB_EXPORT __declspec(dllexport)
#else
#define SMILIB_EXPORT
#endif
#endif

/* include the public project specific macros */
#include "config_smi_default.h"

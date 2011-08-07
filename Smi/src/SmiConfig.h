/* Copyright (C) 2011
 * All Rights Reserved.
 * This code is published under the Eclipse Public License.
 *
 * $Id: OsiConfig.h 1735 2011-06-08 17:00:08Z stefan $
 *
 * Include file for the configuration of Smi.
 *
 * On systems where the code is configured with the configure script
 * (i.e., compilation is always done with HAVE_CONFIG_H defined), this
 * header file includes the automatically generated header file.
 *
 * On systems that are compiled in other ways (e.g., with the
 * Developer Studio), a default header files is included.
 */

#ifndef __SMICONFIG_H__
#define __SMICONFIG_H__

#ifdef HAVE_CONFIG_H
#include "config_smi.h"
#else
#include "config_smi_default.h"
#endif

#endif /*__SMICONFIG_H__*/

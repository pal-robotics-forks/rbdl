/*
 * RBDL - Rigid Body Dynamics Library
 * Copyright (c) 2011-2012 Martin Felis <martin.felis@iwr.uni-heidelberg.de>
 *
 * Licensed under the zlib license. See LICENSE for more details.
 */

#ifndef _RBDLCONFIG_H
#define _RBDLCONFIG_H

#define RBDL_API_VERSION 0x020000

/* #undef RBDL_USE_SIMPLE_MATH */
/* #undef RBDL_ENABLE_LOGGING */
#define RBDL_BUILD_REVISION "unknown"
#define RBDL_BUILD_TYPE "unknown"
#define RBDL_BUILD_BRANCH "unknown"
/* #undef BUILD_ADDON_LUAMODEL */

/* compatibility defines */
#ifdef _WIN32
	#define __func__ __FUNCTION__
#endif

#endif

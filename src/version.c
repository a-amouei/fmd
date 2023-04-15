/*
  version.c: This file is part of Free Molecular Dynamics

  Copyright (C) 2023 Arham Amouye Foumani

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#include "types.h"
#include <stdio.h>
#include <string.h>
#include "general.h"

int fmd_version_getMajor()
{
    return 0;
}

int fmd_version_getMinor()
{
    return 1;
}

int fmd_version_getRevision()
{
    return 0;
}

/* 'r' for release
   'd' for development */
char fmd_version_getType()
{
    return 'r';
}

fmd_string_t fmd_version_getString()
{
    char str[16];

    sprintf(str, "%d.%d.%d%c", fmd_version_getMajor(),
                               fmd_version_getMinor(),
                               fmd_version_getRevision(),
                               fmd_version_getType());

    fmd_string_t out = m_alloc(strlen(str) + 1);

    strcpy(out, str);

    return out;
}

/*
  timer.h: This file is part of Free Molecular Dynamics

  Copyright (C) 2019 Arham Amouye Foumani

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

#ifndef TIMER_H
#define TIMER_H

#include "config.h"
#include "types.h"

typedef enum
{
    TIMER_SIMPLE
} timercat_t;

typedef struct _fmd_timer
{
    timercat_t cat;
    fmd_bool_t enabled;
    double start;
    double stop;
    double interval;
} fmd_timer_t;

typedef struct _fmd fmd_t;

void fmd_timer_free(fmd_t *md);
void fmd_timer_sendTimerTickEvents(fmd_t *md);

#endif /* TIMER_H */

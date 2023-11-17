/*
  events.h: This file is part of Free Molecular Dynamics

  Copyright (C) 2020 Arham Amouye Foumani

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

#ifndef EVENTS_H
#define EVENTS_H

#include "config.h"
#include "types.h"

typedef enum
{
    FMD_EVENT_TIMER_TICK,
    FMD_EVENT_FIELD_UPDATE
} fmd_event_t;

typedef void (*fmd_EventHandler_t)(fmd_t *md, fmd_event_t event, void *usp, fmd_params_t *params);

typedef struct
{
    fmd_handle_t timer;
} fmd_event_params_timer_tick_t;

typedef struct
{
    fmd_handle_t turi;
    fmd_handle_t field;
} fmd_event_params_field_update_t;

#endif /* EVENTS_H */

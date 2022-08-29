/*
  timer.c: This file is part of Free Molecular Dynamics

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

#include "timer.h"
#include "base.h"
#include "events.h"
#include "general.h"

void fmd_timer_free(fmd_t *md)
{
    if (md->timers != NULL)
    {
        free(md->timers);
        md->timers = NULL;
        md->timers_num = 0;
    }
}

fmd_handle_t fmd_timer_makeSimple(fmd_t *md, fmd_real_t start, fmd_real_t interval, fmd_real_t stop)
{
    int i = md->timers_num;

    md->timers = (fmd_timer_t *)re_alloc(md->timers, (i+1) * sizeof(fmd_timer_t));

    md->timers[i].enabled = FMD_TRUE;
    md->timers[i].cat = TIMER_SIMPLE;
    md->timers[i].start = start;
    md->timers[i].interval = interval;
    md->timers[i].stop = stop;
    md->timers_num++;

    return i;
}

void _fmd_timer_sendTimerTickEvents(fmd_t *md)
{
    for (int i=0; i < md->timers_num; i++)
    {
        if ( md->timers[i].enabled && md->time >= md->timers[i].start &&
             !(md->time > md->timers[i].stop && md->timers[i].stop >= md->timers[i].start) )
        {
            if (md->timers[i].cat == TIMER_SIMPLE)
            {
                if (_fmd_timer_is_its_time(md->time, md->timestep/2.0, md->timers[i].start, md->timers[i].interval))
                {
                    fmd_event_params_timer_tick_t params;

                    params.timer = i;
                    md->EventHandler(md, FMD_EVENT_TIMER_TICK, (fmd_params_t *)&params);
                }
            }
        }
    }
}

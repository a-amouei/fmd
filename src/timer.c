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

unsigned fmd_timer_makeSimple(fmd_t *md, double start, double interval, double stop)
{
    int i = md->timers_num;

    md->timers = (fmd_timer_t *)realloc(md->timers, (i+1) * sizeof(fmd_timer_t));
    // TO-DO: handle memory error
    assert(md->timers != NULL);
    md->timers[i].enabled = 1;
    md->timers[i].kind = TIMERKIND_SIMPLE;
    md->timers[i].start = start;
    md->timers[i].interval = interval;
    md->timers[i].stop = stop;
    md->timers_num++;

    return i;
}

void fmd_timer_sendTimerTickEvents(fmd_t *md)
{
    for (unsigned i=0; i < md->timers_num; i++)
    {
        if ( md->timers[i].enabled && md->mdTime >= md->timers[i].start &&
             !(md->mdTime > md->timers[i].stop && md->timers[i].stop >= md->timers[i].start) )
        {
            if (md->timers[i].kind == TIMERKIND_SIMPLE)
            {
                if (fmod( fabs(md->mdTime - md->timers[i].start), md->timers[i].interval ) < md->delta_t)
                    md->eventHandler(md, FMD_EVENT_TIMERTICK, i);
            }
        }
    }
}

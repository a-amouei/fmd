/*
  list.c: This file is part of Free Molecular Dynamics

  Copyright (C) 2019 Arham Amouye Foumani, Hossein Ghorbanfekr

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

/* this module owes some ideas to the doubly-linked list of GLib */

#include <stdlib.h>
#include "list.h"
#include "general.h"

list_t *fmd_list_prepend(list_t *list, void *data)
{
    list_t *item = malloc(sizeof(list_t));

    if (item != NULL)
    {
        item->data = data;
        item->next = list;
        if (list == NULL)
            item->prev = NULL;
        else
        {
            item->prev = list->prev;
            if (list->prev != NULL)
                list->prev->next = item;
            list->prev = item;
        }
    }

    return item;
}

unsigned fmd_list_length(list_t *list)
{
    unsigned length = 0;

    while (list != NULL)
    {
        length++;
        list = list->next;
    }

    return length;
}

/* returns NULL if not found */
list_t *fmd_list_find_custom(list_t *list, const void *data, comparefunc_t func)
{
    while (list != NULL)
    {
        if (func(list->data, data) == 0)
            return list;
        list = list->next;
    }

    return NULL;
}

// not needed at the moment
/*list_t *fmd_list_find(list_t *list, const void *data)
{
    while (list != NULL)
    {
        if (list->data == data)
            return list;
        list = list->next;
    }

    return NULL;
}*/

void fmd_list_free(list_t *list)
{
    while (list != NULL)
    {
        list_t *next = list->next;
        free(list);
        list = next;
    }
}

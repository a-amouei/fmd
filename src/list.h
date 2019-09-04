/*
  list.h: This file is part of Free Molecular Dynamics

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

#ifndef LIST_H
#define LIST_H

#include "config.h"

typedef int (*comparefunc_t)(const void *a, const void *b);

typedef struct _list list_t;

struct _list
{
    void *data;
    list_t *next;
    list_t *prev;
};

list_t *fmd_list_prepend(list_t *list, void *data);
unsigned fmd_list_length(list_t *list);
list_t *fmd_list_find_custom(list_t *list, const void *data, comparefunc_t func);
void fmd_list_free(list_t *list);

#endif /* LIST_H */

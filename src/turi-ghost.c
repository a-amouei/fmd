/*
  turi_ghost.c: This file is part of Free Molecular Dynamics

  Copyright (C) 2021 Arham Amouye Foumani

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

#include "turi-ghost.h"
#include "fmd-private.h"
#include "types.h"
#include "misc.h"
#include "general.h"
#include "turi.h"

void _fmd_turi_update_ghosts_1d(fmd_t *md, turi_t *t, int d, tghost_pack_t *p)
{
    MPI_Status status;
    int count;

    /* sending to lower process, receiving from upper process */

    if (t->rank_of_lower_owner[d] != MPI_PROC_NULL)
        p->pack._1D(md, t, false, p->sendbuf, &count);
    else
        count = 0;

    MPI_Sendrecv(p->sendbuf, count, MPI_PACKED, t->rank_of_lower_owner[d], 85101,
                 p->recvbuf, p->bufsize, MPI_PACKED, t->rank_of_upper_owner[d], 85101,
                 md->MD_comm, &status);

    if (t->rank_of_upper_owner[d] != MPI_PROC_NULL)
        p->unpack._1D(md, t, false, p->recvbuf);

    /* sending to upper process, receiving from lower process */

    if (t->rank_of_upper_owner[d] != MPI_PROC_NULL)
        p->pack._1D(md, t, true, p->sendbuf, &count);
    else
        count = 0;

    MPI_Sendrecv(p->sendbuf, count, MPI_PACKED, t->rank_of_upper_owner[d], 85103,
                 p->recvbuf, p->bufsize, MPI_PACKED, t->rank_of_lower_owner[d], 85103,
                 md->MD_comm, &status);

    if (t->rank_of_lower_owner[d] != MPI_PROC_NULL)
        p->unpack._1D(md, t, true, p->recvbuf);
}

static void communicate_in_direction_d(
    fmd_t *md, turi_t *t, int d, fmd_ituple_t vitc_start_send_lower,
    fmd_ituple_t vitc_stop_send_lower, fmd_ituple_t vitc_start_receive_lower,
    fmd_ituple_t vitc_stop_receive_lower, fmd_ituple_t vitc_start_send_upper,
    fmd_ituple_t vitc_stop_send_upper, fmd_ituple_t vitc_start_receive_upper,
    fmd_ituple_t vitc_stop_receive_upper, fields_packer_t pack,
    fields_unpacker_t unpack)
{
    MPI_Status status;
    size_t packsize;
    void *data_receive, *data_send;

    /* sending to lower process, receiving from upper process */

    if (t->rank_of_lower_owner[d] != MPI_PROC_NULL)
        pack(md, t, vitc_start_send_lower, vitc_stop_send_lower, &packsize, &data_send);
    else
        pack(md, t, vitc_start_send_lower, vitc_stop_send_lower, &packsize, NULL);  /* if no lower process, only calculate packsize */

    data_receive = m_alloc(packsize);

    MPI_Sendrecv(data_send, packsize, MPI_PACKED, t->rank_of_lower_owner[d], 85101,
                 data_receive, packsize, MPI_PACKED, t->rank_of_upper_owner[d], 85101,
                 md->MD_comm, &status);

    if (t->rank_of_lower_owner[d] != MPI_PROC_NULL) free(data_send);

    if (t->rank_of_upper_owner[d] != MPI_PROC_NULL)
        unpack(md, t, vitc_start_receive_upper, vitc_stop_receive_upper, data_receive);

    /* sending to upper process, receiving from lower process */

    if (t->rank_of_upper_owner[d] != MPI_PROC_NULL)
        pack(md, t, vitc_start_send_upper, vitc_stop_send_upper, &packsize, &data_send);

    MPI_Sendrecv(data_send, packsize, MPI_PACKED, t->rank_of_upper_owner[d], 85103,
                 data_receive, packsize, MPI_PACKED, t->rank_of_lower_owner[d], 85103,
                 md->MD_comm, &status);

    if (t->rank_of_upper_owner[d] != MPI_PROC_NULL) free(data_send);

    if (t->rank_of_lower_owner[d] != MPI_PROC_NULL)
        unpack(md, t, vitc_start_receive_lower, vitc_stop_receive_lower, data_receive);

    free(data_receive);
}

/* vitc stands for "virtual index of turi-cell" */
static void prepare_communication_in_direction_d(
    turi_t *t, int d, fmd_ituple_t vitc_start_send_lower,
    fmd_ituple_t vitc_stop_send_lower, fmd_ituple_t vitc_start_receive_lower,
    fmd_ituple_t vitc_stop_receive_lower, fmd_ituple_t vitc_start_send_upper,
    fmd_ituple_t vitc_stop_send_upper, fmd_ituple_t vitc_start_receive_upper,
    fmd_ituple_t vitc_stop_receive_upper)
{
    int dd;

    for (dd=0; dd<3; dd++)
    {
        if (dd == d) /* only ghost */
        {
            vitc_start_send_lower[dd] = t->itc_start_owned[dd];
            vitc_stop_send_lower[dd] = vitc_start_send_lower[dd] + 1;
            vitc_stop_receive_lower[dd] = t->itc_start_owned[dd];
            vitc_start_receive_lower[dd] = vitc_stop_receive_lower[dd] - 1;
            vitc_stop_send_upper[dd] = t->itc_stop[dd];
            vitc_start_send_upper[dd] = vitc_stop_send_upper[dd] - 1;
            vitc_start_receive_upper[dd] = t->itc_stop[dd];
            vitc_stop_receive_upper[dd] = vitc_start_receive_upper[dd] + 1;
        }
        else if (dd > d) /* including ghost */
        {
            int m = (t->has_upper_lower_owner_procs[dd] ? 1 : 0);

            vitc_start_receive_lower[dd] = vitc_start_send_lower[dd]
                                         = vitc_start_receive_upper[dd]
                                         = vitc_start_send_upper[dd]
                                         = t->itc_start_owned[dd] - m;
            vitc_stop_receive_lower[dd]  = vitc_stop_send_lower[dd]
                                         = vitc_stop_receive_upper[dd]
                                         = vitc_stop_send_upper[dd]
                                         = t->itc_stop[dd] + m;
        }
        else /* excluding ghost */
        {
            vitc_start_receive_lower[dd] = vitc_start_send_lower[dd]
                                         = vitc_start_receive_upper[dd]
                                         = vitc_start_send_upper[dd]
                                         = t->itc_start_owned[dd];
            vitc_stop_receive_lower[dd]  = vitc_stop_send_lower[dd]
                                         = vitc_stop_receive_upper[dd]
                                         = vitc_stop_send_upper[dd]
                                         = t->itc_stop[dd];
        }
    }
}

void _fmd_turi_update_ghosts(fmd_t *md, turi_t *t, fields_packer_t pack, fields_unpacker_t unpack)
{
    int d;
    fmd_ituple_t vitc_start_send_lower, vitc_stop_send_lower;
    fmd_ituple_t vitc_start_send_upper, vitc_stop_send_upper;
    fmd_ituple_t vitc_start_receive_lower, vitc_stop_receive_lower;
    fmd_ituple_t vitc_start_receive_upper, vitc_stop_receive_upper;

    for (d = 3-1; d >= 0; d--)
    {
        if (t->has_upper_lower_owner_procs[d])    /* if this isn't the only owner
                                                     subdomain in direction d */
        {
            prepare_communication_in_direction_d(
                t, d, vitc_start_send_lower, vitc_stop_send_lower,
                vitc_start_receive_lower, vitc_stop_receive_lower,
                vitc_start_send_upper, vitc_stop_send_upper,
                vitc_start_receive_upper, vitc_stop_receive_upper);
            communicate_in_direction_d(
                md, t, d, vitc_start_send_lower, vitc_stop_send_lower,
                vitc_start_receive_lower, vitc_stop_receive_lower,
                vitc_start_send_upper, vitc_stop_send_upper, vitc_start_receive_upper,
                vitc_stop_receive_upper, pack, unpack);
        }
    }
}

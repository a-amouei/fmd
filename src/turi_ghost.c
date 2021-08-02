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

#include "turi_ghost.h"
#include "types.h"
#include "base.h"
#include "general.h"
#include "turi.h"

static void communicate_in_direction_d(
    fmd_t *md, turi_t *t, int d, fmd_ituple_t vitc_start_send_lower,
    fmd_ituple_t vitc_stop_send_lower, fmd_ituple_t vitc_start_receive_lower,
    fmd_ituple_t vitc_stop_receive_lower, fmd_ituple_t vitc_start_send_upper,
    fmd_ituple_t vitc_stop_send_upper, fmd_ituple_t vitc_start_receive_upper,
    fmd_ituple_t vitc_stop_receive_upper, fields_packer_t pack,
    fields_unpacker_t unpack)
{
    MPI_Status status;
    MPI_Request request;
    size_t packsize;
    void *data_receive, *data_send;

    if (t->rank_of_lower_owner[d] >= md->SubDomain.myrank && !md->PBC[d])
    {
        /* we need pack() for sending data; it's called here because we need packsize for receiving data */
        pack(md, t, vitc_start_send_upper, vitc_stop_send_upper, &packsize, &data_send);

        /* receiving from upper process */

        data_receive = malloc(packsize);
        assert(data_receive != NULL); /* TO-DO: handle memory error */

        MPI_Recv(data_receive, packsize, MPI_BYTE, t->rank_of_upper_owner[d], 85101, md->MD_comm, &status);

        unpack(md, t, vitc_start_receive_upper, vitc_stop_receive_upper, data_receive);

        free(data_receive);

        /* sending to upper process */

        MPI_Send(data_send, packsize, MPI_BYTE, t->rank_of_upper_owner[d], 85103, md->MD_comm);

        free(data_send);
    }
    else
        if (t->rank_of_upper_owner[d] <= md->SubDomain.myrank && !md->PBC[d])
        {
            /* sending to lower process */

            pack(md, t, vitc_start_send_lower, vitc_stop_send_lower, &packsize, &data_send);

            MPI_Send(data_send, packsize, MPI_BYTE, t->rank_of_lower_owner[d], 85101, md->MD_comm);

            free(data_send);

            /* receiving from lower process */

            data_receive = malloc(packsize);
            assert(data_receive != NULL); /* TO-DO: handle memory error */

            MPI_Recv(data_receive, packsize, MPI_BYTE, t->rank_of_lower_owner[d], 85103, md->MD_comm, &status);

            unpack(md, t, vitc_start_receive_lower, vitc_stop_receive_lower, data_receive);

            free(data_receive);
        }
        else
        {
            /* sending to lower process, receiving from upper process */

            pack(md, t, vitc_start_send_lower, vitc_stop_send_lower, &packsize, &data_send);

            data_receive = malloc(packsize);
            assert(data_receive != NULL); /* TO-DO: handle memory error */

            MPI_Isend(data_send, packsize, MPI_BYTE, t->rank_of_lower_owner[d], 85101, md->MD_comm, &request);
            MPI_Recv(data_receive, packsize, MPI_BYTE, t->rank_of_upper_owner[d], 85101, md->MD_comm, &status);
            MPI_Wait(&request, &status);

            free(data_send);

            unpack(md, t, vitc_start_receive_upper, vitc_stop_receive_upper, data_receive);

            /* sending to upper process, receiving from lower process */

            pack(md, t, vitc_start_send_upper, vitc_stop_send_upper, &packsize, &data_send);

            MPI_Isend(data_send, packsize, MPI_BYTE, t->rank_of_upper_owner[d], 85103, md->MD_comm, &request);
            MPI_Recv(data_receive, packsize, MPI_BYTE, t->rank_of_lower_owner[d], 85103, md->MD_comm, &status);
            MPI_Wait(&request, &status);

            free(data_send);

            unpack(md, t, vitc_start_receive_lower, vitc_stop_receive_lower, data_receive);

            free(data_receive);
        }
}

/* vitc stands for "virtual index of turi-cell" */
static void prepare_communication_in_direction_d(
    fmd_t *md, turi_t *t, int d, fmd_ituple_t vitc_start_send_lower,
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
            int m = (t->rank_of_lower_owner[dd] == md->SubDomain.myrank ? 0 : 1);

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
        if (t->rank_of_lower_owner[d] != md->SubDomain.myrank) /* if this isn't the only owner
                                                                  subdomain in direction d */
        {
            prepare_communication_in_direction_d(
                md, t, d, vitc_start_send_lower, vitc_stop_send_lower,
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

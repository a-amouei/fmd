/*
  md_ghost.c: This file is part of Free Molecular Dynamics

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

/* functions for working with ghost cells of MD subdomains, including
   communication of data from or to them */


#include "base.h"
#include "md_ghost.h"
#include "types.h"

typedef struct
{
    fmd_rtuple_t x;
    unsigned atomkind;
    int GroupID;
} position_struct_t;

static void particles_migrate_in_direction_d(
    fmd_t *md, int d, fmd_ituple_t ic_start_send_lower,
    fmd_ituple_t ic_stop_send_lower, fmd_ituple_t ic_start_receive_lower,
    fmd_ituple_t ic_stop_receive_lower, fmd_ituple_t ic_start_send_upper,
    fmd_ituple_t ic_stop_send_upper, fmd_ituple_t ic_start_receive_upper,
    fmd_ituple_t ic_stop_receive_upper)
{
    MPI_Status status;
    MPI_Request request;
    int sum_length_send, sum_length_receive;
    int k, kreceive, cells_num;
    int *cells_length_send, *cells_length_receive;
    fmd_ituple_t ic;
    particle_core_t *pcs_send, *pcs_receive;
    int dd;

    if ( ((md->SubDomain.is[d] == 0) || (md->SubDomain.is[d] == md->ns[d]-1)) && !md->PBC[d] )
    {
        if (md->SubDomain.is[d] == 0)
        {
            /* receiving from upper process */
            cells_num = 1;
            for (dd=0; dd<3; dd++)
                cells_num *= ic_stop_receive_upper[dd] - ic_start_receive_upper[dd];
            cells_length_receive = (int *)malloc((cells_num+1) * sizeof(int));
            MPI_Recv(cells_length_receive, cells_num+1, MPI_INT, md->SubDomain.rank_of_upper_subd[d],
                     1, md->MD_comm, &status);
            md->SubDomain.NumberOfParticles += cells_length_receive[cells_num];
            sum_length_receive = cells_length_receive[cells_num] * sizeof(particle_core_t);
            pcs_receive = (particle_core_t *)malloc(sum_length_receive);
            MPI_Recv(pcs_receive, sum_length_receive, MPI_CHAR,
                     md->SubDomain.rank_of_upper_subd[d], 2, md->MD_comm, &status);

            kreceive = k = 0;
            ITERATE(ic, ic_start_receive_upper, ic_stop_receive_upper)
            {
                for (int cc=0; cc < cells_length_receive[kreceive]; cc++)
                    INSERT_PART_CORE_IN_CELL(pcs_receive[k++], md->SubDomain.grid[ic[0]][ic[1]][ic[2]]);
                kreceive++;
            }

            free(cells_length_receive);
            free(pcs_receive);
            /* */
            ITERATE(ic, ic_start_send_lower, ic_stop_send_lower)
            {
                md->SubDomain.NumberOfParticles -= md->SubDomain.grid[ic[0]][ic[1]][ic[2]].parts_num;
                MINIMIZE_CELL(md->SubDomain.grid[ic[0]][ic[1]][ic[2]]);
            }

            /* sending to upper process */
            cells_length_send = (int *)malloc((cells_num+1) * sizeof(int));
            k = sum_length_send = 0;
            ITERATE(ic, ic_start_send_upper, ic_stop_send_upper)
            {
                cells_length_send[k] = md->SubDomain.grid[ic[0]][ic[1]][ic[2]].parts_num;
                sum_length_send += cells_length_send[k++];
            }
            cells_length_send[cells_num] = sum_length_send;
            MPI_Send(cells_length_send, cells_num+1, MPI_INT, md->SubDomain.rank_of_upper_subd[d], 3,
                     md->MD_comm);
            free(cells_length_send);
            md->SubDomain.NumberOfParticles -= sum_length_send;
            sum_length_send *= sizeof(particle_core_t);
            pcs_send = (particle_core_t *)malloc(sum_length_send);

            k=0;
            ITERATE(ic, ic_start_send_upper, ic_stop_send_upper)
            {
                for (int cc=0; cc < md->SubDomain.grid[ic[0]][ic[1]][ic[2]].parts_num; cc++)
                    pcs_send[k++] = md->SubDomain.grid[ic[0]][ic[1]][ic[2]].parts[cc].core;

                MINIMIZE_CELL(md->SubDomain.grid[ic[0]][ic[1]][ic[2]]);
            }

            MPI_Send(pcs_send, sum_length_send, MPI_CHAR,
                     md->SubDomain.rank_of_upper_subd[d], 4, md->MD_comm);
            free(pcs_send);
        }
        else /* if (md->SubDomain.is[d] == md->ns[d]-1) */
        {
            /* sending to lower process */
            cells_num = 1;
            for (dd=0; dd<3; dd++)
                cells_num *= ic_stop_send_lower[dd] - ic_start_send_lower[dd];
            cells_length_send = (int *)malloc((cells_num+1) * sizeof(int));
            k = sum_length_send = 0;
            ITERATE(ic, ic_start_send_lower, ic_stop_send_lower)
            {
                cells_length_send[k] = md->SubDomain.grid[ic[0]][ic[1]][ic[2]].parts_num;
                sum_length_send += cells_length_send[k++];
            }
            cells_length_send[cells_num] = sum_length_send;
            MPI_Send(cells_length_send, cells_num+1, MPI_INT, md->SubDomain.rank_of_lower_subd[d], 1,
                     md->MD_comm);
            free(cells_length_send);
            md->SubDomain.NumberOfParticles -= sum_length_send;
            sum_length_send *= sizeof(particle_core_t);
            pcs_send = (particle_core_t *)malloc(sum_length_send);

            k=0;
            ITERATE(ic, ic_start_send_lower, ic_stop_send_lower)
            {
                for (int cc=0; cc < md->SubDomain.grid[ic[0]][ic[1]][ic[2]].parts_num; cc++)
                    pcs_send[k++] = md->SubDomain.grid[ic[0]][ic[1]][ic[2]].parts[cc].core;

                MINIMIZE_CELL(md->SubDomain.grid[ic[0]][ic[1]][ic[2]]);
            }

            MPI_Send(pcs_send, sum_length_send, MPI_CHAR,
                     md->SubDomain.rank_of_lower_subd[d], 2, md->MD_comm);
            free(pcs_send);

            /* receiving from lower process */
            cells_length_receive = (int *)malloc((cells_num+1) * sizeof(int));
            MPI_Recv(cells_length_receive, cells_num+1, MPI_INT, md->SubDomain.rank_of_lower_subd[d], 3,
                     md->MD_comm, &status);
            md->SubDomain.NumberOfParticles += cells_length_receive[cells_num];
            sum_length_receive = cells_length_receive[cells_num] * sizeof(particle_core_t);
            pcs_receive = (particle_core_t *)malloc(sum_length_receive);
            MPI_Recv(pcs_receive, sum_length_receive, MPI_CHAR,
                     md->SubDomain.rank_of_lower_subd[d], 4, md->MD_comm, &status);
            kreceive = k = 0;
            ITERATE(ic, ic_start_receive_lower, ic_stop_receive_lower)
            {
                for (int cc=0; cc<cells_length_receive[kreceive]; cc++)
                    INSERT_PART_CORE_IN_CELL(pcs_receive[k++], md->SubDomain.grid[ic[0]][ic[1]][ic[2]]);
                kreceive++;
            }
            free(cells_length_receive);
            free(pcs_receive);
            /* */
            ITERATE(ic, ic_start_send_upper, ic_stop_send_upper)
            {
                md->SubDomain.NumberOfParticles -= md->SubDomain.grid[ic[0]][ic[1]][ic[2]].parts_num;
                MINIMIZE_CELL(md->SubDomain.grid[ic[0]][ic[1]][ic[2]]);
            }
        }
    }
    else
    {
        /* sending to lower process, receiving from upper process */
        cells_num = 1;
        for (dd=0; dd<3; dd++)
            cells_num *= ic_stop_send_lower[dd] - ic_start_send_lower[dd];
        cells_length_send = (int *)malloc((cells_num+1) * sizeof(int));
        cells_length_receive = (int *)malloc((cells_num+1) * sizeof(int));
        k = sum_length_send = 0;
        ITERATE(ic, ic_start_send_lower, ic_stop_send_lower)
        {
            cells_length_send[k] = md->SubDomain.grid[ic[0]][ic[1]][ic[2]].parts_num;
            sum_length_send += cells_length_send[k++];
        }
        cells_length_send[cells_num] = sum_length_send;
        MPI_Isend(cells_length_send, cells_num+1, MPI_INT, md->SubDomain.rank_of_lower_subd[d], 1,
                  md->MD_comm, &request);
        MPI_Recv(cells_length_receive, cells_num+1, MPI_INT, md->SubDomain.rank_of_upper_subd[d], 1,
                 md->MD_comm, &status);
        MPI_Wait(&request, &status);
        md->SubDomain.NumberOfParticles += cells_length_receive[cells_num] - sum_length_send;
        sum_length_send *= sizeof(particle_core_t);
        pcs_send = (particle_core_t *)malloc(sum_length_send);
        sum_length_receive = cells_length_receive[cells_num] * sizeof(particle_core_t);
        pcs_receive = (particle_core_t *)malloc(sum_length_receive);

        k=0;
        ITERATE(ic, ic_start_send_lower, ic_stop_send_lower)
        {
            for (int cc=0; cc < md->SubDomain.grid[ic[0]][ic[1]][ic[2]].parts_num; cc++)
                pcs_send[k++] = md->SubDomain.grid[ic[0]][ic[1]][ic[2]].parts[cc].core;

            MINIMIZE_CELL(md->SubDomain.grid[ic[0]][ic[1]][ic[2]]);
        }

        MPI_Isend(pcs_send, sum_length_send, MPI_CHAR,
                  md->SubDomain.rank_of_lower_subd[d], 2, md->MD_comm, &request);
        MPI_Recv(pcs_receive, sum_length_receive, MPI_CHAR,
                 md->SubDomain.rank_of_upper_subd[d], 2, md->MD_comm, &status);
        MPI_Wait(&request, &status);
        free(pcs_send);

        kreceive = k = 0;
        ITERATE(ic, ic_start_receive_upper, ic_stop_receive_upper)
        {
            for (int cc=0; cc < cells_length_receive[kreceive]; cc++)
            {
                if (md->SubDomain.is[d] == md->ns[d] - 1)
                    pcs_receive[k].x[d] += md->l[d];
                INSERT_PART_CORE_IN_CELL(pcs_receive[k++], md->SubDomain.grid[ic[0]][ic[1]][ic[2]]);
            }
            kreceive++;
        }

        free(pcs_receive);

        /* sending to upper process, receiving from lower process */
        k = sum_length_send = 0;
        ITERATE(ic, ic_start_send_upper, ic_stop_send_upper)
        {
            cells_length_send[k] = md->SubDomain.grid[ic[0]][ic[1]][ic[2]].parts_num;
            sum_length_send += cells_length_send[k++];
        }
        cells_length_send[cells_num] = sum_length_send;
        MPI_Isend(cells_length_send, cells_num+1, MPI_INT, md->SubDomain.rank_of_upper_subd[d], 3,
                  md->MD_comm, &request);
        MPI_Recv(cells_length_receive, cells_num+1, MPI_INT, md->SubDomain.rank_of_lower_subd[d], 3,
                 md->MD_comm, &status);
        MPI_Wait(&request, &status);
        free(cells_length_send);
        md->SubDomain.NumberOfParticles += cells_length_receive[cells_num] - sum_length_send;
        sum_length_send *= sizeof(particle_core_t);
        pcs_send = (particle_core_t *)malloc(sum_length_send);
        sum_length_receive = cells_length_receive[cells_num] * sizeof(particle_core_t);
        pcs_receive = (particle_core_t *)malloc(sum_length_receive);

        k=0;
        ITERATE(ic, ic_start_send_upper, ic_stop_send_upper)
        {
            for (int cc=0; cc < md->SubDomain.grid[ic[0]][ic[1]][ic[2]].parts_num; cc++)
                pcs_send[k++] = md->SubDomain.grid[ic[0]][ic[1]][ic[2]].parts[cc].core;

            MINIMIZE_CELL(md->SubDomain.grid[ic[0]][ic[1]][ic[2]]);
        }

        MPI_Isend(pcs_send, sum_length_send, MPI_CHAR,
                  md->SubDomain.rank_of_upper_subd[d], 4, md->MD_comm, &request);
        MPI_Recv(pcs_receive, sum_length_receive, MPI_CHAR,
                 md->SubDomain.rank_of_lower_subd[d], 4, md->MD_comm, &status);
        MPI_Wait(&request, &status);
        free(pcs_send);

        kreceive = k = 0;
        ITERATE(ic, ic_start_receive_lower, ic_stop_receive_lower)
        {
            for (int cc=0; cc < cells_length_receive[kreceive]; cc++)
            {
                if (md->SubDomain.is[d] == 0)
                    pcs_receive[k].x[d] -= md->l[d];
                INSERT_PART_CORE_IN_CELL(pcs_receive[k++], md->SubDomain.grid[ic[0]][ic[1]][ic[2]]);
            }
            kreceive++;
        }

        free(cells_length_receive);
        free(pcs_receive);
    }
}

static void ghostparticles_init_in_direction_d(
    fmd_t *md, int d, fmd_ituple_t ic_start_send_lower, fmd_ituple_t ic_stop_send_lower,
    fmd_ituple_t ic_start_receive_lower, fmd_ituple_t ic_stop_receive_lower,
    fmd_ituple_t ic_start_send_upper, fmd_ituple_t ic_stop_send_upper,
    fmd_ituple_t ic_start_receive_upper, fmd_ituple_t ic_stop_receive_upper)
{
    MPI_Status status;
    MPI_Request request;
    int sum_length_send, sum_length_receive;
    int k, kreceive, cells_num;
    int *cells_length_send, *cells_length_receive;
    fmd_ituple_t ic;
    position_struct_t *data_send, *data_receive;
    int dd;

    if ( ((md->SubDomain.is[d] == 0) || (md->SubDomain.is[d] == md->ns[d]-1)) && !md->PBC[d] )
    {
        if (md->SubDomain.is[d] == 0)
        {
            /* receiving from upper process */
            cells_num = 1;
            for (dd=0; dd<3; dd++)
                cells_num *= ic_stop_receive_upper[dd] - ic_start_receive_upper[dd];
            cells_length_receive = (int *)malloc((cells_num+1) * sizeof(int));
            MPI_Recv(cells_length_receive, cells_num+1, MPI_INT, md->SubDomain.rank_of_upper_subd[d], 5,
                     md->MD_comm, &status);
            sum_length_receive = sizeof(position_struct_t) * cells_length_receive[cells_num];
            data_receive = (position_struct_t *)malloc(sum_length_receive);
            MPI_Recv(data_receive, sum_length_receive, MPI_CHAR,
                     md->SubDomain.rank_of_upper_subd[d], 6, md->MD_comm, &status);

            kreceive = k = 0;
            ITERATE(ic, ic_start_receive_upper, ic_stop_receive_upper)
            {
                for (int cc=0; cc < cells_length_receive[kreceive]; cc++)
                {
                    new_particle(&md->SubDomain.grid[ic[0]][ic[1]][ic[2]]);
                    for (dd=0; dd<3; dd++)
                        md->SubDomain.grid[ic[0]][ic[1]][ic[2]].parts[cc].core.x[dd] = data_receive[k].x[dd];
                    md->SubDomain.grid[ic[0]][ic[1]][ic[2]].parts[cc].core.atomkind = data_receive[k].atomkind;
                    md->SubDomain.grid[ic[0]][ic[1]][ic[2]].parts[cc].core.GroupID = data_receive[k].GroupID;
                    k++;
                }
                kreceive++;
            }

            free(cells_length_receive);
            free(data_receive);

            /* sending to upper process */
            cells_length_send = (int *)malloc((cells_num+1) * sizeof(int));
            k = sum_length_send = 0;
            ITERATE(ic, ic_start_send_upper, ic_stop_send_upper)
            {
                cells_length_send[k] = md->SubDomain.grid[ic[0]][ic[1]][ic[2]].parts_num;
                sum_length_send += cells_length_send[k++];
            }
            cells_length_send[cells_num] = sum_length_send;
            MPI_Send(cells_length_send, cells_num+1, MPI_INT, md->SubDomain.rank_of_upper_subd[d], 7,
                     md->MD_comm);
            free(cells_length_send);
            sum_length_send *= sizeof(position_struct_t);
            data_send = (position_struct_t *)malloc(sum_length_send);

            k=0;
            ITERATE(ic, ic_start_send_upper, ic_stop_send_upper)
            {
                for (int cc=0; cc < md->SubDomain.grid[ic[0]][ic[1]][ic[2]].parts_num; cc++)
                {
                    for (dd=0; dd<3; dd++)
                        data_send[k].x[dd] = md->SubDomain.grid[ic[0]][ic[1]][ic[2]].parts[cc].core.x[dd];
                    data_send[k].atomkind = md->SubDomain.grid[ic[0]][ic[1]][ic[2]].parts[cc].core.atomkind;
                    data_send[k].GroupID = md->SubDomain.grid[ic[0]][ic[1]][ic[2]].parts[cc].core.GroupID;
                    k++;
                }
            }

            MPI_Send(data_send, sum_length_send, MPI_CHAR, md->SubDomain.rank_of_upper_subd[d],
                     8, md->MD_comm);
            free(data_send);
        }
        else /* if (md->SubDomain.is[d] == md->ns[d]-1) */
        {
            /* sending to lower process */
            cells_num = 1;
            for (dd=0; dd<3; dd++)
                cells_num *= ic_stop_send_lower[dd] - ic_start_send_lower[dd];
            cells_length_send = (int *)malloc((cells_num+1) * sizeof(int));
            k = sum_length_send = 0;
            ITERATE(ic, ic_start_send_lower, ic_stop_send_lower)
            {
                cells_length_send[k] = md->SubDomain.grid[ic[0]][ic[1]][ic[2]].parts_num;
                sum_length_send += cells_length_send[k++];
            }
            cells_length_send[cells_num] = sum_length_send;
            MPI_Send(cells_length_send, cells_num+1, MPI_INT, md->SubDomain.rank_of_lower_subd[d], 5,
                     md->MD_comm);
            free(cells_length_send);
            sum_length_send *= sizeof(position_struct_t);
            data_send = (position_struct_t *)malloc(sum_length_send);

            k=0;
            ITERATE(ic, ic_start_send_lower, ic_stop_send_lower)
            {
                for (int cc=0; cc < md->SubDomain.grid[ic[0]][ic[1]][ic[2]].parts_num; cc++)
                {
                    for (dd=0; dd<3; dd++)
                        data_send[k].x[dd] = md->SubDomain.grid[ic[0]][ic[1]][ic[2]].parts[cc].core.x[dd];
                    data_send[k].atomkind = md->SubDomain.grid[ic[0]][ic[1]][ic[2]].parts[cc].core.atomkind;
                    data_send[k].GroupID = md->SubDomain.grid[ic[0]][ic[1]][ic[2]].parts[cc].core.GroupID;
                    k++;
                }
            }

            MPI_Send(data_send, sum_length_send, MPI_CHAR, md->SubDomain.rank_of_lower_subd[d],
                     6, md->MD_comm);
            free(data_send);

            /* receiving from lower process */
            cells_length_receive = (int *)malloc((cells_num+1) * sizeof(int));
            MPI_Recv(cells_length_receive, cells_num+1, MPI_INT, md->SubDomain.rank_of_lower_subd[d], 7,
                     md->MD_comm, &status);
            sum_length_receive = sizeof(position_struct_t) * cells_length_receive[cells_num];
            data_receive = (position_struct_t *)malloc(sum_length_receive);
            MPI_Recv(data_receive, sum_length_receive, MPI_CHAR,
                     md->SubDomain.rank_of_lower_subd[d], 8, md->MD_comm, &status);

            kreceive = k = 0;
            ITERATE(ic, ic_start_receive_lower, ic_stop_receive_lower)
            {
                for (int cc=0; cc < cells_length_receive[kreceive]; cc++)
                {
                    new_particle(&md->SubDomain.grid[ic[0]][ic[1]][ic[2]]);
                    for (dd=0; dd<3; dd++)
                        md->SubDomain.grid[ic[0]][ic[1]][ic[2]].parts[cc].core.x[dd] = data_receive[k].x[dd];
                    md->SubDomain.grid[ic[0]][ic[1]][ic[2]].parts[cc].core.atomkind = data_receive[k].atomkind;
                    md->SubDomain.grid[ic[0]][ic[1]][ic[2]].parts[cc].core.GroupID = data_receive[k].GroupID;
                    k++;
                }
                kreceive++;
            }

            free(cells_length_receive);
            free(data_receive);
        }
    }
    else
    {
        /* sending to lower process, receiving from upper process */
        cells_num = 1;
        for (dd=0; dd<3; dd++)
            cells_num *= ic_stop_send_lower[dd] - ic_start_send_lower[dd];
        cells_length_send = (int *)malloc((cells_num+1) * sizeof(int));
        cells_length_receive = (int *)malloc((cells_num+1) * sizeof(int));
        k = sum_length_send = 0;
        ITERATE(ic, ic_start_send_lower, ic_stop_send_lower)
        {
            cells_length_send[k] = md->SubDomain.grid[ic[0]][ic[1]][ic[2]].parts_num;
            sum_length_send += cells_length_send[k++];
        }
        cells_length_send[cells_num] = sum_length_send;
        MPI_Isend(cells_length_send, cells_num+1, MPI_INT, md->SubDomain.rank_of_lower_subd[d], 5,
                  md->MD_comm, &request);
        MPI_Recv(cells_length_receive, cells_num+1, MPI_INT, md->SubDomain.rank_of_upper_subd[d], 5,
                 md->MD_comm, &status);
        MPI_Wait(&request, &status);
        sum_length_send *= sizeof(position_struct_t);
        data_send = (position_struct_t *)malloc(sum_length_send);
        sum_length_receive = sizeof(position_struct_t) * cells_length_receive[cells_num];
        data_receive = (position_struct_t *)malloc(sum_length_receive);

        k=0;
        ITERATE(ic, ic_start_send_lower, ic_stop_send_lower)
        {
            for (int cc=0; cc < md->SubDomain.grid[ic[0]][ic[1]][ic[2]].parts_num; cc++)
            {
                for (dd=0; dd<3; dd++)
                    data_send[k].x[dd] = md->SubDomain.grid[ic[0]][ic[1]][ic[2]].parts[cc].core.x[dd];
                data_send[k].atomkind = md->SubDomain.grid[ic[0]][ic[1]][ic[2]].parts[cc].core.atomkind;
                data_send[k].GroupID = md->SubDomain.grid[ic[0]][ic[1]][ic[2]].parts[cc].core.GroupID;
                k++;
            }
        }

        MPI_Isend(data_send, sum_length_send, MPI_CHAR, md->SubDomain.rank_of_lower_subd[d], 6,
                  md->MD_comm, &request);
        MPI_Recv(data_receive, sum_length_receive, MPI_CHAR,
                 md->SubDomain.rank_of_upper_subd[d], 6, md->MD_comm, &status);
        MPI_Wait(&request, &status);
        free(data_send);

        kreceive = k = 0;
        ITERATE(ic, ic_start_receive_upper, ic_stop_receive_upper)
        {
            for (int cc=0; cc < cells_length_receive[kreceive]; cc++)
            {
                new_particle(&md->SubDomain.grid[ic[0]][ic[1]][ic[2]]);
                if (md->SubDomain.is[d] == md->ns[d] - 1)
                    data_receive[k].x[d] += md->l[d];
                for (dd=0; dd<3; dd++)
                    md->SubDomain.grid[ic[0]][ic[1]][ic[2]].parts[cc].core.x[dd] = data_receive[k].x[dd];
                md->SubDomain.grid[ic[0]][ic[1]][ic[2]].parts[cc].core.atomkind = data_receive[k].atomkind;
                md->SubDomain.grid[ic[0]][ic[1]][ic[2]].parts[cc].core.GroupID = data_receive[k].GroupID;
                k++;
            }
            kreceive++;
        }

        free(data_receive);

        /* sending to upper process, receiving from lower process */
        k = sum_length_send = 0;
        ITERATE(ic, ic_start_send_upper, ic_stop_send_upper)
        {
            cells_length_send[k] = md->SubDomain.grid[ic[0]][ic[1]][ic[2]].parts_num;
            sum_length_send += cells_length_send[k++];
        }
        cells_length_send[cells_num] = sum_length_send;
        MPI_Isend(cells_length_send, cells_num+1, MPI_INT, md->SubDomain.rank_of_upper_subd[d], 7,
                  md->MD_comm, &request);
        MPI_Recv(cells_length_receive, cells_num+1, MPI_INT, md->SubDomain.rank_of_lower_subd[d], 7,
                 md->MD_comm, &status);
        MPI_Wait(&request, &status);
        free(cells_length_send);
        sum_length_send *= sizeof(position_struct_t);
        data_send = (position_struct_t *)malloc(sum_length_send);
        sum_length_receive = sizeof(position_struct_t) * cells_length_receive[cells_num];
        data_receive = (position_struct_t *)malloc(sum_length_receive);

        k=0;
        ITERATE(ic, ic_start_send_upper, ic_stop_send_upper)
        {
            for (int cc=0; cc < md->SubDomain.grid[ic[0]][ic[1]][ic[2]].parts_num; cc++)
            {
                for (dd=0; dd<3; dd++)
                    data_send[k].x[dd] = md->SubDomain.grid[ic[0]][ic[1]][ic[2]].parts[cc].core.x[dd];
                data_send[k].atomkind = md->SubDomain.grid[ic[0]][ic[1]][ic[2]].parts[cc].core.atomkind;
                data_send[k].GroupID = md->SubDomain.grid[ic[0]][ic[1]][ic[2]].parts[cc].core.GroupID;
                k++;
            }
        }

        MPI_Isend(data_send, sum_length_send, MPI_CHAR, md->SubDomain.rank_of_upper_subd[d], 8,
                  md->MD_comm, &request);
        MPI_Recv(data_receive, sum_length_receive, MPI_CHAR,
                 md->SubDomain.rank_of_lower_subd[d], 8, md->MD_comm, &status);
        MPI_Wait(&request, &status);
        free(data_send);

        kreceive = k = 0;
        ITERATE(ic, ic_start_receive_lower, ic_stop_receive_lower)
        {
            for (int cc=0; cc < cells_length_receive[kreceive]; cc++)
            {
                new_particle(&md->SubDomain.grid[ic[0]][ic[1]][ic[2]]);
                if (md->SubDomain.is[d] == 0)
                    data_receive[k].x[d] -= md->l[d];
                for (dd=0; dd<3; dd++)
                    md->SubDomain.grid[ic[0]][ic[1]][ic[2]].parts[cc].core.x[dd] = data_receive[k].x[dd];
                md->SubDomain.grid[ic[0]][ic[1]][ic[2]].parts[cc].core.atomkind = data_receive[k].atomkind;
                md->SubDomain.grid[ic[0]][ic[1]][ic[2]].parts[cc].core.GroupID = data_receive[k].GroupID;
                k++;
            }
            kreceive++;
        }

        free(cells_length_receive);
        free(data_receive);
    }
}

static void ghostparticles_update_Fprime_in_direction_d(
    fmd_t *md, int d, fmd_ituple_t ic_start_send_lower,
    fmd_ituple_t ic_stop_send_lower, fmd_ituple_t ic_start_receive_lower,
    fmd_ituple_t ic_stop_receive_lower, fmd_ituple_t ic_start_send_upper,
    fmd_ituple_t ic_stop_send_upper, fmd_ituple_t ic_start_receive_upper,
    fmd_ituple_t ic_stop_receive_upper)
{
    MPI_Status status;
    MPI_Request request;
    int sum_length_send, sum_length_receive;
    fmd_ituple_t ic;
    fmd_real_t *data_receive, *data_send;
    int k;

    if ( ((md->SubDomain.is[d] == 0) || (md->SubDomain.is[d] == md->ns[d]-1)) && !md->PBC[d] )
    {
        if (md->SubDomain.is[d] == 0)
        {
            /* receiving from upper process */
            MPI_Recv(&sum_length_receive, 1, MPI_INT, md->SubDomain.rank_of_upper_subd[d], 100,
                     md->MD_comm, &status);
            data_receive = (fmd_real_t *)malloc(sum_length_receive * sizeof(fmd_real_t));
            MPI_Recv(data_receive, sum_length_receive, FMD_MPI_REAL,
                     md->SubDomain.rank_of_upper_subd[d], 101, md->MD_comm, &status);

            k = 0;
            ITERATE(ic, ic_start_receive_upper, ic_stop_receive_upper)
                for (int cc=0; cc < md->SubDomain.grid[ic[0]][ic[1]][ic[2]].parts_num; cc++)
                    md->SubDomain.grid[ic[0]][ic[1]][ic[2]].parts[cc].FembPrime = data_receive[k++];

            free(data_receive);

            /* sending to upper process */
            sum_length_send = 0;
            ITERATE(ic, ic_start_send_upper, ic_stop_send_upper)
                sum_length_send += md->SubDomain.grid[ic[0]][ic[1]][ic[2]].parts_num;

            MPI_Send(&sum_length_send, 1, MPI_INT, md->SubDomain.rank_of_upper_subd[d], 102,
                     md->MD_comm);
            data_send = (fmd_real_t *)malloc(sum_length_send * sizeof(fmd_real_t));

            k = 0;
            ITERATE(ic, ic_start_send_upper, ic_stop_send_upper)
                for (int cc=0; cc < md->SubDomain.grid[ic[0]][ic[1]][ic[2]].parts_num; cc++)
                    data_send[k++] = md->SubDomain.grid[ic[0]][ic[1]][ic[2]].parts[cc].FembPrime;

            MPI_Send(data_send, sum_length_send, FMD_MPI_REAL, md->SubDomain.rank_of_upper_subd[d],
                     103, md->MD_comm);
            free(data_send);
        }
        else  /* if (md->SubDomain.is[d] == md->ns[d]-1) */
        {
            /* sending to lower process */
            sum_length_send = 0;
            ITERATE(ic, ic_start_send_lower, ic_stop_send_lower)
                sum_length_send += md->SubDomain.grid[ic[0]][ic[1]][ic[2]].parts_num;
            MPI_Send(&sum_length_send, 1, MPI_INT, md->SubDomain.rank_of_lower_subd[d], 100,
                     md->MD_comm);
            data_send = (fmd_real_t *)malloc(sum_length_send * sizeof(fmd_real_t));

            k = 0;
            ITERATE(ic, ic_start_send_lower, ic_stop_send_lower)
                for (int cc=0; cc < md->SubDomain.grid[ic[0]][ic[1]][ic[2]].parts_num; cc++)
                    data_send[k++] = md->SubDomain.grid[ic[0]][ic[1]][ic[2]].parts[cc].FembPrime;

            MPI_Send(data_send, sum_length_send, FMD_MPI_REAL, md->SubDomain.rank_of_lower_subd[d],
                     101, md->MD_comm);
            free(data_send);

            /* receiving from lower process */
            MPI_Recv(&sum_length_receive, 1, MPI_INT, md->SubDomain.rank_of_lower_subd[d], 102,
                     md->MD_comm, &status);
            data_receive = (fmd_real_t *)malloc(sum_length_receive * sizeof(fmd_real_t));
            MPI_Recv(data_receive, sum_length_receive, FMD_MPI_REAL,
                     md->SubDomain.rank_of_lower_subd[d], 103, md->MD_comm, &status);

            k = 0;
            ITERATE(ic, ic_start_receive_lower, ic_stop_receive_lower)
                for (int cc=0; cc < md->SubDomain.grid[ic[0]][ic[1]][ic[2]].parts_num; cc++)
                    md->SubDomain.grid[ic[0]][ic[1]][ic[2]].parts[cc].FembPrime = data_receive[k++];

            free(data_receive);
        }
    }
    else
    {
        /* sending to lower process, receiving from upper process */
        sum_length_send = 0;
        ITERATE(ic, ic_start_send_lower, ic_stop_send_lower)
            sum_length_send += md->SubDomain.grid[ic[0]][ic[1]][ic[2]].parts_num;
        MPI_Isend(&sum_length_send, 1, MPI_INT, md->SubDomain.rank_of_lower_subd[d], 100,
                  md->MD_comm, &request);
        MPI_Recv(&sum_length_receive, 1, MPI_INT, md->SubDomain.rank_of_upper_subd[d], 100,
                 md->MD_comm, &status);
        MPI_Wait(&request, &status);
        data_send = (fmd_real_t *)malloc(sum_length_send * sizeof(fmd_real_t));
        data_receive = (fmd_real_t *)malloc(sum_length_receive * sizeof(fmd_real_t));

        k = 0;
        ITERATE(ic, ic_start_send_lower, ic_stop_send_lower)
            for (int cc=0; cc < md->SubDomain.grid[ic[0]][ic[1]][ic[2]].parts_num; cc++)
                data_send[k++] = md->SubDomain.grid[ic[0]][ic[1]][ic[2]].parts[cc].FembPrime;

        MPI_Isend(data_send, sum_length_send, FMD_MPI_REAL, md->SubDomain.rank_of_lower_subd[d], 101,
                  md->MD_comm, &request);
        MPI_Recv(data_receive, sum_length_receive, FMD_MPI_REAL, md->SubDomain.rank_of_upper_subd[d], 101,
                 md->MD_comm, &status);
        MPI_Wait(&request, &status);
        free(data_send);

        k = 0;
        ITERATE(ic, ic_start_receive_upper, ic_stop_receive_upper)
            for (int cc=0; cc < md->SubDomain.grid[ic[0]][ic[1]][ic[2]].parts_num; cc++)
                md->SubDomain.grid[ic[0]][ic[1]][ic[2]].parts[cc].FembPrime = data_receive[k++];

        free(data_receive);

        /* sending to upper process, receiving from lower process */
        sum_length_send = 0;
        ITERATE(ic, ic_start_send_upper, ic_stop_send_upper)
            sum_length_send += md->SubDomain.grid[ic[0]][ic[1]][ic[2]].parts_num;
        MPI_Isend(&sum_length_send, 1, MPI_INT, md->SubDomain.rank_of_upper_subd[d], 102,
                  md->MD_comm, &request);
        MPI_Recv(&sum_length_receive, 1, MPI_INT, md->SubDomain.rank_of_lower_subd[d], 102,
                 md->MD_comm, &status);
        MPI_Wait(&request, &status);
        data_send = (fmd_real_t *)malloc(sum_length_send * sizeof(fmd_real_t));
        data_receive = (fmd_real_t *)malloc(sum_length_receive * sizeof(fmd_real_t));

        k = 0;
        ITERATE(ic, ic_start_send_upper, ic_stop_send_upper)
            for (int cc=0; cc < md->SubDomain.grid[ic[0]][ic[1]][ic[2]].parts_num; cc++)
                data_send[k++] = md->SubDomain.grid[ic[0]][ic[1]][ic[2]].parts[cc].FembPrime;

        MPI_Isend(data_send, sum_length_send, FMD_MPI_REAL, md->SubDomain.rank_of_upper_subd[d], 103,
                  md->MD_comm, &request);
        MPI_Recv(data_receive, sum_length_receive, FMD_MPI_REAL, md->SubDomain.rank_of_lower_subd[d], 103,
                 md->MD_comm, &status);
        MPI_Wait(&request, &status);
        free(data_send);

        k = 0;
        ITERATE(ic, ic_start_receive_lower, ic_stop_receive_lower)
            for (int cc=0; cc < md->SubDomain.grid[ic[0]][ic[1]][ic[2]].parts_num; cc++)
                md->SubDomain.grid[ic[0]][ic[1]][ic[2]].parts[cc].FembPrime = data_receive[k++];

        free(data_receive);
    }
}

static void ghostparticles_update_LocOrdParam_in_direction_d(
    fmd_t *md, int d, fmd_ituple_t ic_start_send_lower,
    fmd_ituple_t ic_stop_send_lower, fmd_ituple_t ic_start_receive_lower,
    fmd_ituple_t ic_stop_receive_lower, fmd_ituple_t ic_start_send_upper,
    fmd_ituple_t ic_stop_send_upper, fmd_ituple_t ic_start_receive_upper,
    fmd_ituple_t ic_stop_receive_upper)
{
    MPI_Status status;
    MPI_Request request;
    int sum_length_send, sum_length_receive;
    fmd_ituple_t ic;
    float *data_receive, *data_send;
    int k;

    if ( ((md->SubDomain.is[d] == 0) || (md->SubDomain.is[d] == md->ns[d]-1)) && !md->PBC[d] )
    {
        if (md->SubDomain.is[d] == 0)
        {
            /* receiving from upper process */
            MPI_Recv(&sum_length_receive, 1, MPI_INT, md->SubDomain.rank_of_upper_subd[d], 130,
                     md->MD_comm, &status);
            data_receive = (float *)malloc(sum_length_receive * sizeof(float));
            MPI_Recv(data_receive, sum_length_receive, MPI_FLOAT,
                     md->SubDomain.rank_of_upper_subd[d], 131, md->MD_comm, &status);

            k = 0;
            ITERATE(ic, ic_start_receive_upper, ic_stop_receive_upper)
                for (int cc=0; cc < md->SubDomain.grid[ic[0]][ic[1]][ic[2]].parts_num; cc++)
                    md->SubDomain.grid[ic[0]][ic[1]][ic[2]].parts[cc].core.LocOrdParam = data_receive[k++];

            free(data_receive);

            /* sending to upper process */
            sum_length_send = 0;
            ITERATE(ic, ic_start_send_upper, ic_stop_send_upper)
                sum_length_send += md->SubDomain.grid[ic[0]][ic[1]][ic[2]].parts_num;
            MPI_Send(&sum_length_send, 1, MPI_INT, md->SubDomain.rank_of_upper_subd[d], 132,
                     md->MD_comm);
            data_send = (float *)malloc(sum_length_send * sizeof(float));

            k = 0;
            ITERATE(ic, ic_start_send_upper, ic_stop_send_upper)
                for (int cc=0; cc < md->SubDomain.grid[ic[0]][ic[1]][ic[2]].parts_num; cc++)
                    data_send[k++] = md->SubDomain.grid[ic[0]][ic[1]][ic[2]].parts[cc].core.LocOrdParam;

            MPI_Send(data_send, sum_length_send, MPI_FLOAT, md->SubDomain.rank_of_upper_subd[d],
                     133, md->MD_comm);
            free(data_send);
        }
        else  /* if (md->SubDomain.is[d] == md->ns[d]-1) */
        {
            /* sending to lower process */
            sum_length_send = 0;
            ITERATE(ic, ic_start_send_lower, ic_stop_send_lower)
                sum_length_send += md->SubDomain.grid[ic[0]][ic[1]][ic[2]].parts_num;
            MPI_Send(&sum_length_send, 1, MPI_INT, md->SubDomain.rank_of_lower_subd[d], 130,
                     md->MD_comm);
            data_send = (float *)malloc(sum_length_send * sizeof(float));

            k = 0;
            ITERATE(ic, ic_start_send_lower, ic_stop_send_lower)
                for (int cc=0; cc < md->SubDomain.grid[ic[0]][ic[1]][ic[2]].parts_num; cc++)
                    data_send[k++] = md->SubDomain.grid[ic[0]][ic[1]][ic[2]].parts[cc].core.LocOrdParam;

            MPI_Send(data_send, sum_length_send, MPI_FLOAT, md->SubDomain.rank_of_lower_subd[d],
                     131, md->MD_comm);
            free(data_send);

            /* receiving from lower process */
            MPI_Recv(&sum_length_receive, 1, MPI_INT, md->SubDomain.rank_of_lower_subd[d], 132,
                     md->MD_comm, &status);
            data_receive = (float *)malloc(sum_length_receive * sizeof(float));
            MPI_Recv(data_receive, sum_length_receive, MPI_FLOAT,
                     md->SubDomain.rank_of_lower_subd[d], 133, md->MD_comm, &status);

            k = 0;
            ITERATE(ic, ic_start_receive_lower, ic_stop_receive_lower)
                for (int cc=0; cc < md->SubDomain.grid[ic[0]][ic[1]][ic[2]].parts_num; cc++)
                    md->SubDomain.grid[ic[0]][ic[1]][ic[2]].parts[cc].core.LocOrdParam = data_receive[k++];

            free(data_receive);
        }
    }
    else
    {
        /* sending to lower process, receiving from upper process */
        sum_length_send = 0;
        ITERATE(ic, ic_start_send_lower, ic_stop_send_lower)
            sum_length_send += md->SubDomain.grid[ic[0]][ic[1]][ic[2]].parts_num;
        MPI_Isend(&sum_length_send, 1, MPI_INT, md->SubDomain.rank_of_lower_subd[d], 130,
                  md->MD_comm, &request);
        MPI_Recv(&sum_length_receive, 1, MPI_INT, md->SubDomain.rank_of_upper_subd[d], 130,
                 md->MD_comm, &status);
        MPI_Wait(&request, &status);
        data_send = (float *)malloc(sum_length_send * sizeof(float));
        data_receive = (float *)malloc(sum_length_receive * sizeof(float));

        k = 0;
        ITERATE(ic, ic_start_send_lower, ic_stop_send_lower)
            for (int cc=0; cc < md->SubDomain.grid[ic[0]][ic[1]][ic[2]].parts_num; cc++)
                data_send[k++] = md->SubDomain.grid[ic[0]][ic[1]][ic[2]].parts[cc].core.LocOrdParam;

        MPI_Isend(data_send, sum_length_send, MPI_FLOAT, md->SubDomain.rank_of_lower_subd[d], 131,
                  md->MD_comm, &request);
        MPI_Recv(data_receive, sum_length_receive, MPI_FLOAT, md->SubDomain.rank_of_upper_subd[d], 131,
                 md->MD_comm, &status);
        MPI_Wait(&request, &status);
        free(data_send);

        k = 0;
        ITERATE(ic, ic_start_receive_upper, ic_stop_receive_upper)
            for (int cc=0; cc < md->SubDomain.grid[ic[0]][ic[1]][ic[2]].parts_num; cc++)
                md->SubDomain.grid[ic[0]][ic[1]][ic[2]].parts[cc].core.LocOrdParam = data_receive[k++];

        free(data_receive);

        /* sending to upper process, receiving from lower process */
        sum_length_send = 0;
        ITERATE(ic, ic_start_send_upper, ic_stop_send_upper)
            sum_length_send += md->SubDomain.grid[ic[0]][ic[1]][ic[2]].parts_num;
        MPI_Isend(&sum_length_send, 1, MPI_INT, md->SubDomain.rank_of_upper_subd[d], 132,
                  md->MD_comm, &request);
        MPI_Recv(&sum_length_receive, 1, MPI_INT, md->SubDomain.rank_of_lower_subd[d], 132,
                 md->MD_comm, &status);
        MPI_Wait(&request, &status);
        data_send = (float *)malloc(sum_length_send * sizeof(float));
        data_receive = (float *)malloc(sum_length_receive * sizeof(float));

        k = 0;
        ITERATE(ic, ic_start_send_upper, ic_stop_send_upper)
            for (int cc=0; cc < md->SubDomain.grid[ic[0]][ic[1]][ic[2]].parts_num; cc++)
                data_send[k++] = md->SubDomain.grid[ic[0]][ic[1]][ic[2]].parts[cc].core.LocOrdParam;

        MPI_Isend(data_send, sum_length_send, MPI_FLOAT, md->SubDomain.rank_of_upper_subd[d], 133,
                  md->MD_comm, &request);
        MPI_Recv(data_receive, sum_length_receive, MPI_FLOAT, md->SubDomain.rank_of_lower_subd[d], 133,
                 md->MD_comm, &status);
        MPI_Wait(&request, &status);
        free(data_send);

        k = 0;
        ITERATE(ic, ic_start_receive_lower, ic_stop_receive_lower)
            for (int cc=0; cc < md->SubDomain.grid[ic[0]][ic[1]][ic[2]].parts_num; cc++)
                md->SubDomain.grid[ic[0]][ic[1]][ic[2]].parts[cc].core.LocOrdParam = data_receive[k++];

        free(data_receive);
    }
}

static void particles_prepare_migration_in_direction_d(
    SubDomain_t *s_p, int d, fmd_ituple_t ic_start_send_lower,
    fmd_ituple_t ic_stop_send_lower, fmd_ituple_t ic_start_receive_lower,
    fmd_ituple_t ic_stop_receive_lower, fmd_ituple_t ic_start_send_upper,
    fmd_ituple_t ic_stop_send_upper, fmd_ituple_t ic_start_receive_upper,
    fmd_ituple_t ic_stop_receive_upper)
{
    int dd;

    for (dd=0; dd<3; dd++)
    {
        if (dd == d) // only ghost
        {
            ic_start_send_lower[dd] = s_p->ic_start[dd];
            ic_stop_send_lower[dd] = ic_start_send_lower[dd] + s_p->ic_start[dd];
            ic_start_receive_lower[dd] = 0;
            ic_stop_receive_lower[dd] = s_p->ic_start[dd];
            ic_stop_send_upper[dd] = s_p->ic_stop[dd];
            ic_start_send_upper[dd] = ic_stop_send_upper[dd] - s_p->ic_start[dd];
            ic_start_receive_upper[dd] = s_p->ic_stop[dd];
            ic_stop_receive_upper[dd] = s_p->cell_num[dd];
        }
        else if (dd > d) // including ghost
        {
            ic_start_receive_lower[dd] = ic_start_send_lower[dd]
                                       = ic_start_receive_upper[dd]
                                       = ic_start_send_upper[dd] = 0;
            ic_stop_receive_lower[dd]  = ic_stop_send_lower[dd] = ic_stop_receive_upper[dd]
                                       = ic_stop_send_upper[dd] = s_p->cell_num[dd];
        }
        else // excluding ghost
        {
            ic_start_receive_lower[dd] = ic_start_send_lower[dd]
                                       = ic_start_receive_upper[dd]
                                       = ic_start_send_upper[dd] = s_p->ic_start[dd];
            ic_stop_receive_lower[dd]  = ic_stop_send_lower[dd] = ic_stop_receive_upper[dd]
                                       = ic_stop_send_upper[dd] = s_p->ic_stop[dd];
        }
    }
}

static void ghostparticles_prepare_init_update_in_direction_d(
    fmd_t *md, int d, fmd_ituple_t ic_start_send_lower,
    fmd_ituple_t ic_stop_send_lower, fmd_ituple_t ic_start_receive_lower,
    fmd_ituple_t ic_stop_receive_lower, fmd_ituple_t ic_start_send_upper,
    fmd_ituple_t ic_stop_send_upper, fmd_ituple_t ic_start_receive_upper,
    fmd_ituple_t ic_stop_receive_upper)
{
    int dd;

    for (dd=0; dd<3; dd++)
    {
        if (dd == d) // only ghost
        {
            ic_start_send_lower[dd] = md->SubDomain.ic_start[dd];
            ic_stop_send_lower[dd] = ic_start_send_lower[dd] + 1;
            ic_stop_receive_lower[dd] = md->SubDomain.ic_start[dd];
            ic_start_receive_lower[dd] = ic_stop_receive_lower[dd] - 1;
            ic_stop_send_upper[dd] = md->SubDomain.ic_stop[dd];
            ic_start_send_upper[dd] = ic_stop_send_upper[dd] - 1;
            ic_start_receive_upper[dd] = md->SubDomain.ic_stop[dd];
            ic_stop_receive_upper[dd] = ic_start_receive_upper[dd] + 1;
        }
        else if (dd > d) // including ghost
        {
            ic_start_receive_lower[dd] = ic_start_send_lower[dd]
                                       = ic_start_receive_upper[dd]
                                       = ic_start_send_upper[dd]
                                       = md->SubDomain.ic_start[dd] - (md->ns[dd] == 1 ? 0 : 1);
            ic_stop_receive_lower[dd]  = ic_stop_send_lower[dd]
                                       = ic_stop_receive_upper[dd]
                                       = ic_stop_send_upper[dd]
                                       = md->SubDomain.ic_stop[dd] + (md->ns[dd] == 1 ? 0 : 1);
        }
        else // excluding ghost
        {
            ic_start_receive_lower[dd] = ic_start_send_lower[dd]
                                       = ic_start_receive_upper[dd]
                                       = ic_start_send_upper[dd]
                                       = md->SubDomain.ic_start[dd];
            ic_stop_receive_lower[dd]  = ic_stop_send_lower[dd]
                                       = ic_stop_receive_upper[dd]
                                       = ic_stop_send_upper[dd]
                                       = md->SubDomain.ic_stop[dd];
        }
    }
}

void fmd_ghostparticles_update_LocOrdParam(fmd_t *md)
{
    int d;
    fmd_ituple_t ic_start_send_lower, ic_stop_send_lower;
    fmd_ituple_t ic_start_send_upper, ic_stop_send_upper;
    fmd_ituple_t ic_start_receive_lower, ic_stop_receive_lower;
    fmd_ituple_t ic_start_receive_upper, ic_stop_receive_upper;

    for (d = 3-1; d >= 0; d--)
    {
        if (md->ns[d] != 1)
        {
            ghostparticles_prepare_init_update_in_direction_d(
                md, d, ic_start_send_lower, ic_stop_send_lower,
                ic_start_receive_lower, ic_stop_receive_lower,
                ic_start_send_upper, ic_stop_send_upper, ic_start_receive_upper,
                ic_stop_receive_upper);
            ghostparticles_update_LocOrdParam_in_direction_d(
                md, d, ic_start_send_lower, ic_stop_send_lower,
                ic_start_receive_lower, ic_stop_receive_lower,
                ic_start_send_upper, ic_stop_send_upper, ic_start_receive_upper,
                ic_stop_receive_upper);
        }
    }
}

void fmd_ghostparticles_update_Femb(fmd_t *md)
{
    int d;
    fmd_ituple_t ic_start_send_lower, ic_stop_send_lower;
    fmd_ituple_t ic_start_send_upper, ic_stop_send_upper;
    fmd_ituple_t ic_start_receive_lower, ic_stop_receive_lower;
    fmd_ituple_t ic_start_receive_upper, ic_stop_receive_upper;

    for (d = 3-1; d >= 0; d--)
    {
        if (md->ns[d] != 1)
        {
            ghostparticles_prepare_init_update_in_direction_d(
                md, d, ic_start_send_lower, ic_stop_send_lower,
                ic_start_receive_lower, ic_stop_receive_lower,
                ic_start_send_upper, ic_stop_send_upper, ic_start_receive_upper,
                ic_stop_receive_upper);
            ghostparticles_update_Fprime_in_direction_d(
                md, d, ic_start_send_lower, ic_stop_send_lower,
                ic_start_receive_lower, ic_stop_receive_lower,
                ic_start_send_upper, ic_stop_send_upper, ic_start_receive_upper,
                ic_stop_receive_upper);
        }
    }
}

void fmd_ghostparticles_init(fmd_t *md)
{
    int d;
    fmd_ituple_t ic_start_send_lower, ic_stop_send_lower;
    fmd_ituple_t ic_start_send_upper, ic_stop_send_upper;
    fmd_ituple_t ic_start_receive_lower, ic_stop_receive_lower;
    fmd_ituple_t ic_start_receive_upper, ic_stop_receive_upper;

    for (d = 3-1; d >= 0; d--)
    {
        if (md->ns[d] != 1)
        {
            ghostparticles_prepare_init_update_in_direction_d(
                md, d, ic_start_send_lower, ic_stop_send_lower,
                ic_start_receive_lower, ic_stop_receive_lower,
                ic_start_send_upper, ic_stop_send_upper, ic_start_receive_upper,
                ic_stop_receive_upper);
            ghostparticles_init_in_direction_d(
                md, d, ic_start_send_lower, ic_stop_send_lower,
                ic_start_receive_lower, ic_stop_receive_lower,
                ic_start_send_upper, ic_stop_send_upper, ic_start_receive_upper,
                ic_stop_receive_upper);
        }
    }
}

void fmd_ghostparticles_delete(fmd_t *md)
{
    int d;
    fmd_ituple_t ic_from, ic_to;
    fmd_ituple_t jc;

    for (d=0; d<3; d++)
    {
        ic_from[d] = md->SubDomain.ic_start[d] - (md->ns[d] == 1 ? 0 : 1);
        jc[d] = ic_to[d] = md->SubDomain.ic_stop[d] + (md->ns[d] == 1 ? 0 : 1);
    }

    for (d=0; d<3; d++)
    {
        jc[d] = md->SubDomain.ic_start[d];
        _fmd_cleanGridSegment(md->SubDomain.grid, ic_from, jc);
        ic_from[d] = md->SubDomain.ic_stop[d];
        _fmd_cleanGridSegment(md->SubDomain.grid, ic_from, ic_to);
        ic_from[d] = md->SubDomain.ic_start[d];
        jc[d] = ic_to[d] = md->SubDomain.ic_stop[d];
    }
}

void fmd_particles_migrate(fmd_t *md)
{
    fmd_ituple_t ic_start_send_lower, ic_stop_send_lower;
    fmd_ituple_t ic_start_send_upper, ic_stop_send_upper;
    fmd_ituple_t ic_start_receive_lower, ic_stop_receive_lower;
    fmd_ituple_t ic_start_receive_upper, ic_stop_receive_upper;
    int d;

    for (d = 0; d < 3; d++)
    {
        if (md->ns[d] != 1)
        {
            particles_prepare_migration_in_direction_d(
                &md->SubDomain, d,
                ic_start_receive_lower, ic_stop_receive_lower,
                ic_start_send_lower, ic_stop_send_lower, ic_start_receive_upper,
                ic_stop_receive_upper, ic_start_send_upper, ic_stop_send_upper);
            particles_migrate_in_direction_d(
                md, d, ic_start_send_lower, ic_stop_send_lower,
                ic_start_receive_lower, ic_stop_receive_lower,
                ic_start_send_upper, ic_stop_send_upper, ic_start_receive_upper,
                ic_stop_receive_upper);
        }
    }
}

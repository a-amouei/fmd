/*
  h5.c: This file is part of Free Molecular Dynamics

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

#include "h5.h"
#include "types.h"
#include "array.h"
#include "turi.h"
#include "general.h"
#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include "base.h"

void _fmd_h5_ds_init(h5_dataspaces_t *ds)
{
    hsize_t sz;

    ds->ds_scalar = H5Screate(H5S_SCALAR);
    assert(ds->ds_scalar >= 0);

    sz = 3;
    ds->ds_simple_1_3 = H5Screate_simple(1, &sz, NULL);
    assert(ds->ds_simple_1_3 >= 0);
}

void _fmd_h5_ds_free(h5_dataspaces_t *ds)
{
    herr_t status;

    status = H5Sclose(ds->ds_scalar);
    assert(status >= 0);

    status = H5Sclose(ds->ds_simple_1_3);
    assert(status >= 0);
}

static void set_attr_ftuple(h5_dataspaces_t *ds, hid_t obj, const char *name, float *f)
{
    hid_t attr;
    herr_t status;

    attr = H5Acreate(obj, name, H5T_IEEE_F32LE, ds->ds_simple_1_3, H5P_DEFAULT, H5P_DEFAULT);
    assert(attr >= 0);
    status = H5Awrite(attr, H5T_NATIVE_FLOAT, f);
    assert(status >= 0);

    status = H5Aclose(attr);
    assert(status >= 0);
}

static void set_attr_ituple(h5_dataspaces_t *ds, hid_t obj, const char *name, fmd_ituple_t i)
{
    hid_t attr;
    herr_t status;

    attr = H5Acreate(obj, name, H5T_STD_I32LE, ds->ds_simple_1_3, H5P_DEFAULT, H5P_DEFAULT);
    assert(attr >= 0);
    status = H5Awrite(attr, H5T_NATIVE_INT, i);
    assert(status >= 0);

    status = H5Aclose(attr);
    assert(status >= 0);
}

static void set_attr_string(h5_dataspaces_t *ds, hid_t obj, const char *name, const char *value)
{
    hid_t str_type_mem, str_type_file, attr;
    herr_t status;

    hsize_t len = strlen(value);

    str_type_mem = H5Tcopy(H5T_C_S1);
    assert(str_type_mem >= 0);
    status = H5Tset_size(str_type_mem, len+1);
    assert(status >= 0);

    str_type_file = H5Tcopy(H5T_FORTRAN_S1);
    assert(str_type_file >= 0);
    status = H5Tset_size(str_type_file, len);
    assert(status >= 0);

    attr = H5Acreate(obj, name, str_type_file, ds->ds_scalar, H5P_DEFAULT, H5P_DEFAULT);
    assert(attr >= 0);
    status = H5Awrite(attr, str_type_mem, value);
    assert(status >= 0);

    status = H5Aclose(attr);
    assert(status >= 0);
    status = H5Tclose(str_type_mem);
    assert(status >= 0);
    status = H5Tclose(str_type_file);
    assert(status >= 0);
}

static void create_mesh(h5_dataspaces_t *ds, hid_t parent, fmd_ituple_t dims, fmd_ftriple_t UpperBounds)
{
    hid_t group_id;
    herr_t status;

    group_id = H5Gcreate(parent, "meshA", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    set_attr_string(ds, group_id, "vsType", "mesh");
    set_attr_string(ds, group_id, "vsKind", "uniform");
    set_attr_ituple(ds, group_id, "vsStartCell", _fmd_ThreeZeros_int);
    set_attr_ituple(ds, group_id, "vsNumCells", dims);
    set_attr_ftuple(ds, group_id, "vsLowerBounds", _fmd_ThreeZeros_float);
    set_attr_ftuple(ds, group_id, "vsUpperBounds", UpperBounds);

    status = H5Gclose(group_id);
    assert(status >= 0);
}

void _fmd_h5_save_scalar_field_float(fmd_t *md, fmd_string_t fieldname, turi_t *t, fmd_string_t path, fmd_array3D_t *arr)
{
    hid_t file_id, dataset_id, dataspace_id;
    herr_t status;

    file_id = H5Fcreate(path, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    assert(file_id >= 0);

    fmd_ftriple_t ubound = {md->l[0], md->l[1], md->l[2]};
    create_mesh(&md->h5_dataspaces, file_id, t->tdims_global, ubound);

    hsize_t s[3];
    s[0] = t->tdims_global[0];
    s[1] = t->tdims_global[1];
    s[2] = t->tdims_global[2];

    dataspace_id = H5Screate_simple(3, s, NULL);
    assert(dataspace_id >= 0);

    dataset_id = H5Dcreate1(file_id, fieldname, H5T_IEEE_F32LE, dataspace_id, H5P_DEFAULT);
    assert(dataset_id >= 0);

    set_attr_string(&md->h5_dataspaces, dataset_id, "vsType", "variable");
    set_attr_string(&md->h5_dataspaces, dataset_id, "vsMesh", "meshA");
    set_attr_string(&md->h5_dataspaces, dataset_id, "vsCentering", "zonal");
    set_attr_string(&md->h5_dataspaces, dataset_id, "vsIndexOrder", "compMinorC");

    float *data = _fmd_array_convert_numerical_scalar_3d_to_flat_float(arr);
    assert(data != NULL);

    status = H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
    assert(status >= 0);

    free(data);

    status = H5Dclose(dataset_id);
    assert(status >= 0);

    status = H5Sclose(dataspace_id);
    assert(status >= 0);

    status = H5Fclose(file_id);
    assert(status >= 0);
}

void _fmd_h5_save_tuple_field_float(fmd_t *md, fmd_string_t fieldname, turi_t *t, fmd_string_t path, fmd_array3D_t *arr)
{
    hid_t file_id, dataset_id, dataspace_id;
    herr_t status;

    file_id = H5Fcreate(path, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    assert(file_id >= 0);

    fmd_ftriple_t ubound = {md->l[0], md->l[1], md->l[2]};
    create_mesh(&md->h5_dataspaces, file_id, t->tdims_global, ubound);

    hsize_t s[4];
    s[0] = t->tdims_global[0];
    s[1] = t->tdims_global[1];
    s[2] = t->tdims_global[2];
    s[3] = 3;

    dataspace_id = H5Screate_simple(4, s, NULL);
    assert(dataspace_id >= 0);

    dataset_id = H5Dcreate1(file_id, fieldname, H5T_IEEE_F32LE, dataspace_id, H5P_DEFAULT);
    assert(dataset_id >= 0);

    set_attr_string(&md->h5_dataspaces, dataset_id, "vsType", "variable");
    set_attr_string(&md->h5_dataspaces, dataset_id, "vsMesh", "meshA");
    set_attr_string(&md->h5_dataspaces, dataset_id, "vsCentering", "zonal");
    set_attr_string(&md->h5_dataspaces, dataset_id, "vsIndexOrder", "compMinorC");

    float *data = _fmd_array_convert_numerical_tuple_3d_to_flat_float(arr);
    assert(data != NULL);

    status = H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
    assert(status >= 0);

    free(data);

    status = H5Dclose(dataset_id);
    assert(status >= 0);

    status = H5Sclose(dataspace_id);
    assert(status >= 0);

    status = H5Fclose(file_id);
    assert(status >= 0);
}

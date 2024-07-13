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
#include "fmd-private.h"
#include "types.h"
#include "array.h"
#include "turi.h"
#include "general.h"
#include <string.h>
#include <stdlib.h>
#include "misc.h"

void _fmd_h5_ds_init(fmd_t *md, h5_dataspaces_t *ds)
{
    hsize_t sz;

    ds->ds_scalar = H5Screate(H5S_SCALAR);
    if (ds->ds_scalar < 0) _fmd_error_unsuccessful_hdf5(md, false, __FILE__, (fmd_string_t)__func__, __LINE__);

    sz = 3;
    ds->ds_simple_1_3 = H5Screate_simple(1, &sz, NULL);
    if (ds->ds_simple_1_3 < 0)
        _fmd_error_unsuccessful_hdf5(md, false, __FILE__, (fmd_string_t)__func__, __LINE__);
}

void _fmd_h5_ds_free(fmd_t *md, h5_dataspaces_t *ds)
{
    herr_t status;

    status = H5Sclose(ds->ds_scalar);
    if (status < 0) _fmd_error_unsuccessful_hdf5(md, false, __FILE__, (fmd_string_t)__func__, __LINE__);

    status = H5Sclose(ds->ds_simple_1_3);
    if (status < 0) _fmd_error_unsuccessful_hdf5(md, false, __FILE__, (fmd_string_t)__func__, __LINE__);
}

static void set_attr_ftuple(fmd_t *md, h5_dataspaces_t *ds, hid_t obj, const char *name, const float *f)
{
    hid_t attr;
    herr_t status;

    attr = H5Acreate(obj, name, H5T_IEEE_F32LE, ds->ds_simple_1_3, H5P_DEFAULT, H5P_DEFAULT);
    if (attr < 0) _fmd_error_unsuccessful_hdf5(md, false, __FILE__, (fmd_string_t)__func__, __LINE__);

    status = H5Awrite(attr, H5T_NATIVE_FLOAT, f);
    if (status < 0) _fmd_error_unsuccessful_hdf5(md, false, __FILE__, (fmd_string_t)__func__, __LINE__);

    status = H5Aclose(attr);
    if (status < 0) _fmd_error_unsuccessful_hdf5(md, false, __FILE__, (fmd_string_t)__func__, __LINE__);
}

static void set_attr_ituple(fmd_t *md, h5_dataspaces_t *ds, hid_t obj, const char *name, const fmd_ituple_t i)
{
    hid_t attr;
    herr_t status;

    attr = H5Acreate(obj, name, H5T_STD_I32LE, ds->ds_simple_1_3, H5P_DEFAULT, H5P_DEFAULT);
    if (attr < 0) _fmd_error_unsuccessful_hdf5(md, false, __FILE__, (fmd_string_t)__func__, __LINE__);

    status = H5Awrite(attr, H5T_NATIVE_INT, i);
    if (status < 0) _fmd_error_unsuccessful_hdf5(md, false, __FILE__, (fmd_string_t)__func__, __LINE__);

    status = H5Aclose(attr);
    if (status < 0) _fmd_error_unsuccessful_hdf5(md, false, __FILE__, (fmd_string_t)__func__, __LINE__);
}

static void set_attr_utuple(fmd_t *md, h5_dataspaces_t *ds, hid_t obj, const char *name, const fmd_utuple_t u)
{
    hid_t attr;
    herr_t status;

    attr = H5Acreate(obj, name, H5T_STD_I32LE, ds->ds_simple_1_3, H5P_DEFAULT, H5P_DEFAULT);
    if (attr < 0) _fmd_error_unsuccessful_hdf5(md, false, __FILE__, (fmd_string_t)__func__, __LINE__);

    status = H5Awrite(attr, H5T_NATIVE_UINT, u);
    if (status < 0) _fmd_error_unsuccessful_hdf5(md, false, __FILE__, (fmd_string_t)__func__, __LINE__);

    status = H5Aclose(attr);
    if (status < 0) _fmd_error_unsuccessful_hdf5(md, false, __FILE__, (fmd_string_t)__func__, __LINE__);
}

static void set_attr_string(fmd_t *md, h5_dataspaces_t *ds, hid_t obj, const char *name, const char *value)
{
    hid_t str_type_mem, str_type_file, attr;
    herr_t status;

    hsize_t len = strlen(value);

    str_type_mem = H5Tcopy(H5T_C_S1);
    if (str_type_mem < 0) _fmd_error_unsuccessful_hdf5(md, false, __FILE__, (fmd_string_t)__func__, __LINE__);

    status = H5Tset_size(str_type_mem, len+1);
    if (status < 0) _fmd_error_unsuccessful_hdf5(md, false, __FILE__, (fmd_string_t)__func__, __LINE__);

    str_type_file = H5Tcopy(H5T_FORTRAN_S1);
    if (str_type_file < 0) _fmd_error_unsuccessful_hdf5(md, false, __FILE__, (fmd_string_t)__func__, __LINE__);

    status = H5Tset_size(str_type_file, len);
    if (status < 0) _fmd_error_unsuccessful_hdf5(md, false, __FILE__, (fmd_string_t)__func__, __LINE__);

    attr = H5Acreate(obj, name, str_type_file, ds->ds_scalar, H5P_DEFAULT, H5P_DEFAULT);
    if (attr < 0) _fmd_error_unsuccessful_hdf5(md, false, __FILE__, (fmd_string_t)__func__, __LINE__);

    status = H5Awrite(attr, str_type_mem, value);
    if (status < 0) _fmd_error_unsuccessful_hdf5(md, false, __FILE__, (fmd_string_t)__func__, __LINE__);

    status = H5Aclose(attr);
    if (status < 0) _fmd_error_unsuccessful_hdf5(md, false, __FILE__, (fmd_string_t)__func__, __LINE__);

    status = H5Tclose(str_type_mem);
    if (status < 0) _fmd_error_unsuccessful_hdf5(md, false, __FILE__, (fmd_string_t)__func__, __LINE__);

    status = H5Tclose(str_type_file);
    if (status < 0) _fmd_error_unsuccessful_hdf5(md, false, __FILE__, (fmd_string_t)__func__, __LINE__);
}

static void create_mesh(fmd_t *md, h5_dataspaces_t *ds, hid_t parent, fmd_utuple_t dims, fmd_ftriple_t UpperBounds)
{
    hid_t group_id;
    herr_t status;

    group_id = H5Gcreate(parent, "meshA", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    set_attr_string(md, ds, group_id, "vsType", "mesh");
    set_attr_string(md, ds, group_id, "vsKind", "uniform");
    set_attr_ituple(md, ds, group_id, "vsStartCell", _fmd_ThreeZeros_int);
    set_attr_utuple(md, ds, group_id, "vsNumCells", dims);
    set_attr_ftuple(md, ds, group_id, "vsLowerBounds", _fmd_ThreeZeros_float);
    set_attr_ftuple(md, ds, group_id, "vsUpperBounds", UpperBounds);

    status = H5Gclose(group_id);
    if (status < 0) _fmd_error_unsuccessful_hdf5(md, false, __FILE__, (fmd_string_t)__func__, __LINE__);
}

void _fmd_h5_save_scalar_field_float(fmd_t *md, fmd_string_t fieldname, turi_t *t, fmd_string_t path, fmd_array3s_t *arr)
{
    hid_t file_id, dataset_id, dataspace_id;
    herr_t status;

    file_id = H5Fcreate(path, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (file_id < 0) _fmd_error_unsuccessful_hdf5(md, false, __FILE__, (fmd_string_t)__func__, __LINE__);

    fmd_ftriple_t ubound = {md->l[0], md->l[1], md->l[2]};
    create_mesh(md, &md->h5_dataspaces, file_id, t->tdims_global, ubound);

    hsize_t s[3];
    s[0] = t->tdims_global[0];
    s[1] = t->tdims_global[1];
    s[2] = t->tdims_global[2];

    dataspace_id = H5Screate_simple(3, s, NULL);
    if (dataspace_id < 0) _fmd_error_unsuccessful_hdf5(md, false, __FILE__, (fmd_string_t)__func__, __LINE__);

    dataset_id = H5Dcreate1(file_id, fieldname, H5T_IEEE_F32LE, dataspace_id, H5P_DEFAULT);
    if (dataset_id < 0) _fmd_error_unsuccessful_hdf5(md, false, __FILE__, (fmd_string_t)__func__, __LINE__);

    set_attr_string(md, &md->h5_dataspaces, dataset_id, "vsType", "variable");
    set_attr_string(md, &md->h5_dataspaces, dataset_id, "vsMesh", "meshA");
    set_attr_string(md, &md->h5_dataspaces, dataset_id, "vsCentering", "zonal");
    set_attr_string(md, &md->h5_dataspaces, dataset_id, "vsIndexOrder", "compMinorC");

    float *data = _fmd_array_convert_numerical_scalar_3d_to_flat_float(md, arr);
    if (data == NULL)
        _fmd_error_function_failed(md, false, __FILE__, (fmd_string_t)__func__, __LINE__,
                                   "_fmd_array_convert_numerical_scalar_3d_to_flat_float");

    status = H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
    if (status < 0) _fmd_error_unsuccessful_hdf5(md, false, __FILE__, (fmd_string_t)__func__, __LINE__);

    free(data);

    status = H5Dclose(dataset_id);
    if (status < 0) _fmd_error_unsuccessful_hdf5(md, false, __FILE__, (fmd_string_t)__func__, __LINE__);

    status = H5Sclose(dataspace_id);
    if (status < 0) _fmd_error_unsuccessful_hdf5(md, false, __FILE__, (fmd_string_t)__func__, __LINE__);

    status = H5Fclose(file_id);
    if (status < 0) _fmd_error_unsuccessful_hdf5(md, false, __FILE__, (fmd_string_t)__func__, __LINE__);
}

void _fmd_h5_save_tuple_field_float(fmd_t *md, fmd_string_t fieldname, turi_t *t, fmd_string_t path, fmd_array3s_t *arr)
{
    hid_t file_id, dataset_id, dataspace_id;
    herr_t status;

    file_id = H5Fcreate(path, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (file_id < 0) _fmd_error_unsuccessful_hdf5(md, false, __FILE__, (fmd_string_t)__func__, __LINE__);

    fmd_ftriple_t ubound = {md->l[0], md->l[1], md->l[2]};
    create_mesh(md, &md->h5_dataspaces, file_id, t->tdims_global, ubound);

    hsize_t s[4];
    s[0] = t->tdims_global[0];
    s[1] = t->tdims_global[1];
    s[2] = t->tdims_global[2];
    s[3] = 3;

    dataspace_id = H5Screate_simple(4, s, NULL);
    if (dataspace_id < 0) _fmd_error_unsuccessful_hdf5(md, false, __FILE__, (fmd_string_t)__func__, __LINE__);

    dataset_id = H5Dcreate1(file_id, fieldname, H5T_IEEE_F32LE, dataspace_id, H5P_DEFAULT);
    if (dataset_id < 0) _fmd_error_unsuccessful_hdf5(md, false, __FILE__, (fmd_string_t)__func__, __LINE__);

    set_attr_string(md, &md->h5_dataspaces, dataset_id, "vsType", "variable");
    set_attr_string(md, &md->h5_dataspaces, dataset_id, "vsMesh", "meshA");
    set_attr_string(md, &md->h5_dataspaces, dataset_id, "vsCentering", "zonal");
    set_attr_string(md, &md->h5_dataspaces, dataset_id, "vsIndexOrder", "compMinorC");

    float *data = _fmd_array_convert_numerical_tuple_3d_to_flat_float(md, arr);
    if (data == NULL)
        _fmd_error_function_failed(md, false, __FILE__, (fmd_string_t)__func__, __LINE__,
                                   "_fmd_array_convert_numerical_tuple_3d_to_flat_float");

    status = H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
    if (status < 0) _fmd_error_unsuccessful_hdf5(md, false, __FILE__, (fmd_string_t)__func__, __LINE__);

    free(data);

    status = H5Dclose(dataset_id);
    if (status < 0) _fmd_error_unsuccessful_hdf5(md, false, __FILE__, (fmd_string_t)__func__, __LINE__);

    status = H5Sclose(dataspace_id);
    if (status < 0) _fmd_error_unsuccessful_hdf5(md, false, __FILE__, (fmd_string_t)__func__, __LINE__);

    status = H5Fclose(file_id);
    if (status < 0) _fmd_error_unsuccessful_hdf5(md, false, __FILE__, (fmd_string_t)__func__, __LINE__);
}

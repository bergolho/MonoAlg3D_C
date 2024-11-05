//
// Created by bergolho and jennyhelyanwe on 03/06/24.
//

#include "domain_helpers.h"

#include "../3dparty/sds/sds.h"
#include "../3dparty/stb_ds.h"
#include "../config/domain_config.h"
#include "../config_helpers/config_helpers.h"
#include "../libraries_common/common_data_structures.h"
#include "../logger/logger.h"
#include "../utils/stop_watch.h"
#include "../utils/utils.h"
#include <assert.h>
#include <time.h>

#include "mesh_info_data.h"
#include "custom_mesh_info_data.h"

#include <float.h>
#include <unistd.h>

SET_CUSTOM_DATA_FOR_MESH(set_custom_data_for_dti_mesh_with_dense_and_sparse_regions_twave) {

    INITIALIZE_DTI_MESH_INFO(cell);

    // Ruben mesh labels (transmuralityLabels): 1=ENDO, 2=MID, 3=EPI
    // Dense/Sparse regions (DenseSparseRegion): 1=NORMAL, 2=DENSE, 3=SPARSE
    // Cellular model labels: 0=ENDO, 1=MID, 2=EPI, 3=DENSE, 4=SPARSE
    DTI_MESH_TRANSMURALITY(cell) = custom_data[0];
    DTI_MESH_TRANSMURALITY_LABELS(cell) = custom_data[1]-1;
    DTI_MESH_BASE_APEX_HETEROGENEITY(cell) = custom_data[2];
    DTI_MESH_FAST_ENDO(cell) = (int)custom_data[3];
    
    cell->sigma.fibers.f[0] = custom_data[5];
    cell->sigma.fibers.f[1] = custom_data[6];
    cell->sigma.fibers.f[2] = custom_data[7];

    cell->sigma.fibers.s[0] = custom_data[8];
    cell->sigma.fibers.s[1] = custom_data[9];
    cell->sigma.fibers.s[2] = custom_data[10];

    cell->sigma.fibers.n[0] = custom_data[11];
    cell->sigma.fibers.n[1] = custom_data[12];
    cell->sigma.fibers.n[2] = custom_data[13];

    DTI_MESH_SCALE_FACTOR_IKS(cell) = custom_data[14];
}

// Domain function to load the ALG file
// Works with all DTI meshes (T-wave personalisation) - Fiber orientation, transmurality, dense/sparse regions 
// and scale factor for IKs (resolution = 500um)
SET_SPATIAL_DOMAIN(initialize_grid_with_dti_mesh_with_dense_and_sparse_regions_twave) {

    char *mesh_file = NULL;
    GET_PARAMETER_STRING_VALUE_OR_REPORT_ERROR(mesh_file, config, "mesh_file");

    real_cpu start_h = 500.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real_cpu, start_h, config, "original_discretization");

    real_cpu desired_h = 500.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real_cpu, desired_h, config, "desired_discretization");

    assert(the_grid);

    log_info("DTI Mesh with discretization: %lf\n", desired_h);

    uint32_t num_volumes = 0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(uint32_t, num_volumes, config, "num_volumes");

    // ALG mesh column config:
    // center_x, center_y, center_z, dx, dy, dz, transmurality, transmurality_label, base_to_apex, dense_sparse_region_tag, fibrotic, fx, fy, fz, sx, sy, sz, nx, ny, nz, sf_IKs  
    uint32_t num_extra_fields = 15;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(uint32_t, num_extra_fields, config, "num_extra_fields");

    int num_loaded = set_custom_mesh_from_file(the_grid, mesh_file, num_volumes, start_h, num_extra_fields, set_custom_data_for_dti_mesh_with_dense_and_sparse_regions_twave);

    log_info("Read %d volumes from file (expected %d): %s\n", num_loaded, num_volumes, mesh_file);

    int num_refs_remaining = (int)(start_h / desired_h) - 1;

    if(num_refs_remaining > 0) {
        refine_grid(the_grid, num_refs_remaining);
    }

    the_grid->start_discretization = SAME_POINT3D(desired_h);

    free(mesh_file);

    return num_loaded;
}

SET_SPATIAL_DOMAIN(initialize_grid_with_dti_mesh_with_dense_and_sparse_regions_twave_mesh_number) {

	char *mesh_directory = NULL;
	GET_PARAMETER_STRING_VALUE_OR_REPORT_ERROR(mesh_directory, config, "mesh_directory");
	int mesh_number = 0;
	GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(int, mesh_number, config, "mesh_number");

    real_cpu start_h = 500.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real_cpu, start_h, config, "original_discretization");

    real_cpu desired_h = 500.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real_cpu, desired_h, config, "desired_discretization");

    assert(the_grid);

    log_info("DTI Mesh with discretization: %lf\n", desired_h);

	char mesh_file[500];
	sprintf(mesh_file,"%s/%d.alg", mesh_directory, mesh_number);
    uint32_t num_volumes = 0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(uint32_t, num_volumes, config, "num_volumes");

    // ALG mesh column config:
    // center_x, center_y, center_z, dx, dy, dz, transmurality, transmurality_label, base_to_apex, dense_sparse_region_tag, fibrotic, fx, fy, fz, sx, sy, sz, nx, ny, nz, sf_IKs  
    uint32_t num_extra_fields = 15;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(uint32_t, num_extra_fields, config, "num_extra_fields");

    int num_loaded = set_custom_mesh_from_file(the_grid, mesh_file, num_volumes, start_h, num_extra_fields, set_custom_data_for_dti_mesh_with_dense_and_sparse_regions_twave);

    log_info("Read %d volumes from file (expected %d): %s\n", num_loaded, num_volumes, mesh_file);

    int num_refs_remaining = (int)(start_h / desired_h) - 1;

    if(num_refs_remaining > 0) {
        refine_grid(the_grid, num_refs_remaining);
    }

    the_grid->start_discretization = SAME_POINT3D(desired_h);

    free(mesh_file);

    return num_loaded;
}

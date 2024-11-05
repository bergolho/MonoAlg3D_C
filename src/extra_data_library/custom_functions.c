//
// Created by bergolho and jennyhelyanwe on 03/06/24.
//

#include <unistd.h>

#include "../config/extra_data_config.h"
#include "../config_helpers/config_helpers.h"
#include "../libraries_common/common_data_structures.h"
#include "../domains_library/mesh_info_data.h"
#include "../domains_library/custom_mesh_info_data.h"
#include "../logger/logger.h"
#include "helper_functions.h"

// Extra data function to configure the transmurality and IKs scaling factor
// Works with all DTI meshes - Extra data function to load the ToRORd cellular model with modifications
// Function used for the T-wave personalisation paper.
SET_EXTRA_DATA(set_extra_data_for_dti_mesh_twave_with_torord_gksgkrtjca) {

	uint32_t num_active_cells = the_grid->num_active_cells;
	struct cell_node ** ac = the_grid->active_cells;

	struct extra_data_for_torord_gksgkrtjca_twave *extra_data = NULL;
	extra_data = set_common_torord_gksgkrtjca_twave_data(config, num_active_cells);		// ToRORd used for T-wave personalisation

	// Ruben mesh labels (transmuralityLabels): 1=ENDO, 2=MID, 3=EPI
    // Dense/Sparse regions (DenseSparseRegion): 1=NORMAL, 2=DENSE, 3=SPARSE
    // Cellular model labels: 0=ENDO, 1=MID, 2=EPI, 3=DENSE, 4=SPARSE
	OMP(parallel for)
	for (int i = 0; i < num_active_cells; i++) {
		real tag = DTI_MESH_TRANSMURALITY_LABELS(ac[i]);
        real sf_IKs = DTI_MESH_SCALE_FACTOR_IKS(ac[i]);
		extra_data->transmurality[i] = tag;
        extra_data->sf_IKs[i] = sf_IKs;
	}

	SET_EXTRA_DATA_SIZE(sizeof(struct extra_data_for_torord_gksgkrtjca_twave));

	return (void*)extra_data;
}

//
// Created by sachetto on 13/10/17.
//

#include <stdbool.h>
#include <stdlib.h>

#include "../3dparty/sds/sds.h"
#include "../alg/grid/grid.h"
#include "../config/save_mesh_config.h"
#include "../utils/utils.h"

#include "../vtk_utils/vtk_unstructured_grid.h"
#include "../vtk_utils/vtk_polydata_grid.h"
#include "../libraries_common/common_data_structures.h"

#ifdef COMPILE_CUDA
#include "../gpu_utils/gpu_utils.h"
#endif

static char *file_prefix;
static char *file_prefix_purkinje;
static bool binary = false;
static bool clip_with_plain = false;
static bool clip_with_bounds = false;
static bool save_pvd = true;
static bool compress = false;
static int compression_level = 3;
char *output_dir;

static bool initialized = false;

static struct vtk_polydata_grid *vtk_polydata = NULL;

void add_file_to_pvd(real_cpu current_t, const char *output_dir, const char *base_name, bool first_save_call);

static sds create_base_name(char *f_prefix, int iteration_count, char *extension) {
    return sdscatprintf(sdsempty(), "%s_it_%d.%s", f_prefix, iteration_count, extension);
}

SAVE_MESH(save_as_adjacency_list) {

    int iteration_count = time_info->iteration;

    if(!initialized) {
        GET_PARAMETER_VALUE_CHAR_OR_REPORT_ERROR(output_dir, config->config_data, "output_dir");
        GET_PARAMETER_VALUE_CHAR_OR_REPORT_ERROR(file_prefix, config->config_data, "file_prefix");
        initialized = true;
    }

    sds tmp = sdsnew(output_dir);
    tmp = sdscat(tmp, "/");

    sds base_name = NULL;
    if(binary) {
        base_name = create_base_name(file_prefix, iteration_count, "bin");
    }
    else {
        base_name = create_base_name(file_prefix, iteration_count, "txt");
    }

    tmp = sdscat(tmp, base_name);

    FILE *output_file = fopen(tmp, "w");

    sdsfree(base_name);
    sdsfree(tmp);

    struct cell_node *neighbour;

    FOR_EACH_CELL(the_grid) {

        if(cell->active) {

            fprintf(output_file, "%d ", cell->grid_position);

            neighbour = get_cell_neighbour(cell, cell->z_front);
            if(neighbour) {
                fprintf(output_file, "%d ", neighbour->grid_position);
            }
            neighbour = get_cell_neighbour(cell, cell->z_back);
            if(neighbour) {
                fprintf(output_file, "%d ", neighbour->grid_position);
            }

            neighbour = get_cell_neighbour(cell, cell->y_down);
            if(neighbour) {
                fprintf(output_file, "%d ", neighbour->grid_position);
            }

            neighbour = get_cell_neighbour(cell, cell->y_top);
            if(neighbour) {
                fprintf(output_file, "%d ", neighbour->grid_position);
            }

            neighbour = get_cell_neighbour(cell, cell->x_right);
            if(neighbour) {
                fprintf(output_file, "%d ", neighbour->grid_position);
            }

            neighbour = get_cell_neighbour(cell, cell->x_left);
            if(neighbour) {
                fprintf(output_file, "%d", neighbour->grid_position);
            }

            fprintf(output_file, "\n");
        }
    }

    fclose(output_file);


}

struct save_one_cell_state_variables_persistent_data {
    FILE *file;
    char *file_name;
    real_cpu cell_center_x;
    real_cpu cell_center_y;
    real_cpu cell_center_z;
    uint32_t cell_sv_position;

};

INIT_SAVE_MESH(init_save_one_cell_state_variables) {
    config->persistent_data = malloc(sizeof(struct save_one_cell_state_variables_persistent_data));
    GET_PARAMETER_VALUE_CHAR_OR_REPORT_ERROR( ((struct save_one_cell_state_variables_persistent_data *) config->persistent_data)->file_name, config->config_data, "file_name");
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, ((struct save_one_cell_state_variables_persistent_data *) config->persistent_data)->cell_center_x, config->config_data, "cell_center_x");
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, ((struct save_one_cell_state_variables_persistent_data *) config->persistent_data)->cell_center_y, config->config_data, "cell_center_y");
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, ((struct save_one_cell_state_variables_persistent_data *) config->persistent_data)->cell_center_z, config->config_data, "cell_center_z");

    ((struct save_one_cell_state_variables_persistent_data *) config->persistent_data)->file = fopen(((struct save_one_cell_state_variables_persistent_data *) config->persistent_data)->file_name, "w");
    ((struct save_one_cell_state_variables_persistent_data *) config->persistent_data)->cell_sv_position = -1;

}

SAVE_MESH(save_one_cell_state_variables) {

    struct save_one_cell_state_variables_persistent_data *params = ((struct save_one_cell_state_variables_persistent_data *) config->persistent_data);

    if(params->cell_sv_position == -1) {
        if(!the_grid->adaptive) {
            FOR_EACH_CELL(the_grid) {
                if(cell->center.x == params->cell_center_x && cell->center.y == params->cell_center_y && cell->center.z == params->cell_center_z) {
                    params->cell_sv_position = cell->sv_position;
                    printf("%d\n", params->cell_sv_position);
                    break;
                }
            }
        }
    }

    if(ode_solver->gpu) {
            real *cell_sv;

            // cell_sv = (real*)malloc(sizeof(real)*ode_solver->model_data.number_of_ode_equations*the_grid->num_active_cells);
            cell_sv = (real *)malloc(sizeof(real) * ode_solver->model_data.number_of_ode_equations);

            // check_cuda_errors(cudaMemcpy2D(cell_sv, the_grid->num_active_cells*sizeof(real), ode_solver->sv, ode_solver->pitch, the_grid->num_active_cells*sizeof(real), ode_solver->model_data.number_of_ode_equations, cudaMemcpyDeviceToHost));
            check_cuda_errors(cudaMemcpy2D(cell_sv, sizeof(real), ode_solver->sv + params->cell_sv_position,
                                           ode_solver->pitch, sizeof(real),
                                           ode_solver->model_data.number_of_ode_equations, cudaMemcpyDeviceToHost));

            fprintf(params->file, "%lf %lf ", time_info->current_t, cell_sv[0]);
            for(int i = 2; i < 11; i++) {
                fprintf(params->file, " %lf ", cell_sv[i]);
            }

            fprintf(params->file, " %lf ", cell_sv[13]);
            fprintf(params->file, " %lf ", cell_sv[11]);
            fprintf(params->file, " %lf ", cell_sv[12]);

            for(int i = 14; i < 30; i++) {
                fprintf(params->file, " %lf ", cell_sv[i]);
            }

            fprintf(params->file, " %lf ", cell_sv[32]);
            fprintf(params->file, " %lf ", cell_sv[33]);

            fprintf(params->file, " %lf ", cell_sv[30]);
            fprintf(params->file, " %lf ", cell_sv[31]);

            fprintf(params->file, " %lf ", cell_sv[39]);
            fprintf(params->file, " %lf ", cell_sv[40]);
            fprintf(params->file, " %lf ", cell_sv[41]);
            fprintf(params->file, " %lf ", cell_sv[1]);

            fprintf(params->file, " %lf ", cell_sv[34]);
            fprintf(params->file, " %lf ", cell_sv[35]);
            fprintf(params->file, " %lf ", cell_sv[36]);
            fprintf(params->file, " %lf ", cell_sv[37]);
            fprintf(params->file, " %lf ", cell_sv[38]);
            fprintf(params->file, " %lf\n", cell_sv[42]);
            free(cell_sv);

    }
    else {

        real *cell_sv =  &ode_solver->sv[params->cell_sv_position * ode_solver->model_data.number_of_ode_equations];

        fprintf(params->file, "%lf %lf ", time_info->current_t, cell_sv[0]);
        for(int i = 2; i < 11; i++) {
            fprintf(params->file, " %lf ", cell_sv[i]);
        }

        fprintf(params->file, " %lf ", cell_sv[13]);
        fprintf(params->file, " %lf ", cell_sv[11]);
        fprintf(params->file, " %lf ", cell_sv[12]);

        for(int i = 14; i < 30; i++) {
            fprintf(params->file, " %lf ", cell_sv[i]);
        }

        fprintf(params->file, " %lf ", cell_sv[32]);
        fprintf(params->file, " %lf ", cell_sv[33]);

        fprintf(params->file, " %lf ", cell_sv[30]);
        fprintf(params->file, " %lf ", cell_sv[31]);

        fprintf(params->file, " %lf ", cell_sv[39]);
        fprintf(params->file, " %lf ", cell_sv[40]);
        fprintf(params->file, " %lf ", cell_sv[41]);
        fprintf(params->file, " %lf ", cell_sv[1]);

        fprintf(params->file, " %lf ", cell_sv[34]);
        fprintf(params->file, " %lf ", cell_sv[35]);
        fprintf(params->file, " %lf ", cell_sv[36]);
        fprintf(params->file, " %lf ", cell_sv[37]);
        fprintf(params->file, " %lf ", cell_sv[38]);
        fprintf(params->file, " %lf\n", cell_sv[42]);


    }

}

END_SAVE_MESH(end_save_one_cell_state_variables) {
    free(((struct save_one_cell_state_variables_persistent_data *) config->persistent_data)->file_name);
    fclose(((struct save_one_cell_state_variables_persistent_data *) config->persistent_data)->file);
    free(config->persistent_data);
}

SAVE_MESH(save_as_text_or_binary) {

    int iteration_count = time_info->iteration;

//    if(!initialized) {
        GET_PARAMETER_VALUE_CHAR_OR_REPORT_ERROR(output_dir, config->config_data, "output_dir");
        GET_PARAMETER_VALUE_CHAR_OR_REPORT_ERROR(file_prefix, config->config_data, "file_prefix");
        GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(binary, config->config_data, "binary");
        GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(clip_with_plain, config->config_data, "clip_with_plain");
        GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(clip_with_bounds, config->config_data, "clip_with_bounds");
        initialized = true;
//    }

    real_cpu min_x = 0.0;
    real_cpu min_y = 0.0;
    real_cpu min_z = 0.0;
    real_cpu max_x = 0.0;
    real_cpu max_y = 0.0;
    real_cpu max_z = 0.0;

    float p0[3] = {0, 0, 0};
    float n[3] = {0, 0, 0};

    if(clip_with_plain) {
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, n[0], config->config_data, "normal_x");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, n[1], config->config_data, "normal_y");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, n[2], config->config_data, "normal_z");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, p0[0], config->config_data, "origin_x");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, p0[1], config->config_data, "origin_y");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, p0[2], config->config_data, "origin_z");
    }

    if(clip_with_bounds) {
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, min_x, config->config_data, "min_x");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, min_y, config->config_data, "min_y");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, min_z, config->config_data, "min_z");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, max_x, config->config_data, "max_x");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, max_y, config->config_data, "max_y");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, max_z, config->config_data, "max_z");
    }

    real_cpu l = sqrtf(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
    real_cpu A = n[0] / l;
    real_cpu B = n[1] / l;
    real_cpu C = n[2] / l;
    real_cpu D = -(n[0] * p0[0] + n[1] * p0[1] + n[2] * p0[2]);

    real_cpu side;

    sds tmp = sdsnew(output_dir);
    tmp = sdscat(tmp, "/");

    sds base_name = NULL;
    if(binary) {
        base_name = create_base_name(file_prefix, iteration_count, "bin");
    }
    else {
        base_name = create_base_name(file_prefix, iteration_count, "txt");
    }

    tmp = sdscat(tmp, base_name);

    FILE *output_file = fopen(tmp, "w");

    sdsfree(base_name);
    sdsfree(tmp);

    struct cell_node *grid_cell = the_grid->first_cell;

    float center_x, center_y, center_z, dx, dy, dz;
    float v;

    while(grid_cell != 0) {

        if(grid_cell->active) {

            center_x = grid_cell->center.x;
            center_y = grid_cell->center.y;
            center_z = grid_cell->center.z;

            if(clip_with_plain) {
                side = A * center_x + B * center_y + C * center_z + D;
                if(side < 0) {
                    grid_cell = grid_cell->next;
                    continue;
                }
            }

            if(clip_with_bounds) {
                bool ignore_cell = center_x < min_x || center_x > max_x || center_y < min_y || center_y > max_y ||
                                   center_z < min_z || center_z > max_z;

                if(ignore_cell) {
                    grid_cell = grid_cell->next;
                    continue;
                }
            }

            v = grid_cell->v;
            dx = grid_cell->discretization.x/2.0;
            dy = grid_cell->discretization.y/2.0;
            dz = grid_cell->discretization.z/2.0;

            if(binary) {
                fwrite(&center_x, sizeof(center_x), 1, output_file);
                fwrite(&center_y, sizeof(center_y), 1, output_file);
                fwrite(&center_z, sizeof(center_z), 1, output_file);
                fwrite(&dx, sizeof(dx), 1, output_file);
                fwrite(&dy, sizeof(dy), 1, output_file);
                fwrite(&dz, sizeof(dz), 1, output_file);
                fwrite(&v, sizeof(v), 1, output_file);
            } else {
                fprintf(output_file, "%g,%g,%g,%g,%g,%g,%g\n", center_x, center_y, center_z, dx, dy, dz, v);
            }
        }
        grid_cell = grid_cell->next;
    }

    fclose(output_file);
}

struct save_as_vtk_or_vtu_persistent_data {
    struct vtk_unstructured_grid * grid;
    bool first_save_call;
};

INIT_SAVE_MESH(init_save_as_vtk_or_vtu) {
    config->persistent_data = malloc(sizeof(struct save_as_vtk_or_vtu_persistent_data));
    ((struct save_as_vtk_or_vtu_persistent_data *) config->persistent_data)->grid = NULL;
    ((struct save_as_vtk_or_vtu_persistent_data *) config->persistent_data)->first_save_call = true;
}

END_SAVE_MESH(end_save_as_vtk_or_vtu) {
    free_vtk_unstructured_grid(((struct save_as_vtk_or_vtu_persistent_data *) config->persistent_data)->grid);
    free(config->persistent_data);
}

SAVE_MESH(save_as_vtk) {

    int iteration_count = time_info->iteration;

    if(((struct save_as_vtk_or_vtu_persistent_data *) config->persistent_data)->first_save_call) {
        GET_PARAMETER_VALUE_CHAR_OR_REPORT_ERROR(output_dir, config->config_data, "output_dir");
        GET_PARAMETER_VALUE_CHAR_OR_REPORT_ERROR(file_prefix, config->config_data, "file_prefix");
        GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(clip_with_plain, config->config_data, "clip_with_plain");
        GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(clip_with_bounds, config->config_data, "clip_with_bounds");
        GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(binary, config->config_data, "binary");

        ((struct save_as_vtk_or_vtu_persistent_data *) config->persistent_data)->first_save_call = false;

    }
    float plain_coords[6] = {0, 0, 0, 0, 0, 0};
    float bounds[6] = {0, 0, 0, 0, 0, 0};

    if(clip_with_plain) {
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, plain_coords[0], config->config_data, "origin_x");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, plain_coords[1], config->config_data, "origin_y");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, plain_coords[2], config->config_data, "origin_z");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, plain_coords[3], config->config_data, "normal_x");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, plain_coords[4], config->config_data, "normal_y");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, plain_coords[5], config->config_data, "normal_z");
    }

    if(clip_with_bounds) {
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, bounds[0], config->config_data, "min_x");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, bounds[1], config->config_data, "min_y");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, bounds[2], config->config_data, "min_z");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, bounds[3], config->config_data, "max_x");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, bounds[4], config->config_data, "max_y");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, bounds[5], config->config_data, "max_z");
    }

    sds output_dir_with_file = sdsnew(output_dir);
    output_dir_with_file = sdscat(output_dir_with_file, "/");
    sds base_name = create_base_name(file_prefix, iteration_count, "vtk");

    real_cpu current_t = time_info->current_t;

    //TODO: change this. We dont need the current_t here
    output_dir_with_file = sdscatprintf(output_dir_with_file, base_name, current_t);

    bool read_only_data = ((struct save_as_vtk_or_vtu_persistent_data *) config->persistent_data)->grid != NULL;
    new_vtk_unstructured_grid_from_alg_grid(&(((struct save_as_vtk_or_vtu_persistent_data *) config->persistent_data)->grid), the_grid, clip_with_plain, plain_coords, clip_with_bounds, bounds, read_only_data);

    save_vtk_unstructured_grid_as_legacy_vtk(((struct save_as_vtk_or_vtu_persistent_data *) config->persistent_data)->grid, output_dir_with_file, binary);

    if(the_grid->adaptive) {
        free_vtk_unstructured_grid(((struct save_as_vtk_or_vtu_persistent_data *) config->persistent_data)->grid);
        ((struct save_as_vtk_or_vtu_persistent_data *) config->persistent_data)->grid = NULL;
    }

    sdsfree(output_dir_with_file);
    sdsfree(base_name);
}

static inline void write_pvd_header(FILE *pvd_file) {
    fprintf(pvd_file, "<VTKFile type=\"Collection\" version=\"0.1\" compressor=\"vtkZLibDataCompressor\">\n");
    fprintf(pvd_file, "\t<Collection>\n");
    fprintf(pvd_file, "\t</Collection>\n");
    fprintf(pvd_file, "</VTKFile>");
}

void add_file_to_pvd(real_cpu current_t, const char *output_dir, const char *base_name, bool first_call) {

    sds pvd_name = sdsnew(output_dir);
    pvd_name = sdscat(pvd_name, "/simulation_result.pvd");

    static FILE *pvd_file = NULL;
    pvd_file = fopen(pvd_name, "r+");

    if(!pvd_file) {
        pvd_file = fopen(pvd_name, "w");
        write_pvd_header(pvd_file);
    }
    else {
        if(first_call) {
            fclose(pvd_file);
            pvd_file = fopen(pvd_name, "w");
            write_pvd_header(pvd_file);
        }
    }

    sdsfree(pvd_name);

    fseek(pvd_file, -26, SEEK_END);

    fprintf(pvd_file, "\n\t\t<DataSet timestep=\"%lf\" group=\"\" part=\"0\" file=\"%s\"/>\n", current_t, base_name);
    fprintf(pvd_file, "\t</Collection>\n");
    fprintf(pvd_file, "</VTKFile>");
    fclose(pvd_file);
}

SAVE_MESH(save_as_vtu) {

    int iteration_count = time_info->iteration;

    if(((struct save_as_vtk_or_vtu_persistent_data *) config->persistent_data)->first_save_call) {
        GET_PARAMETER_VALUE_CHAR_OR_REPORT_ERROR(output_dir, config->config_data, "output_dir");
        GET_PARAMETER_VALUE_CHAR_OR_REPORT_ERROR(file_prefix, config->config_data, "file_prefix");
        GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(clip_with_plain, config->config_data, "clip_with_plain");
        GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(clip_with_bounds, config->config_data, "clip_with_bounds");
        GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(binary, config->config_data, "binary");
        GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(save_pvd, config->config_data, "save_pvd");
        GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(compress, config->config_data, "compress");
        GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(int, compression_level, config->config_data, "compression_level");

        if(compress) binary = true;

        if(!save_pvd) {
            ((struct save_as_vtk_or_vtu_persistent_data *) config->persistent_data)->first_save_call = false;
        }
    }

    float plain_coords[6] = {0, 0, 0, 0, 0, 0};
    float bounds[6] = {0, 0, 0, 0, 0, 0};

    if(clip_with_plain) {
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, plain_coords[0], config->config_data, "origin_x");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, plain_coords[1], config->config_data, "origin_y");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, plain_coords[2], config->config_data, "origin_z");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, plain_coords[3], config->config_data, "normal_x");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, plain_coords[4], config->config_data, "normal_y");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, plain_coords[5], config->config_data, "normal_z");
    }

    if(clip_with_bounds) {
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, bounds[0], config->config_data, "min_x");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, bounds[1], config->config_data, "min_y");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, bounds[2], config->config_data, "min_z");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, bounds[3], config->config_data, "max_x");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, bounds[4], config->config_data, "max_y");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, bounds[5], config->config_data, "max_z");
    }


    sds output_dir_with_file = sdsnew(output_dir);
    output_dir_with_file = sdscat(output_dir_with_file, "/");
    sds base_name = create_base_name(file_prefix, iteration_count, "vtu");

    real_cpu current_t = time_info->current_t;

    output_dir_with_file = sdscatprintf(output_dir_with_file, base_name, current_t);

    if(save_pvd) {
        add_file_to_pvd(current_t, output_dir, base_name, ((struct save_as_vtk_or_vtu_persistent_data *) config->persistent_data)->first_save_call);
        ((struct save_as_vtk_or_vtu_persistent_data *) config->persistent_data)->first_save_call = false;
    }

    bool read_only_data = ((struct save_as_vtk_or_vtu_persistent_data *) config->persistent_data)->grid != NULL;
    new_vtk_unstructured_grid_from_alg_grid(&((struct save_as_vtk_or_vtu_persistent_data *) config->persistent_data)->grid, the_grid, clip_with_plain, plain_coords, clip_with_bounds, bounds, read_only_data);

    if(compress) {
        save_vtk_unstructured_grid_as_vtu_compressed(((struct save_as_vtk_or_vtu_persistent_data *) config->persistent_data)->grid, output_dir_with_file, compression_level);
    }
    else {
        save_vtk_unstructured_grid_as_vtu(((struct save_as_vtk_or_vtu_persistent_data *) config->persistent_data)->grid, output_dir_with_file, binary);
    }

    //TODO: I do not know if we should to this here or call the end and init save functions on the adaptivity step.....
    if(the_grid->adaptive) {
        free_vtk_unstructured_grid(((struct save_as_vtk_or_vtu_persistent_data *) config->persistent_data)->grid);
        ((struct save_as_vtk_or_vtu_persistent_data *) config->persistent_data)->grid = NULL;
    }

}

// SAVE_MESH(save_as_vtk_purkinje) {

//     char *output_dir;
//     GET_PARAMETER_VALUE_CHAR_OR_REPORT_ERROR(output_dir, config->config_data, "output_dir");


//     if(!initialized)
//     {
//         GET_PARAMETER_VALUE_CHAR_OR_REPORT_ERROR(file_prefix, config->config_data, "file_prefix");
//         GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(clip_with_plain, config->config_data, "clip_with_plain");
//         GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(clip_with_bounds, config->config_data, "clip_with_bounds");
//         GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(binary, config->config_data, "binary");
//         initialized = true;
//     }
//     float plain_coords[6] = {0, 0, 0, 0, 0, 0};
//     float bounds[6] = {0, 0, 0, 0, 0, 0};

//     if(clip_with_plain) {
//         GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, plain_coords[0], config->config_data, "origin_x");
//         GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, plain_coords[1], config->config_data, "origin_y");
//         GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, plain_coords[2], config->config_data, "origin_z");
//         GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, plain_coords[3], config->config_data, "normal_x");
//         GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, plain_coords[3], config->config_data, "normal_x");
//         GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, plain_coords[4], config->config_data, "normal_y");
//         GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, plain_coords[5], config->config_data, "normal_z");
//     }

//     if(clip_with_bounds) {
//         GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, bounds[0], config->config_data, "min_x");
//         GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, bounds[1], config->config_data, "min_y");
//         GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, bounds[2], config->config_data, "min_z");
//         GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, bounds[3], config->config_data, "max_x");
//         GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, bounds[4], config->config_data, "max_y");
//         GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, bounds[5], config->config_data, "max_z");
//     }

//     sds output_dir_with_file = sdsnew(output_dir);
//     output_dir_with_file = sdscat(output_dir_with_file, "/");
//     sds base_name = create_base_name(file_prefix, iteration_count, "vtk");
//     output_dir_with_file = sdscatprintf(output_dir_with_file, base_name, current_t);

//     new_vtk_polydata_grid_from_purkinje_grid(&vtk_polydata, the_grid, clip_with_plain, plain_coords, clip_with_bounds, bounds, !the_grid->adaptive,'v');
//     save_vtk_polydata_grid_as_legacy_vtk(vtk_polydata, output_dir_with_file, binary);

//     if(the_grid->adaptive)
//         free_vtk_polydata_grid(vtk_polydata);

//     sdsfree(output_dir_with_file);
//     sdsfree(base_name);

// }

void write_transmembrane_potential_vtp (struct vtk_polydata_grid **vtk_polydata, struct grid *the_grid,
                                        char *output_dir, char *file_prefix, int iteration_count, real_cpu current_t,
                                        bool save_pvd, bool compress, int compression_level, bool binary,
                                        bool clip_with_plain, float *plain_coords,
                                        bool clip_with_bounds, float *bounds)
{
    assert(the_grid->purkinje);

    sds output_dir_with_file = sdsnew(output_dir);
    output_dir_with_file = sdscat(output_dir_with_file, "/");
    sds base_name = create_base_name(file_prefix, iteration_count, "vtp");

    output_dir_with_file = sdscatprintf(output_dir_with_file, base_name, current_t);

    if(save_pvd)
    {
        //TODO: we need an INIT en END function for this
        static bool first_call = true;
        add_file_to_pvd(current_t, output_dir, base_name, first_call);
        first_call = false;
    }

    new_vtk_polydata_grid_from_purkinje_grid(vtk_polydata, the_grid->purkinje, clip_with_plain, plain_coords, clip_with_bounds, bounds, *vtk_polydata!=NULL,'v');

    if(compress)
    {
        save_vtk_polydata_grid_as_vtp_compressed(*vtk_polydata, output_dir_with_file, compression_level);
    }
    else
    {
        save_vtk_polydata_grid_as_vtp(*vtk_polydata, output_dir_with_file, binary);
    }

    if(the_grid->adaptive)
        free_vtk_polydata_grid(*vtk_polydata);

    sdsfree(output_dir_with_file);
    sdsfree(base_name);
}

void write_transmembrane_potential_vtu (struct vtk_unstructured_grid **vtk_grid, struct grid *the_grid,
                                        char *output_dir, char *file_prefix, int iteration_count, real_cpu current_t,
                                        bool save_pvd, bool compress, int compression_level, bool binary,
                                        bool clip_with_plain, float *plain_coords,
                                        bool clip_with_bounds, float *bounds)
{
    sds output_dir_with_file = sdsnew(output_dir);
    output_dir_with_file = sdscat(output_dir_with_file, "/");
    sds base_name = create_base_name(file_prefix, iteration_count, "vtu");

    output_dir_with_file = sdscatprintf(output_dir_with_file, base_name, current_t);

    if(save_pvd) {
        static bool first_call = true;
        add_file_to_pvd(current_t, output_dir, base_name, first_call);
        first_call = false;
    }

    new_vtk_unstructured_grid_from_alg_grid(vtk_grid, the_grid, clip_with_plain, plain_coords, clip_with_bounds, bounds, *vtk_grid!=NULL);

    if(compress) {
        save_vtk_unstructured_grid_as_vtu_compressed(*vtk_grid, output_dir_with_file, compression_level);
    }
    else {
        save_vtk_unstructured_grid_as_vtu(*vtk_grid, output_dir_with_file, binary);
    }

    sdsfree(output_dir_with_file);
    sdsfree(base_name);
}


SAVE_MESH(save_as_vtp_purkinje) {

    char *output_dir;
    GET_PARAMETER_VALUE_CHAR_OR_REPORT_ERROR(output_dir, config->config_data, "output_dir");

    if(!initialized) {
        GET_PARAMETER_VALUE_CHAR_OR_REPORT_ERROR(file_prefix, config->config_data, "file_prefix");
        GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(clip_with_plain, config->config_data, "clip_with_plain");
        GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(clip_with_bounds, config->config_data, "clip_with_bounds");
        GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(binary, config->config_data, "binary");
        GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(save_pvd, config->config_data, "save_pvd");
        GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(compress, config->config_data, "compress");
        GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(int, compression_level, config->config_data, "compression_level");

        if(compress) binary = true;

        initialized = true;
    }
    float plain_coords[6] = {0, 0, 0, 0, 0, 0};
    float bounds[6] = {0, 0, 0, 0, 0, 0};

    if(clip_with_plain) {
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, plain_coords[0], config->config_data, "origin_x");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, plain_coords[1], config->config_data, "origin_y");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, plain_coords[2], config->config_data, "origin_z");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, plain_coords[3], config->config_data, "normal_x");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, plain_coords[4], config->config_data, "normal_y");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, plain_coords[5], config->config_data, "normal_z");
    }

    if(clip_with_bounds) {
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, bounds[0], config->config_data, "min_x");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, bounds[1], config->config_data, "min_y");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, bounds[2], config->config_data, "min_z");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, bounds[3], config->config_data, "max_x");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, bounds[4], config->config_data, "max_y");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, bounds[5], config->config_data, "max_z");
    }

    write_transmembrane_potential_vtp(&vtk_polydata, the_grid, output_dir, file_prefix, time_info->iteration, time_info->current_t,\
                                             save_pvd, compress, compression_level, binary, clip_with_plain, plain_coords, clip_with_bounds, bounds);
}

static struct vtk_unstructured_grid *vtk_grid = NULL;
static char *file_prefix;
static char *file_prefix_purkinje;

SAVE_MESH(save_as_vtu_tissue_coupled_vtp_purkinje) {

    char *output_dir;
    GET_PARAMETER_VALUE_CHAR_OR_REPORT_ERROR(output_dir, config->config_data, "output_dir");


    if(!initialized)
    {
        GET_PARAMETER_VALUE_CHAR_OR_REPORT_ERROR(file_prefix, config->config_data, "file_prefix");
        GET_PARAMETER_VALUE_CHAR_OR_REPORT_ERROR(file_prefix_purkinje, config->config_data, "file_prefix_purkinje");
        GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(clip_with_plain, config->config_data, "clip_with_plain");
        GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(clip_with_bounds, config->config_data, "clip_with_bounds");
        GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(binary, config->config_data, "binary");
        GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(save_pvd, config->config_data, "save_pvd");
        GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(compress, config->config_data, "compress");
        GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(int, compression_level, config->config_data, "compression_level");

        if(compress) binary = true;

        initialized = true;
    }
    float plain_coords[6] = {0, 0, 0, 0, 0, 0};
    float bounds[6] = {0, 0, 0, 0, 0, 0};

    if(clip_with_plain)
    {
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, plain_coords[0], config->config_data, "origin_x");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, plain_coords[1], config->config_data, "origin_y");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, plain_coords[2], config->config_data, "origin_z");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, plain_coords[3], config->config_data, "normal_x");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, plain_coords[4], config->config_data, "normal_y");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, plain_coords[5], config->config_data, "normal_z");
    }

    if(clip_with_bounds)
    {
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, bounds[0], config->config_data, "min_x");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, bounds[1], config->config_data, "min_y");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, bounds[2], config->config_data, "min_z");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, bounds[3], config->config_data, "max_x");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, bounds[4], config->config_data, "max_y");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, bounds[5], config->config_data, "max_z");
    }


    write_transmembrane_potential_vtu(&vtk_grid, the_grid, output_dir, file_prefix, time_info->iteration, time_info->current_t,\
                                            save_pvd, compress, compression_level, binary, clip_with_plain, plain_coords, clip_with_bounds, bounds);
    write_transmembrane_potential_vtp(&vtk_polydata, the_grid, output_dir, file_prefix_purkinje, time_info->iteration, time_info->current_t,\
                                            save_pvd, compress, compression_level, binary, clip_with_plain, plain_coords, clip_with_bounds, bounds);


    if(the_grid->adaptive)
        free_vtk_unstructured_grid(vtk_grid);

}

struct save_with_activation_times_persistent_data {
    struct point_hash_entry *last_time_v;
    struct point_hash_entry *num_activations;
    struct point_hash_entry *cell_was_active;
    struct point_voidp_hash_entry *activation_times;
    struct point_voidp_hash_entry *apds;
    bool first_save_call;

};

INIT_SAVE_MESH(init_save_with_activation_times) {

    config->persistent_data = calloc(1, sizeof(struct save_with_activation_times_persistent_data));
    hmdefault(((struct save_with_activation_times_persistent_data*)config->persistent_data)->cell_was_active, 0.0);
    hmdefault(((struct save_with_activation_times_persistent_data*)config->persistent_data)->last_time_v, -100.0);
    hmdefault(((struct save_with_activation_times_persistent_data*)config->persistent_data)->num_activations, 0);
    hmdefault(((struct save_with_activation_times_persistent_data*)config->persistent_data)->activation_times, NULL);
    hmdefault(((struct save_with_activation_times_persistent_data*)config->persistent_data)->apds, NULL);
    ((struct save_with_activation_times_persistent_data*)config->persistent_data)->first_save_call = true;

}

END_SAVE_MESH(end_save_with_activation_times) {
    free(config->persistent_data);
}

SAVE_MESH(save_with_activation_times) {

    int mesh_output_pr = 0;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(int, mesh_output_pr, config->config_data, "mesh_print_rate");

    int iteration_count = time_info->iteration;

    if(mesh_output_pr) {
        if (iteration_count % mesh_output_pr == 0)
            save_as_text_or_binary(time_info, config, the_grid, ode_solver);
    }

    float time_threshold = 10.0f;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(float, time_threshold, config->config_data, "time_threshold");

    GET_PARAMETER_VALUE_CHAR_OR_REPORT_ERROR(output_dir, config->config_data, "output_dir");

    float activation_threshold = -30.0f;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(float, activation_threshold, config->config_data, "activation_threshold");

    float apd_threshold = -83.0f;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(float, apd_threshold, config->config_data, "apd_threshold");

    real_cpu current_t = time_info->current_t;
    real_cpu last_t = time_info->final_t;
    real_cpu dt = time_info->dt;

    sds output_dir_with_file = sdsnew(output_dir);
    output_dir_with_file = sdscat(output_dir_with_file, "/");
    sds base_name = create_base_name("activation_info", 0, "acm");
    output_dir_with_file = sdscatprintf(output_dir_with_file, base_name, current_t);

    struct save_with_activation_times_persistent_data *persistent_data =
            (struct save_with_activation_times_persistent_data*)config->persistent_data;

    struct cell_node *grid_cell = the_grid->first_cell;

    real_cpu center_x, center_y, center_z, dx, dy, dz;
    real_cpu v;

    FILE *act_file = fopen(output_dir_with_file, "w");

    fprintf(act_file, "%d\n", (last_t-current_t) <= dt ); //rounding errors

    while(grid_cell != 0) {

        if( grid_cell->active || ( grid_cell->mesh_extra_info && ( FIBROTIC(grid_cell) || BORDER_ZONE(grid_cell) ) ) ) {

            center_x = grid_cell->center.x;
            center_y = grid_cell->center.y;
            center_z = grid_cell->center.z;

            v = grid_cell->v;

            struct point_3d cell_coordinates;
            cell_coordinates.x = center_x;
            cell_coordinates.y = center_y;
            cell_coordinates.z = center_z;

            dx = grid_cell->discretization.x / 2.0;
            dy = grid_cell->discretization.y / 2.0;
            dz = grid_cell->discretization.z / 2.0;

            fprintf(act_file, "%g,%g,%g,%g,%g,%g,%d,%d,%d ", center_x, center_y, center_z, dx, dy, dz, grid_cell->active, FIBROTIC(grid_cell), BORDER_ZONE(grid_cell));

            int n_activations = 0;
            float *apds_array = NULL;
            float *activation_times_array = NULL;

            if(grid_cell->active) {

                float last_v = hmget(persistent_data->last_time_v, cell_coordinates);

                n_activations = (int) hmget(persistent_data->num_activations, cell_coordinates);
                activation_times_array = (float *) hmget(persistent_data->activation_times, cell_coordinates);
                apds_array = (float *) hmget(persistent_data->apds, cell_coordinates);

                int act_times_len = arrlen(activation_times_array);

                if (current_t == 0.0f) {
                    hmput(persistent_data->last_time_v, cell_coordinates, v);
                } else {
                    if ((last_v < activation_threshold) && (v >= activation_threshold)) {

                        if (act_times_len == 0) {
                            n_activations++;
                            hmput(persistent_data->num_activations, cell_coordinates, n_activations);
                                    arrput(activation_times_array, current_t);
                            float tmp = hmget(persistent_data->cell_was_active, cell_coordinates);
                            hmput(persistent_data->cell_was_active, cell_coordinates, tmp + 1);
                            hmput(persistent_data->activation_times, cell_coordinates, activation_times_array);
                        } else { //This is to avoid spikes in the middle of an Action Potential
                            float last_act_time = activation_times_array[act_times_len - 1];
                            if (current_t - last_act_time > time_threshold) {
                                n_activations++;
                                hmput(persistent_data->num_activations, cell_coordinates, n_activations);
                                        arrput(activation_times_array, current_t);
                                float tmp = hmget(persistent_data->cell_was_active, cell_coordinates);
                                hmput(persistent_data->cell_was_active, cell_coordinates, tmp + 1);
                                hmput(persistent_data->activation_times, cell_coordinates, activation_times_array);
                            }
                        }
                    }

                    //CHECK APD
                    bool was_active = (hmget(persistent_data->cell_was_active, cell_coordinates) != 0.0);
                    if (was_active) {
                        if (v <= apd_threshold || (hmget(persistent_data->cell_was_active, cell_coordinates) == 2.0) || (last_t-current_t) <= dt) {

                            int tmp = (int)hmget(persistent_data->cell_was_active, cell_coordinates);
                            int act_time_array_len = arrlen(activation_times_array);
                            //if this in being calculated because we had a new activation before the cell achieved the rest potential,
                            // we need to get the activation before this one
                            real_cpu last_act_time = activation_times_array[act_time_array_len  - tmp];
                            real_cpu apd = current_t - last_act_time;
                                    arrput(apds_array, apd);
                            hmput(persistent_data->apds, cell_coordinates, apds_array);
                            hmput(persistent_data->cell_was_active, cell_coordinates, tmp - 1);
                        }
                    }

                    hmput(persistent_data->last_time_v, cell_coordinates, v);
                }
            }

            fprintf(act_file, "%d [ ", n_activations);

            for (unsigned long i = 0; i < n_activations; i++) {
                fprintf(act_file, "%lf ", activation_times_array[i]);
            }
            fprintf(act_file, "] ");

            fprintf(act_file, "[ ");

            for (unsigned long i = 0; i < arrlen(apds_array); i++) {
                fprintf(act_file, "%lf ", apds_array[i]);
            }
            fprintf(act_file, "]\n");
        }

        grid_cell = grid_cell->next;
    }

    fclose(act_file);


}

SAVE_MESH(no_save) {
    //Nop
}

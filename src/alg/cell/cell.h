//
// Created by sachetto on 29/09/17.
//

#ifndef MONOALG3D_CELL_H
#define MONOALG3D_CELL_H

#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "../../common_types/common_types.h"
#include "../../monodomain/constants.h"

#if defined(_OPENMP)
#include <omp.h>
#endif

#define CELL_NODE_TYPE 'b'
#define TRANSITION_NODE_TYPE 'w'
#define CELL_CENTER_EQUALS(cell, x1, y1, z1) (cell->center.x == x1 && cell->center.y == y1 && cell->center.z == z1)


struct element {
#ifdef ENABLE_DDM
    char direction;     // NEW parameter !!!
#endif
    real_cpu value;
    uint32_t column; // Column of the matrix to which this element belongs.
    struct cell_node *cell;
};

struct basic_cell_data {
    char type;
    uint8_t level; // This should be enough for the refinement levels
};

struct cell_node {

    struct basic_cell_data cell_data; // DO NOT CHANGE THIS MEMBER POSITION

    uint64_t bunch_number; // Bunch identifier

    struct point_3d center;

    void *z_front; // Points to cell node or transition node above this cell. Z front
    void *z_back; // Points to cell node or transition node below this cell. Z back
    void *y_top;  // Points to cell node or transition node rightward this cell.Y top
    void *y_down;  // Points to cell node or transition node leftward this cell. Y down
    void *x_right; // Points to cell node or transition node in front of this cell. X right
    void *x_left;  // Points to cell node or transition node behind this cell. X left

    struct cell_node *previous; // Previous cell in the Hilbert curve ordering.
    struct cell_node *next;     // Next cell of in the Hilbert curve ordering.

    // Indicates position of cell on grid according to ordering provided by
    // the modified Hilbert curve.
    uint32_t grid_position;

    // Variable used to storage the form of the Hilbert curve in the  bunch
    // created when this cell is refined.
    uint8_t hilbert_shape_number;

    // Fluxes used to decide if a cell should be refined or if a bunch
    // should be derefined.
    real_cpu z_front_flux, // Flux coming from north direction.
        z_back_flux,   // Flux coming from south direction.
        y_top_flux,    // Flux coming from east direction.
        y_down_flux,    // Flux coming from west direction.
        x_right_flux,   // Flux coming from front direction.
        x_left_flux;    // Flux coming from back direction.

    /* The matrix row. The elements[0] corresponds to the diagonal element of the row. */
    element_array elements;

    struct point_3d discretization;

    bool active;
    bool can_change;
    bool visited;

    //______________________________________________________________________________
    /* Variables used in solving the discretized system Ax = b through the conjugate gradient
   method.
   The grid discretization matrix and its resolution are directly implemented on the grid,
   which improves performance. There is no independent linear algebra package. */
    real_cpu Ax; /* Element of vector Ax = b associated to this cell. Also plays the role of Ap.*/
    real_cpu b;  /* In Ax = b, corresponds to the element in vector b associated to this cell. */

    void *linear_system_solver_extra_info;
    size_t linear_system_solver_extra_info_size;

    uint32_t sv_position;

    void *mesh_extra_info;
    size_t mesh_extra_info_size;

    // Variables used by some applications of partial differential equations.
    real_cpu v;

    struct point_3d sigma;

#ifdef ENABLE_DDM
    struct point_3d kappa;
#endif


#if defined(_OPENMP)
    omp_lock_t updating;
#endif

};

// ----------------------------------------------------

struct transition_node {
    struct basic_cell_data cell_data; // DO NOT CHANGE THIS STRUCT POSITION

    void *single_connector;
    void *quadruple_connector1;
    void *quadruple_connector2;
    void *quadruple_connector3;
    void *quadruple_connector4;

    /* Directions that a transition node may assume:
     * 'e': east
     * 'w': west
     * 'n': north
     * 's': south
     * 'f': front
     * 'b': back
     */
    char direction;
};

// ----------------------------------------------------
struct terminal {
    bool active;

    struct node *purkinje_cell;

    struct cell_node **tissue_cells;

};
// ----------------------------------------------------

void init_basic_cell_data_with_type(struct basic_cell_data *data, char type);

struct cell_node* new_cell_node();

void init_cell_node (struct cell_node *cell_node);

void free_cell_node (struct cell_node *cell_node);

void lock_cell_node (struct cell_node *cell_node);

void unlock_cell_node (struct cell_node *cell_node);

struct transition_node* new_transition_node();

void init_transition_node (struct transition_node *transition_node);

void set_transition_node_data (struct transition_node *the_transition_node, uint16_t level,
                               char direction, void *single_connector, void *quadruple_connector1,
                               void *quadruple_connector2, void *quadruple_connector3,
                               void *quadruple_connector4);

void set_cell_node_data(struct cell_node *the_cell, struct point_3d discretization,
                        uint64_t bunch_number, void *east, void *north, void *west, void *south,
                        void *front, void *back, void *previous, void *next,
                        uint32_t grid_position, uint8_t hilbert_shape_number, struct point_3d center);

void set_cell_flux (struct cell_node *the_cell, char direction);
real_cpu get_cell_maximum_flux (struct cell_node *the_cell);

void set_refined_cell_data (struct cell_node *the_cell, struct cell_node *other_cell,
                            struct point_3d discretization, struct point_3d center,
                            uint64_t bunch_number, ui32_array free_sv_positions, ui32_array *refined_this_step);

void set_refined_transition_node_data (struct transition_node *the_node,
                                       struct cell_node *other_node, char direction);

void simplify_refinement (struct transition_node *transition_node);
void refine_cell (struct cell_node *cell, ui32_array free_sv_positions,
                  ui32_array *refined_this_step);

bool cell_needs_derefinement (struct cell_node *grid_cell, real_cpu derefinement_bound);
struct cell_node *get_front_northeast_cell (struct cell_node *first_bunch_cell);
uint8_t get_father_bunch_number (struct cell_node *first_bunch_cell);
void simplify_derefinement(struct transition_node *transition_node);

void derefine_cell_bunch (struct cell_node *first_bunch_cell, ui32_array *free_sv_positions);

struct cell_node * get_cell_neighbour(struct cell_node *grid_cell, void *neighbour_grid_cell);
bool cell_has_neighbour(struct cell_node *grid_cell, void *neighbour_grid_cell);
void * get_cell_neighbour_as_void(struct cell_node *grid_cell, void *neighbour_grid_cell, char *type);


#endif // MONOALG3D_CELL_H

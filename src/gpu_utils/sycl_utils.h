//
// Created by bergolho on 19/02/25.
//

#ifndef SYCL_UTILS_H
#define SYCL_UTILS_H

#include <stdlib.h>
#include <stdio.h>
#include <stddef.h>
#include <stdint.h>

#define __SYCLCC__ 2

#include "../common_types/common_types.h"

#ifdef COMPILE_SYCL 
#include <sycl/sycl.hpp>
#include <dpct/dpct.hpp>
//#include <dpct/sparse_utils.hpp>
//#include <dpct/blas_utils.hpp>
#endif

void hello_sycl ();

enum copy_direction {
    HOST_TO_DEVICE,
    DEVICE_TO_HOST,
};

enum sparse_index_type {
    INDEX_INT32,
    INDEX_INT64,
};

enum sparse_index_base {
    INDEX_BASE_ZERO,
};

enum data_type {
    REAL_FLOAT,
    REAL_DOUBLE,
};

//#ifdef __cplusplus
//extern "C" {
//#endif
//struct sycl_manager {
//    sycl::queue queue;
//    dpct::device_ext *dev_ct1;
//
//    void set_queue (const bool use_gpu) {
//        int device_id;
//        // GPU
//        if (use_gpu) {
//            device_id = 0;
//        }
//        // CPU
//        else {
//            device_id = 1;
//        }   
//        dev_ct1 = &dpct::get_device(device_id);
//        queue = dev_ct1->in_order_queue();
//    }
//};
//#ifdef __cplusplus
//}
//#endif


void malloc_device_sycl(void **ptr, size_t n);
void free_device_sycl(void *ptr);
void memcpy_device(void *dest, const void *src, size_t n, enum copy_direction kind);
/*
void create_sparse_handle(void **handle);
void create_blas_handle(void **handle);
void sparse_create_scr(void *mat, int64_t rows, int64_t cols, int64_t nnz,
                       void* row_ptr,
                       void* cols_ind_ptr,
                       void* vals_ptr,
                       enum sparse_index_type csr_row_offsets_type,
                       enum sparse_index_type csr_col_ind_type,
                       enum sparse_index_base idx_base,
                       enum data_type         value_type);

void create_dense_vector(void** descr, int64_t size, void *values, enum data_type valueType);
*/

/*
#ifdef COMPILE_SYCL
void sparse_spmv(dpct::sparse::descriptor_ptr handle, oneapi::mkl::transpose opA,
                            const void *alpha, dpct::sparse::sparse_matrix_desc_t matA,
                            std::shared_ptr<dpct::sparse::dense_vector_desc> vecX, const void *beta,
                            std::shared_ptr<dpct::sparse::dense_vector_desc> vecY,
                            enum data_type computeType, void *externalBuffer);
#endif

void blas_dot(void *handle, int n, real *x, int incx, real *y, int incy, real *result);
*/

#endif //SYCL_UTILS_H

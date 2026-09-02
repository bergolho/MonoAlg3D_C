//
// Created by sachetto on 21/08/2020.
//

struct gpu_persistent_data {
    int *d_col, *d_row;
    float *d_val, *d_x;

    float *d_r, *d_p, *d_Ax, *d_rw, *d_v, *d_t, *d_omega;

    int N, nz;

    /* Get handle to the CUBLAS context */
    cublasHandle_t cublasHandle;

    /* Get handle to the CUSPARSE context */
    cusparseHandle_t cusparseHandle;

    cusparseSpMatDescr_t matA;
    cusparseDnVecDescr_t vecx;
    cusparseDnVecDescr_t vecp;
    cusparseDnVecDescr_t vecv;
    cusparseDnVecDescr_t vecr;
    cusparseDnVecDescr_t vect;
    cusparseDnVecDescr_t vecAx;
    cusparseDnVecDescr_t vecomega;

    float *d_valsILU0;
    float *d_zm1, *d_zm2, *d_rm2;
    float *d_y;
    
    int                 bufferSizeLU;
    size_t              bufferSizeMV, bufferSizeL, bufferSizeU;
    void*               d_bufferLU, *d_bufferMV,  *d_bufferL, *d_bufferU;
    cusparseSpSVDescr_t spsvDescrL, spsvDescrU;
    cusparseMatDescr_t   matLU;
    csrilu02Info_t      infoILU;

    cusparseSpMatDescr_t matM_lower, matM_upper;  

    /* Double-precision data used by the GPU BiCGSTAB solver. */
    double *d_val_bicg;
    double *d_x_bicg, *d_r_bicg, *d_rw_bicg;
    double *d_p_bicg, *d_v_bicg, *d_t_bicg;
    double *d_Ax_bicg;

    cusparseSpMatDescr_t matA_bicg;
    cusparseDnVecDescr_t vecx_bicg;
    cusparseDnVecDescr_t vecp_bicg;
    cusparseDnVecDescr_t vecv_bicg;
    cusparseDnVecDescr_t vecr_bicg;
    cusparseDnVecDescr_t vect_bicg;
    cusparseDnVecDescr_t vecAx_bicg;

    size_t bufferSizeMV_bicg;
    void *d_bufferMV_bicg;

};

const float floatone = 1.0f;
const float floatzero = 0.0f;

INIT_LINEAR_SYSTEM(init_gpu_conjugate_gradient) {

    if(the_grid->adaptive) {
        log_error_and_exit("The gpu conjugate gradient does not support mesh adaptivity. Aborting!\n");
    }

    struct gpu_persistent_data *persistent_data = CALLOC_ONE_TYPE(struct gpu_persistent_data);

    int_array I = NULL, J = NULL;
    f32_array val = NULL;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real_cpu, tol, config, "tolerance");
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(int, max_its, config, "max_iterations");

    check_cublas_error(cublasCreate(&(persistent_data->cublasHandle)));
    check_cublas_error(cusparseCreate(&(persistent_data->cusparseHandle)));

    uint32_t num_active_cells;
    struct cell_node **active_cells = NULL;

    if (is_purkinje) {
        grid_to_csr(the_grid, &val, &I, &J, true);
        num_active_cells = the_grid->purkinje->num_active_purkinje_cells;
        active_cells = the_grid->purkinje->purkinje_cells;
    } else {
        grid_to_csr(the_grid, &val, &I, &J, false);
        num_active_cells = the_grid->num_active_cells;
        active_cells = the_grid->active_cells;
    }

    int nz = arrlen(val);
    int N = num_active_cells;

    persistent_data->nz = nz;
    persistent_data->N = N;

    check_cuda_error(cudaMalloc((void **)&(persistent_data->d_col), nz * sizeof(int)));
    check_cuda_error(cudaMalloc((void **)&(persistent_data->d_row), (N + 1) * sizeof(int)));
    check_cuda_error(cudaMalloc((void **)&(persistent_data->d_val), nz * sizeof(float)));
    check_cuda_error(cudaMalloc((void **)&(persistent_data->d_x), N * sizeof(float)));
    check_cuda_error(cudaMalloc((void **)&(persistent_data->d_r), N * sizeof(float)));
    check_cuda_error(cudaMalloc((void **)&(persistent_data->d_p), N * sizeof(float)));
    check_cuda_error(cudaMalloc((void **)&(persistent_data->d_Ax), N * sizeof(float)));

    /* Wrap raw data into cuSPARSE generic API objects */
    check_cuda_error(cusparseCreateCsr(&(persistent_data->matA), N, N, nz, persistent_data->d_row, persistent_data->d_col, persistent_data->d_val, CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I, CUSPARSE_INDEX_BASE_ZERO, CUDA_R_32F));
    check_cuda_error(cusparseCreateDnVec(&(persistent_data->vecx), N,  persistent_data->d_x, CUDA_R_32F));
    check_cuda_error(cusparseCreateDnVec(&(persistent_data->vecp), N,  persistent_data->d_p, CUDA_R_32F));
    check_cuda_error(cusparseCreateDnVec(&(persistent_data->vecAx), N, persistent_data->d_Ax, CUDA_R_32F));

    cudaMemcpy(persistent_data->d_col, J, nz * sizeof(int), cudaMemcpyHostToDevice);      // JA
    cudaMemcpy(persistent_data->d_row, I, (N + 1) * sizeof(int), cudaMemcpyHostToDevice); // IA
    cudaMemcpy(persistent_data->d_val, val, nz * sizeof(float), cudaMemcpyHostToDevice);  // A
    float *rhs = (float *)malloc(sizeof(float) * num_active_cells);

    OMP(parallel for)
    for(uint32_t i = 0; i < num_active_cells; i++) {
        rhs[i] = active_cells[i]->b;
    }

    check_cuda_error(cudaMemcpy(persistent_data->d_x, rhs, N * sizeof(float), cudaMemcpyHostToDevice)); // Result

    float alpha = 1.0f;
    float beta = 0.0f;
    
    check_cuda_error(cusparseSpMV_bufferSize(persistent_data->cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, &alpha, persistent_data->matA, persistent_data->vecx, &beta, persistent_data->vecAx,
                     CUDA_R_32F, CUSPARSE_SPMV_ALG_DEFAULT, &(persistent_data->bufferSizeMV)));
    check_cuda_error(cudaMalloc(&(persistent_data->d_bufferMV), persistent_data->bufferSizeMV));

    config->persistent_data = persistent_data;
    free(rhs);
    arrfree(I);
    arrfree(J);
    arrfree(val);
}

END_LINEAR_SYSTEM(end_gpu_conjugate_gradient) {

    struct gpu_persistent_data *persistent_data = (struct gpu_persistent_data *)config->persistent_data;

    if(!persistent_data) return;

    check_cuda_error((cudaError_t)cusparseDestroy(persistent_data->cusparseHandle));
    check_cuda_error((cudaError_t)cublasDestroy(persistent_data->cublasHandle));

    if (persistent_data->matA)  { check_cuda_error(cusparseDestroySpMat(persistent_data->matA)); }
    if (persistent_data->vecx)  { check_cuda_error(cusparseDestroyDnVec(persistent_data->vecx)); }
    if (persistent_data->vecAx) { check_cuda_error(cusparseDestroyDnVec(persistent_data->vecAx)); }
    if (persistent_data->vecp)  { check_cuda_error(cusparseDestroyDnVec(persistent_data->vecp)); }

    check_cuda_error(cudaFree(persistent_data->d_col));
    check_cuda_error(cudaFree(persistent_data->d_row));
    check_cuda_error(cudaFree(persistent_data->d_val));
    check_cuda_error(cudaFree(persistent_data->d_x));
    check_cuda_error(cudaFree(persistent_data->d_r));
    check_cuda_error(cudaFree(persistent_data->d_p));
    check_cuda_error(cudaFree(persistent_data->d_Ax));

    free(persistent_data);
}

SOLVE_LINEAR_SYSTEM(gpu_conjugate_gradient) {

    /* Conjugate gradient.
       ------------------------------------------
       Follows the description by Golub & Van Loan, "Matrix Computations 3rd ed.", Section 10.2.6
    */
    struct gpu_persistent_data *persistent_data = (struct gpu_persistent_data *)config->persistent_data;

    if(!persistent_data) {
        log_error_and_exit("The gpu_conjugate_gradient solver needs to be initialized before being called. Add a init_function in the [linear_system_solver] section of the .ini file!\n");
    }

    float dot;

    float a, b, na, r0, r1;

    int k;
    float alpha, alpham1, beta;

    float *rhs; // Vector B
    rhs = (float *)malloc(sizeof(float) * num_active_cells);

    OMP(parallel for)
    for(uint32_t i = 0; i < num_active_cells; i++) {
        rhs[i] = active_cells[i]->b;
    }

    int N = persistent_data->N;    

    cudaMemcpy(persistent_data->d_r, rhs, N * sizeof(float), cudaMemcpyHostToDevice); // B

    alpha = 1.0f;
    alpham1 = -1.0f;
    beta = 0.0f;
    r0 = 0.0f;    

    check_cuda_error(cusparseSpMV(persistent_data->cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, &alpha, persistent_data->matA, persistent_data->vecx, &beta, persistent_data->vecAx, CUDA_R_32F,
                                  CUSPARSE_SPMV_ALG_DEFAULT, persistent_data->d_bufferMV));

    cublasSaxpy(persistent_data->cublasHandle, N, &alpham1, persistent_data->d_Ax, 1, persistent_data->d_r, 1);
    cublasSdot(persistent_data->cublasHandle, N, persistent_data->d_r, 1, persistent_data->d_r, 1, &r1);

    k = 1;

    while(r1 >= tol && k <= max_its) {

        if(k > 1) {
            b = r1 / r0;
            cublasSscal(persistent_data->cublasHandle, N, &b, persistent_data->d_p, 1);
            cublasSaxpy(persistent_data->cublasHandle, N, &alpha, persistent_data->d_r, 1, persistent_data->d_p, 1);
 
        } else {
            cublasScopy(persistent_data->cublasHandle, N, persistent_data->d_r, 1, persistent_data->d_p, 1);
        }

        check_cuda_error(cusparseSpMV(persistent_data->cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, &alpha, persistent_data->matA, persistent_data->vecp, &beta, persistent_data->vecAx, CUDA_R_32F,
                                      CUSPARSE_SPMV_ALG_DEFAULT, persistent_data->d_bufferMV));
        
        cublasSdot(persistent_data->cublasHandle, N, persistent_data->d_p, 1, persistent_data->d_Ax, 1, &dot);

        a = r1 / dot;

        cublasSaxpy(persistent_data->cublasHandle, N, &a, persistent_data->d_p, 1, persistent_data->d_x, 1);
        na = -a;
        cublasSaxpy(persistent_data->cublasHandle, N, &na, persistent_data->d_Ax, 1, persistent_data->d_r, 1);

        r0 = r1;
        cublasSdot(persistent_data->cublasHandle, N, persistent_data->d_r, 1, persistent_data->d_r, 1, &r1);


        cudaDeviceSynchronize();
        k++;
    }

    cudaMemcpy(rhs, persistent_data->d_x, N * sizeof(float), cudaMemcpyDeviceToHost);

    *number_of_iterations = k - 1;
    *error = r1;

    OMP(parallel for)
    for(uint32_t i = 0; i < num_active_cells; i++) {
        active_cells[i]->v = rhs[i];
    }

    free(rhs);
}


INIT_LINEAR_SYSTEM(init_gpu_biconjugate_gradient) {

    struct gpu_persistent_data *persistent_data = CALLOC_ONE_TYPE(struct gpu_persistent_data);

    int_array I = NULL, J = NULL;
    f32_array val = NULL;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real_cpu, tol, config, "tolerance");
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(int, max_its, config, "max_iterations");

    check_cuda_error((cudaError_t)cublasCreate(&(persistent_data->cublasHandle)));
    check_cuda_error((cudaError_t)cusparseCreate(&(persistent_data->cusparseHandle)));

    uint32_t num_active_cells;
    struct cell_node **active_cells = NULL;

    if (is_purkinje) {
        grid_to_csr(the_grid, &val, &I, &J, true);
        num_active_cells = the_grid->purkinje->num_active_purkinje_cells;
        active_cells = the_grid->purkinje->purkinje_cells;
    } else {
        grid_to_csr(the_grid, &val, &I, &J, false);
        num_active_cells = the_grid->num_active_cells;
        active_cells = the_grid->active_cells;
    }

    int nz = arrlen(val);
    int N = num_active_cells;

    persistent_data->N = N;
    persistent_data->nz = nz;

    /* CSR structure is shared with the other GPU solvers. */
    check_cuda_error(cudaMalloc((void **)&(persistent_data->d_col), nz * sizeof(int)));
    check_cuda_error(cudaMalloc((void **)&(persistent_data->d_row), (N + 1) * sizeof(int)));

    /* BiCGSTAB uses a fully double-precision matrix/vector path. */
    check_cuda_error(cudaMalloc((void **)&(persistent_data->d_val_bicg), nz * sizeof(double)));
    check_cuda_error(cudaMalloc((void **)&(persistent_data->d_x_bicg), N * sizeof(double)));
    check_cuda_error(cudaMalloc((void **)&(persistent_data->d_r_bicg), N * sizeof(double)));
    check_cuda_error(cudaMalloc((void **)&(persistent_data->d_rw_bicg), N * sizeof(double)));
    check_cuda_error(cudaMalloc((void **)&(persistent_data->d_p_bicg), N * sizeof(double)));
    check_cuda_error(cudaMalloc((void **)&(persistent_data->d_v_bicg), N * sizeof(double)));
    check_cuda_error(cudaMalloc((void **)&(persistent_data->d_t_bicg), N * sizeof(double)));
    check_cuda_error(cudaMalloc((void **)&(persistent_data->d_Ax_bicg), N * sizeof(double)));

    check_cuda_error(cudaMemcpy(persistent_data->d_col, J, nz * sizeof(int), cudaMemcpyHostToDevice));
    check_cuda_error(cudaMemcpy(persistent_data->d_row, I, (N + 1) * sizeof(int), cudaMemcpyHostToDevice));

    /* Convert the float CSR values produced by grid_to_csr() to double once at initialization. */
    double *val_double = (double *)malloc(sizeof(double) * nz);
    if (!val_double) {
        log_error_and_exit("Unable to allocate host memory for the double-precision BiCGSTAB matrix.\n");
    }

    OMP(parallel for)
    for (int i = 0; i < nz; i++) {
        val_double[i] = (double)val[i];
    }

    check_cuda_error(cudaMemcpy(persistent_data->d_val_bicg, val_double,
                                nz * sizeof(double), cudaMemcpyHostToDevice));

    free(val_double);

    check_cuda_error(cudaMemset(persistent_data->d_r_bicg, 0, N * sizeof(double)));
    check_cuda_error(cudaMemset(persistent_data->d_rw_bicg, 0, N * sizeof(double)));
    check_cuda_error(cudaMemset(persistent_data->d_p_bicg, 0, N * sizeof(double)));
    check_cuda_error(cudaMemset(persistent_data->d_t_bicg, 0, N * sizeof(double)));
    check_cuda_error(cudaMemset(persistent_data->d_v_bicg, 0, N * sizeof(double)));

    double *rhs = (double *)malloc(sizeof(double) * num_active_cells);
    if (!rhs) {
        log_error_and_exit("Unable to allocate host RHS for the double-precision BiCGSTAB solver.\n");
    }

    OMP(parallel for)
    for (uint32_t i = 0; i < num_active_cells; i++) {
        rhs[i] = (double)active_cells[i]->b;
    }

    /* Preserve the original initial guess x0 = b used by this solver. */
    check_cuda_error(cudaMemcpy(persistent_data->d_x_bicg, rhs,
                                N * sizeof(double), cudaMemcpyHostToDevice));

    double alpha = 1.0;
    double beta = 0.0;

    check_cublas_error(cusparseCreateCsr(&(persistent_data->matA_bicg),
                                         N, N, nz,
                                         persistent_data->d_row,
                                         persistent_data->d_col,
                                         persistent_data->d_val_bicg,
                                         CUSPARSE_INDEX_32I,
                                         CUSPARSE_INDEX_32I,
                                         CUSPARSE_INDEX_BASE_ZERO,
                                         CUDA_R_64F));

    check_cublas_error(cusparseCreateDnVec(&(persistent_data->vecx_bicg),
                                           N, persistent_data->d_x_bicg, CUDA_R_64F));
    check_cublas_error(cusparseCreateDnVec(&(persistent_data->vecp_bicg),
                                           N, persistent_data->d_p_bicg, CUDA_R_64F));
    check_cublas_error(cusparseCreateDnVec(&(persistent_data->vecv_bicg),
                                           N, persistent_data->d_v_bicg, CUDA_R_64F));
    check_cublas_error(cusparseCreateDnVec(&(persistent_data->vecr_bicg),
                                           N, persistent_data->d_r_bicg, CUDA_R_64F));
    check_cublas_error(cusparseCreateDnVec(&(persistent_data->vect_bicg),
                                           N, persistent_data->d_t_bicg, CUDA_R_64F));
    check_cublas_error(cusparseCreateDnVec(&(persistent_data->vecAx_bicg),
                                           N, persistent_data->d_Ax_bicg, CUDA_R_64F));

    check_cublas_error(cusparseSpMV_bufferSize(persistent_data->cusparseHandle,
                                               CUSPARSE_OPERATION_NON_TRANSPOSE,
                                               &alpha,
                                               persistent_data->matA_bicg,
                                               persistent_data->vecx_bicg,
                                               &beta,
                                               persistent_data->vecAx_bicg,
                                               CUDA_R_64F,
                                               CUSPARSE_SPMV_ALG_DEFAULT,
                                               &(persistent_data->bufferSizeMV_bicg)));

    check_cuda_error(cudaMalloc(&(persistent_data->d_bufferMV_bicg),
                                persistent_data->bufferSizeMV_bicg));

    arrfree(I);
    arrfree(J);
    arrfree(val);

    config->persistent_data = persistent_data;
}

SOLVE_LINEAR_SYSTEM(gpu_biconjugate_gradient) {

    struct gpu_persistent_data *persistent_data = (struct gpu_persistent_data *)config->persistent_data;

    if(!persistent_data) {
        log_error_and_exit("The gpu_biconjugate_gradient solver needs to be initialized before being called. Add a init_function in the [linear_system_solver] section of the .ini file\n");
    }

    double rho = 0.0;
    double rhop = 0.0;
    double beta = 0.0;
    double alpha = 1.0;
    double negalpha = 0.0;
    double omega = 1.0;
    double negomega = 0.0;
    double temp = 0.0;
    double temp2 = 0.0;
    double nrmr = 0.0;
    double nrmr0 = 0.0;

    const double zero = 0.0;
    const double one = 1.0;
    const double mone = -1.0;

    double *rhs = (double *)malloc(sizeof(double) * num_active_cells);
    if (!rhs) {
        log_error_and_exit("Unable to allocate host RHS for the double-precision BiCGSTAB solver.\n");
    }

    OMP(parallel for)
    for(uint32_t i = 0; i < num_active_cells; i++) {
        rhs[i] = (double)active_cells[i]->b;
    }

    int N = persistent_data->N;

    check_cuda_error(cudaMemcpy(persistent_data->d_r_bicg, rhs,
                                N * sizeof(double), cudaMemcpyHostToDevice));

    /* Compute initial residual r0 = b - A*x0. */
    check_cublas_error(cusparseSpMV(persistent_data->cusparseHandle,
                                    CUSPARSE_OPERATION_NON_TRANSPOSE,
                                    &one,
                                    persistent_data->matA_bicg,
                                    persistent_data->vecx_bicg,
                                    &zero,
                                    persistent_data->vecAx_bicg,
                                    CUDA_R_64F,
                                    CUSPARSE_SPMV_ALG_DEFAULT,
                                    persistent_data->d_bufferMV_bicg));

    check_cublas_error(cublasDaxpy(persistent_data->cublasHandle,
                                   N, &mone,
                                   persistent_data->d_Ax_bicg, 1,
                                   persistent_data->d_r_bicg, 1));

    /* r_hat = r0 and p0 = r0. */
    check_cublas_error(cublasDcopy(persistent_data->cublasHandle,
                                   N, persistent_data->d_r_bicg, 1,
                                   persistent_data->d_rw_bicg, 1));
    check_cublas_error(cublasDcopy(persistent_data->cublasHandle,
                                   N, persistent_data->d_r_bicg, 1,
                                   persistent_data->d_p_bicg, 1));
    check_cublas_error(cublasDnrm2(persistent_data->cublasHandle,
                                   N, persistent_data->d_r_bicg, 1, &nrmr0));

    /* Check if the initial solution is already good enough. */
    if (nrmr0 < tol) {
        *number_of_iterations = 0;
        *error = nrmr0;

        check_cuda_error(cudaMemcpy(rhs,
                                    persistent_data->d_x_bicg,
                                    N * sizeof(double),
                                    cudaMemcpyDeviceToHost));

        OMP(parallel for)
        for(uint32_t j = 0; j < num_active_cells; j++) {
            active_cells[j]->v = (real_cpu)rhs[j];
        }

        free(rhs);
        return;
    }

    int i;
    for(i = 0; i < max_its;) {
        rhop = rho;

        check_cublas_error(cublasDdot(persistent_data->cublasHandle,
                                      N,
                                      persistent_data->d_rw_bicg, 1,
                                      persistent_data->d_r_bicg, 1,
                                      &rho));

        if(i > 0) {
            beta = (rho / rhop) * (alpha / omega);

            negomega = -omega;
            check_cublas_error(cublasDaxpy(persistent_data->cublasHandle,
                                           N, &negomega,
                                           persistent_data->d_v_bicg, 1,
                                           persistent_data->d_p_bicg, 1));
            check_cublas_error(cublasDscal(persistent_data->cublasHandle,
                                           N, &beta,
                                           persistent_data->d_p_bicg, 1));
            check_cublas_error(cublasDaxpy(persistent_data->cublasHandle,
                                           N, &one,
                                           persistent_data->d_r_bicg, 1,
                                           persistent_data->d_p_bicg, 1));
        }

        /* v = A*p */
        check_cublas_error(cusparseSpMV(persistent_data->cusparseHandle,
                                        CUSPARSE_OPERATION_NON_TRANSPOSE,
                                        &one,
                                        persistent_data->matA_bicg,
                                        persistent_data->vecp_bicg,
                                        &zero,
                                        persistent_data->vecv_bicg,
                                        CUDA_R_64F,
                                        CUSPARSE_SPMV_ALG_DEFAULT,
                                        persistent_data->d_bufferMV_bicg));

        check_cublas_error(cublasDdot(persistent_data->cublasHandle,
                                      N,
                                      persistent_data->d_rw_bicg, 1,
                                      persistent_data->d_v_bicg, 1,
                                      &temp));

        alpha = rho / temp;
        negalpha = -alpha;

        check_cublas_error(cublasDaxpy(persistent_data->cublasHandle,
                                       N, &negalpha,
                                       persistent_data->d_v_bicg, 1,
                                       persistent_data->d_r_bicg, 1));
        check_cublas_error(cublasDaxpy(persistent_data->cublasHandle,
                                       N, &alpha,
                                       persistent_data->d_p_bicg, 1,
                                       persistent_data->d_x_bicg, 1));
        check_cublas_error(cublasDnrm2(persistent_data->cublasHandle,
                                       N, persistent_data->d_r_bicg, 1, &nrmr));

        if(nrmr < tol * nrmr0) {
            break;
        }

        /* t = A*r */
        check_cublas_error(cusparseSpMV(persistent_data->cusparseHandle,
                                        CUSPARSE_OPERATION_NON_TRANSPOSE,
                                        &one,
                                        persistent_data->matA_bicg,
                                        persistent_data->vecr_bicg,
                                        &zero,
                                        persistent_data->vect_bicg,
                                        CUDA_R_64F,
                                        CUSPARSE_SPMV_ALG_DEFAULT,
                                        persistent_data->d_bufferMV_bicg));

        check_cublas_error(cublasDdot(persistent_data->cublasHandle,
                                      N,
                                      persistent_data->d_t_bicg, 1,
                                      persistent_data->d_r_bicg, 1,
                                      &temp));
        check_cublas_error(cublasDdot(persistent_data->cublasHandle,
                                      N,
                                      persistent_data->d_t_bicg, 1,
                                      persistent_data->d_t_bicg, 1,
                                      &temp2));

        omega = temp / temp2;
        negomega = -omega;

        check_cublas_error(cublasDaxpy(persistent_data->cublasHandle,
                                       N, &omega,
                                       persistent_data->d_r_bicg, 1,
                                       persistent_data->d_x_bicg, 1));
        check_cublas_error(cublasDaxpy(persistent_data->cublasHandle,
                                       N, &negomega,
                                       persistent_data->d_t_bicg, 1,
                                       persistent_data->d_r_bicg, 1));

        check_cublas_error(cublasDnrm2(persistent_data->cublasHandle,
                                       N, persistent_data->d_r_bicg, 1, &nrmr));

        if(nrmr < tol * nrmr0) {
            i++;
            break;
        }
        i++;
    }

    check_cuda_error(cudaMemcpy(rhs,
                                persistent_data->d_x_bicg,
                                N * sizeof(double),
                                cudaMemcpyDeviceToHost));

    *number_of_iterations = i - 1;
    *error = nrmr;

    OMP(parallel for)
    for(uint32_t j = 0; j < num_active_cells; j++) {
        active_cells[j]->v = (real_cpu)rhs[j];
    }

    free(rhs);
}

END_LINEAR_SYSTEM(end_gpu_biconjugate_gradient) {

    struct gpu_persistent_data *persistent_data = (struct gpu_persistent_data *)config->persistent_data;

    if(!persistent_data) return;

    check_cuda_error((cudaError_t)cusparseDestroy(persistent_data->cusparseHandle));
    check_cuda_error((cudaError_t)cublasDestroy(persistent_data->cublasHandle));

    if (persistent_data->matA_bicg) { check_cuda_error((cudaError_t)cusparseDestroySpMat(persistent_data->matA_bicg)); }
    if (persistent_data->vecx_bicg) { check_cuda_error((cudaError_t)cusparseDestroyDnVec(persistent_data->vecx_bicg)); }
    if (persistent_data->vecAx_bicg) { check_cuda_error((cudaError_t)cusparseDestroyDnVec(persistent_data->vecAx_bicg)); }
    if (persistent_data->vecp_bicg) { check_cuda_error((cudaError_t)cusparseDestroyDnVec(persistent_data->vecp_bicg)); }
    if (persistent_data->vecv_bicg) { check_cuda_error((cudaError_t)cusparseDestroyDnVec(persistent_data->vecv_bicg)); }
    if (persistent_data->vecr_bicg) { check_cuda_error((cudaError_t)cusparseDestroyDnVec(persistent_data->vecr_bicg)); }
    if (persistent_data->vect_bicg) { check_cuda_error((cudaError_t)cusparseDestroyDnVec(persistent_data->vect_bicg)); }

    check_cuda_error(cudaFree(persistent_data->d_col));
    check_cuda_error(cudaFree(persistent_data->d_row));
    check_cuda_error(cudaFree(persistent_data->d_val_bicg));
    check_cuda_error(cudaFree(persistent_data->d_x_bicg));
    check_cuda_error(cudaFree(persistent_data->d_r_bicg));
    check_cuda_error(cudaFree(persistent_data->d_rw_bicg));
    check_cuda_error(cudaFree(persistent_data->d_p_bicg));
    check_cuda_error(cudaFree(persistent_data->d_t_bicg));
    check_cuda_error(cudaFree(persistent_data->d_v_bicg));
    check_cuda_error(cudaFree(persistent_data->d_Ax_bicg));
    check_cuda_error(cudaFree(persistent_data->d_bufferMV_bicg));

    free(persistent_data);
}


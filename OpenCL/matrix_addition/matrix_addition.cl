__kernel void matrix_add(__global float* matAA, __global float* matBB, __global float* matCC) {
    
    // index of current element to be processed
    int i = get_global_id(0);
    
    // do operation
    matCC[i] = matAA[i] + matBB[i];
}

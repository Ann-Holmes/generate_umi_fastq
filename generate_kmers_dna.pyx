from libc.stdlib cimport rand, malloc, free

cpdef str generate_DNA_cy(int nt):
    cdef char *bases = b"ATCG"
    cdef char *result = <char*> malloc(nt + 1)  # 注意：现在导入了malloc
    cdef int i
    for i in range(nt):
        result[i] = bases[rand() % 4]
    result[nt] = '\0'
    
    dna = result[:nt].decode('utf-8')
    
    free(result)
    
    return dna

def generate_kmers(str seq, int read_size):
    cdef int i, n = len(seq) - read_size + 1
    return [seq[i:i + read_size] for i in range(n)]

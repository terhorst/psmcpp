from libcpp.vector cimport vector
from libcpp.pair cimport pair
from libcpp.map cimport map
from libcpp cimport bool
from libcpp.memory cimport unique_ptr
from libcpp.string cimport string

cdef extern from "common.h":
    ctypedef vector[vector[adouble]] ParameterVector
    cdef cppclass Vector[T]:
        T& operator()(int)
        int size()
    adouble double_vec_to_adouble(double x, vector[double] v)
    cdef cppclass adouble:
        adouble()
        adouble(double)
        double value()
        Vector[double] derivatives()
    cdef cppclass Matrix[T]:
        int rows()
        int cols()
        T& operator()(int, int)
    cdef double toDouble(const adouble &)
    void init_eigen()
    void init_logger_cb(void(*)(const char*, const char*, const char*))
    void fill_jacobian(const adouble &, double*)
    void store_matrix(const Matrix[double]&, double*)
    void store_matrix(const Matrix[adouble]&, double*)
    void store_matrix(const Matrix[adouble]&, double*, double*)

cdef extern from "block_key.h":
    cdef cppclass block_key:
        int size() const
        int operator()(int) const

ctypedef Matrix[double]* pMatrixD
ctypedef Matrix[adouble]* pMatrixAd
ctypedef map[block_key, Vector[double]]* pBlockMap

cdef extern from "inference_manager.h":
    cdef cppclass InferenceManager nogil:
        void setTheta(const double)
        void setRho(const double)
        void Estep(bool)
        void setParams(const ParameterVector &) except +
        vector[double] loglik()
        vector[adouble] Q()
        bool debug
        bool saveGamma
        bool folded
        vector[double] hidden_states
        vector[pMatrixD] getGammas()
        vector[pMatrixD] getXisums()
        vector[pBlockMap] getGammaSums()
        Matrix[adouble]& getPi()
        Matrix[adouble]& getTransition()
        Matrix[adouble]& getEmission()
        map[block_key, Vector[adouble]]& getEmissionProbs()
    cdef cppclass OnePopInferenceManager(InferenceManager) nogil:
        OnePopInferenceManager(const int, const vector[int],
                const vector[int*], const vector[double], 
                double* const, const vector[int],
                const bool) except +
    cdef cppclass TwoPopInferenceManager(InferenceManager) nogil:
        TwoPopInferenceManager(const int, const int, const int, const int,
                const vector[int], const vector[int*], const vector[double], 
                double* const, const vector[int], const bool) except +
        void setParams(const ParameterVector&, const ParameterVector&, const ParameterVector&, const double)
    Matrix[adouble] sfs_cython(const int, const ParameterVector, const double, const double, bool) nogil


cdef extern from "piecewise_constant_rate_function.h":
    cdef cppclass PiecewiseConstantRateFunction[T] nogil:
        PiecewiseConstantRateFunction(const ParameterVector, const vector[double])
        T R(T x)
        T random_time(const double, const double, const long long)

# This code is only used for testing purposes
cdef extern from "jcsfs.h":
    cdef cppclass JointCSFS[T] nogil:
        JointCSFS(int, int, int, int, vector[double], int)
        void pre_compute(const ParameterVector&, const ParameterVector&, const double)
        vector[Matrix[T]] compute(const PiecewiseConstantRateFunction[T]&) const

cdef extern from "matrix_cache.h":
    void init_cache(const string)

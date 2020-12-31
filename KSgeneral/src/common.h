#ifndef __common_h__
#define __common_h__

#include <vector>

inline void convolve_same_size(int size, const double* src0, const double* src1, double* dest)
{
    for (int j = 0; j < size; ++j) {
        double convolution_at_j = 0.0;
        for (int k = 0; k <= j; ++k) {
            convolution_at_j += src0[k] * src1[j-k];
        }
        dest[j] = convolution_at_j;
    }
}

template<class T>
class DoubleBuffer {
    public:
        DoubleBuffer(int n, T value);
        std::vector<T>& get_src();
        std::vector<T>& get_dest();
        void flip();

    private:
        std::vector<T> buf0;
        std::vector<T> buf1;
        bool buf0_is_src;
};

template<class T>
DoubleBuffer<T>::DoubleBuffer(int n, T value) :
    buf0(n, value), buf1(n, value), buf0_is_src(true)
{
}

template<class T>
inline std::vector<T>& DoubleBuffer<T>::get_src()
{
    return buf0_is_src ? buf0 : buf1;
}

template<class T>
inline std::vector<T>& DoubleBuffer<T>::get_dest()
{
    return buf0_is_src ? buf1 : buf0;
}

template<class T>
inline void DoubleBuffer<T>::flip()
{
    buf0_is_src = !buf0_is_src;
}

#endif

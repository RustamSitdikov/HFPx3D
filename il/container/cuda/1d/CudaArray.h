//==============================================================================
//
//                                  InsideLoop
//
// This file is distributed under the University of Illinois Open Source
// License. See LICENSE.txt for details.
//
//==============================================================================

#ifndef IL_CUDAARRAY_H
#define IL_CUDAARRAY_H

#include <cuda.h>
#include <cuda_runtime_api.h>

// <initializer_list> is needed for std::initializer_list<T>
#include <initializer_list>

#include <il/base.h>

namespace il {

template <typename T>
class CudaArray {
 private:
  T* data_;
  il::int_t size_;
  il::int_t capacity_;

 public:
  /* \brief Construct an array of n elements
   */
  explicit CudaArray(il::int_t n);
  CudaArray(const il::CudaArray<T>& x);
  CudaArray(CudaArray<T>&& x);

  /* \brief The destructor
   */
  ~CudaArray();

  /* \brief Get the size of the il::CudaArray<T>
   */
  il::int_t size() const;

  /* \brief Get a pointer to the first element of the array
   */
  const T* data() const;

  /* \brief Get a pointer to the first element of the array
   */
  T* data();
};

template <typename T>
CudaArray<T>::CudaArray(il::int_t n) {
  IL_EXPECT_FAST(n >= 0);

  cudaMalloc(&data_, n * sizeof(T));
  size_ = n;
  capacity_ = n;
}

template <typename T>
CudaArray<T>::CudaArray(const il::CudaArray<T>& x) {
  const il::int_t n = x.size();

  cudaError_t error;
  error = cudaMalloc(&data_, n * sizeof(T));
  IL_EXPECT_FAST(error == 0);
  error = cudaMemcpy(data_, x.data_, n * sizeof(T), cudaMemcpyDeviceToDevice);
  IL_EXPECT_FAST(error == 0);

  size_ = n;
  capacity_ = n;
}

template <typename T>
CudaArray<T>::CudaArray(il::CudaArray<T>&& x) {
  data_ = x.data_;
  size_ = x.size_;
  capacity_ = x.capacity_;
  x.data_ = nullptr;
  x.size_ = 0;
  x.capacity_ = 0;
}

template <typename T>
CudaArray<T>::~CudaArray() {
  if (data_ != nullptr) {
    cudaFree(data_);
  }
}

template <typename T>
il::int_t CudaArray<T>::size() const {
  return size_;
}

template <typename T>
const T* CudaArray<T>::data() const {
  return data_;
}

template <typename T>
T* CudaArray<T>::data() {
  return data_;
}
}  // namespace il
#endif  // IL_CUDAARRAY_H

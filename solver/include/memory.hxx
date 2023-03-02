/**
 * \file    memory.hxx
 * \author  akirby
 */

#ifndef MEMORY_HXX
#define MEMORY_HXX

/* system header files */
#include <memory>

namespace dg {

template<typename T>
class memory {
  template <typename U> friend class memory;

  private:
    using size_t = std::size_t;
    using ptrdiff_t = std::ptrdiff_t;

    std::shared_ptr<T[]> shrdPtr;
    size_t lngth;
    size_t offset;

  public:
    memory() :
        shrdPtr{nullptr},
        lngth{0},
        offset{0} {}

    memory(const size_t lngth_) :
        shrdPtr(new T[lngth_]),
        lngth{lngth_},
        offset{0} {}

    memory(const size_t lngth_,
           const T val) :
        shrdPtr(new T[lngth_]),
        lngth{lngth_},
        offset{0} {
        for (size_t i=0;i<lngth;++i) {
          shrdPtr[i] = val;
        }
    }

    memory(const memory<T> &m)=default;
    memory& operator = (const memory<T> &m)=default;
   ~memory()=default;

    void malloc(const size_t lngth_) {
        *this = memory<T>(lngth_);
    }

    void malloc(const size_t lngth_, const T val) {
        *this = memory<T>(lngth_, val);
    }

    void calloc(const size_t lngth_) {
        *this = memory<T>(lngth_, T{0});
    }

    void realloc(const size_t lngth_) {
        memory<T> m(lngth_);
        const ptrdiff_t cnt = std::min(lngth, lngth_);
        m.copyFrom(*this, cnt);
        *this = m;
    }

    memory& swap(memory<T> &m) {
        std::swap(shrdPtr, m.shrdPtr);
        std::swap(lngth, m.lngth);
        std::swap(offset, m.offset);
        return *this;
    }

    T* ptr() {
        if (shrdPtr != nullptr) {
            return shrdPtr.get()+offset;
        } else {
            return nullptr;
        }
    }
    const T* ptr() const {
        if (shrdPtr != nullptr) {
            return shrdPtr.get()+offset;
        } else {
            return nullptr;
        }
    }

    T* begin() {return ptr();}
    T* end() {return ptr() + length();}

    size_t length() const {
        return lngth;
    }

    size_t size() const {
        return lngth*sizeof(T);
    }

    size_t use_count() const {
        return shrdPtr.use_count();
    }

    T& operator[](const ptrdiff_t idx) const {
        return shrdPtr[idx+offset];
    }

    bool operator == (const memory<T> &other) const {
        return (shrdPtr==other.shrdPtr && offset==other.offset);
    }

    bool operator != (const memory<T> &other) const {
        return (shrdPtr!=other.shrdPtr || offset!=other.offset);
    }

    memory<T> operator + (const ptrdiff_t offset_) const {
        return slice(offset_);
    }

    memory<T>& operator += (const ptrdiff_t offset_) {
        *this = slice(offset_);
        return *this;
    }

    memory<T> slice(const ptrdiff_t offset_,
                    const ptrdiff_t count = -1) const {
        memory<T> m(*this);
        m.offset = offset + offset_;
        m.lngth = (count==-1)
                ? (lngth - offset_)
                : count;
        return m;
    }

    /* copy from raw pointer */
    void copyFrom(const T* src,
                  const ptrdiff_t count = -1,
                  const ptrdiff_t offset_ = 0) {

        const ptrdiff_t cnt = (count==-1) ? lngth : count;

        std::copy(src,
                  src+cnt,
                  ptr()+offset_);
    }

    /* copy from memory */
    void copyFrom(const memory<T> src,
                  const ptrdiff_t count = -1,
                  const ptrdiff_t offset_ = 0) {
        const ptrdiff_t cnt = (count==-1) ? lngth : count;

        std::copy(src.ptr(),
                  src.ptr()+cnt,
                  ptr()+offset_);
    }

    /* copy to raw pointer */
    void copyTo(T *dest,
                const ptrdiff_t count = -1,
                const ptrdiff_t offset_ = 0) const {
        const ptrdiff_t cnt = (count==-1) ? lngth : count;

        std::copy(ptr()+offset_,
                  ptr()+offset_+cnt,
                  dest);
    }

    /* copy to memory */
    void copyTo(memory<T> dest,
                const ptrdiff_t count = -1,
                const ptrdiff_t offset_ = 0) const {
        const ptrdiff_t cnt = (count==-1) ? lngth : count;

        std::copy(ptr()+offset_,
                  ptr()+offset_+cnt,
                  dest.ptr());
    }

    memory<T> clone() const {
        memory<T> m(lngth);
        m.copyFrom(*this);
        return m;
    }

    void free() {
        shrdPtr = nullptr;
        lngth=0;
        offset=0;
    }
};

/* extern declare common instantiations for faster compilation*/
template class memory<int>;
template class memory<float>;
template class memory<double>;

} //namespace dg

#endif /* MEMORY_HXX */

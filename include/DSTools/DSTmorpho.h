//
//  DSTmorpho.h
//
//  File created by GILLARD William
//  Centre de Physic des Particules de Marseille
//  Licensed under CC BY-NC 4.0
//  You may share and adapt this code with attribution, 
//  but not for commercial purposes.
//  Licence text: https://creativecommons.org/licenses/by-nc/4.0/

#ifndef __DST_morpho__
#define __DST_morpho__

#include <stdexcept>
#include <type_traits>
#include <limits>
#include <cmath>
#include <thread>                 // added for std::thread
#include <future>                 // added for std::async / std::future

#include <DSTfits/FITSdata.h>
#include <DSTfits/FITSexception.h>
#include <DSTfits/FITSimg.h>
#include <DSTfits/FITShdu.h>

#include <DSTfits/DSF_version.h>

#include "DST_version.h"

using namespace DSL;

namespace DST
{

    namespace Morphology
    {
#pragma region - kernel class
        enum ktype
        {
            MIN,
            MAX,
            MEAN,
            MEDIAN,
            SUM,
            RMS,
            RMSE,
            DILATE,
            ERODE,
            OPEN,
            CLOSE,
            WHITE_TOPHAT,
            BLACK_TOPHAT,
            CORRELATE,
            nCORRELATE
        };

        template <typename T>
        class kernel
        {

#pragma region • members
            private:
                std::valarray<T> fkernel;
                std::vector<size_t> fsize;
                ktype ftype;

#pragma endregion

            /**
             * @brief Validity check of the kernel dimensions
             */
            void validityCheck() const
            {
                size_t total_size = 1;
                for(size_t i=0; i < fsize.size(); i++)
                    total_size *= fsize[i];

                if(fkernel.size() != total_size)
                    throw std::invalid_argument("Kernel size missmatche.["+std::to_string(__LINE__)+"]");
            }

            public:
#pragma region • constructors / destructor

                /**
                 * @brief Construct a new kernel object     
                 * 
                 * @param k kernel array
                 * @param s kernel size vector
                 */
                kernel(const std::valarray<T>& k, const std::vector<size_t>& s, ktype kt = ktype::SUM):fkernel(k), fsize(s),ftype(kt)
                {
                    if(is_allowed_storage_type_v<T> == false)
                        throw std::invalid_argument("Kernel type not supported.["+std::to_string(__LINE__)+"]");

                    validityCheck();
                }

                /**
                 * @brief Construct a new kernel object
                 * 
                 * @param k kernel array
                 * @param s kernel size vector
                 */
                kernel(const std::valarray<T>& k, const std::initializer_list<size_t>& s, ktype kt = ktype::SUM):kernel(k, std::vector<size_t>(s),kt)
                {
                    if(is_allowed_storage_type_v<T> == false)
                        throw std::invalid_argument("Kernel type not supported.["+std::to_string(__LINE__)+"]");

                    validityCheck();
                }


                /**
                 * @brief Copy constructor
                 * 
                 * @param k Kernel to copy to \c this
                 */
                kernel(const kernel<T>& k):fkernel(k.fkernel), fsize(k.fsize),ftype(k.ftype)
                {
                    validityCheck();
                }

                /**
                 * @brief destructor
                 */
                ~kernel(){fsize.clear();}
#pragma endregion

#pragma region • accessors
                const ktype& type() const { return ftype;}
#pragma endregion

#pragma region • modifiers
                void setType(const ktype& t) { ftype = t; }

                inline void setWeightAt(const size_t& i, const T& v)
                {
                    if(i >= fkernel.size())
                        throw std::out_of_range("Index out of range.[" +std::to_string(__LINE__)+"]");

                    fkernel[i] = v;
                }
#pragma endregion


#pragma region • methods
                /**
                 * @brief Get the kernel data array
                 */
                inline const std::valarray<T>& data() const {return fkernel;};

                /**
                 * @brief Get kernel dimension along specified axis
                 * 
                 * @param idx axis index
                 * @return integer giving the diemsnion of the kernel along \c idx axis. 
                 */
                inline const size_t& size(const size_t& idx) const { if(idx >= fsize.size()) throw std::out_of_range("Index out of range.[" +std::to_string(__LINE__)+"]"); return fsize[idx]; };

                /**
                 * @brief Get the kernel size vector
                 */
                inline const std::vector<size_t>& sizes() const {return fsize;};

                /**
                 * @brief Get the total number of elements in the kernel
                 */
                inline size_t size() const { return fkernel.size(); }
#pragma endregion

#pragma region • type utilities (promotion)
            private:
                // Determine if numeric range of A includes B (arithmetic types only)
                template <typename A, typename B, bool Aok = std::is_arithmetic_v<A>, bool Bok = std::is_arithmetic_v<B>>
                struct range_superset_impl : std::false_type {};

                template <typename A, typename B>
                struct range_superset_impl<A,B,true,true>
                {
                    static constexpr long double a_min = static_cast<long double>(std::numeric_limits<A>::lowest());
                    static constexpr long double a_max = static_cast<long double>(std::numeric_limits<A>::max());
                    static constexpr long double b_min = static_cast<long double>(std::numeric_limits<B>::lowest());
                    static constexpr long double b_max = static_cast<long double>(std::numeric_limits<B>::max());
                    static constexpr bool value = (a_min <= b_min) && (a_max >= b_max);
                };

            public:
                /**
                 * @brief Promoted type when combining this kernel<T> with an Other pixel type.
                 *        If T includes Other's numeric range, yields T; otherwise yields Other.
                 */
                template <typename Other>
                using promote_with = std::conditional_t<range_superset_impl<T,Other>::value, T, Other>;
#pragma endregion

#pragma region • operators

                /**
                 * @brief assignement operator
                 */
                kernel<T>& operator=(const kernel<T>& k)
                {
                    kernel<T> temp(k.fkernel,k.fsize);
                    std::swap(fkernel,temp.fkernel);
                    std::swap(fsize,temp.fsize);

                    validityCheck();
                    return *this;
                }

                /**
                 * @brief Addition by an other kernel
                 */
                kernel<T>& operator+=(const kernel<T>& k)
                {
                    if(fsize != k.fsize)
                        throw std::invalid_argument("Kernel size missmatche.["+std::to_string(__LINE__)+"]");

                    fkernel += k.fkernel;
                    return *this;
                }

                /**
                 * @brief Subtraction by an other kernel
                 */
                kernel<T>& operator-=(const kernel<T>& k)
                {
                    if(fsize != k.fsize)
                        throw std::invalid_argument("Kernel size missmatche.["+std::to_string(__LINE__)+"]");

                    fkernel -= k.fkernel;
                    return *this;
                }

                /**
                 * @brief Multiply by an other kernel
                 */
                kernel<T>& operator*=(const kernel<T>& k)
                {
                    if(fsize != k.fsize)
                        throw std::invalid_argument("Kernel size missmatche.["+std::to_string(__LINE__)+"]");

                    fkernel *= k.fkernel;
                    return *this;
                }

                /**
                 * @brief Addition with an other kernel
                 */
                kernel<T>& operator/=(const kernel<T>& k)
                {
                    if(fsize != k.fsize)
                        throw std::invalid_argument("Kernel size missmatche.["+std::to_string(__LINE__)+"]");

                    fkernel /= k.fkernel;
                    return *this;
                }

                /**
                 * @brief Addition by an other array
                 */
                kernel<T>& operator+=(const std::valarray<T>& a)
                {
                    if(is_allowed_storage_type_v<T> == false)
                        throw std::invalid_argument("Kernel type not supported.["+std::to_string(__LINE__)+"]");

                    if(fkernel.size() != a.size())
                        throw std::invalid_argument("Kernel size missmatche.["+std::to_string(__LINE__)+"]");

                    fkernel += a;
                    return *this;
                }

                /**
                 * @brief Subtract by an other array
                 */
                kernel<T>& operator-=(const std::valarray<T>& a)
                {
                    if(is_allowed_storage_type_v<T> == false)
                        throw std::invalid_argument("Kernel type not supported.["+std::to_string(__LINE__)+"]");

                    if(fkernel.size() != a.size())
                        throw std::invalid_argument("Kernel size missmatche.["+std::to_string(__LINE__)+"]");

                    fkernel -= a;
                    return *this;
                }

                /**
                 * @brief Multiply by an other array
                 */
                kernel<T>& operator*=(const std::valarray<T>& a)
                {
                    if(is_allowed_storage_type_v<T> == false)
                        throw std::invalid_argument("Kernel type not supported.["+std::to_string(__LINE__)+"]");

                    if(fkernel.size() != a.size())
                        throw std::invalid_argument("Kernel size missmatche.["+std::to_string(__LINE__)+"]");

                    fkernel *= a;
                    return *this;
                }

                /**
                 * @brief Divide by an other array
                 */
                kernel<T>& operator/=(const std::valarray<T>& a)
                {
                    if(is_allowed_storage_type_v<T> == false)
                        throw std::invalid_argument("Kernel type not supported.["+std::to_string(__LINE__)+"]");

                    if(fkernel.size() != a.size())
                        throw std::invalid_argument("Kernel size missmatche.["+std::to_string(__LINE__)+"]");

                    fkernel /= a;
                    return *this;
                }

                /**
                 * @brief Addition by a scalar
                 */
                kernel<T>& operator+=(const T& a)
                {
                    if(is_allowed_storage_type_v<T> == false)
                        throw std::invalid_argument("Kernel type not supported.["+std::to_string(__LINE__)+"]");

                    fkernel += a;
                    return *this;
                }

                /**
                 * @brief Substract by a scalar
                 */
                kernel<T>& operator-=(const T& a)
                {
                    if(is_allowed_storage_type_v<T> == false)
                        throw std::invalid_argument("Kernel type not supported.["+std::to_string(__LINE__)+"]");

                    fkernel -= a;
                    return *this;
                }

                /**
                 * @brief Multiply by a scalar
                 */
                kernel<T>& operator*=(const T& a)
                {
                    if(is_allowed_storage_type_v<T> == false)
                        throw std::invalid_argument("Kernel type not supported.["+std::to_string(__LINE__)+"]");

                    fkernel *= a;
                    return *this;
                }

                /**
                 * @brief Divide by a scalar
                 */
                kernel<T>& operator/=(const T& a)
                {
                    if(is_allowed_storage_type_v<T> == false)
                        throw std::invalid_argument("Kernel type not supported.["+std::to_string(__LINE__)+"]");

                    fkernel /= a;
                    return *this;
                }
#pragma endregion
#pragma region • predifined kernel

                /**
                 * @brief Create kernel of 1
                 */
                static kernel<T> one(const std::vector<size_t>& s)
                {
                    size_t total_size = 1;
                    for(size_t i=0; i < s.size(); i++)
                        total_size *= s[i];

                    return kernel<T>(std::valarray<T>(T(1),total_size),s);
                };

                static kernel<T> one(const std::initializer_list<size_t>& s)
                {
                    return one(std::vector<size_t>(s));
                };

                /**
                 * @brief Create kernel of 0
                 */
                static kernel<T> zero(const std::vector<size_t>& s)
                {
                    size_t total_size = 1;
                    for(size_t i=0; i < s.size(); i++)
                        total_size *= s[i];

                    return kernel<T>(std::valarray<T>(T(0),total_size),s);
                };

                static kernel<T> zero(const std::initializer_list<size_t>& s)
                {
                    return zero(std::vector<size_t>(s));
                };

#pragma endregion
        };
#pragma endregion

#pragma region - folding function

        /**
         * @brief Fold an image with a given kernel
         * 
         * @param in Input image
         * @param k Kernel element
         * @return Pointer to a \c DSL::FITScube image representation of the folded image
         */
        template <typename T, typename K>
        std::shared_ptr<FITScube> folding(const std::shared_ptr< FITSimg<K> >& in, const kernel<T>& k)
        {
            return folding<T,K>(in,k,1);
        }

        /**
         * @brief Fold an image with a given kernel
         * 
         * @param in Input image
         * @param k Kernel element
         * @param nThreads Number of threads to use for the operation
         * @return Pointer to a \c DSL::FITScube image representation of the folded image
         * @throw FITSexception if input image pointer is null, kernel size is zero, or kernel/image dimension mismatch
         */
        template <typename T, typename K>
        std::shared_ptr<FITScube> folding(const std::shared_ptr< FITSimg<K> >& in, const kernel<T>& k, size_t nThreads)
        {

            if(k.type() == ktype::DILATE)
                return dilation<T>(in,k,nThreads);
            else if(k.type() == ktype::ERODE)
                return erosion<T>(in,k,nThreads);
            else if( k.type() == ktype::OPEN )
                return opening<T>(in,k,nThreads);
            else if( k.type() == ktype::CLOSE )
                return closing<T>(in,k,nThreads);
            else if( k.type() == ktype::WHITE_TOPHAT )
                return wTopHat<T>(in,k,nThreads);
            else if( k.type() == ktype::BLACK_TOPHAT )
                return bTopHat<T>(in,k,nThreads);
            else if( k.type() == ktype::CORRELATE )
                return correlate<T>(in,k,nThreads);
            else if( k.type() == ktype::nCORRELATE )
                return correlate<T>(in,k,nThreads, true);

            if(in == nullptr)
                throw FITSexception(BAD_DIMEN,"DST","folding","Input FITSimg<K> pointer is null.["+std::to_string(__LINE__)+"]");

            if(k.size() == 0)
                throw FITSexception(BAD_DIMEN,"DST","folding","Kernel size is zero.["+std::to_string(__LINE__)+"]");

            if(k.sizes().size() != in->GetDimension() || k.sizes().size() < 1)
                throw FITSexception(BAD_DIMEN,"DST","folding","Kernel missmatch.["+std::to_string(__LINE__)+"]");

            static_assert(is_allowed_storage_type_v<T>, "Kernel type not supported.");
            static_assert(is_allowed_storage_type_v<K>, "Image pixel type not supported.");

            using U = typename kernel<T>::template promote_with<K>;
            static_assert(is_allowed_storage_type_v<U>, "Promoted pixel type not supported.");

            // Construct your output data as std::valarray<U>.
            std::valarray<U> out_data(U(0),in->Nelements());
            const std::valarray<K>* in_data = in->template GetData<K>();
            const std::valarray<bool> msk = in->GetMask();

            // Prepare kernel data in promoted type U
            std::valarray<U> tmpKernel(U(0),k.data().size());

            if(tmpKernel.size() != k.data().size())
                throw std::invalid_argument("Kernel size missmatche.["+std::to_string(__LINE__)+"]");

            for(size_t i = 0; i < k.data().size(); ++i)
                tmpKernel[i] += static_cast<U>( k.data()[i] );

            // Prepare kernel/image dimension helpers used by both single- and multi-threaded code
            const std::vector<size_t> kSizes = k.sizes();
            const size_t kDim = kSizes.size();
            const size_t imgDim = in->GetDimension();

            // signed kernel centers (allow negative offsets when computing target coords)
            std::vector<ptrdiff_t> kCenter(kDim);
            for(size_t a = 0; a < kDim; ++a)
                kCenter[a] = static_cast<ptrdiff_t>(kSizes[a] / 2);
        
            if(nThreads <= 1)
            {
                // Single-threaded implementation
                for(size_t i = 0; i < in_data->size(); i++)
                {
                    // First step is to get slice of the image which correspond to the kernel position. When the kernel is near the border  of the image, only portion of the kernel is used and the slice of the image should be aligned with the kernel    center (i.e. tmp array should be null when outside image bounds)

                    // Get pixel coordinates of the current central pixel (0-based)
                    const std::vector<size_t> pxCoords = in->PixelCoordinates(i);

                    // Temporary accumulator (promoted type)
                    std::valarray<U> tmp(U(0), k.data().size());

                    // Map kernel axes to image axes:
                    // If kernel has fewer dims than image, align kernel to the last kDim axes of the image
                    // image axis index for kernel axis 'ka' is: imgAxis = imgDim - kDim + ka
                    // (this allows 2D kernel applied to last two axes of N-D image)
                    for(size_t kj = 0; kj < k.data().size(); ++kj)
                    {
                        // decode linear kernel index into multi-index (row-major: axis 0 fastest)
                        size_t rem = kj;
                        std::vector<size_t> kIdx(kDim);
                        for(size_t ka = 0; ka < kDim; ++ka)
                        {
                            kIdx[ka] = rem % kSizes[ka];
                            rem /= kSizes[ka];
                        }

                        // compute corresponding image coordinates for this kernel element
                        bool outside = false;
                        std::vector<size_t> tgtCoords = pxCoords; // copy original coords for axes not touched by kernel

                        for(size_t ka = 0; ka < kDim; ++ka)
                        {
                            const ptrdiff_t offset = static_cast<ptrdiff_t>(kIdx[ka]) - kCenter[ka];
                            const size_t imgAxis = imgDim - kDim + ka; // map kernel axis to image axis

                            // safe check (should not happen): if mapping yields invalid axis, mark outside
                            if(imgAxis >= imgDim)
                            {
                                outside = true;
                                break;
                            }

                            // compute target coordinate (may be negative)
                            ptrdiff_t tgt = static_cast<ptrdiff_t>(pxCoords[imgAxis]) + offset;
                            if(tgt < 0)
                            {
                                outside = true;
                                break;
                            }

                            const size_t tgt_sz = in->Size(imgAxis+1); // Size is 1-based
                            if(static_cast<size_t>(tgt) >= tgt_sz)
                            {
                                outside = true;
                                break;
                            }

                            tgtCoords[imgAxis] = static_cast<size_t>(tgt);
                        }

                        if(outside)
                        {
                            tmp[kj] = U(0);
                            continue;
                        }

                        // get linear pixel index corresponding to tgtCoords
                        size_t tgtIndex = in->PixelIndex(tgtCoords);

                        // fetch value from input data and promote to U
                        if(msk[tgtIndex])
                            tmp[kj] = std::numeric_limits<U>::quiet_NaN();
                        else
                            tmp[kj] = static_cast<U>((*in_data)[tgtIndex]);
                    }

                    // perform folding operation
                    std::valarray<U> ttmp = tmp*tmpKernel;
                    switch(k.type())
                    {
                        case ktype::MIN:
                            out_data[i] = ((std::valarray<U>) ttmp[ttmp==ttmp]).min();
                            break;
                        
                        case ktype::MAX:
                            out_data[i] = ((std::valarray<U>) ttmp[ttmp==ttmp]).max();
                            break;
                        
                        case ktype::MEAN:
                            out_data[i] = (((std::valarray<U>) ttmp[ttmp==ttmp]).size() > 0)?((std::valarray<U>) ttmp[ttmp==ttmp]).sum()/static_cast<U>(((std::valarray<U>) ttmp[ttmp==ttmp]).size()):U(0);
                            break;
                        
                        case ktype::MEDIAN:
                            {
                                std::valarray<U> prod = ttmp[ttmp==ttmp];
                                if(prod.size() == 0)
                                    out_data[i] = U(0);
                                else
                                {
                                    std::vector<U> vec(std::begin(prod), std::end(prod));
                                    std::sort(vec.begin(), vec.end());
                                    size_t n = vec.size();
                                    
                                    if(n % 2 == 0)
                                        out_data[i] = (vec[n/2 - 1] + vec[n/2]) / U(2);
                                    else
                                        out_data[i] = vec[n/2];
                                }
                            }
                            break;
                        
                        case ktype::RMS:
                            {
                                std::valarray<U> prod = ttmp[ttmp==ttmp];
                                if(prod.size() == 0)
                                    out_data[i] = U(0);
                                else
                                {
                                    U sumsq = (prod*prod).sum();
                                    out_data[i] = std::sqrt( sumsq / static_cast<U>(prod.size()) );
                                }
                            }
                            break;

                        case ktype::RMSE:
                            {
                                std::valarray<U> prod = ttmp[ttmp==ttmp];
                                if(prod.size() == 0)
                                    out_data[i] = U(0);
                                else
                                {
                                    U meanP = prod.sum()/static_cast<U>(prod.size());
                                    U sumsq = (prod*prod-meanP*meanP).sum();
                                    out_data[i] = std::sqrt( sumsq );
                                }
                            }
                            break;

                        case ktype::SUM:
                        default:
                            out_data[i] = ((std::valarray<U>) ttmp[ttmp==ttmp]).sum();
                            break;
                    }
                }
            }
            else
            {
                // Multi-threaded implementation: partition output indices between threads.
                const size_t N = in_data->size();
                const size_t nT = std::min(nThreads, N); // don't spawn more threads than pixels

                std::vector<std::thread> workers;
                workers.reserve(nT);

                const size_t base = N / nT;
                const size_t rem = N % nT;
            
                // Capture stable copies for threads:
                const std::valarray<K>* in_data_ptr = in_data; // copy pointer by value
                auto in_shared = in;                            // copy shared_ptr (keeps object alive)
            
                for(size_t t = 0; t < nT; ++t)
                {
                    const size_t start = t * base + std::min(t, rem);
                    const size_t len = base + (t < rem ? 1 : 0);
                    const size_t end = start + len;
                
                    workers.emplace_back([start, end, &out_data, in_data_ptr, &k, &tmpKernel, &kSizes, kDim, imgDim, in_shared, &kCenter, & msk]()
                     {
                         for(size_t i = start; i < end; ++i)
                         {
                             const std::vector<size_t> pxCoords = in_shared->PixelCoordinates(i);

                            std::valarray<U> tmp(U(0),k.data().size() );

                            for(size_t kj = 0; kj < k.data().size(); ++kj)
                            {
                                size_t remk = kj;
                                std::vector<size_t> kIdx(kDim);
                                for(size_t ka = 0; ka < kDim; ++ka)
                                {
                                    kIdx[ka] = remk % kSizes[ka];
                                    remk /= kSizes[ka];
                                }

                                bool outside = false;
                                std::vector<size_t> tgtCoords = pxCoords;

                                for(size_t ka = 0; ka < kDim; ++ka)
                                {
                                    const ptrdiff_t offset = static_cast<ptrdiff_t>(kIdx[ka]) - kCenter[ka];
                                    const size_t imgAxis = imgDim - kDim + ka;

                                    if(imgAxis >= imgDim)
                                    {
                                        outside = true;
                                        break;
                                    }

                                    ptrdiff_t tgt = static_cast<ptrdiff_t>(pxCoords[imgAxis]) + offset;
                                    if(tgt < 0)
                                    {
                                        outside = true;
                                        break;
                                    }

                                    const size_t tgt_sz = in_shared->Size(imgAxis+1);
                                    if(static_cast<size_t>(tgt) >= tgt_sz)
                                    {
                                        outside = true;
                                        break;
                                    }

                                    tgtCoords[imgAxis] = static_cast<size_t>(tgt);
                                }

                                if(outside)
                                {
                                    tmp[kj] = U(0);
                                    continue;
                                }

                                size_t tgtIndex = in_shared->PixelIndex(tgtCoords);
                                if(msk[tgtIndex])
                                    tmp[kj] = std::numeric_limits<U>::quiet_NaN();
                                else
                                    tmp[kj] = static_cast<U>((*in_data_ptr)[tgtIndex]);
                            }

                            std::valarray<U> ttmp = tmp*tmpKernel;

                            switch(k.type())
                            {
                                case ktype::MIN:
                                    out_data[i] = ((std::valarray<U>) ttmp[ttmp==ttmp]).min();
                                    break;
                                
                                case ktype::MAX:
                                    out_data[i] = ((std::valarray<U>) ttmp[ttmp==ttmp]).max();
                                    break;
                                
                                case ktype::MEAN:
                                    out_data[i] = (((std::valarray<U>) ttmp[ttmp==ttmp]).size() > 0)?((std::valarray<U>) ttmp[ttmp==ttmp]).sum()/       static_cast<U>(((std::valarray<U>) ttmp[ttmp==ttmp]).size()):U(0);
                                    break;
                                
                                case ktype::MEDIAN:
                                    {
                                        std::valarray<U> prod = ttmp[ttmp==ttmp];
                                        if(prod.size() == 0)
                                            out_data[i] = U(0);
                                        else
                                        {
                                            std::vector<U> vec(std::begin(prod), std::end(prod));
                                            std::sort(vec.begin(), vec.end());
                                            size_t n = vec.size();
                                            
                                            if(n % 2 == 0)
                                                out_data[i] = (vec[n/2 - 1] + vec[n/2]) / U(2);
                                            else
                                                out_data[i] = vec[n/2];
                                        }
                                    }
                                    break;
                                
                                case ktype::RMS:
                                    {
                                        std::valarray<U> prod = ttmp[ttmp==ttmp];
                                        if(prod.size() == 0)
                                            out_data[i] = U(0);
                                        else
                                        {
                                            U sumsq = (prod*prod).sum();
                                            out_data[i] = std::sqrt( sumsq / static_cast<U>(prod.size()) );
                                        }
                                    }
                                    break;
                                
                                case ktype::RMSE:
                                    {
                                        std::valarray<U> prod = ttmp[ttmp==ttmp];
                                        if(prod.size() == 0)
                                            out_data[i] = U(0);
                                        else
                                        {
                                            U meanP = prod.sum()/static_cast<U>(prod.size());
                                            U sumsq = (prod*prod-meanP*meanP).sum();
                                            out_data[i] = std::sqrt( sumsq );
                                        }
                                    }
                                    break;
                                
                                case ktype::SUM:
                                default:
                                    out_data[i] = ((std::valarray<U>) ttmp[ttmp==ttmp]).sum();
                                    break;
                            }
                        }
                    });
                }

                for(auto &th : workers) th.join();
            }

            std::vector<size_t> out_data_shape;
            for(size_t d = 0; d < in->GetDimension(); d++)
                out_data_shape.push_back( in->Size(d+1) );

            FITSimg<U>* outImg = new FITSimg<U>(out_data_shape);
            *(outImg->template GetData<U>()) += out_data;

            for(FITSDictionary::const_iterator it = in->HDU().begin(); it != in->HDU().end(); ++it)
            {
                if(!outImg->HDU().Exists(it->first))
                {
                    outImg->HDU().ValueForKey( it->first, it->second.value(), it->second.type(), it->second.comment() );
                }
            }
            outImg->HDU().ValueForKey("DST_PROC","FOLDING",key_type::fChar,"Processing history keyword");
            outImg->HDU().ValueForKey("DST_VER",DST::gGIT::version_short(),key_type::fChar,"DST framwork version library");
            outImg->HDU().ValueForKey("DST_SHA",DST::gGIT::sha1(),key_type::fChar,"DST framwork GIT identifier");

            return std::shared_ptr<FITScube>(outImg);
        }

        /**
         * @brief Fold an image with a given kernel
         * 
         * @param in Input image
         * @param k Kernel element
         * @return Pointer to a \c DSL::FITScube image representation of the folded image
         */
        template <typename T>
        std::shared_ptr<FITScube> folding(const std::shared_ptr<FITScube>& in, const kernel<T>& k)
        {
            return folding<T>(in,k,1);
        }

        /**
         * @brief Fold an image with a given kernel
         * 
         * @param in Input image
         * @param k Kernel element
         * @param nThreads Number of threads to use for the operation
         * @return Pointer to a \c DSL::FITScube image representation of the folded image
         * @throw FITSexception if input image pointer is null, kernel size is zero, or kernel/image dimension mismatch
         */
        template <typename T>
        std::shared_ptr<FITScube> folding(const std::shared_ptr<FITScube>& in, const kernel<T>& k, size_t nThreads)
        {
            if(in == nullptr)
                throw FITSexception(BAD_DIMEN,"DST","folding","Input FITScube pointer is null.["+std::to_string(__LINE__)+"]");

            if(k.size() == 0)
                throw FITSexception(BAD_DIMEN,"DST","folding","Kernel size is zero.["+std::to_string(__LINE__)+"]");

            switch(in->GetEqBitPerPixel())
            {
                case SBYTE_IMG:
                    return folding<T, int8_t> (std::dynamic_pointer_cast< FITSimg< int8_t> >(in), k, nThreads);
                case BYTE_IMG:
                    return folding<T, uint8_t>(std::dynamic_pointer_cast< FITSimg< uint8_t> >(in), k, nThreads);
                case SHORT_IMG:
                    return folding<T,int16_t> (std::dynamic_pointer_cast< FITSimg<int16_t> >(in), k, nThreads);
                case USHORT_IMG:
                    return folding<T,uint16_t>(std::dynamic_pointer_cast< FITSimg<uint16_t> >(in), k, nThreads);
                case LONG_IMG:
                    return folding<T,int32_t> (std::dynamic_pointer_cast< FITSimg<int32_t> >(in), k, nThreads);
                case ULONG_IMG:
                    return folding<T,uint32_t>(std::dynamic_pointer_cast< FITSimg<uint32_t> >(in), k, nThreads);
                case LONGLONG_IMG:
                    return folding<T,int64_t> (std::dynamic_pointer_cast< FITSimg<int64_t> >(in), k, nThreads);
                case ULONGLONG_IMG:
                    return folding<T,uint64_t>(std::dynamic_pointer_cast< FITSimg<uint64_t> >(in), k, nThreads);
                case FLOAT_IMG:
                    return folding<T,float>   (std::dynamic_pointer_cast< FITSimg<float> >(in), k, nThreads);
                case DOUBLE_IMG:
                    return folding<T,double>  (std::dynamic_pointer_cast< FITSimg<double> >(in), k, nThreads);
                default:
                    break;
            }

            throw FITSexception(BAD_DIMEN,"DST","folding","Unsupported pixel type.["+std::to_string(__LINE__)+"]");
        }
#pragma endregion

#pragma region - correlate function

        /**
         * @brief Correlate an image with a floating kernel (weighted sum).
         *
         * Computes, for each pixel p (kernel centered on p):
         *   out(p) = sum_x in(p+x) * k(x)
         *
         * If normalize=true, divides by the sum of absolute kernel weights actually used
         * (i.e., excluding masked pixels and pixels outside image bounds).
         *
         * @param in Input image
         * @param k Kernel weights and shape
         * @param nThreads Number of threads to use
         * @param normalize If true, normalizes by sum(|w|) over valid samples
         */
        template <typename T, typename K>
        std::shared_ptr<FITScube> correlate(const std::shared_ptr< FITSimg<K> >& in, const kernel<T>& k, size_t nThreads = 1, bool normalize = false)
        {
            if(in == nullptr)
                throw FITSexception(BAD_DIMEN,"DST","correlate","Input FITSimg<K> pointer is null.["+std::to_string(__LINE__)+"]");

            if(k.size() == 0)
                throw FITSexception(BAD_DIMEN,"DST","correlate","Kernel size is zero.["+std::to_string(__LINE__)+"]");

            if(k.sizes().size() != in->GetDimension() || k.sizes().size() < 1)
                throw FITSexception(BAD_DIMEN,"DST","correlate","Kernel missmatch.["+std::to_string(__LINE__)+"]");

            static_assert(is_allowed_storage_type_v<T>, "Kernel type not supported.");
            static_assert(is_allowed_storage_type_v<K>, "Image pixel type not supported.");

            using U = typename kernel<T>::template promote_with<K>;
            static_assert(is_allowed_storage_type_v<U>, "Promoted pixel type not supported.");

            std::valarray<U> out_data(U(0), in->Nelements());
            const std::valarray<K>* in_data = in->template GetData<K>();
            const std::valarray<bool> msk = in->GetMask();

            // Promote kernel weights into U
            std::valarray<U> tmpKernel(U(0), k.data().size());
            if(tmpKernel.size() != k.data().size())
                throw std::invalid_argument("Kernel size missmatche.["+std::to_string(__LINE__)+"]");
            for(size_t i = 0; i < k.data().size(); ++i)
                tmpKernel[i] = static_cast<U>(k.data()[i]);

            const std::vector<size_t> kSizes = k.sizes();
            const size_t kDim = kSizes.size();
            const size_t imgDim = in->GetDimension();

            std::vector<ptrdiff_t> kCenter(kDim);
            for(size_t a = 0; a < kDim; ++a)
                kCenter[a] = static_cast<ptrdiff_t>(kSizes[a] / 2);

            auto absU = [](U v) -> U
            {
                if(v < U(0))
                    return U(-v);
                return v;
            };

            auto processRange = [&](size_t start, size_t end)
            {
                for(size_t i = start; i < end; ++i)
                {
                    const std::vector<size_t> pxCoords = in->PixelCoordinates(i);
                    U acc = U(0);
                    U wsum = U(0);

                    for(size_t kj = 0; kj < k.data().size(); ++kj)
                    {
                        size_t remk = kj;
                        std::vector<size_t> kIdx(kDim);
                        for(size_t ka = 0; ka < kDim; ++ka)
                        {
                            kIdx[ka] = remk % kSizes[ka];
                            remk /= kSizes[ka];
                        }

                        bool outside = false;
                        std::vector<size_t> tgtCoords = pxCoords;

                        for(size_t ka = 0; ka < kDim; ++ka)
                        {
                            const ptrdiff_t offset = static_cast<ptrdiff_t>(kIdx[ka]) - kCenter[ka];
                            const size_t imgAxis = imgDim - kDim + ka;

                            if(imgAxis >= imgDim)
                            {
                                outside = true;
                                break;
                            }

                            ptrdiff_t tgt = static_cast<ptrdiff_t>(pxCoords[imgAxis]) + offset;
                            if(tgt < 0)
                            {
                                outside = true;
                                break;
                            }

                            const size_t tgt_sz = in->Size(imgAxis+1);
                            if(static_cast<size_t>(tgt) >= tgt_sz)
                            {
                                outside = true;
                                break;
                            }

                            tgtCoords[imgAxis] = static_cast<size_t>(tgt);
                        }

                        if(outside)
                            continue;

                        const size_t tgtIndex = in->PixelIndex(tgtCoords);
                        if(msk[tgtIndex])
                            continue;

                        const U w = tmpKernel[kj];
                        const U v = static_cast<U>((*in_data)[tgtIndex]);
                        acc += v * w;

                        if(normalize)
                        {
                            const U aw = absU(w);
                            if(aw == aw) // ignore NaN weights
                                wsum += aw;
                        }
                    }

                    if(normalize)
                    {
                        if(wsum != U(0))
                            acc /= wsum;
                        else
                            acc = U(0);
                    }

                    out_data[i] = acc;
                }
            };

            if(nThreads <= 1)
            {
                processRange(0, in_data->size());
            }
            else
            {
                const size_t N = in_data->size();
                const size_t nT = std::min(nThreads, N);

                std::vector<std::thread> workers;
                workers.reserve(nT);

                const size_t base = N / nT;
                const size_t rem = N % nT;

                for(size_t t = 0; t < nT; ++t)
                {
                    const size_t start = t * base + std::min(t, rem);
                    const size_t len = base + (t < rem ? 1 : 0);
                    const size_t end = start + len;
                    workers.emplace_back([&, start, end]() { processRange(start, end); });
                }
                for(auto& th : workers) th.join();
            }

            std::vector<size_t> out_data_shape;
            for(size_t d = 0; d < in->GetDimension(); d++)
                out_data_shape.push_back(in->Size(d+1));

            FITSimg<U>* outImg = new FITSimg<U>(out_data_shape);
            *(outImg->template GetData<U>()) += out_data;

            for(FITSDictionary::const_iterator it = in->HDU().begin(); it != in->HDU().end(); ++it)
            {
                if(!outImg->HDU().Exists(it->first))
                    outImg->HDU().ValueForKey(it->first, it->second.value(), it->second.type(), it->second.comment());
            }

            outImg->HDU().ValueForKey("DST_PROC","CORRELATE",key_type::fChar,"Processing history keyword");
            outImg->HDU().ValueForKey("DST_NORM", normalize ? "T" : "F", key_type::fChar, "Correlation normalization flag");
            outImg->HDU().ValueForKey("DST_VER",DST::gGIT::version_short(),key_type::fChar,"DST framwork version library");
            outImg->HDU().ValueForKey("DST_SHA",DST::gGIT::sha1(),key_type::fChar,"DST framwork GIT identifier");

            return std::shared_ptr<FITScube>(outImg);
        }

        /**
         * @brief Correlate a dynamically-typed FITScube with a floating kernel.
         */
        template <typename T>
        std::shared_ptr<FITScube> correlate(const std::shared_ptr<FITScube>& in, const kernel<T>& k, size_t nThreads = 1, bool normalize = false)
        {
            if(in == nullptr)
                throw FITSexception(BAD_DIMEN,"DST","correlate","Input FITScube pointer is null.["+std::to_string(__LINE__)+"]");

            if(k.size() == 0)
                throw FITSexception(BAD_DIMEN,"DST","correlate","Kernel size is zero.["+std::to_string(__LINE__)+"]");

            switch(in->GetEqBitPerPixel())
            {
                case SBYTE_IMG:
                    return correlate<T, int8_t> (std::dynamic_pointer_cast< FITSimg< int8_t> >(in), k, nThreads, normalize);
                case BYTE_IMG:
                    return correlate<T, uint8_t>(std::dynamic_pointer_cast< FITSimg< uint8_t> >(in), k, nThreads, normalize);
                case SHORT_IMG:
                    return correlate<T,int16_t> (std::dynamic_pointer_cast< FITSimg<int16_t> >(in), k, nThreads, normalize);
                case USHORT_IMG:
                    return correlate<T,uint16_t>(std::dynamic_pointer_cast< FITSimg<uint16_t> >(in), k, nThreads, normalize);
                case LONG_IMG:
                    return correlate<T,int32_t> (std::dynamic_pointer_cast< FITSimg<int32_t> >(in), k, nThreads, normalize);
                case ULONG_IMG:
                    return correlate<T,uint32_t>(std::dynamic_pointer_cast< FITSimg<uint32_t> >(in), k, nThreads, normalize);
                case LONGLONG_IMG:
                    return correlate<T,int64_t> (std::dynamic_pointer_cast< FITSimg<int64_t> >(in), k, nThreads, normalize);
                case ULONGLONG_IMG:
                    return correlate<T,uint64_t>(std::dynamic_pointer_cast< FITSimg<uint64_t> >(in), k, nThreads, normalize);
                case FLOAT_IMG:
                    return correlate<T,float>   (std::dynamic_pointer_cast< FITSimg<float> >(in), k, nThreads, normalize);
                case DOUBLE_IMG:
                    return correlate<T,double>  (std::dynamic_pointer_cast< FITSimg<double> >(in), k, nThreads, normalize);
                default:
                    break;
            }

            throw FITSexception(BAD_DIMEN,"DST","correlate","Unsupported pixel type.["+std::to_string(__LINE__)+"]");
        }

#pragma endregion

#pragma region - dilation function

        /**
         * @brief Dilate an image with a given structure element kernel
         * 
         * @param in Input image
         * @param k Structure Kernel element
         * @return Pointer to a \c DSL::FITScube image representation of the folded image
         * @throw FITSexception if input image pointer is null, kernel size is zero, or kernel/image dimension mismatch
         */
        template <typename T, typename K>
        std::shared_ptr<FITScube> dilation(const std::shared_ptr< FITSimg<K> >& in, const kernel<T>& k)
        {
            return dilation<T,K>(in,k,1);
        }

        /**
         * @brief Dilate an image with a given structure element kernel
         * 
         * @param in Input image
         * @param k Structure Kernel element
         * @param nThreads Number of threads to use for the operation
         * @return Pointer to a \c DSL::FITScube image representation of the folded image
         * @throw FITSexception if input image pointer is null, kernel size is zero, or kernel/image dimension mismatch
         */
        template <typename T, typename K>
        std::shared_ptr<FITScube> dilation(const std::shared_ptr< FITSimg<K> >& in, const kernel<T>& k, size_t nThreads)
        {
            if(in == nullptr)
                throw FITSexception(BAD_DIMEN,"DST","dilation","Input FITSimg<K> pointer is null.["+std::to_string(__LINE__)+"]");

            if(k.size() == 0)
                throw FITSexception(BAD_DIMEN,"DST","dilation","Kernel size is zero.["+std::to_string(__LINE__)+"]");

            if(k.sizes().size() != in->GetDimension() || k.sizes().size() < 1)
                throw FITSexception(BAD_DIMEN,"DST","dilation","Kernel missmatch.["+std::to_string(__LINE__)+"]");

            static_assert(is_allowed_storage_type_v<T>, "Kernel type not supported.");
            static_assert(is_allowed_storage_type_v<K>, "Image pixel type not supported.");

            using U = typename kernel<T>::template promote_with<K>;
            static_assert(is_allowed_storage_type_v<U>, "Promoted pixel type not supported.");

            // Construct your output data as std::valarray<U>.
            std::valarray<U> out_data(U(0),in->Nelements());
            const std::valarray<K>* in_data = in->template GetData<K>();
            const std::valarray<bool> msk = in->GetMask();

            // Prepare kernel data in promoted type U
            std::valarray<U> tmpKernel(U(0),k.data().size());

            if(tmpKernel.size() != k.data().size())
                throw std::invalid_argument("Kernel size missmatche.["+std::to_string(__LINE__)+"]");

            for(size_t i = 0; i < k.data().size(); ++i)
                tmpKernel[i] += static_cast<U>( k.data()[i] );

            // Prepare kernel/image dimension helpers used by both single- and multi-threaded code
            const std::vector<size_t> kSizes = k.sizes();
            const size_t kDim = kSizes.size();
            const size_t imgDim = in->GetDimension();

            // signed kernel centers (allow negative offsets when computing target coords)
            std::vector<ptrdiff_t> kCenter(kDim);
            for(size_t a = 0; a < kDim; ++a)
                kCenter[a] = static_cast<ptrdiff_t>(kSizes[a] / 2);
        
            if(nThreads <= 1)
            {
                // Single-threaded implementation
                for(size_t i = 0; i < in_data->size(); i++)
                {
                    // First step is to get slice of the image which correspond to the kernel position. When the kernel is near the border  of the image, only portion of the kernel is used and the slice of the image should be aligned with the kernel    center (i.e. tmp array should be null when outside image bounds)

                    // Get pixel coordinates of the current central pixel (0-based)
                    const std::vector<size_t> pxCoords = in->PixelCoordinates(i);

                    // Temporary accumulator (promoted type)
                    std::valarray<U> tmp(U(0), k.data().size());

                    // Map kernel axes to image axes:
                    // If kernel has fewer dims than image, align kernel to the last kDim axes of the image
                    // image axis index for kernel axis 'ka' is: imgAxis = imgDim - kDim + ka
                    // (this allows 2D kernel applied to last two axes of N-D image)
                    for(size_t kj = 0; kj < k.data().size(); ++kj)
                    {
                        // decode linear kernel index into multi-index (row-major: axis 0 fastest)
                        size_t rem = kj;
                        std::vector<size_t> kIdx(kDim);
                        for(size_t ka = 0; ka < kDim; ++ka)
                        {
                            kIdx[ka] = rem % kSizes[ka];
                            rem /= kSizes[ka];
                        }

                        // compute corresponding image coordinates for this kernel element
                        bool outside = false;
                        std::vector<size_t> tgtCoords = pxCoords; // copy original coords for axes not touched by kernel

                        for(size_t ka = 0; ka < kDim; ++ka)
                        {
                            const ptrdiff_t offset = static_cast<ptrdiff_t>(kIdx[ka]) - kCenter[ka];
                            const size_t imgAxis = imgDim - kDim + ka; // map kernel axis to image axis

                            // safe check (should not happen): if mapping yields invalid axis, mark outside
                            if(imgAxis >= imgDim)
                            {
                                outside = true;
                                break;
                            }

                            // compute target coordinate (may be negative)
                            ptrdiff_t tgt = static_cast<ptrdiff_t>(pxCoords[imgAxis]) + offset;
                            if(tgt < 0)
                            {
                                outside = true;
                                break;
                            }

                            const size_t tgt_sz = in->Size(imgAxis+1); // Size is 1-based
                            if(static_cast<size_t>(tgt) >= tgt_sz)
                            {
                                outside = true;
                                break;
                            }

                            tgtCoords[imgAxis] = static_cast<size_t>(tgt);
                        }

                        if(outside)
                        {
                            tmp[kj] = std::numeric_limits<U>::lowest();
                            continue;
                        }

                        // get linear pixel index corresponding to tgtCoords
                        size_t tgtIndex = in->PixelIndex(tgtCoords);

                        // fetch value from input data and promote to U
                        if(msk[tgtIndex])
                            tmp[kj] = std::numeric_limits<U>::lowest();
                        else
                            tmp[kj] = static_cast<U>((*in_data)[tgtIndex]);
                    }

                    // flat dilation: max( tmp [ kernel ]) over footprint, ignoring NaNs
                    const std::valarray<U> vals = tmp[ tmpKernel > U(0) ];
                    out_data[i] = (vals.size() > 0) ? vals.max() : U(0);
                }
            }
            else
            {
                // Multi-threaded implementation: partition output indices between threads.
                const size_t N = in_data->size();
                const size_t nT = std::min(nThreads, N); // don't spawn more threads than pixels

                std::vector<std::thread> workers;
                workers.reserve(nT);

                const size_t base = N / nT;
                const size_t rem = N % nT;
            
                // Capture stable copies for threads:
                const std::valarray<K>* in_data_ptr = in_data; // copy pointer by value
                auto in_shared = in;                            // copy shared_ptr (keeps object alive)
            
                for(size_t t = 0; t < nT; ++t)
                {
                    const size_t start = t * base + std::min(t, rem);
                    const size_t len = base + (t < rem ? 1 : 0);
                    const size_t end = start + len;
                
                    workers.emplace_back([start, end, &out_data, in_data_ptr, &k, &tmpKernel, &kSizes, kDim, imgDim, in_shared, &kCenter,&  msk]()
                     {
                         for(size_t i = start; i < end; ++i)
                         {
                             const std::vector<size_t> pxCoords = in_shared->PixelCoordinates(i);

                            std::valarray<U> tmp(U(0),k.data().size() );

                            for(size_t kj = 0; kj < k.data().size(); ++kj)
                            {
                                size_t remk = kj;
                                std::vector<size_t> kIdx(kDim);
                                for(size_t ka = 0; ka < kDim; ++ka)
                                {
                                    kIdx[ka] = remk % kSizes[ka];
                                    remk /= kSizes[ka];
                                }

                                bool outside = false;
                                std::vector<size_t> tgtCoords = pxCoords;

                                for(size_t ka = 0; ka < kDim; ++ka)
                                {
                                    const ptrdiff_t offset = static_cast<ptrdiff_t>(kIdx[ka]) - kCenter[ka];
                                    const size_t imgAxis = imgDim - kDim + ka;

                                    if(imgAxis >= imgDim)
                                    {
                                        outside = true;
                                        break;
                                    }

                                    ptrdiff_t tgt = static_cast<ptrdiff_t>(pxCoords[imgAxis]) + offset;
                                    if(tgt < 0)
                                    {
                                        outside = true;
                                        break;
                                    }

                                    const size_t tgt_sz = in_shared->Size(imgAxis+1);
                                    if(static_cast<size_t>(tgt) >= tgt_sz)
                                    {
                                        outside = true;
                                        break;
                                    }

                                    tgtCoords[imgAxis] = static_cast<size_t>(tgt);
                                }

                                if(outside)
                                {
                                    tmp[kj] = std::numeric_limits<U>::lowest();
                                    continue;
                                }

                                size_t tgtIndex = in_shared->PixelIndex(tgtCoords);

                                if(msk[tgtIndex])
                                    tmp[kj] = std::numeric_limits<U>::lowest();
                                else
                                    tmp[kj] = static_cast<U>((*in_data_ptr)[tgtIndex]);
                            }

                            {
                                // non-flat dilation: max( tmp  ) over footprint, ignoring NaNs
                                const std::valarray<U> vals = tmp[ tmpKernel > U(0) ];
                                out_data[i] = (vals.size() > 0) ? vals.max() : U(0);
                            }
                        }
                    });
                }

                for(auto &th : workers) th.join();
            }

            std::vector<size_t> out_data_shape;
            for(size_t d = 0; d < in->GetDimension(); d++)
                out_data_shape.push_back( in->Size(d+1) );

            FITSimg<U>* outImg = new FITSimg<U>(out_data_shape);
            *(outImg->template GetData<U>()) += out_data;

            for(FITSDictionary::const_iterator it = in->HDU().begin(); it != in->HDU().end(); ++it)
            {
                if(!outImg->HDU().Exists(it->first))
                {
                    outImg->HDU().ValueForKey( it->first, it->second.value(), it->second.type(), it->second.comment() );
                }
            }
            outImg->HDU().ValueForKey("DST_PROC","DILATION",key_type::fChar,"Processing history keyword");
            outImg->HDU().ValueForKey("DST_VER",DST::gGIT::version_short(),key_type::fChar,"DST framwork version library");
            outImg->HDU().ValueForKey("DST_SHA",DST::gGIT::sha1(),key_type::fChar,"DST framwork GIT identifier");

            return std::shared_ptr<FITScube>(outImg);
        }


        /**
         * @brief Dilate an image with a given kernel
         * 
         * @param in Input image
         * @param k Kernel element
         * @return Pointer to a \c DSL::FITScube image representation of the dilated image
         */
        template <typename T>
        std::shared_ptr<FITScube> dilation(const std::shared_ptr<FITScube>& in, const kernel<T>& k)
        {
            return dilation<T>(in,k,1);
        }

        /**
         * @brief Dilate an image with a given kernel
         * 
         * @param in Input image
         * @param k Kernel element
         * @param nThreads Number of threads to use for the operation
         * @return Pointer to a \c DSL::FITScube image representation of the dilated image
         * @throw FITSexception if input image pointer is null, kernel size is zero, or kernel/image dimension mismatch
         */
        template <typename T>
        std::shared_ptr<FITScube> dilation(const std::shared_ptr<FITScube>& in, const kernel<T>& k, size_t nThreads)
        {
            if(in == nullptr)
                throw FITSexception(BAD_DIMEN,"DST","dilation","Input FITScube pointer is null.["+std::to_string(__LINE__)+"]");

            if(k.size() == 0)
                throw FITSexception(BAD_DIMEN,"DST","dilation","Kernel size is zero.["+std::to_string(__LINE__)+"]");

            switch(in->GetEqBitPerPixel())
            {
                case SBYTE_IMG:
                    return dilation<T, int8_t> (std::dynamic_pointer_cast< FITSimg< int8_t> >(in), k, nThreads);
                case BYTE_IMG:
                    return dilation<T, uint8_t>(std::dynamic_pointer_cast< FITSimg< uint8_t> >(in), k, nThreads);
                case SHORT_IMG:
                    return dilation<T,int16_t> (std::dynamic_pointer_cast< FITSimg<int16_t> >(in), k, nThreads);
                case USHORT_IMG:
                    return dilation<T,uint16_t>(std::dynamic_pointer_cast< FITSimg<uint16_t> >(in), k, nThreads);
                case LONG_IMG:
                    return dilation<T,int32_t> (std::dynamic_pointer_cast< FITSimg<int32_t> >(in), k, nThreads);
                case ULONG_IMG:
                    return dilation<T,uint32_t>(std::dynamic_pointer_cast< FITSimg<uint32_t> >(in), k, nThreads);
                case LONGLONG_IMG:
                    return dilation<T,int64_t> (std::dynamic_pointer_cast< FITSimg<int64_t> >(in), k, nThreads);
                case ULONGLONG_IMG:
                    return dilation<T,uint64_t>(std::dynamic_pointer_cast< FITSimg<uint64_t> >(in), k, nThreads);
                case FLOAT_IMG:
                    return dilation<T,float>   (std::dynamic_pointer_cast< FITSimg<float> >(in), k, nThreads);
                case DOUBLE_IMG:
                    return dilation<T,double>  (std::dynamic_pointer_cast< FITSimg<double> >(in), k, nThreads);
                default:
                    break;
            }

            throw FITSexception(BAD_DIMEN,"DST","dilation","Unsupported pixel type.["+std::to_string(__LINE__)+"]");
        }
#pragma endregion
#pragma region - erosion function

        /**
         * @brief Erode an image with a given structure element kernel
         * 
         * @param in Input image
         * @param k Structure Kernel element
         * @return Pointer to a \c DSL::FITScube image representation of the eroded image
         * @throw FITSexception if input image pointer is null, kernel size is zero, or kernel/image dimension mismatch
         */
        template <typename T, typename K>
        std::shared_ptr<FITScube> erosion(const std::shared_ptr< FITSimg<K> >& in, const kernel<T>& k)
        {
            return erosion<T,K>(in,k,1);
        }

        /**
         * @brief Erode an image with a given structure element kernel
         * 
         * @param in Input image
         * @param k Structure Kernel element
         * @param nThreads Number of threads to use for the operation
         * @return Pointer to a \c DSL::FITScube image representation of the eroded image
         * @throw FITSexception if input image pointer is null, kernel size is zero, or kernel/image dimension mismatch
         */
        template <typename T, typename K>
        std::shared_ptr<FITScube> erosion(const std::shared_ptr< FITSimg<K> >& in, const kernel<T>& k, size_t nThreads)
        {
            if(in == nullptr)
                throw FITSexception(BAD_DIMEN,"DST","erosion","Input FITSimg<K> pointer is null.["+std::to_string(__LINE__)+"]");

            if(k.size() == 0)
                throw FITSexception(BAD_DIMEN,"DST","erosion","Kernel size is zero.["+std::to_string(__LINE__)+"]");

            if(k.sizes().size() != in->GetDimension() || k.sizes().size() < 1)
                throw FITSexception(BAD_DIMEN,"DST","erosion","Kernel missmatch.["+std::to_string(__LINE__)+"]");

            static_assert(is_allowed_storage_type_v<T>, "Kernel type not supported.");
            static_assert(is_allowed_storage_type_v<K>, "Image pixel type not supported.");

            using U = typename kernel<T>::template promote_with<K>;
            static_assert(is_allowed_storage_type_v<U>, "Promoted pixel type not supported.");

            // Construct your output data as std::valarray<U>.
            std::valarray<U> out_data(U(0),in->Nelements());
            const std::valarray<K>* in_data = in->template GetData<K>();
            const std::valarray<bool> msk = in->GetMask();

            // Prepare kernel data in promoted type U
            std::valarray<U> tmpKernel(U(0),k.data().size());

            if(tmpKernel.size() != k.data().size())
                throw std::invalid_argument("Kernel size missmatche.["+std::to_string(__LINE__)+"]");

            for(size_t i = 0; i < k.data().size(); ++i)
                tmpKernel[i] += static_cast<U>( k.data()[i] );

            // Prepare kernel/image dimension helpers used by both single- and multi-threaded code
            const std::vector<size_t> kSizes = k.sizes();
            const size_t kDim = kSizes.size();
            const size_t imgDim = in->GetDimension();

            // signed kernel centers (allow negative offsets when computing target coords)
            std::vector<ptrdiff_t> kCenter(kDim);
            for(size_t a = 0; a < kDim; ++a)
                kCenter[a] = static_cast<ptrdiff_t>(kSizes[a] / 2);
        
            if(nThreads <= 1)
            {
                // Single-threaded implementation
                for(size_t i = 0; i < in_data->size(); i++)
                {
                    // First step is to get slice of the image which correspond to the kernel position. When the kernel is near the border  of the image, only portion of the kernel is used and the slice of the image should be aligned with the kernel    center (i.e. tmp array should be null when outside image bounds)

                    // Get pixel coordinates of the current central pixel (0-based)
                    const std::vector<size_t> pxCoords = in->PixelCoordinates(i);

                    // Temporary accumulator (promoted type)
                    std::valarray<U> tmp(U(0), k.data().size());

                    // Map kernel axes to image axes:
                    // If kernel has fewer dims than image, align kernel to the last kDim axes of the image
                    // image axis index for kernel axis 'ka' is: imgAxis = imgDim - kDim + ka
                    // (this allows 2D kernel applied to last two axes of N-D image)
                    for(size_t kj = 0; kj < k.data().size(); ++kj)
                    {
                        // decode linear kernel index into multi-index (kernel coordinates)
                        size_t rem = kj;
                        std::vector<size_t> kIdx(kDim);
                        for(ptrdiff_t ka = static_cast<ptrdiff_t>(kDim) - 1; ka >= 0; --ka)
                        {
                            kIdx[ka] = rem % kSizes[ka];
                            rem /= kSizes[ka];
                        }

                        // compute corresponding image coordinates for this kernel element
                        bool outside = false;
                        std::vector<size_t> tgtCoords = pxCoords; // copy original coords for axes not touched by kernel

                        for(size_t ka = 0; ka < kDim; ++ka)
                        {
                            const ptrdiff_t offset = static_cast<ptrdiff_t>(kIdx[ka]) - kCenter[ka];
                            const size_t imgAxis = imgDim - kDim + ka; // map kernel axis to image axis

                            // safe check (should not happen): if mapping yields invalid axis, mark outside
                            if(imgAxis >= imgDim)
                            {
                                outside = true;
                                break;
                            }

                            // compute target coordinate (may be negative)
                            ptrdiff_t tgt = static_cast<ptrdiff_t>(pxCoords[imgAxis]) + offset;
                            if(tgt < 0)
                            {
                                outside = true;
                                break;
                            }

                            const size_t tgt_sz = in->Size(imgAxis+1); // Size is 1-based
                            if(static_cast<size_t>(tgt) >= tgt_sz)
                            {
                                outside = true;
                                break;
                            }

                            tgtCoords[imgAxis] = static_cast<size_t>(tgt);
                        }

                        if(outside)
                        {
                            tmp[kj] = std::numeric_limits<U>::max();
                            continue;
                        }

                        // get linear pixel index corresponding to tgtCoords
                        size_t tgtIndex = in->PixelIndex(tgtCoords);

                        // fetch value from input data and promote to U
                        if(msk[tgtIndex])
                            tmp[kj] = std::numeric_limits<U>::max();
                        else
                            tmp[kj] = static_cast<U>((*in_data)[tgtIndex]);
                    }

                    {
                        // perform folding operation
                        const std::valarray<U> vals = tmp[ tmpKernel > U(0) ];  // footprint  
                        out_data[i] = (vals.size() > 0) ? vals.min() : U(0);
                    }
                }
            }
            else
            {
                // Multi-threaded implementation: partition output indices between threads.
                const size_t N = in_data->size();
                const size_t nT = std::min(nThreads, N); // don't spawn more threads than pixels

                std::vector<std::thread> workers;
                workers.reserve(nT);

                const size_t base = N / nT;
                const size_t rem = N % nT;
            
                // Capture stable copies for threads:
                const std::valarray<K>* in_data_ptr = in_data; // copy pointer by value
                auto in_shared = in;                            // copy shared_ptr (keeps object alive)
            
                for(size_t t = 0; t < nT; ++t)
                {
                    const size_t start = t * base + std::min(t, rem);
                    const size_t len = base + (t < rem ? 1 : 0);
                    const size_t end = start + len;
                
                    workers.emplace_back([start, end, &out_data, in_data_ptr, &k, &tmpKernel, &kSizes, kDim, imgDim, in_shared, &kCenter,&  msk]()
                     {
                         for(size_t i = start; i < end; ++i)
                         {
                             const std::vector<size_t> pxCoords = in_shared->PixelCoordinates(i);

                            std::valarray<U> tmp(U(0),k.data().size() );

                            for(size_t kj = 0; kj < k.data().size(); ++kj)
                            {
                                size_t remk = kj;
                                std::vector<size_t> kIdx(kDim);
                                for(ptrdiff_t ka = static_cast<ptrdiff_t>(kDim) - 1; ka >= 0; --ka)
                                {
                                    kIdx[ka] = remk % kSizes[ka];
                                    remk /= kSizes[ka];
                                }

                                bool outside = false;
                                std::vector<size_t> tgtCoords = pxCoords;

                                for(size_t ka = 0; ka < kDim; ++ka)
                                {
                                    const ptrdiff_t offset = static_cast<ptrdiff_t>(kIdx[ka]) - kCenter[ka];
                                    const size_t imgAxis = imgDim - kDim + ka;

                                    if(imgAxis >= imgDim)
                                    {
                                        outside = true;
                                        break;
                                    }

                                    ptrdiff_t tgt = static_cast<ptrdiff_t>(pxCoords[imgAxis]) + offset;
                                    if(tgt < 0)
                                    {
                                        outside = true;
                                        break;
                                    }

                                    const size_t tgt_sz = in_shared->Size(imgAxis+1);
                                    if(static_cast<size_t>(tgt) >= tgt_sz)
                                    {
                                        outside = true;
                                        break;
                                    }

                                    tgtCoords[imgAxis] = static_cast<size_t>(tgt);
                                }

                                if(outside)
                                {
                                    tmp[kj] = std::numeric_limits<U>::max();
                                    continue;
                                }

                                size_t tgtIndex = in_shared->PixelIndex(tgtCoords);

                                if(msk[tgtIndex])
                                    tmp[kj] = std::numeric_limits<U>::max();
                                else
                                    tmp[kj] = static_cast<U>((*in_data_ptr)[tgtIndex]);
                            }

                            {
                                const std::valarray<U> vals = tmp[ tmpKernel > U(0) ];  // footprint
                                out_data[i] = (vals.size() > 0) ? vals.min() : U(0);
                            }
                        }
                    });
                }

                for(auto &th : workers) th.join();
            }

            std::vector<size_t> out_data_shape;
            for(size_t d = 0; d < in->GetDimension(); d++)
                out_data_shape.push_back( in->Size(d+1) );

            FITSimg<U>* outImg = new FITSimg<U>(out_data_shape);
            *(outImg->template GetData<U>()) += out_data;

            for(FITSDictionary::const_iterator it = in->HDU().begin(); it != in->HDU().end(); ++it)
            {
                if(!outImg->HDU().Exists(it->first))
                {
                    outImg->HDU().ValueForKey( it->first, it->second.value(), it->second.type(), it->second.comment() );
                }
            }
            outImg->HDU().ValueForKey("DST_PROC","EROSION",key_type::fChar,"Processing history keyword");
            outImg->HDU().ValueForKey("DST_VER",DST::gGIT::version_short(),key_type::fChar,"DST framwork version library");
            outImg->HDU().ValueForKey("DST_SHA",DST::gGIT::sha1(),key_type::fChar,"DST framwork GIT identifier");

            return std::shared_ptr<FITScube>(outImg);
        }

        /**
         * @brief Erode an image with a given kernel
         * 
         * @param in Input image
         * @param k Kernel element
         * @return Pointer to a \c DSL::FITScube image representation of the dilated image
         */
        template <typename T>
        std::shared_ptr<FITScube> erosion(const std::shared_ptr<FITScube>& in, const kernel<T>& k)
        {
            return erosion<T>(in,k,1);
        }

        /**
         * @brief Erode an image with a given kernel
         * 
         * @param in Input image
         * @param k Kernel element
         * @param nThreads Number of threads to use for the operation
         * @return Pointer to a \c DSL::FITScube image representation of the dilated image
         * @throw FITSexception if input image pointer is null, kernel size is zero, or kernel/image dimension mismatch
         */
        template <typename T>
        std::shared_ptr<FITScube> erosion(const std::shared_ptr<FITScube>& in, const kernel<T>& k, size_t nThreads)
        {
            if(in == nullptr)
                throw FITSexception(BAD_DIMEN,"DST","erosion","Input FITScube pointer is null.["+std::to_string(__LINE__)+"]");

            if(k.size() == 0)
                throw FITSexception(BAD_DIMEN,"DST","erosion","Kernel size is zero.["+std::to_string(__LINE__)+"]");

            switch(in->GetEqBitPerPixel())
            {
                case SBYTE_IMG:
                    return erosion<T, int8_t> (std::dynamic_pointer_cast< FITSimg< int8_t> >(in), k, nThreads);
                case BYTE_IMG:
                    return erosion<T, uint8_t>(std::dynamic_pointer_cast< FITSimg< uint8_t> >(in), k, nThreads);
                case SHORT_IMG:
                    return erosion<T,int16_t> (std::dynamic_pointer_cast< FITSimg<int16_t> >(in), k, nThreads);
                case USHORT_IMG:
                    return erosion<T,uint16_t>(std::dynamic_pointer_cast< FITSimg<uint16_t> >(in), k, nThreads);
                case LONG_IMG:
                    return erosion<T,int32_t> (std::dynamic_pointer_cast< FITSimg<int32_t> >(in), k, nThreads);
                case ULONG_IMG:
                    return erosion<T,uint32_t>(std::dynamic_pointer_cast< FITSimg<uint32_t> >(in), k, nThreads);
                case LONGLONG_IMG:
                    return erosion<T,int64_t> (std::dynamic_pointer_cast< FITSimg<int64_t> >(in), k, nThreads);
                case ULONGLONG_IMG:
                    return erosion<T,uint64_t>(std::dynamic_pointer_cast< FITSimg<uint64_t> >(in), k, nThreads);
                case FLOAT_IMG:
                    return erosion<T,float>   (std::dynamic_pointer_cast< FITSimg<float> >(in), k, nThreads);
                case DOUBLE_IMG:
                    return erosion<T,double>  (std::dynamic_pointer_cast< FITSimg<double> >(in), k, nThreads);
                default:
                    break;
            }

            throw FITSexception(BAD_DIMEN,"DST","erosion","Unsupported pixel type.["+std::to_string(__LINE__)+"]");
        }
#pragma endregion
#pragma region - opening function

        /**
         * @brief Open an image with a given structure element kernel
         * 
         * @param in Input image
         * @param k Structure Kernel element
         * @return Pointer to a \c DSL::FITScube image representation of the eroded image
         * @throw FITSexception if input image pointer is null, kernel size is zero, or kernel/image dimension mismatch
         */
        template <typename T, typename K>
        std::shared_ptr<FITScube> opening(const std::shared_ptr< FITSimg<K> >& in, const kernel<T>& k)
        {
            return opening<T,K>(in,k,1);
        }

        /**
         * @brief Open an image with a given structure element kernel
         * 
         * @param in Input image
         * @param k Structure Kernel element
         * @param nThreads Number of threads to use for the operation
         * @return Pointer to a \c DSL::FITScube image representation of the eroded image
         * @throw FITSexception if input image pointer is null, kernel size is zero, or kernel/image dimension mismatch
         */
        template <typename T, typename K>
        std::shared_ptr<FITScube> opening(const std::shared_ptr< FITSimg<K> >& in, const kernel<T>& k, size_t nThreads)
        {
            std::shared_ptr<FITScube> eroded = erosion<T,K>(in,k,nThreads);
            std::shared_ptr<FITScube> opened = dilation<T,K>(std::dynamic_pointer_cast< FITSimg<K> >(eroded),k,nThreads);

            opened->HDU().ValueForKey("DST_PROC","OPENING",key_type::fChar,"Processing history keyword");
            opened->HDU().ValueForKey("DST_VER",DST::gGIT::version_short(),"DST framwork version library");
            opened->HDU().ValueForKey("DST_SHA",DST::gGIT::sha1(),"DST framwork GIT identifier");

            eroded.reset();

            return opened;
        }

        /**
         * @brief Open an image with a given kernel
         * 
         * @param in Input image
         * @param k Kernel element
         * @return Pointer to a \c DSL::FITScube image representation of the dilated image
         */
        template <typename T>
        std::shared_ptr<FITScube> opening(const std::shared_ptr<FITScube>& in, const kernel<T>& k)
        {
            return opening<T>(in,k,1);
        }

        /**
         * @brief Open an image with a given kernel
         * 
         * @param in Input image
         * @param k Kernel element
         * @param nThreads Number of threads to use for the operation
         * @return Pointer to a \c DSL::FITScube image representation of the dilated image
         * @throw FITSexception if input image pointer is null, kernel size is zero, or kernel/image dimension mismatch
         */
        template <typename T>
        std::shared_ptr<FITScube> opening(const std::shared_ptr<FITScube>& in, const kernel<T>& k, size_t nThreads)
        {
            if(in == nullptr)
                throw FITSexception(BAD_DIMEN,"DST","opening","Input FITScube pointer is null.["+std::to_string(__LINE__)+"]");

            if(k.size() == 0)
                throw FITSexception(BAD_DIMEN,"DST","opening","Kernel size is zero.["+std::to_string(__LINE__)+"]");

            switch(in->GetEqBitPerPixel())
            {
                case SBYTE_IMG:
                    return opening<T, int8_t> (std::dynamic_pointer_cast< FITSimg< int8_t> >(in), k, nThreads);
                case BYTE_IMG:
                    return opening<T, uint8_t>(std::dynamic_pointer_cast< FITSimg< uint8_t> >(in), k, nThreads);
                case SHORT_IMG:
                    return opening<T,int16_t> (std::dynamic_pointer_cast< FITSimg<int16_t> >(in), k, nThreads);
                case USHORT_IMG:
                    return opening<T,uint16_t>(std::dynamic_pointer_cast< FITSimg<uint16_t> >(in), k, nThreads);
                case LONG_IMG:
                    return opening<T,int32_t> (std::dynamic_pointer_cast< FITSimg<int32_t> >(in), k, nThreads);
                case ULONG_IMG:
                    return opening<T,uint32_t>(std::dynamic_pointer_cast< FITSimg<uint32_t> >(in), k, nThreads);
                case LONGLONG_IMG:
                    return opening<T,int64_t> (std::dynamic_pointer_cast< FITSimg<int64_t> >(in), k, nThreads);
                case ULONGLONG_IMG:
                    return opening<T,uint64_t>(std::dynamic_pointer_cast< FITSimg<uint64_t> >(in), k, nThreads);
                case FLOAT_IMG:
                    return opening<T,float>   (std::dynamic_pointer_cast< FITSimg<float> >(in), k, nThreads);
                case DOUBLE_IMG:
                    return opening<T,double>  (std::dynamic_pointer_cast< FITSimg<double> >(in), k, nThreads);
                default:
                    break;
            }

            throw FITSexception(BAD_DIMEN,"DST","opening","Unsupported pixel type.["+std::to_string(__LINE__)+"]");
        }
#pragma endregion
#pragma region - closing function

        /**
         * @brief Close an image with a given structure element kernel
         * 
         * @param in Input image
         * @param k Structure Kernel element
         * @return Pointer to a \c DSL::FITScube image representation of the eroded image
         * @throw FITSexception if input image pointer is null, kernel size is zero, or kernel/image dimension mismatch
         */
        template <typename T, typename K>
        std::shared_ptr<FITScube> closing(const std::shared_ptr< FITSimg<K> >& in, const kernel<T>& k)
        {
            return closing<T,K>(in,k,1);
        }

        /**
         * @brief Close an image with a given structure element kernel
         * 
         * @param in Input image
         * @param k Structure Kernel element
         * @param nThreads Number of threads to use for the operation
         * @return Pointer to a \c DSL::FITScube image representation of the eroded image
         * @throw FITSexception if input image pointer is null, kernel size is zero, or kernel/image dimension mismatch
         */
        template <typename T, typename K>
        std::shared_ptr<FITScube> closing(const std::shared_ptr< FITSimg<K> >& in, const kernel<T>& k, size_t nThreads)
        {
            std::shared_ptr<FITScube> dilated = dilation<T,K>(in,k,nThreads);
            std::shared_ptr<FITScube> closed = erosion<T,K>(std::dynamic_pointer_cast< FITSimg<K> >(dilated),k,nThreads);

            closed->HDU().ValueForKey("DST_PROC","CLOSING",key_type::fChar,"Processing history keyword");
            closed->HDU().ValueForKey("DST_VER",DST::gGIT::version_short(),key_type::fChar,"DST framwork version library");
            closed->HDU().ValueForKey("DST_SHA",DST::gGIT::sha1(),key_type::fChar,"DST framwork GIT identifier");

            dilated.reset();

            return closed;
        }

        /**
         * @brief Close an image with a given kernel
         * 
         * @param in Input image
         * @param k Kernel element
         * @return Pointer to a \c DSL::FITScube image representation of the dilated image
         */
        template <typename T>
        std::shared_ptr<FITScube> closing(const std::shared_ptr<FITScube>& in, const kernel<T>& k)
        {
            return closing<T>(in,k,1);
        }

        /**
         * @brief Close an image with a given kernel
         * 
         * @param in Input image
         * @param k Kernel element
         * @param nThreads Number of threads to use for the operation
         * @return Pointer to a \c DSL::FITScube image representation of the dilated image
         * @throw FITSexception if input image pointer is null, kernel size is zero, or kernel/image dimension mismatch
         */
        template <typename T>
        std::shared_ptr<FITScube> closing(const std::shared_ptr<FITScube>& in, const kernel<T>& k, size_t nThreads)
        {
            if(in == nullptr)
                throw FITSexception(BAD_DIMEN,"DST","closing","Input FITScube pointer is null.["+std::to_string(__LINE__)+"]");

            if(k.size() == 0)
                throw FITSexception(BAD_DIMEN,"DST","closing","Kernel size is zero.["+std::to_string(__LINE__)+"]");

            switch(in->GetEqBitPerPixel())
            {
                case SBYTE_IMG:
                    return closing<T, int8_t> (std::dynamic_pointer_cast< FITSimg< int8_t> >(in), k, nThreads);
                case BYTE_IMG:
                    return closing<T, uint8_t>(std::dynamic_pointer_cast< FITSimg< uint8_t> >(in), k, nThreads);
                case SHORT_IMG:
                    return closing<T,int16_t> (std::dynamic_pointer_cast< FITSimg<int16_t> >(in), k, nThreads);
                case USHORT_IMG:
                    return closing<T,uint16_t>(std::dynamic_pointer_cast< FITSimg<uint16_t> >(in), k, nThreads);
                case LONG_IMG:
                    return closing<T,int32_t> (std::dynamic_pointer_cast< FITSimg<int32_t> >(in), k, nThreads);
                case ULONG_IMG:
                    return closing<T,uint32_t>(std::dynamic_pointer_cast< FITSimg<uint32_t> >(in), k, nThreads);
                case LONGLONG_IMG:
                    return closing<T,int64_t> (std::dynamic_pointer_cast< FITSimg<int64_t> >(in), k, nThreads);
                case ULONGLONG_IMG:
                    return closing<T,uint64_t>(std::dynamic_pointer_cast< FITSimg<uint64_t> >(in), k, nThreads);
                case FLOAT_IMG:
                    return closing<T,float>   (std::dynamic_pointer_cast< FITSimg<float> >(in), k, nThreads);
                case DOUBLE_IMG:
                    return closing<T,double>  (std::dynamic_pointer_cast< FITSimg<double> >(in), k, nThreads);
                default:
                    break;
            }

            throw FITSexception(BAD_DIMEN,"DST","closing","Unsupported pixel type.["+std::to_string(__LINE__)+"]");
        }
#pragma endregion
#pragma region - white top-hat transformation

        /**
         * @brief Perform a white top-hat transformation on an image with a given structure element kernel
         * 
         * @param in Input image
         * @param k Structure Kernel element
         * @return Pointer to a \c DSL::FITScube image representation of the eroded image
         * @throw FITSexception if input image pointer is null, kernel size is zero, or kernel/image dimension mismatch
         */
        template <typename T, typename K>
        std::shared_ptr<FITScube> wTopHat(const std::shared_ptr< FITSimg<K> >& in, const kernel<T>& k)
        {
            return wTopHat<T,K>(in,k,1);
        }

        /**
         * @brief Perform a white top-hat transformation on an image with a given structure element kernel
         * 
         * @param in Input image
         * @param k Structure Kernel element
         * @param nThreads Number of threads to use for the operation
         * @return Pointer to a \c DSL::FITScube image representation of the eroded image
         * @throw FITSexception if input image pointer is null, kernel size is zero, or kernel/image dimension mismatch
         */
        template <typename T, typename K>
        std::shared_ptr<FITScube> wTopHat(const std::shared_ptr< FITSimg<K> >& in, const kernel<T>& k, size_t nThreads)
        {
            // Copy the original image
            std::shared_ptr<FITScube> transformed = std::shared_ptr<FITScube>( new FITSimg<K>(*(in.get())) );
            std::shared_ptr<FITScube> opened = opening<T,K>(in,k,nThreads);

            *(transformed.get()) -= *(opened.get());

            transformed->HDU().ValueForKey("DST_PROC","White Top-Hat",key_type::fChar,"Processing history keyword");
            transformed->HDU().ValueForKey("DST_VER",DST::gGIT::version_short(),key_type::fChar,"DST framwork version library");
            transformed->HDU().ValueForKey("DST_SHA",DST::gGIT::sha1(),key_type::fChar,"DST framwork GIT identifier");

            opened.reset();

            return transformed;
        }

        /**
         * @brief Perform a white top-hat transformation on an image with a given structure element kernel
         * 
         * @param in Input image
         * @param k Kernel element
         * @return Pointer to a \c DSL::FITScube image representation of the dilated image
         */
        template <typename T>
        std::shared_ptr<FITScube> wTopHat(const std::shared_ptr<FITScube>& in, const kernel<T>& k)
        {
            return wTopHat<T>(in,k,1);
        }

        /**
         * @brief Perform a white top-hat transformation on an image with a given structure element kernel
         * 
         * @param in Input image
         * @param k Kernel element
         * @param nThreads Number of threads to use for the operation
         * @return Pointer to a \c DSL::FITScube image representation of the dilated image
         * @throw FITSexception if input image pointer is null, kernel size is zero, or kernel/image dimension mismatch
         */
        template <typename T>
        std::shared_ptr<FITScube> wTopHat(const std::shared_ptr<FITScube>& in, const kernel<T>& k, size_t nThreads)
        {
            if(in == nullptr)
                throw FITSexception(BAD_DIMEN,"DST","wTopHat","Input FITScube pointer is null.["+std::to_string(__LINE__)+"]");

            if(k.size() == 0)
                throw FITSexception(BAD_DIMEN,"DST","wTopHat","Kernel size is zero.["+std::to_string(__LINE__)+"]");

            switch(in->GetEqBitPerPixel())
            {
                case SBYTE_IMG:
                    return wTopHat<T, int8_t> (std::dynamic_pointer_cast< FITSimg< int8_t> >(in), k, nThreads);
                case BYTE_IMG:
                    return wTopHat<T, uint8_t>(std::dynamic_pointer_cast< FITSimg< uint8_t> >(in), k, nThreads);
                case SHORT_IMG:
                    return wTopHat<T,int16_t> (std::dynamic_pointer_cast< FITSimg<int16_t> >(in), k, nThreads);
                case USHORT_IMG:
                    return wTopHat<T,uint16_t>(std::dynamic_pointer_cast< FITSimg<uint16_t> >(in), k, nThreads);
                case LONG_IMG:
                    return wTopHat<T,int32_t> (std::dynamic_pointer_cast< FITSimg<int32_t> >(in), k, nThreads);
                case ULONG_IMG:
                    return wTopHat<T,uint32_t>(std::dynamic_pointer_cast< FITSimg<uint32_t> >(in), k, nThreads);
                case LONGLONG_IMG:
                    return wTopHat<T,int64_t> (std::dynamic_pointer_cast< FITSimg<int64_t> >(in), k, nThreads);
                case ULONGLONG_IMG:
                    return wTopHat<T,uint64_t>(std::dynamic_pointer_cast< FITSimg<uint64_t> >(in), k, nThreads);
                case FLOAT_IMG:
                    return wTopHat<T,float>   (std::dynamic_pointer_cast< FITSimg<float> >(in), k, nThreads);
                case DOUBLE_IMG:
                    return wTopHat<T,double>  (std::dynamic_pointer_cast< FITSimg<double> >(in), k, nThreads);
                default:
                    break;
            }

            throw FITSexception(BAD_DIMEN,"DST","wTopHat","Unsupported pixel type.["+std::to_string(__LINE__)+"]");
        }
#pragma endregion
#pragma region - black top-hat transformation

        /**
         * @brief Perform a black top-hat transformation on an image with a given structure element kernel
         * 
         * @param in Input image
         * @param k Structure Kernel element
         * @return Pointer to a \c DSL::FITScube image representation of the eroded image
         * @throw FITSexception if input image pointer is null, kernel size is zero, or kernel/image dimension mismatch
         */
        template <typename T, typename K>
        std::shared_ptr<FITScube> bTopHat(const std::shared_ptr< FITSimg<K> >& in, const kernel<T>& k)
        {
            return bTopHat<T,K>(in,k,1);
        }

        /**
         * @brief Perform a black top-hat transformation on an image with a given structure element kernel
         * 
         * @param in Input image
         * @param k Structure Kernel element
         * @param nThreads Number of threads to use for the operation
         * @return Pointer to a \c DSL::FITScube image representation of the eroded image
         * @throw FITSexception if input image pointer is null, kernel size is zero, or kernel/image dimension mismatch
         */
        template <typename T, typename K>
        std::shared_ptr<FITScube> bTopHat(const std::shared_ptr< FITSimg<K> >& in, const kernel<T>& k, size_t nThreads)
        {
            // Copy the original image
            std::shared_ptr<FITScube> transformed = closing<T,K>(in,k,nThreads);
            *(transformed.get()) -= *(in.get());

            transformed->HDU().ValueForKey("DST_PROC","black Top-Hat",key_type::fChar,"Processing history keyword");
            transformed->HDU().ValueForKey("DST_VER",DST::gGIT::version_short(),key_type::fChar,"DST framwork version library");
            transformed->HDU().ValueForKey("DST_SHA",DST::gGIT::sha1(),key_type::fChar,"DST framwork GIT identifier");

            return transformed;
        }

        /**
         * @brief Perform a black top-hat transformation on an image with a given structure element kernel
         * 
         * @param in Input image
         * @param k Kernel element
         * @return Pointer to a \c DSL::FITScube image representation of the dilated image
         */
        template <typename T>
        std::shared_ptr<FITScube> bTopHat(const std::shared_ptr<FITScube>& in, const kernel<T>& k)
        {
            return bTopHat<T>(in,k,1);
        }

        /**
         * @brief Perform a black top-hat transformation on an image with a given structure element kernel
         * 
         * @param in Input image
         * @param k Kernel element
         * @param nThreads Number of threads to use for the operation
         * @return Pointer to a \c DSL::FITScube image representation of the dilated image
         * @throw FITSexception if input image pointer is null, kernel size is zero, or kernel/image dimension mismatch
         */
        template <typename T>
        std::shared_ptr<FITScube> bTopHat(const std::shared_ptr<FITScube>& in, const kernel<T>& k, size_t nThreads)
        {
            if(in == nullptr)
                throw FITSexception(BAD_DIMEN,"DST","bTopHat","Input FITScube pointer is null.["+std::to_string(__LINE__)+"]");

            if(k.size() == 0)
                throw FITSexception(BAD_DIMEN,"DST","bTopHat","Kernel size is zero.["+std::to_string(__LINE__)+"]");

            switch(in->GetEqBitPerPixel())
            {
                case SBYTE_IMG:
                    return bTopHat<T, int8_t> (std::dynamic_pointer_cast< FITSimg< int8_t> >(in), k, nThreads);
                case BYTE_IMG:
                    return bTopHat<T, uint8_t>(std::dynamic_pointer_cast< FITSimg< uint8_t> >(in), k, nThreads);
                case SHORT_IMG:
                    return bTopHat<T,int16_t> (std::dynamic_pointer_cast< FITSimg<int16_t> >(in), k, nThreads);
                case USHORT_IMG:
                    return bTopHat<T,uint16_t>(std::dynamic_pointer_cast< FITSimg<uint16_t> >(in), k, nThreads);
                case LONG_IMG:
                    return bTopHat<T,int32_t> (std::dynamic_pointer_cast< FITSimg<int32_t> >(in), k, nThreads);
                case ULONG_IMG:
                    return bTopHat<T,uint32_t>(std::dynamic_pointer_cast< FITSimg<uint32_t> >(in), k, nThreads);
                case LONGLONG_IMG:
                    return bTopHat<T,int64_t> (std::dynamic_pointer_cast< FITSimg<int64_t> >(in), k, nThreads);
                case ULONGLONG_IMG:
                    return bTopHat<T,uint64_t>(std::dynamic_pointer_cast< FITSimg<uint64_t> >(in), k, nThreads);
                case FLOAT_IMG:
                    return bTopHat<T,float>   (std::dynamic_pointer_cast< FITSimg<float> >(in), k, nThreads);
                case DOUBLE_IMG:
                    return bTopHat<T,double>  (std::dynamic_pointer_cast< FITSimg<double> >(in), k, nThreads);
                default:
                    break;
            }

            throw FITSexception(BAD_DIMEN,"DST","bTopHat","Unsupported pixel type.["+std::to_string(__LINE__)+"]");
        }
#pragma endregion
#pragma region - gradiant

        /**
         * @brief Perform gradiant transformation on an image with a given structure element kernel
         * 
         * @param in Input image
         * @param k Structure Kernel element
         * @return Pointer to a \c DSL::FITScube image representation of the eroded image
         * @throw FITSexception if input image pointer is null, kernel size is zero, or kernel/image dimension mismatch
         */
        template <typename T, typename K>
        std::shared_ptr<FITScube> gradiant(const std::shared_ptr< FITSimg<K> >& in, const kernel<T>& k)
        {
            return gradiant<T,K>(in,k,1);
        }

        /**
         * @brief Perform gradiant transformation on an image with a given structure element kernel
         * 
         * @param in Input image
         * @param k Structure Kernel element
         * @param nThreads Number of threads to use for the operation
         * @return Pointer to a \c DSL::FITScube image representation of the eroded image
         * @throw FITSexception if input image pointer is null, kernel size is zero, or kernel/image dimension mismatch
         */
        template <typename T, typename K>
        std::shared_ptr<FITScube> gradiant(const std::shared_ptr< FITSimg<K> >& in, const kernel<T>& k, size_t nThreads)
        {
            // Launch dilation and erosion in parallel to utilize multiple cores.
            // Each operation still receives the requested nThreads; adjust if you
            // want to avoid oversubscription (e.g. pass nThreads/2 to each).
            auto fut_dil = std::async(std::launch::async, [&](){ return dilation<T,K>(in, k, nThreads); });
            auto fut_ero = std::async(std::launch::async, [&](){ return erosion<T,K>(in, k, nThreads); });

            // Wait for results (exceptions from tasks will propagate here)
            std::shared_ptr<FITScube> dilated = fut_dil.get();
            std::shared_ptr<FITScube> eroded  = fut_ero.get();

            if(!dilated)
                throw FITSexception(BAD_DIMEN,"DST","gradiant","Dilation returned null.["+std::to_string(__LINE__)+"]");
            if(!eroded)
                throw FITSexception(BAD_DIMEN,"DST","gradiant","Erosion returned null.["+std::to_string(__LINE__)+"]");

            // Compute gradient = dilated - eroded
            *(dilated.get()) -= *(eroded.get());

            dilated->HDU().ValueForKey("DST_PROC","GRADIANT",key_type::fChar,"Processing history keyword");
            dilated->HDU().ValueForKey("DST_VER",DST::gGIT::version_short(),key_type::fChar,"DST framwork version library");
            dilated->HDU().ValueForKey("DST_SHA",DST::gGIT::sha1(),key_type::fChar,"DST framwork GIT identifier");

            eroded.reset();
            
            return dilated;
        }

        /**
         * @brief Perform gradiant transformation on an image with a given structure element kernel
         * 
         * @param in Input image
         * @param k Kernel element
         * @return Pointer to a \c DSL::FITScube image representation of the dilated image
         */
        template <typename T>
        std::shared_ptr<FITScube> gradiant(const std::shared_ptr<FITScube>& in, const kernel<T>& k)
        {
            return gradiant<T>(in,k,1);
        }

        /**
         * @brief Perform gradiant transformation on an image with a given structure element kernel
         * 
         * @param in Input image
         * @param k Kernel element
         * @param nThreads Number of threads to use for the operation
         * @return Pointer to a \c DSL::FITScube image representation of the dilated image
         * @throw FITSexception if input image pointer is null, kernel size is zero, or kernel/image dimension mismatch
         */
        template <typename T>
        std::shared_ptr<FITScube> gradiant(const std::shared_ptr<FITScube>& in, const kernel<T>& k, size_t nThreads)
        {
            if(in == nullptr)
                throw FITSexception(BAD_DIMEN,"DST","gradiant","Input FITScube pointer is null.["+std::to_string(__LINE__)+"]");

            if(k.size() == 0)
                throw FITSexception(BAD_DIMEN,"DST","gradiant","Kernel size is zero.["+std::to_string(__LINE__)+"]");

            switch(in->GetEqBitPerPixel())
            {
                case SBYTE_IMG:
                    return gradiant<T, int8_t> (std::dynamic_pointer_cast< FITSimg< int8_t> >(in), k, nThreads);
                case BYTE_IMG:
                    return gradiant<T, uint8_t>(std::dynamic_pointer_cast< FITSimg< uint8_t> >(in), k, nThreads);
                case SHORT_IMG:
                    return gradiant<T,int16_t> (std::dynamic_pointer_cast< FITSimg<int16_t> >(in), k, nThreads);
                case USHORT_IMG:
                    return gradiant<T,uint16_t>(std::dynamic_pointer_cast< FITSimg<uint16_t> >(in), k, nThreads);
                case LONG_IMG:
                    return gradiant<T,int32_t> (std::dynamic_pointer_cast< FITSimg<int32_t> >(in), k, nThreads);
                case ULONG_IMG:
                    return gradiant<T,uint32_t>(std::dynamic_pointer_cast< FITSimg<uint32_t> >(in), k, nThreads);
                case LONGLONG_IMG:
                    return gradiant<T,int64_t> (std::dynamic_pointer_cast< FITSimg<int64_t> >(in), k, nThreads);
                case ULONGLONG_IMG:
                    return gradiant<T,uint64_t>(std::dynamic_pointer_cast< FITSimg<uint64_t> >(in), k, nThreads);
                case FLOAT_IMG:
                    return gradiant<T,float>   (std::dynamic_pointer_cast< FITSimg<float> >(in), k, nThreads);
                case DOUBLE_IMG:
                    return gradiant<T,double>  (std::dynamic_pointer_cast< FITSimg<double> >(in), k, nThreads);
                default:
                    break;
            }

            throw FITSexception(BAD_DIMEN,"DST","gradiant","Unsupported pixel type.["+std::to_string(__LINE__)+"]");
        }
#pragma endregion
    }

}

#endif /* defined(__DST_morpho__) */
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
#include <thread>                 // added for std::thread

#include <DSTfits/FITSdata.h>
#include <DSTfits/FITSexception.h>
#include <DSTfits/FITSimg.h>
#include <DSTfits/FITShdu.h>

#include <DSTfits/DSF_version.h>

using namespace DSL;

namespace DST
{

#pragma region - kernel class
    template <typename T>
    class kernel
    {
#pragma region • members
        private:
            std::valarray<T> fkernel;
            std::vector<size_t> fsize;
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
            kernel(const std::valarray<T>& k, const std::vector<size_t>& s):fkernel(k), fsize(s)
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
            kernel(const std::valarray<T>& k, const std::initializer_list<size_t>& s):kernel(k, std::vector<size_t>(s))
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
            kernel(const kernel<T>& k):fkernel(k.fkernel), fsize(k.fsize)
            {
                validityCheck();
            }

            /**
             * @brief destructor
             */
            ~kernel(){fsize.clear();}
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
            inline const size_t& size(const size_t& idx) const { if(idx >= fsize.size()) throw std::out_of_range("Index out of range.["+std::to_string(__LINE__)+"]"); return fsize[idx]; };

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

    // Strongly typed overload (preferred): K is known from FITSimg<K>
    // Strongly typed overload (preferred): K is known from FITSimg<K>
    template <typename T, typename K>
    std::shared_ptr<FITScube> folding(const std::shared_ptr< FITSimg<K> >& in, const kernel<T>& k)
    {
        return folding<T,K>(in,k,1);
    }

    template <typename T, typename K>
    std::shared_ptr<FITScube> folding(const std::shared_ptr< FITSimg<K> >& in, const kernel<T>& k, size_t nThreads)
    {
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
                // First step is to get slice of the image which correspond to the kernel position. When the kernel is near the border of the image, only portion of the kernel is used and the slice of the image should be aligned with the kernel center (i.e. tmp array should be null when outside image bounds)
                
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
                        tmp[kj] = U(0);
                        continue;
                    }

                    // get linear pixel index corresponding to tgtCoords
                    size_t tgtIndex = in->PixelIndex(tgtCoords);

                    // fetch value from input data and promote to U
                    tmp[kj] = static_cast<U>((*in_data)[tgtIndex]);
                }

                // perform folding operation
                out_data[i] = (tmp * tmpKernel).sum();
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
 
                workers.emplace_back([start, end, &out_data, in_data_ptr, &k, &tmpKernel, &kSizes, kDim, imgDim, in_shared, &kCenter]()
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
                                tmp[kj] = U(0);
                                continue;
                            }

                            size_t tgtIndex = in_shared->PixelIndex(tgtCoords);
                            tmp[kj] = static_cast<U>((*in_data_ptr)[tgtIndex]);
                        }

                        out_data[i] = (tmp * tmpKernel).sum();
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
        
        return std::shared_ptr<FITScube>(outImg);
    }

    // Keep the original signature; fix the kernel size check. You can make this a wrapper.
    template <typename T>
    std::shared_ptr<FITScube> folding(const std::shared_ptr<FITScube>& in, const kernel<T>& k)
    {
        return folding<T>(in,k,1);
    }

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

    template <typename T>
    std::shared_ptr<FITScube> unfolding(const std::shared_ptr<FITScube>& in, const kernel<T>& k)
    {
        if(in == nullptr)
            throw FITSexception(BAD_DIMEN,"DST","folding","Input FITScube pointer is null.["+std::to_string(__LINE__)+"]");

        if(k.size() == 0)
            throw FITSexception(BAD_DIMEN,"DST","folding","Kernel size is zero.["+std::to_string(__LINE__)+"]");

        return nullptr;
    }
}

#endif /* defined(__DST_morpho__) */